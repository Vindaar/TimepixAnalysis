import std / [sequtils, algorithm]
import pkg / seqmath

proc toEdf*(x: seq[float], bins: seq[float]): seq[float] =
  ## Computes the EDF of the input data `x` given some `bins`.
  ##
  ## The bins are used as boundaries to count elements below each bin edge. The
  ## result contains `bins.len` elements, all in [0, 1]
  let xS = x.sorted
  var i = 0
  result = newSeq[float](bins.len)
  for j, b in bins:
    while i < xS.high and xS[i] <= b:
      inc i
    result[j] = i.float / x.len.float
  #for x in mitems(result):
  #  x = x / result[^1]
  doAssert result.allIt(it <= 1.0)

proc kolmogorovSmirnovBinned*(x1, x2: seq[float]): float =
  ## Compute the Kolmogorov-Smirnov test for two datasets.
  ##
  ## The data is binned first to min & max of the combined data range and based on the
  ## associated EDF the KS test is performed.
  ##
  ## ``KS(x) = sup_x | EDF₁(x) - EDF₂(x) |``
  ##
  ## or in ~code
  ##
  ## ``KS(x) = max(abs(EDF₁ -. EDF₂(x)))``
  let range = (min: min(x1.min, x2.min), max: max(x1.max, x2.max))
  let bins = linspace(range[0], range[1], 200)
  let x1Edf = x1.toEdf(bins)
  let x2Edf = x2.toEdf(bins)
  result = 0.0
  for i, b in bins:
    result = max( result, abs(x1Edf[i] - x2Edf[i]) )

proc kolmogorovSmirnovPVal*(n1, n2: int, ksDiff: float): float =
  ## Computes the p-Value of the Kolmogorov-Smirnov test, assuming large samples.
  ## For large samples the KS statistic converges to the Kolmogorov distribution:
  ##
  ## `P(K \leq x) = 1 - 2 \sum_{k=1}^∞ (-1)^{k-1} e^{-2k² x²}`
  ##
  ## with the p-Value:
  ##
  ## `p = 2 \sum_{k=1}^∞ (-1)^{k-1} e^{-2k² λ²}`
  ## where `λ = \sqrt{\frac{n m}{n+m}} D_{\text{observed}}`.
  ##
  ## We truncate the sum after 100 terms or when the contribution becomes small.
  result = 0.0
  let λ = sqrt((n1 * n2) / (n1 + n2)) * ksDiff
  for k in 1 .. 100:
    let kf = k.float
    let term = pow(-1.0, kf - 1.0) * exp(-2.0 * kf * kf * λ * λ)
    result += term
    if abs(term) < 1e-10:  # Convergence check
      break
  result *= 2.0
  doAssert result >= 0.0 and result <= 1.0, "Calculated p-Value outside [0, 1]: " & $result

proc kolmogorovSmirnov*(x1, x2: seq[float]): tuple[statistic: float, pValue: float] =
  ## Compute the Kolmogorov-Smirnov test for two datasets.
  ##
  ## The KS test is performed on the unbinned data using an unbinned EDF.
  ##
  ## ``KS(x) = sup_x | EDF₁(x) - EDF₂(x) |``
  ##
  ## or in ~code
  ##
  ## ``KS(x) = max(abs(EDF₁ -. EDF₂(x)))``
  ##
  ## Returns the KS statistic and the p-value.
  let all = concat(x1, x2).sorted
  let x1 = x1.sorted
  let x2 = x2.sorted
  let m = x1.len.float
  let n = x2.len.float
  var KS = 0.0
  for x in all:
    let c1 = (x1.upperBound(x)).float / m # first idx > x, hence # of elements less or equal
    let c2 = (x2.upperBound(x)).float / n
    let diff = abs(c1 - c2)
    KS = max(KS, diff)
  result = (KS, kolmogorovSmirnovPVal(m.int, n.int, KS))

proc cramerVonMises*(xIn, yIn: seq[float]): float =
  ## Computes the Cramér-von Mises criterion for two samples x and y. The
  ## two samples do not need to be of the same size.
  ##
  ## NOTE: This is a Bing Chat implementation. It directly computes the
  ## T value by buildin the empirical cumulative distribution functions.
  ## (Yes, really. I was lazy :D)
  let x = xIn.sorted
  let y = yIn.sorted
  let m = x.len # size of first sample
  let n = y.len # size of second sample
  let N = m + n # size of combined sample
  let dHN = 1.0 / N.float # change in combined ECDF
  var i = 0 # index for x
  var j = 0 # index for y
  var Fm = 0.0 # empirical CDF of x
  var Gn = 0.0 # empirical CDF of y
  result = 0.0 # Cramér-von Mises criterion
  while i < m and j < n: # loop until one sample is exhausted
    if x[i] <= y[j]: # next value is from x
      Fm += 1.0 / m.float # update Fm
      inc i # increment i
    else: # next value is from y
      Gn += 1.0 / n.float # update Gn
      inc j # increment j
    result += (Fm - Gn) ^ 2 * dHN # update result using change in HN
  for k in i ..< m: # loop until x is exhausted
    Fm += 1.0 / m.float # update Fm
    result += (Fm - Gn) ^ 2 * dHN # update result using change in HN
  for k in j ..< n: # loop until y is exhausted
    Gn += 1.0 / n.float # update Gn
    result += (Fm - Gn) ^ 2 * dHN # update result using change in HN
  result *= m * n / N # scale result by sample sizes and combined size

when isMainModule:
  import random, strformat

  proc testKS() =
    # Test with samples from same distribution (should have high p-value)
    randomize(42)
    var sample1: seq[float] = @[]
    var sample2: seq[float] = @[]

    for i in 0..<100:
      sample1.add(rand(1.0))  # Uniform [0,1]
      sample2.add(rand(1.0))  # Uniform [0,1]

    let result1 = kolmogorovSmirnov(sample1, sample2)
    echo &"Same distribution test:"
    echo &"  KS statistic: {result1.statistic:.4f}"
    echo &"  p-value: {result1.pValue:.4f}"
    echo ""

    # Test with samples from different distributions (should have low p-value)
    var sample3: seq[float] = @[]
    for i in 0..<100:
      sample3.add(rand(1.0) + 0.5)  # Uniform [0.5, 1.5]

    let result2 = kolmogorovSmirnov(sample1, sample3)
    echo &"Different distribution test:"
    echo &"  KS statistic: {result2.statistic:.4f}"
    echo &"  p-value: {result2.pValue:.4f}"
    echo ""

    # Test with different sample sizes
    var smallSample: seq[float] = @[]
    for i in 0..<20:
      smallSample.add(rand(1.0))

    let result3 = kolmogorovSmirnov(smallSample, sample1)
    echo &"Different sample sizes test (n1=20, n2=100):"
    echo &"  KS statistic: {result3.statistic:.4f}"
    echo &"  p-value: {result3.pValue:.4f}"

  testKS()
