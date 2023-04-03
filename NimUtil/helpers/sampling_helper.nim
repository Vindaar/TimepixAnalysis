import std / [random, sequtils, algorithm]
import pkg / seqmath

type
  Sampler* = object
    xs*: seq[float]
    ys*: seq[float]
    edf*: seq[float]

template toEDF*(data: seq[float], isCumSum = false): untyped =
  ## Computes the EDF of binned data
  var dataCdf = data
  if not isCumSum:
    seqmath.cumsum(dataCdf)
  let integral = dataCdf[^1]
  ## XXX: why min?
  let baseline = min(data) # 0.0
  dataCdf.mapIt((it - baseline) / (integral - baseline))

proc sample*(cdf: seq[float], xs: seq[float]): float =
  let point = rand(1.0)
  let idx = cdf.lowerBound(point)
  if idx < cdf.len:
    result = xs[idx]
  else:
    result = Inf

proc sample*(s: Sampler): float =
  result = sample(s.edf, s.xs)

proc expFn(x: float, λ: float): float =
  result = 1.0 / λ * exp(- x / λ)

proc sampler*(
  fn: proc(x: float): float,
  low, high: float,
  num = 1000
            ): Sampler =
  ## Note: it may be useful to hand a closure with wrapped arguments!
  let xs = linspace(low, high, num)
  let ys = xs.mapIt( fn(it) )

  # now sample 100,000 points
  let cdf = ys.toEdf()
  result = Sampler(edf: cdf, xs: xs, ys: ys)

proc sampleFrom*(
  fn: proc(x: float): float,
  low, high: float,
  num = 1000,
  samples = 1_000_000
               ): seq[float] =

  let sampler = sampler(fn, low, high, num)
  result = newSeq[float](samples)
  for i in 0 ..< samples:
    result[i] = sample(sampler)

when isMainModule:
  ## Mini test: Compare with plot output from /tmp/test_sample.nim! (see `StatusAndProgress.org`)
  import pkg / ggplotnim
  let Upper = 3.0
  let λ = 2.0

  let fnSample = (proc(x: float): float =
    result = expFn(x, λ)
  )

  let ySampled = sampleFrom(fnSample, 0.0, 3.0)
  ggplot(toDf(ySampled), aes("ySampled")) +
    geom_histogram(bins = 100, density = true, alpha = 0.5, hdKind = hdOutline, fillColor = "red") +
    ggsave("/t/test_sample_from.pdf")
