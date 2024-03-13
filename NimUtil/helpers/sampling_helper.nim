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
  # Start from 0! We calc from PDF which means our 0 point is at _beginning_. We don't
  # want samples at x = 0
  # Essentially the difference between starting from *binned data* and a pdf. In case of
  # binned data cumSum of 0 contains entries in first bin. But PDF is more equal to having
  # the _left bin edge_ without any content
  let baseline = dataCdf[0] # 0.0 # min(data) # 0.0
  #echo min(data)
  dataCdf.mapIt((it - baseline) / (integral - baseline))

proc sample*(rnd: var Rand, cdf: seq[float], xs, ys: seq[float]): float =
  let u = rnd.rand(1.0)
  let idx = cdf.lowerBound(u)
  #echo idx, " from u ", u
  if idx == 0:
    result = xs[idx]
  elif idx < cdf.len:
    if u == cdf[idx]:
      result = xs[idx]
    else: # linear interpolation, point between previous and next at `u`
      result = xs[idx-1] + (u - cdf[idx-1]) * (xs[idx] - xs[idx-1]) / (cdf[idx] - cdf[idx-1])
      #echo "Instead of ", xs[idx], " give ", result
  else:
    result = Inf

proc sample*(rnd: var Rand, s: Sampler): float =
  result = rnd.sample(s.edf, s.xs, s.ys)

proc expFn*(x: float, λ: float): float =
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
  rnd: var Rand,
  fn: proc(x: float): float,
  low, high: float,
  num = 1000,
  samples = 1_000_000
               ): seq[float] =

  let sampler = sampler(fn, low, high, num)
  result = newSeq[float](samples)
  for i in 0 ..< samples:
    result[i] = rnd.sample(sampler)

when isMainModule:
  ## Mini test: Compare with plot output from /tmp/test_sample.nim! (see `StatusAndProgress.org`)
  import pkg / ggplotnim
  let Upper = 3.0
  let λ = 2.0

  let fnSample = (proc(x: float): float =
    result = expFn(x, λ)
  )

  var rnd = initRand(1337)
  let ySampled = rnd.sampleFrom(fnSample, 0.0, 3.0)
  ggplot(toDf(ySampled), aes("ySampled")) +
    geom_histogram(bins = 100, density = true, alpha = 0.5, hdKind = hdOutline, fillColor = "red") +
    ggsave("/home/basti/Sync/test_sample_from.pdf")
