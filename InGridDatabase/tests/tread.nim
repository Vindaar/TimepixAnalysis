import typetraits
import arraymancer
import ingridDatabase

when isMainModule:
  const chipName = "E4W66"
  let runPeriod = "Run1"
  let t = getTotCalib(chipName, runPeriod)
  echo t

  let scvs = getScurveSeq(chipName, runPeriod)
  echo scvs

  let thl = getThreshold(chipName, runPeriod)
  echo thl
