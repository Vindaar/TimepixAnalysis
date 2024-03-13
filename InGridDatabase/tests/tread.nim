import typetraits
import arraymancer
import ingridDatabase

when isMainModule:
  const chipName = "E4W66"
  let t = getTotCalib(chipName)
  echo t

  let scvs = getScurveSeq(chipName)
  echo scvs

  let thl = getThreshold(chipName)
  echo thl
  
