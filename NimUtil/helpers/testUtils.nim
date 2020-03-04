import json
from ggplotnim import almostEqual
export json
export almostEqual

proc echoFields(j1, j2: JsonNode) =
  doAssert j1.len >= j2.len
  for kx, vx in pairs(j1):
    echo "Is ", kx, ", with val: ", vx, " in j2? ", kx in j2
    if kx notin j2:
      echo "^^^^^^\n"

proc compareJObjects*(j1, j2: JsonNode): bool =
  template returnOnFalse(c1, c2: untyped): untyped =
    result = c1 == c2
    if not result:
      echo "Didn't match ", astToStr(c1), " == ", astToStr(c2)
      echo "Was ", c1, " and ", c2
      return false
  returnOnFalse(j1.kind, JObject)
  returnOnFalse(j2.kind, JObject)
  if j1.len > j2.len:
    echo "j1 contains more elements than j2!"
    echoFields(j1, j2)
  elif j2.len > j1.len:
    echo "j2 contains more elements than j1!"
    echoFields(j2, j1)
  else:
    # all good, same field number
    discard
  #returnOnFalse(j1.len, j2.len)
  for k, v in pairs(j1):
    returnOnFalse(k in j2, true)
    returnOnFalse(v.kind, j2[k].kind)
    case v.kind
    of JObject:
      returnOnFalse(compareJObjects(v, j2[k]), true)
    of JFloat:
      returnOnFalse(almostEqual(v.getFloat, j2[k].getFloat), true)
    else:
      returnOnFalse(v, j2[k])
