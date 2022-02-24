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

template returnOnFalse(c1, c2: untyped, key = ""): untyped =
  result = c1 == c2
  if not result:
    echo "Didn't match ", astToStr(c1), " == ", astToStr(c2)
    echo "Was ", c1, " and ", c2
    if key.len > 0:
      echo "For keys: ", key
    return false

proc compareJson*(j1, j2: JsonNode, exceptKeys: seq[string] = @[]): bool =
  ## `exceptKeys` can be given to *not* compare the given object keys
  ## This is (currently) used for keys that depend on the compiled version
  ## of a program (e.g. the date it was compiled)
  returnOnFalse(j1.kind, j2.kind)
  case j1.kind
  of JObject:
    returnOnFalse(j1.len, j2.len)
    for k, v in pairs(j1):
      if k in exceptKeys: continue
      when not defined(linux):
        if k == "txtPos" or k == "txtText": # txtText broken for multiline
          echo "INFO: Skipping key ", k, " due to cairo differences in text " &
            "printing on different platforms"
          continue
      returnOnFalse(k in j2, true)
      returnOnFalse(compareJson(v, j2[k]), true, k)
  of JFloat:
    when defined(linux):
      let cmpFloat = almostEqual(j1.getFloat, j2.getFloat, 1e-4)
    else:
      ## TODO: due to some cairo issue related to different platforms we get different
      ## positions on mac/windows. For now we just use a much larger epsilon. Need to
      ## investigate this difference once I have access to a machine running windows again.
      let cmpFloat = abs(j1.getFloat - j2.getFloat) < 0.1 # crude manual...
    if not cmpFloat:
      echo "Float compare failed: ", j1.getFloat, " <-> ", j2.getFloat
    returnOnFalse(cmpFloat, true)
  of JArray:
    returnOnFalse(j1.len, j2.len)
    for i in 0 ..< j1.len:
      returnOnFalse(compareJson(j1[i], j2[i]), true)
  else:
    returnOnFalse(j1, j2)

func almostEq*[T: SomeNumber](x, y: T, ep = 1e-6): bool =
  ## checks very roughly if the values match. Anything beyond
  ## 1e-5 should be of no issue for us
  ## NOTE: We explicitly do *not* use the `stdlib` `almostEqual` as we don't care about
  ## any kind of smart comparison!
  result = abs((x - y).float) < ep
