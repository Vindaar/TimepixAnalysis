import jsffi, jsbind, macros, strutils
include karax / prelude
import karax / kdom
import karax / jjson
import plot_types

when defined(client):
  const UseWs* = true
else:
  const UseWs* = false

# allows to convert `JsObject` to string
#proc toString*(x: JsObject): cstring {.jsimportgWithName: "String".}

proc parseJsonToJs(json: cstring): JsObject {.jsimportgWithName: "JSON.parse".}
proc stringify*(value: JsObject | JsonNode,
                replacer: JsObject,
                space: JsObject): cstring {.jsimportgWithName: "JSON.stringify".}
proc toString*(x: JsObject | JsonNode): cstring =
  result = x.stringify(nil, toJs(2))

proc selectedIndex*(n: Node): cint {.importcpp: "#.selectedIndex".}

proc pretty*(x: JsonNode): cstring =
  result = x.stringify(nil, toJs(2))

proc parseJson*(s: kstring): JsonNode =
  result = % parseJsonToJs(s)

proc parseEnum*[T: enum](s: cstring, default: T): T =
  result = strutils.parseEnum[T]($s, default)

proc parseFloat*(s: kstring): float =
  result = ($s).parseFloat

proc traverseTree(input: NimNode): NimNode =
  # iterate children
  for i in 0 ..< input.len:
    case input[i].kind
    of nnkSym:
      # if we found a symbol, take it
      result = input[i]
    of nnkBracketExpr:
      # has more children, traverse
      result = traverseTree(input[i])
    else:
      error("Unsupported type: " & $input.kind)

macro getInnerType*(TT: typed): untyped =
  ## macro to get the subtype of a nested type by iterating
  ## the AST
  # traverse the AST
  let res = traverseTree(TT.getTypeInst)
  # assign symbol to result
  result = quote do:
    `res`

proc to*(node: JsonNode, dtype: typedesc): dtype =
  ## simplified version to convert JS JsonNode to dtype
  ## currently only supports sequence types!
  type innerType = getInnerType(dtype)
  for x in node:
    when innerType is int:
      result.add x.getInt
    elif innerType is float:
      result.add $(x.getFNum).parseFloat
    else:
      result.add x.getStr

proc getFloat*(x: JsonNode): float =
  result = ($(x.getFNum)).parseFloat

proc getNextStatic*(pState: PlotState): kstring =
  if pState.staticP.idx + 1 < pState.staticP.keys.len:
    kstring(pState.staticP.keys[pState.staticP.idx + 1])
  else:
    kstring""

proc getPrevStatic*(pState: PlotState): kstring =
  if pState.staticP.idx > 0:
    kstring(pState.staticP.keys[pState.staticP.idx - 1])
  else:
    kstring""

func decInRangeStatic*(pState: var PlotState) {.inline.} =
  if pState.staticP.idx > 1:
    dec pState.staticP.idx
  else:
    pState.staticP.idx = 0

func incInRangeStatic*(pState: var PlotState) {.inline.} =
  if pState.staticP.idx < pState.staticP.keys.high:
    inc pState.staticP.idx
  else:
    pState.staticP.idx = pState.staticP.keys.high

when UseWs:
  proc getNextServer*(pState: PlotState): kstring =
    if pState.serverP.idx + 1 < pState.serverP.keys.len:
      kstring(pState.serverP.keys[pState.serverP.idx + 1])
    else:
      kstring""

  proc getPrevServer*(pState: PlotState): kstring =
    if pState.serverP.idx > 0:
      kstring(pState.serverP.keys[pState.serverP.idx - 1])
    else:
      kstring""

  func decInRangeServer*(pState: var PlotState) {.inline.} =
    if pState.serverP.idx > 1:
      dec pState.serverP.idx
    else:
      pState.serverP.idx = 0

  func incInRangeServer*(pState: var PlotState) {.inline.} =
    if pState.serverP.idx < (pState.serverP.nObj - 1):
      inc pState.serverP.idx
    else:
      pState.serverP.idx = max(pState.serverP.nObj - 1, 0)

  proc toggleServer*(conf: var Config) =
    echo "Plot ", conf.plotViaServer
    conf.plotViaServer = not conf.plotViaServer
    echo "Plot now ", conf.plotViaServer

  proc toggleWebsocketConnection*(conf: var Config) =
    ## TODO: clarify with above proc!!!
    echo "websocket active: ", conf.websocketActive
    conf.websocketActive = not conf.websocketActive
    echo "Plot now ", conf.websocketActive


  proc toggleEventDisplay*(conf: var Config) =
    conf.plotEventDisplay = not conf.plotEventDisplay

  proc toggleInteractivePlot*(conf: var Config) =
    conf.interactivePlot = not conf.interactivePlot

  proc fromServerPrev*(pState: var PlotState) {.inline.} =
    pState.decInRangeServer

  proc fromServerNext*(pState: var PlotState) {.inline.} =
    pState.incInRangeServer
