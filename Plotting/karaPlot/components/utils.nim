import jsffi, jsbind
include karax / prelude
import plot_types

when defined(client):
  const UseWs* = true
else:
  const UseWs* = false

# allows to convert `JsObject` to string
proc toString*(x: JsObject): cstring {.jsimportgWithName: "String".}

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

  proc toggleEventDisplay*(conf: var Config) =
    conf.plotEventDisplay = not conf.plotEventDisplay

  proc toggleInteractivePlot*(conf: var Config) =
    conf.interactivePlot = not conf.interactivePlot

  proc fromServerPrev*(pState: var PlotState) {.inline.} =
    pState.decInRangeServer

  proc fromServerNext*(pState: var PlotState) {.inline.} =
    pState.incInRangeServer
