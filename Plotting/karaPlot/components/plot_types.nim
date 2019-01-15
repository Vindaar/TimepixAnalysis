import jsffi
include karax / prelude
import ../common_types
import ingrid / ingrid_types

type
  Config* = object
    useWs*: bool
    connected*: bool
    receiving*: bool
    doneReceiving*: bool
    plotViaServer*: bool
    plotEventDisplay*: bool
    interactivePlot*: bool
    websocketActive*: bool

  #PlotState = object
  #  staticIdx: int    # index for static data from JSON file
  #  evDisplayIdx: int # index for event display / received from server

  StaticPlot* = object
    svg*: JsAssoc[int, JsObject]
    nSvg*: int
    plt*: JsAssoc[int, JsObject]
    nPlt*: int
    keys*: seq[kstring]
    idx*: int

  ServerPlot* = object
    data*: JsAssoc[int, JsObject]
    nObj*: int
    svg*: JsAssoc[int, JsObject]
    nSvg*: int
    plt*: JsAssoc[int, JsObject]
    nPlt*: int
    keys*: seq[kstring]
    idx*: int

  PlotState* = object
    staticP*: StaticPlot
    serverP*: ServerPlot
    pds*: seq[(kstring, PlotDescriptor)]
    activePd*: PlotDescriptor
    fileInfo*: FileInfo
