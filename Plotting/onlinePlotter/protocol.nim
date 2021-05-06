import plotly
import karax / kbase
import strutils

when defined(js):
  import karax / jjson
  import jsffi
  import jsbind
else:
  import json


type
  Messages* {.pure.} = enum
    MessageUnknown = "Error"
    Connected = "Connected"
    Data = "Data"
    Request = "Request"

  DataPacket* = object
    case kind*: Messages
    of Messages.Data:
      payload*: kstring
    else: discard

  JsonPlot* = object
    traces*: JsonNode
    layout*: JsonNode

when not defined(js):
  proc asData*(packet: DataPacket): kstring =
    let node = % packet
    toUgly(result, node)

  proc `%`*(grid: Grid): JsonNode =
    ## serialize a plotly grid
    result = newJObject()
    let pltJson = grid.toPlotJson
    result["Traces"] = pltJson.traces
    result["Layout"] = pltJson.layout

else:
  proc parseJsonToJs(json: cstring): JsObject {.jsimportgWithName: "JSON.parse".}
  proc stringify*(value: JsObject | JsonNode,
                  replacer: JsObject,
                  space: JsObject): cstring {.jsimportgWithName: "JSON.stringify".}
  proc toString*(x: JsObject | JsonNode): cstring =
    result = x.stringify(nil, toJs(2))

  proc pretty*(x: JsonNode): cstring =
    result = x.stringify(nil, toJs(2))

  proc parseJson*(s: kstring): JsonNode =
    result = % parseJsonToJs(s)

  proc parseEnum*[T: enum](s: cstring, default: T): T =
    result = strutils.parseEnum[T]($s, default)

proc initDataPacket*(kind: Messages, payload = kstring""): DataPacket =
  case kind
  of Messages.Data:
    result = DataPacket(kind: Messages.Data,
                        payload: payload)
  else:
    result = DataPacket(kind: kind)

proc `%`*(dp: DataPacket): JsonNode =
  result = newJObject()
  result["kind"] = % $dp.kind
  case dp.kind
  of Messages.Data:
    result["payload"] = % dp.payload
  else: discard

proc `$`*(packet: DataPacket): kstring =
  result = (% packet).pretty

proc parseGrid*(jnode: JsonNode): JsonPlot =
  ## parses a JSON string containing a `JsonPlot` object, which itself contains
  ## a Grid plot
  result.traces = jnode[kstring"Traces"]
  result.layout = jnode[kstring"Layout"]

when defined(js):
  proc parseGrid*(data: kstring): JsonPlot =
    ## parses a JSON string containing a `JsonPlot` object, which itself contains
    ## a Grid plot
    #echo "DAta ", data
    let jnode = data.parseJson
    #echo jnode.len
    result = jnode.parseGrid
    #result = data.parseJsonToJs

proc parseDataPacket*(dp: kstring): DataPacket =
  let pJson = parseJson(dp)
  result.kind = parseEnum[Messages](pJson["kind"].getStr,
                                    Messages.MessageUnknown)
  case result.kind
  of Messages.Data:
    result.payload = pJson["payload"].getStr
  else: discard

when defined(js):
  proc asData*(packet: DataPacket): kstring =
    let node = % packet
    result = node.toString
