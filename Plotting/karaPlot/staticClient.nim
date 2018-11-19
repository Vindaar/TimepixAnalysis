import json
import plotly
import strutils, strformat, tables, sugar, sequtils
import chroma
import jsffi except `&`
import jsbind
#import json
from dom import getElementById
import jswebsockets# except Event
include karax / prelude
import karax / [kdom, vstyles]

import components / [button, plotWindow, figSelect, utils, plot_types]
import frontend / [menu, figDropdown]
import protocol

let plt = newPlotly()

# include the following containing the path of the file we wish
# to read
include resources/data_input
const data = staticRead(fname)

#type
#  Dropdown = object

# for dynamic client:
# add fields for
# - bins  <- to change binning of an existing plot. Entering that will request
#   the same plot with the desired binning
# - ...

# add event option to show events / clusters (field for event number)
# incuding centerX, centerY

when defined(client):
  import components / serverComm
  static:
    echo "USING"
else:
  static:
    echo "NOT USING"

proc initPlotState(jData: JsObject): PlotState =
  result.staticP.svg = newJsAssoc[int, JsObject]()
  result.staticP.plt = newJsAssoc[int, JsObject]()
  let svgPairs = jData["svg"]
  let pltPairs = jData["plotly"]
  let svgKeys = toSeq(keys(svgPairs))
  let pltKeys = toSeq(keys(pltPairs))
  for i, k in svgKeys:
    result.staticP.svg[i] = svgPairs[svgKeys[i]]
    result.staticP.keys.add $k
  result.staticP.nSvg = svgKeys.len
  for i, k in pltKeys:
    result.staticP.plt[i + result.staticP.nSvg] = pltPairs[pltKeys[i]]
    result.staticP.keys.add k
  result.staticP.nSvg = svgKeys.len
  result.staticP.nPlt = pltKeys.len

  when UseWS:
    # server data is empty at the beginning
    result.serverP.data = newJsAssoc[int, JsObject]()
    result.serverP.svg = newJsAssoc[int, JsObject]()
    result.serverP.plt = newJsAssoc[int, JsObject]()

proc main =

  echo "Start parsing..."
  let jData = parseJsonToJs(data) #parseJson(data)
  var plotState = initPlotState(jData)
  var conf = Config(plotViaServer: false, useWs: false)
  when UseWs:
    conf.useWs = UseWs

  when UseWS:
    # start websocket connection
    var socket = newWebSocket("ws://localhost:8080")
    socket.onOpen = proc (e: dom.Event) =
      echo("sent: Connected!")
      socket.send($Messages.Connected)
      conf.connected = true
      conf.receiving = true
    socket.onClose = proc (e: CloseEvent) =
      # TODO: check if proper close!
      conf.connected = false
      conf.receiving = false
      conf.doneReceiving = true
      redraw()

  echo "...parsed"
  var packet = initDataPacket()

  proc getData(pState: PlotState): JsObject =
    ## gets the appropriate data to plot given the current state
    when not UseWS:
      result = pState.staticP.plt[pState.staticP.idx]
    else:
      if conf.plotViaServer and pState.serverP.nObj > 0:
        # TODO: need to be able to selet whether SVG or plotly
        result = pState.serverP.plt[pState.serverP.idx]
      else:
        echo "nooo"
        result = pState.staticP.plt[pState.staticP.idx]

  proc renderPlotly(pState: PlotState, conf: Config) =
    # make a plotly plot
    if pState.staticP.idx >= pState.staticP.nSvg:
      let pltData = getData(pState)
      # TODO: implement a case, store last plot kind. If last plot kind is the same
      # use react, else use newPlot
      if conf.plotViaServer:
        echo "Ddd ", toSeq(keys(pltData))
        echo "Ddd2 ", toSeq(keys(pltData["Layout"]))
      plt.newPlot("plotly", pltData["Traces"], pltData["Layout"])
      let plotlyPlot = kdom.document.getElementById("plotly")
      plotlyPlot.style.visibility = "visible"
      redraw()
    else:
      # hide plot
      let plotlyPlot = kdom.document.getElementById("plotly")
      plotlyPlot.style.visibility = "hidden"

  when UseWS:
    proc finished(packet: DataPacket): bool =
      result = not packet.recvData and packet.done

    proc handleData(packet: var DataPacket, s: kstring) =
      if s == $Messages.DataStart:
        packet.recvData = true
        echo "Receving data now!"
      else:
        if packet.recvData:
          if s == $Messages.DataStop:
            echo "Data done! ", packet.data.len
            packet.done = true
            packet.recvData = false
          else:
            packet.data &= s
            echo "Got ", s.len, " bytes, now total ", packet.data.len

      if packet.finished:
        echo "packet done "

    proc parsePacket(pState: var PlotState, data: kstring) =
      pState.serverP.data[pState.serverP.nObj] = parseJsonToJs(packet.data)
      # parse that packet, add to svg / plotly
      template getPairsKeys(data, nObj: untyped, name: kstring): untyped =
        let pPairs = data[nObj][name]
        let keys = toSeq(keys(pPairs))
        (pPairs, keys)
      let (svgPairs, svgKeys) = getPairsKeys(pState.serverP.data,
                                             pState.serverP.nObj,
                                             "svg")
      let (pltPairs, pltKeys) = getPairsKeys(pState.serverP.data,
                                             pState.serverP.nObj,
                                             "plotly")
      # increase object count
      inc pState.serverP.nObj
      template addPlt(pPlot, pKeys, idx, pPairs, keys: untyped): untyped =
        if keys.len > 0:
          for i, k in keys:
            let nPlots = idx
            pPlot[nPlots] = pPairs[k]
            pKeys.add k
            inc idx
      addPlt(pState.serverP.svg, pState.serverP.keys, pState.serverP.nSvg, svgPairs, svgKeys)
      addPlt(pState.serverP.plt, pState.serverP.keys, pState.serverP.nPlt, pltPairs, pltKeys)

  proc postRender() =
    ## this is called after rendering via karax. First render the plotly plot
    ## then, handle websocket
    when UseWS:
      # request a new frame from server
      if conf.connected:
        socket.fromServer()
      if not conf.plotViaServer or conf.doneReceiving:
        renderPlotly(plotState, conf)
      var din {.global.}: string = ""
      socket.onMessage = proc (e: MessageEvent) =
        packet.handleData(e.data)
        if packet.finished:
          # parse new object
          plotState.parsePacket(packet.data)
          # reset packet
          packet = initDataPacket()
          echo "Obj count now ", plotState.serverP.nObj
          echo "Packet reset"
          renderPlotly(plotState, conf)
    else:
      renderPlotly(plotState, conf)

  proc render(): VNode =

    result = buildHtml(tdiv):
      h1(text "Static karaPlot")
      renderMenu(plotState, conf)
      p:
        br()
        text "Next: " & $plotState.staticP.idx & " " & plotState.getNextStatic
        br()
        text "Previous: " & $plotState.staticP.idx & " " & plotState.getPrevStatic
      renderFigDropdown(plotState)
      p:
        span(text "Number of static plots available: " & $plotState.staticP.keys.len)
      p:
        if plotState.staticP.idx < plotState.staticP.nSvg:
          renderSvgPlot(plotState.staticP.svg[plotState.staticP.idx])
        # create `div` for the plotly Plot
        tdiv(id = "plotly",
             class = "plot-style")

  setRenderer render, "ROOT", postRender
  setForeignNodeId "plotly"

when isMainModule:
  main()
