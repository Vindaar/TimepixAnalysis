#import json
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
import karax / jjson

import components / [button, plotWindow, figSelect, utils, plot_types, dropdownList]
import frontend / [menu, figDropdown, plotSelect, plotDescriptorSelect]
import protocol
import common_types

let plt = newPlotly()

# include the following containing the path of the file we wish
# to read
include resources/data_input
const data = staticRead(fname).kstring

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
  #let jData = parseJson(data)
  echo "done"
  let jData = parseJsonToJs(data) #parseJson(data)
  var plotState = initPlotState(jData)
  var conf = Config(plotViaServer: false, useWs: false,
                    websocketActive: true)

  var fields: seq[(kstring, PlotKind)]
  for f in PlotKind:
    fields.add (kstring($f), f)

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

      # first request FileInfo
      socket.send($initDataPacket(PacketKind.Request, rqFileInfo))
      # then PlotDescriptors
      socket.send($initDataPacket(PacketKind.Request, rqPlotDescriptors))

    socket.onClose = proc (e: CloseEvent) =
      # TODO: check if proper close!
      conf.connected = false
      conf.receiving = false
      conf.doneReceiving = true
      echo "Done receving all data!"
      redraw()

  echo "...parsed"
  var packetStore = initDataPacketStorage()

  proc getData(pState: PlotState): JsObject =
    ## gets the appropriate data to plot given the current state
    when not UseWS:
      result = pState.staticP.plt[pState.staticP.idx]
    else:
      if conf.plotViaServer and pState.serverP.nObj > 0:
        # TODO: need to be able to selet whether SVG or plotly
        result = pState.serverP.plt[pState.serverP.idx]
      else:
        result = pState.staticP.plt[pState.staticP.idx]

  proc renderPlotly(pState: PlotState, conf: Config) =
    # make a plotly plot
    if pState.staticP.idx >= pState.staticP.nSvg:
      let pltData = getData(pState)
      # TODO: implement a case, store last plot kind. If last plot kind is the same
      # use react, else use newPlot
      if conf.plotViaServer:
        echo "data keys ", toSeq(keys(pltData))
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
      result = not packet.header.recvData and packet.header.done

    proc getPacket(packets: DataPacketStorage,
                   pKind: PacketKind): DataPacket =
      case pKind
      of PacketKind.Descriptors: result = packets.descPacket
      of PacketKind.Plots: result = packets.pltPacket
      of PacketKind.FileInfo: result = packets.fiPacket
      of PacketKind.Request: result = packets.reqPacket
      else: echo "WARNING! Cannot get packet for kind: ", pKind

    proc getPacket(packets: DataPacketStorage,
                   sPacket: DataPacket): DataPacket =
      result = packets.getPacket(sPacket.kind)

    proc assignPacket(packets: var DataPacketStorage,
                      sPacket: DataPacket) =
      case sPacket.kind
      of PacketKind.Descriptors: packets.descPacket = sPacket
      of PacketKind.Plots: packets.pltPacket = sPacket
      of PacketKind.FileInfo: packets.fiPacket = sPacket
      of PacketKind.Request: packets.reqPacket = sPacket
      else: echo "WARNING! Cannot assign packet for kind: ", sPacket.kind

    proc resetPacket(packets: var DataPacketStorage,
                     pKind: PacketKind) =
      case pKind
      of PacketKind.Descriptors:
        packets.descPacket = initDataPacket(kind = PacketKind.Descriptors)
      of PacketKind.Plots:
        packets.pltPacket = initDataPacket(kind = PacketKind.Plots)
      of PacketKind.FileInfo:
        packets.fiPacket = initDataPacket(kind = PacketKind.FileInfo)
      of PacketKind.Request:
        packets.reqPacket = initDataPacket(kind = PacketKind.Request)
      else: echo "WARNING! Cannot reset packet for kind: ", pKind

    proc handleData(packets: var DataPacketStorage, s: kstring): PacketKind =
      ## handles a single data packet received and returns the `PacketKind`
      ## to know what packet arrived in the calling scope
      let singlePacket = parseDataPacket(s)
      var packet = packets.getPacket(singlePacket)
      # copy the header

      case singlePacket.kind
      of PacketKind.Descriptors, PacketKind.Plots, PacketKind.FileInfo:
        packet.header = singlePacket.header
        echo "Yes"
        echo "Packet header! ", packet.header
        echo "no"
        case packet.header.msg
        of Messages.DataStart:
          echo "Receving data now!"
          packet.payload = singlePacket.payload
        of Messages.Data, Messages.DataStop:
          packet.payload &= singlePacket.payload
        of Messages.DataSingle:
          echo "Single packet received!"
          packet.payload &= singlePacket.payload
        else:
          echo "WARNING: couldn't parse message kind! " & $packet.header.msg
      of PacketKind.Request: discard
      else:
        echo "WARNING: Unknown data packet!"

      # re assign the variable packet to the parent
      # TODO: achieve by using ref type instead?
      packets.assignPacket(packet)

      result = singlePacket.kind

    proc assignPlotPacket(pState: var PlotState, packet: DataPacket) =
      doAssert packet.kind == PacketKind.Plots
      pState.serverP.data[pState.serverP.nObj] = parseJsonToJs(packet.payload)
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
      template addPlt(pPlot, pKeys, idx, pPairs, keys: untyped): untyped =
        if keys.len > 0:
          for i, k in keys:
            let nPlots = idx
            pPlot[nPlots] = pPairs[k]
            pKeys.add k
            inc idx
      addPlt(pState.serverP.svg, pState.serverP.keys, pState.serverP.nSvg, svgPairs, svgKeys)
      addPlt(pState.serverP.plt, pState.serverP.keys, pState.serverP.nPlt, pltPairs, pltKeys)
      # increase object count
      inc pState.serverP.nObj

    proc assignDescriptorPacket(pState: var PlotState, packet: DataPacket) =
      doAssert packet.kind == PacketKind.Descriptors
      let pJson = parseJson(packet.payload)
      for x in pJson:
        let pd = parsePd(x)
        let outf = buildOutfile(pd)
        pState.pds.add (outf, pd)
      echo "PDs are ", $pState.pds

    proc assignFileInfoPacket(pState: var PlotState, packet: DataPacket) =
      doAssert packet.kind == PacketKind.FileInfo
      let fiJson = parseJson(packet.payload).parseFileInfo
      echo "received fileInfo: ", $fiJson
      pState.fileInfo = fiJson

    proc parsePacket(pState: var PlotState, packet: DataPacket) =
      case packet.kind
      of PacketKind.Plots:
        echo "Is a plot packet!"
        assignPlotPacket(pState, packet)
      of PacketKind.Descriptors:
        echo "Is a descriptors packet!"
        assignDescriptorPacket(pState, packet)
      of PacketKind.FileInfo:
        echo "Is a file info packet!"
        assignFileInfoPacket(pState, packet)
      of PacketKind.Request:
        if packet.reqKind == rqPing:
          # send a pong
          socket.send($(initDataPacket(kind = PacketKind.Request)))
        else:
          echo "WARNING: Requests from server to client only valid if ping!"
      else:
        echo "WARNING: couldn't parse packet kind! " & $packet.kind

  proc postRender() =
    ## this is called after rendering via karax. First render the plotly plot
    ## then, handle websocket
    when UseWS:
      # request a new frame from server
      if not conf.websocketActive:
        socket.close()
      if conf.connected:
        socket.fromServer()
      if not conf.plotViaServer or conf.doneReceiving:
        renderPlotly(plotState, conf)

      socket.onMessage = proc (e: MessageEvent) =
        let pKind = packetStore.handleData(e.data)
        var sPacket = packetStore.getPacket(pKind)
        if sPacket.finished:
          # parse new object
          plotState.parsePacket(sPacket)
          # reset packet
          packetStore.resetPacket(pKind)
          echo "Obj count now ", plotState.serverP.nObj
          renderPlotly(plotState, conf)
    else:
      renderPlotly(plotState, conf)

  proc render(): VNode =
    result = buildHtml(tdiv):
      h1(text "Static karaPlot")
      renderMenu(plotState, conf)

      if conf.interactivePlot:
        when UseWs:
          renderPlotSelect(plotState, socket, fields)
        else: discard
      p:
        br()
        text "Next: " & $plotState.staticP.idx & " " & plotState.getNextStatic
        br()
        text "Previous: " & $plotState.staticP.idx & " " & plotState.getPrevStatic
      renderFigDropdown(plotState)

      renderPlotDescriptorSelect(plotState)

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
