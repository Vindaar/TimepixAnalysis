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

import components / [button, plotWindow, figSelect, utils]
import protocol

let plt = newPlotly()

# include the following containing the path of the file we wish
# to read
include resources/data_input
#const fname = "eventDisplay.json" #"calibration_cfNoInGrid_cfNoOccupancy_cfNoPolya_cfNoFeSpectrum.json" #"calibration_cfNoFadc_cfNoPolya.json"
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

var i = 0

when defined(client):
  const UseWS = true
  static:
    echo "USING"
else:
  const UseWs = false
  static:
    echo "NOT USING"

proc main =

  echo "Start parsing..."
  let jData = parseJsonToJs(data) #parseJson(data)

  var socket: WebSocket

  if UseWS:
    # start websocket connection
    socket = newWebSocket("ws://localhost:8080")
    socket.onOpen = proc (e: dom.Event) =
      echo("sent: Connected!")
      socket.send($Messages.Connected)
    #socket.onclose = proc (e:CloseEvent) =
    #  echo("closing: ",e.reason)

  echo toString(jData)
  echo "...parsed"
  let svgPairs = jData["svg"] #.getFields
  let pltPairs = jData["plotly"] #.getFields
  let svgKeys = toSeq(keys(svgPairs))
  let pltKeys = toSeq(keys(pltPairs))
  let allKeys = concat(svgKeys, pltKeys)
  let nSvg = svgKeys.len
  let nPlots = svgKeys.len + pltKeys.len
  var packet = initDataPacket()
  for k in keys(svgPairs):
     echo k

  for k in keys(pltPairs):
     echo k
  var i = 0

  template getNext(idx: int): kstring =
    if idx + 1 < allKeys.len:
      kstring(allKeys[idx + 1])
    else:
      kstring""

  template getPrev(idx: int): kstring =
    if idx > 0:
      kstring(allKeys[idx - 1])
    else:
      kstring""

  func decInRange(idx: var int) {.inline.} =
    if idx > 1:
      dec idx
    else:
      idx = 0

  func incInRange(idx: var int) {.inline.} =
    if idx < allKeys.high:
      inc idx
    else:
      idx = allKeys.high

  #func recv() =


  proc fromServer(idx: int) =
    ## request next plot from server
    debugecho "Sent: request"
    socket.send($Messages.Request)
    #socket.onMessage = proc (e: MessageEvent) =
    #  echo "receving"
    #  echo("received: ",e.data)


  func fromServerPrev(idx: var int) {.inline.} =
    decInRange i
    fromServer(i)

  func fromServerNext(idx: var int) {.inline.} =
    incInRange i
    fromServer(i)

  proc renderPlotly() =
    # make a plotly plot
    if i >= nSvg:
      let dd = pltPairs[pltKeys[i - nSvg]]
      # TODO: implement a case, store last plot kind. If last plot kind is the same
      # use react, else use newPlot
      plt.newPlot("plotly", dd["Traces"], dd["Layout"])
      let plotlyPlot = kdom.document.getElementById("plotly")
      plotlyPlot.style.visibility = "visible"
    else:
      # hide plot
      let plotlyPlot = kdom.document.getElementById("plotly")
      plotlyPlot.style.visibility = "hidden"

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
      echo "Is ", packet.data

  proc postRender() =
    ## this is called after rendering via karax. First render the plotly plot
    ## then, handle websocket
    renderPlotly()
    var din {.global.}: string = ""
    socket.onMessage = proc (e: MessageEvent) =
      packet.handleData(e.data)
    if packet.finished:
      packet = initDataPacket()

  #echo pltPairs[pltKeys[0]]
  proc render(): VNode =

    #var svgPplt = fnamesSvg[0]
    #for x in jData["svg"]:
    result = buildHtml(tdiv):
      h1(text "Static karaPlot")
      tdiv(id = "grid"):
        tdiv:
          renderButton("Previous",
                       onClickProc = () => decInRange i)
        tdiv:
          renderButton("Next",
                       onClickProc = () => incInRange i)
        tdiv:
          renderButton("EventDisplay",
                       onClickProc = () => echo "")
        tdiv:
          renderButton("Previous from server",
                       onClickProc = () => fromServerPrev(i))
        tdiv:
          renderButton("Next from server",
                       onClickProc = () => fromServerNext(i))

      p:
        br()
        text "Next: " & $i & " " & getNext(i)
        br()
        text "Previous: " & $i & " " & getPrev(i)
      p:
        tdiv(class = "dropdown")
        renderButton("Dropdown",
                     class = "dropbtn",
                     onClickProc = () => kdom.document.getElementById("myDropdown").classList.toggle("show"))
        tdiv(id = "myDropdown",
             class = "dropdown-content"):
          var idx = 0
          for k in keys(svgPairs):
            p:
              renderFigSelect($k,
                              idx,
                              onClickProc = (event: kdom.Event, node: VNode) => (i = node.id.parseInt))
            inc idx
          for k in keys(pltPairs):
            p:
              renderFigSelect($k,
                              idx,
                              onClickProc = (event: kdom.Event, node: VNode) => (i = node.id.parseInt))
            inc idx
      p:
        span(text $svgKeys)
        span(text $i)
      p:
        if i < nSvg:
          renderSvgPlot(svgPairs[svgKeys[i]])
        # create `div` for the plotly Plot
        tdiv(id = "plotly",
             class = "plot-style")

  setRenderer render, "ROOT", postRender
  setForeignNodeId "plotly"

when isMainModule:
  main()
