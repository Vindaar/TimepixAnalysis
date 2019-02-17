import plotly
import strutils, strformat, tables, sugar, sequtils
import jswebsockets# except Event
import protocol
import jsffi
import karax / kbase
import karax / jjson
import dom

let plt = newPlotly()

proc animateClient(socket: var WebSocket) {.exportc.} =
  ## create Plotly plots and update them as we reiceve more data
  ## via a socket. Uses `setInterval` to loop
  ## This proc is run once the user clicks on the "Start training!" button

  var inited = false

  # send command to server to start the training
  #socket.send($Messages.Request)
  proc doAgain() =
    socket.send(initDataPacket(kind = Messages.Request).asData)
    socket.onMessage = proc (e: MessageEvent) =
      #echo("received: ", e.data)
      # parse the data packet to get new data and layout
      let dp = parseDataPacket(e.data)
      let pltJson = parseGrid(dp.payload)
      # replace data with new data
      if not inited:
        plt.newPlot("ROOT".cstring, pltJson.traces, pltJson.layout)
        inited = true
      else:
        plt.react("ROOT".cstring, pltJson.traces, pltJson.layout)

  discard window.setInterval(doAgain, 100)

proc main() =
  ## main proc of the client (animated plotting using plotly.js). Open a WebSocket,
  ## create plotly `Plot`'s and then wait for data from the socket and update
  ## w/ Plotly.react

  # start websocket connection
  var socket = newWebSocket("ws://localhost:8080")
  socket.onOpen = proc (e: Event) =
    echo("sent: Connected!")
    socket.send(initDataPacket(kind = Messages.Connected).asData)

  animateClient(socket)

  # when done, close...
  socket.onClose = proc (e:CloseEvent) =
    echo("closing: ",e.reason)

when isMainModule:
  main()
