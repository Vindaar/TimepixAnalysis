import websocket, asynchttpserver#, asyncnet, asyncdispatch
import plotData

let server = newAsyncHttpServer()

proc main =

  asyncCheck server.serve(Port(8080, cb))


when isMainModule:
  main()
