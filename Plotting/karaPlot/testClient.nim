import websocket, asyncnet, asyncdispatch, os


let ws = waitFor newAsyncWebsocketClient("localhost",
                                         Port 8080, "/?encoding=text", ssl = false)
echo "connected!"

proc reader() {.async.} =
  sleep(500)
  await ws.sendText("Connected!")
  while true:
    await ws.sendText("request")
    let read = await ws.readData()
    echo "read: ", read

proc ping() {.async.} =
  while true:
    await sleepAsync(6000)
    echo "ping"
    await ws.sendPing(masked = true)

asyncCheck reader()
#asyncCheck ping()
runForever()
