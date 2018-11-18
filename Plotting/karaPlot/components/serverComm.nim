include karax / prelude
import plot_types
import jswebsockets
import ../protocol
import utils

proc fromServer*(socket: WebSocket) =
  ## request next plot from server
  debugecho "Sent: request"
  socket.send($Messages.Request)
  #socket.onMessage = proc (e: MessageEvent) =
  #  echo "receving"
  #  echo("received: ",e.data)
