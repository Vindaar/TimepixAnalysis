include karax / prelude
import plot_types
import jswebsockets
import ../protocol
import utils

proc fromServer*(socket: WebSocket) =
  ## request next plot from server
  debugecho "Sent: request"
  # TODO: replace by proper DataPacket
  socket.send($initDataPacket(kind = PacketKind.Request))
  #socket.onMessage = proc (e: MessageEvent) =
  #  echo "receving"
  #  echo("received: ",e.data)
