import karax / kbase
export kbase

type
  Messages* = enum
    Connected = "Connected"
    DataStart = "DataStart"
    DataStop = "DataStop"
    Request = "Request"

  DataPacket* = object
    recvData*: bool
    done*: bool
    data*: kstring

const FakeFrameSize* = 32760

proc initDataPacket*(): DataPacket =
  result.recvData = false
  result.done = false
  result.data = kstring""
