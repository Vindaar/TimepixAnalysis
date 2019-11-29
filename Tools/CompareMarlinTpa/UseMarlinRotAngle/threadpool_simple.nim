type
  ThreadPool* = object
    discard

  FlowVar*[T] = object
    when T isnot void:
      v: T

proc sync*(tp: ThreadPool) =
  # dummy
  discard

proc newThreadPool*(): ThreadPool =
  result = Threadpool()

template spawn*(tp: ThreadPool, e: typed{nkCall | nkCommand}): untyped =
  when compiles(e isnot void):
    type RetType = type(e)
  else:
    type RetType = void
  FlowVar[RetType](v: e)

proc read*[T](v: FlowVar[T]): T =
  result = v.v

proc `^`*[T](fv: FlowVar[T]): T {.inline.} = fv.read()
