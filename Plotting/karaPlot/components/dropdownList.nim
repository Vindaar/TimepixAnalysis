include karax / prelude
import sequtils
import macros
import sugar

proc onChangeP(ev: Event, n: VNode) =
  discard

proc renderDropdownList*(names: seq[kstring], values: seq[kstring],
                         onChangeProc:
                           proc(ev: Event, n: VNode): void = onChangeP,
                         id = ""): VNode =
  doAssert names.len == values.len
  result = buildHtml:
    select(id = id,
           onChange = onChangeProc):
      for i, xy in zip(names, values):
        option(id = $i,
               value = xy[1]):
          text xy[0]

proc renderDropdownList*[T](tups: seq[(kstring, T)],
                            onChangeProc:
                              proc(ev: Event, n: VNode): void = onChangeP,
                            id = ""): VNode =
  var names: seq[kstring]
  var values: seq[kstring]
  for i, xy in tups:
    names.add $xy[0]
    values.add $xy[1]
  result = renderDropdownList(names, values, onChangeProc, id)

proc renderNamedList*(title: kstring, list: VNode, id = ""): VNode =
  result = buildHtml:
    tdiv(id = id):
      text title
      br()
      list
