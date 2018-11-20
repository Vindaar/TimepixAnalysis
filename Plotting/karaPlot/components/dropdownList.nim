include karax / prelude
import sequtils
import macros

proc renderDropdownList*(names: seq[kstring], values: seq[kstring]): VNode =
  doAssert names.len == values.len
  result = buildHtml(select):
    for xy in zip(names, values):
      option(value = xy[1]):
        text xy[0]

proc renderDropdownList*[T](tups: seq[(kstring, T)]): VNode =
  var names: seq[kstring]
  var values: seq[kstring]
  for xy in tups:
    names.add $xy[0]
    values.add $xy[1]
  expandMacros:
    result = buildHtml(select):
      for i in 0 .. names.high:
        echo xy
      #ption(value = xy[0]):
      #j  text xy[1]

  #result = renderDropdownList(names, values)

proc renderNamedList*(title: kstring, list: VNode): VNode =
  result = buildHtml(tdiv):
    text title
    br()
    list
