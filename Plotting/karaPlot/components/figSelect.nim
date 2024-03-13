import sugar
include karax/prelude

proc renderFigSelect*(caption: string,
                      idx: int,
                      onClickProc: (event: Event, node: VNode) -> void): VNode =
  buildHtml:
    button(class = "clear-completed",
           id = $idx,
           onClick = onClickProc):
      text caption
