import sugar
include karax/prelude

proc renderButton*(caption: string,
                   class = "",
                   onClickProc: () -> void): VNode =
  buildHtml:
    button(class = "clear-completed",
           onClick = onClickProc):
      text caption
