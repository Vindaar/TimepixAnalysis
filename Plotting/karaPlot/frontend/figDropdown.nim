include karax / prelude
import karax / kdom
import sugar
import ../protocol
import ../components/ [plot_types, button, figSelect]

proc renderFigDropdown*(pState: var PlotState): VNode =
  result = buildHtml(p):
    tdiv(class = "dropdown")
    renderButton("Dropdown",
                 class = "dropbtn",
                 onClickProc = () => kdom.document.getElementById("myDropdown").classList.toggle("show"))
    tdiv(id = "myDropdown",
         class = "dropdown-content"):
      var idx = 0
      for k in pState.staticP.keys:
        p:
          renderFigSelect(
            $k,
            idx,
            onClickProc = (event: kdom.Event, node: VNode) => (
              pState.staticP.idx = node.id.parseInt)
          )
        inc idx
