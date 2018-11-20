include karax / prelude
import sequtils
import karax / kdom
import ../components/dropdownList
import ../components/plot_types
import ../protocol

proc selectedIndex*(n: Node): cint {.importcpp: "#.selectedIndex".}

proc renderPlotDescriptorSelect*(pState: var PlotState): VNode =
  proc onChangeHandler(ev: Event, n: VNode) =
    echo "N val ", n.value
    echo "N index ", n.index
    let x = document
              .getElementById("pdselect").selectedOptions[0].text
    let y = document
              .getElementById("pdselect").selectedIndex
    pState.serverP.idx = y
    echo "x ", x
    echo " y ", y
    #pState.activePd = n.value
    #pState.serverP.idx = n.index
    #pState.serverP.idx = n.index

  if pState.pds.len > 0:
    result = buildHtml(tdiv):
      renderNamedList("PlotDescriptors",
                      renderDropdownList(pState.pds,
                                         onChangeProc = onChangeHandler,
                                         id = "pdselect"))

  else:
    result = buildHtml(tdiv())
