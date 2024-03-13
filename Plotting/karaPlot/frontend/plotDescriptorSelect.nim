include karax / prelude
import sequtils
import karax / kdom
import ../components / [dropdownList, plot_types, utils]
import ../protocol

proc renderPlotDescriptorSelect*(pState: var PlotState): VNode =
  ## TODO: selected element from dropdown list first check whether
  ## the corresponding PD is already part of the table of our
  ## PDs. If yes, show the corresponding plot, if not, send
  ## request to server to show it
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
    pState.activePd = n.value.parseJson.parsePd
    #pState.serverP.idx = n.index

  if pState.pds.len > 0:
    result = buildHtml(tdiv):
      renderNamedList("PlotDescriptors",
                      renderDropdownList(pState.pds,
                                         onChangeProc = onChangeHandler,
                                         id = "pdselect"))

  else:
    result = buildHtml(tdiv())
