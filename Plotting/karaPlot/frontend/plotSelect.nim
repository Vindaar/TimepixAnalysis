include karax / prelude
import sequtils
import jswebsockets
import ../components / [dropdownList, utils, plot_types]
import ../protocol
import karax / kdom
import ../common_types

proc getSelectVal(id: kstring): kstring =
  let idx = document.getElementById(id).selectedIndex
  result = document.getElementById(id).options[idx].text

proc renderPlotSelect*(pState: var PlotState,
                       socket: WebSocket,
                       fields: seq[(kstring, PlotKind)]): VNode =
  var curPKind: PlotKind
  proc onChangeHandler(ev: kdom.Event, n: VNode) =
    echo "N val ", n.value
    echo "N index ", n.index
    let x = document
              .getElementById("pdselect").selectedOptions[0].text
    let runIdx = getSelectVal("runSelect")
    let chipIdx = getSelectVal("chipSelect")
    let plotKindIdx = getSelectVal("plotKindSelect")
    #let plotKindIdx = getSelectIdx("plotKindSelect")
    #let pd = buildPd(run
    echo "Run idx ", runIdx
    echo "chip idx ", chipIdx
    echo "plt type idx ", plotKindIdx

  let runs = pState.fileInfo.runs.mapIt((kstring($it), it))
  let chips = pState.fileInfo.chips.mapIt((kstring($it), it))
  result = buildHtml:
    tdiv(id = "grid"):
      renderNamedList("RunNumber",
                      renderDropdownList(runs,
                                         onChangeProc = onChangeHandler,
                                         id = "runSelect"))


      renderNamedList("ChipNumber",
                      renderDropdownList(chips,
                                         onChangeProc = onChangeHandler,
                                         id = "chipSelect"))

      renderNamedList("PlotType",
                      renderDropdownList(fields,
                                         onChangeProc = onChangeHandler,
                                         id = "plotKindSelect"))



      # get active element of plot type and then decide rest of
