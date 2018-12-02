include karax / prelude
import sequtils
import ../components/dropdownList
import ../protocol

proc renderPlotSelect*(fields: seq[kstring],
                       fileInfo: FileInfo): VNode =
  let runs = fileInfo.runs.mapIt(kstring($it))
  let chips = fileInfo.chips.mapIt(kstring($it))
  result = buildHtml:
    tdiv(id = "grid"):
      renderNamedList("RunNumber", renderDropdownList(runs, runs))


      renderNamedList("ChipNumber", renderDropdownList(chips, chips))

      renderNamedList("PlotType", renderDropdownList(fields, fields))
      # get active element of plot type and then decide rest of
