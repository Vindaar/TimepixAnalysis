include karax / prelude
import sequtils
import ../components/dropdownList

proc renderPlotSelect*(fields: seq[kstring],
                       runs: seq[kstring],
                       chips: seq[kstring]): VNode =
  result = buildHtml:
    tdiv(id = "grid"):
      renderNamedList("RunNumber", renderDropdownList(runs, runs))

      renderNamedList("ChipNumber", renderDropdownList(chips, chips))

      renderNamedList("PlotType", renderDropdownList(fields, fields))
      # get active element of plot type and then decide rest of
