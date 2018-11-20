include karax / prelude
import sequtils
import ../components/dropdownList
import ../components/plot_types
import ../protocol

proc renderPlotDescriptorSelect*(pState: PlotState): VNode =
  result = buildHtml(tdiv):
    renderNamedList("PlotDescriptors", renderDropdownList(pState.pds))
