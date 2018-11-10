import sugar
import plotly
include karax/prelude

proc renderPlot*(data: string,
                 isSvg: bool): VNode =
  if isSvg:
    result = buildHtml:
      verbatim(data)
  else:
    # make a plotly plot
    let plt = newPlotly()
    #plt.newPlot(
      
