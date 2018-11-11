import sugar
include karax/prelude
import karax / [kdom, vstyles]
import json
import sequtils
import jsffi
import utils

proc renderSvgPlot*(data: JsObject = nil, disable = false): VNode =
  if not disable:
    result = buildHtml:
      tdiv(id = "plotSvg",
           disabled = "false",
           class = "plot-style"):
           #style = style(StyleAttr.width, kstring"60%")):
                         #(StyleAttr.display, "table".kstring))):
        verbatim(toString(data))
#proc renderPlotly*(plt: PlotlyObj, data: JsonNode): VNode =
#  # make a plotly plot
#  echo toSeq(keys(data.getFields))
#  result = buildHtml():
#    react(plt, "plot0",
#          data["Traces"], data["Layout"], output_type = "div")
  #let el = kdom.document.getElementById("plot0")
  #echo el.repr
  #result = react(plt, el, data["Traces"], data["Layout"], output_type = "div")
