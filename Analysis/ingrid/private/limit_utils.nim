import std / [os, strutils]
import ggplotnim

let UseTeX* = getEnv("USE_TEX", "false").parseBool
let FWIDTH* = getEnv("F_WIDTH", "0.0").parseFloat
let Newline* = if UseTeX: r"\\" else: ""
let Height* = getEnv("Height", "420").parseFloat

proc tripleSide*(): Theme =
  result = Theme(titleFont: some(font(7.0)),
                 labelFont: some(font(7.0)),
                 tickLabelFont: some(font(6.0)),
                 tickLength: some(4.0),
                 tickWidth: some(0.8),
                 gridLineWidth: some(0.8),
                 legendFont: some(font(6.0)),
                 legendTitleFont: some(font(6.0, bold = true)),
                 facetHeaderFont: some(font(6.0, alignKind = taCenter)),
                 baseLabelMargin: some(0.25),
                 annotationFont: some(font(6.0, family = "monospace")),
                 baseScale: some(1.25)) # won't be scaled!
  result = result + margin(left = 3.5, bottom = 3.0, right = 4.5) +
    continuousLegendWidth(0.75) + continuousLegendHeight(3.0)

proc thL*(fWidth: float, width: float,
         baseTheme: (proc(): Theme) = nil,
         height = -1.0, ratio = -1.0,
         textWidth = 458.29268, # 455.24411
         forceWidth = false
        ): Theme =
  if UseTeX:
    let fWidth = if forceWidth: fWidth elif FWIDTH > 0.0: FWIDTH else: fWidth
    let baseTheme = if baseTheme != nil: baseTheme
                    elif fWidth <= 0.5 and fWidth >= 0.4: sideBySide
                    elif fWidth < 0.4: tripleSide
                    else: singlePlot
    result = themeLatex(fWidth, width, baseTheme, height, ratio, textWidth,
                        useTeX = UseTeX)
  else:
    result = Theme()
