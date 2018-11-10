import json
import plotly
import strutils, strformat, tables, sugar, sequtils
import chroma
import jsffi
#import json
#import dom
import jswebsockets# except Event
include karax / prelude
import karax / kdom

import components / [button, plotWindow, figSelect]


const data = staticRead("/home/basti/CastData/ExternCode/TimepixAnalysis/Plotting/karaPlot/calibration_cfNoOccupancy_cfNoPolya.json")

#type
#  Dropdown = object
    
    
proc main =

  echo "Start parsing..."
  let jData = parseJson(data)
  echo "...parsed"  
  let svgPairs = jData["svg"].getFields
  let svgKeys = toSeq(keys(svgPairs))
  for k in keys(svgPairs):
     echo k
  var i = 0
  proc render(): VNode =

    #var svgPlt = fnamesSvg[0]
    #for x in jData["svg"]:
    buildHtml(tdiv):
      h1(text "Static karaPlot")
      p:
        renderButton("TestButton",
                     onClickProc = () => inc i)
      p:
        tdiv(class = "dropdown")
        renderButton("Dropdown",
                     class = "dropbtn",
                     onClickProc = () => kdom.document.getElementById("myDropdown").classList.toggle("show"))
        tdiv(id = "myDropdown",
             class = "dropdown-content"):
          var idx = 0
          for k in keys(svgPairs):
            echo "K is ", k, " idx " , idx
            p:
              renderFigSelect($k,
                              idx,
                              onClickProc = (event: kdom.Event, node: VNode) => (i = node.id.parseInt))
            inc idx
      p:
        span(text $svgKeys)
        span(text $i)
      p:
        renderPlot(svgPairs[svgKeys[i]].getStr, isSvg = true)

  setRenderer render

when isMainModule:
  main()
