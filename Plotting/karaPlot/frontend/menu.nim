include karax / prelude
import jswebsockets
import sugar
import ../components / [plot_types, button, utils, serverComm]

proc renderMenu*(pState: var PlotState, conf: var Config): VNode =
  result = buildHtml:
    tdiv(id = "grid"):
      tdiv:
        renderButton("Previous",
                     onClickProc = () => pState.decInRangeStatic)
      tdiv:
        renderButton("Next",
                     onClickProc = () => pState.incInRangeStatic)
      when UseWs:
        tdiv:
          renderButton("EventDisplay",
                       onClickProc = () => conf.toggleServer())
        tdiv:
          renderButton("Interactive Plotting",
                       onClickProc = () => conf.toggleServer())
        text("Plotting via server: " & $(conf.plotViaServer))

        if conf.plotViaServer:
          tdiv:
            text("Plotting event: " & $pState.serverP.idx)
          tdiv:
            renderButton("Previous from server",
                         onClickProc = () => pState.fromServerPrev())
          tdiv:
            renderButton("Next from server",
                         onClickProc = () => pState.fromServerNext())
        br()
        tdiv:
          if conf.receiving:
            text("Receiving plots from server... " & $pState.serverP.nObj &
              "/" & $pState.pds.len)
          else:
            text("no data receiving")
          if conf.doneReceiving:
            text("All plots received from server... " & $pState.serverP.nObj)
