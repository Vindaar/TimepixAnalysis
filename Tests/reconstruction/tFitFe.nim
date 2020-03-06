import ggplotnim, sequtils, strutils, strformat, ginger
import nimhdf5
import os

import ingrid / calibration / fit_functions
import ingrid / calibration / calib_fitting

const runNumber = 241
const chipNumber = 3
var h5f = H5file(paramStr(1), "r")
var feDset = h5f[("/reconstruction/run_241/chip_3/" & "FeSpectrum").dsetStr]
let feData = feDset[int64]
# call python function with data
let resNim = fitFeSpectrum(feData)
let data = resNim[0].mapIt(feSpectrumFunc(resNim[2], it))
echo data
echo type(data)
let dfNim = seqsToDf({ "hist" : resNim[1],
                       "bins" : resNim[0],
                       "fit" : data})
echo dfNim.pretty(-1)
ggplot(dfNim, aes("bins", "hist")) +
  geom_histogram(stat = "identity") +
  geom_line(aes(y = "fit"), color = some(ggColorHue(4)[3])) +
  ggsave(&"fe_spec_run_{runNumber}_chip_{chipNumber}_nim.pdf")
