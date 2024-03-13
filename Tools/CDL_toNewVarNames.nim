## a simple program which extracts the names of the properties
## in the calibration-cdl.h5 file, in order to relate them to
## the new variable names.

import nimhdf5
import os, strutils
import strtabs
import tables


const pathCDL = "/mnt/Daten/CAST/CDL-reference/calibration-cdl.h5"
const pathRun = "/home/basti/CastData/ExternCode/TimepixAnalysis/Analysis/ingrid/run_file.h5"
# take a random group as the path in order to extract the names
const groupCdlPath = "/calibration-cdl-apr2014-Ag-Ag-6kV"
const groupRunPath = "/reconstruction/run_114/chip_3"

var h5cdl = H5file(pathCDL, "r")
var h5run = H5file(pathRun, "r")

var groupCDl = h5cdl[groupCdlPath.grp_str]
var groupRun = h5run[groupRunPath.grp_str]

for d in groupCDL:
  echo d.name.extractFilename

echo ""

for d in groupRun:
  echo d.name.extractFilename

# the resulting tab is the following
let tab = {"NumberOfPixels"              : "hits",
           "Width"                       : "width",
           "SkewnessTransverse"          : "skewnessTransverse",
           "KurtosisLongitudinal"        : "kurtosisLongitudinal",
           "Length"                      : "length",
           "FractionWithinRmsTransverse" : "fractionInTransverseRms",
           "PositionX"                   : "centerX",
           "RmsLongitudinal"             : "rmsLongitudinal",
           "KurtosisTransverse"          : "kurtosisTransverse",
           "EnergyFromCharge"            : "energyFromCharge",
           "RotationAngle"               : "rotationAngle",
           "SkewnessLongitudinal"        : "skewnessLongitudinal",
           "EventNumber"                 : "eventNumber",
           "RmsTransverse"               : "rmsTransverse",
           "PositionY"                   : "centerY",
           "Excentricity"                : "eccentricity",
           "TotalCharge"                 : "totalCharge" }.toTable()

# using this we will store all information from the CDL file in a
# new version of it

var h5f = H5file("calibration-newTos.h5", "rw")
var recoGroup = h5f.create_group("/reconstrution")
var dataTab = initTable[string, seq[float32]]()
for k in keys(tab):
  dataTab[k] = @[]
for g in h5cdl:
  # iterating over all groups and read data in `tab`, then add to dataTab
  echo "in g ", g
  for key in keys(tab):
    echo "accessing ", g.name / key
    let data = h5cdl[g.name / key, float32]
    # given data add to tab
    dataTab[key].add data

for k, v in dataTab:
  # get key for output
  let outKey = tab[k]
  var dset = h5f.create_dataset(recoGroup.name / outKey, dtype = float32, shape_raw = @[v.len, 1])
  dset[dset.all] = v

discard h5f.close()
