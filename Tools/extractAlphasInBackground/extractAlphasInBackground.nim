import std / [strutils, os, options, sequtils]
import ingrid / tos_helpers
import ggplotnim, cligen, nimhdf5

## Given H5 files containing background data, this script extracts those
## events with energies above a certain threshold (O(1000 keV) maybe) and
## plots those. To see the amount of possibly alpha like particles.

proc readFiles(files: seq[string]): DataFrame =
  ## reads all energy data for all runs in the given file and returns
  ## a DF of them as well as the total run time.
  var totalDuration = 0.0
  for idx, file in files:
    let h5f = H5open(file, "r")
    var df = h5f.readDsets(recoBase(), some((3, @["energyFromCharge"])),
                           commonDsets = @["eventDuration"])
    let fname = file.extractFilename
    df["File"] = constantColumn(fname, df.len)
    echo df
    result.add df

    for num, group in runs(h5f):
      echo "Looking at run number ", num
      let evDurations = h5f[group / "eventDuration", float64]
      totalDuration += evDurations.foldl(a + b, 0.0)

    discard h5f.close()

  echo "Total duration for files ", totalDuration, " s"


proc main(fnames: seq[string]) =
  let df = readFiles(fnames)
  echo "Total run time ", df["eventDuration", float].sum(), " s"
  let dfα = df.filter(f{float: `energyFromCharge` > 1000.0})
  echo "Total number of events > 1000 keV :", dfα

when isMainModule:
  dispatch main
