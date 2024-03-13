import std / [tables, os, strutils, strformat, sugar]
import nimhdf5, cligen
import ingrid / tos_helpers

## Extracts the `CdlSpectrum` and `CdlSpectrumCharge` files
## from a `calibration-*.h5` file and stores it in CSV files

proc main(fname: string) =
  let h5f = H5File(fname, "r")
  h5f.visit_file()
  for f, grp in h5f.groups:
    let cdl = h5f[grp.name / "CdlSpectrum", float]
    let cdlCharge = h5f[grp.name / "CdlSpectrumCharge", float]
    let df = toDf({"Hits" : cdl, "Charge [e⁻]" : cdlCharge})
    let name = f.dup(removePrefix("/"))
    df.writeCsv(&"/tmp/{name}.csv")

  discard h5f.close()


when isMainModule:
  dispatch main
