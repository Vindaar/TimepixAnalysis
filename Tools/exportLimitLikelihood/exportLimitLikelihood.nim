## This is a simple program, which can be used to extract the likelihood
## function obtained from a limit calculation done with
## `mcmc_limit_calculation`.
## It reads the HDF5 file produced by the former and exports it to the desired
## format.

import datamancer, nimhdf5
import nimhdf5 / serialize_tables
import datamancer / serialize

import numericalnim, seqmath, ggplotnim

proc limit(fname: string, outfile = "/tmp/test.csv") =
  var h5f = H5open(fname, "r")
  let df = deserializeH5[DataFrame](h5f, "chainDf", "/", exclude = @["ϑs_s", "ϑs_b", "ϑs_x", "ϑs_y"])
    .filter(f{`gs` < 2.5e-40})
  let (histo, bins) = histogram(df["gs", float].toSeq1D, bins = 1000)

  var binCenters = newSeq[float](histo.len)
  for i in 0 ..< bins.len - 1:
    binCenters[i] = (bins[i+1] + bins[i]) / 2.0
  var dfOut = toDf({ "g4s" : binCenters,
                     "L" : histo})
  echo dfOut
  dfOut = dfOut
    .mutate(f{float: "L" ~ `L` / col("L").max})
  dfOut.writeCsv(outfile, precision = 16)

  #ggplot(df.filter(f{`gs` < 2.5e-40}), aes("gs")) +
  #  geom_histogram(bins = 1000) +
  #  ggsave("/tmp/histo_test.pdf")

import unchained
type
  LimitData = object
    m_a: eV
    bins: int
    upper: float # upper cut value (99.9 percentile of prebinned data)
    gs: seq[float] # coupling constants
    histo: seq[int]

proc massScan(fname: string, outfile = "/tmp/test.csv") =
  var h5f = H5open(fname, "r")
  var lData: Table[eV, LimitData]
  lData = deserializeH5[Table[eV, LimitData]](h5f, "lData", "/")

  var dfs = newSeq[DataFrame]()
  for k, v in lData:
    dfs.add toDf({"m_a [eV]" : k.float, "gs" : v.gs, "L" : v.histo})
  let df = dfs.assignStack()

  df.arrange(["m_a [eV]", "gs"]).writeCsv(outfile, precision = 16)

when isMainModule:
  import cligen
  dispatchMulti([limit], [massScan])
