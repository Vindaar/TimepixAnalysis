import os, strutils, strformat, sequtils, sugar

import ggplotnim, shell

if true: # paramCount() > 1 and paramStr(1) == "--fit":
  let files = toSeq(walkFiles("out/lhood_*.h5"))
  let f2017s = files.filterIt("2017" in it)
  let f2018s = files.filterIt("2018" in it)
  echo f2017s
  for (f2017, f2018) in zip(f2017s, f2018s):
    let eff = f2017.dup(removePrefix("out/lhood_2017_eff_")).dup(removeSuffix(".h5")).parseFloat
    let res = shellVerbose:
      limitCalculation ($f2017) ($f2018) "--axionModel ../../../AxionElectronLimit/axion_gae_1e13_gagamma_1e-12_flux_after_exp_N_1000000.csv --optimizeBy CLs --outfile /tmp/limit_result.csv --eff" ($eff)

let df = readCsvTyped("/tmp/limit_result.csv")
echo df
ggplot(df, aes("eff", "limit")) +
  geom_point() +
  ggtitle("Expected limit vs. signal efficiency") +
  ggsave("out/limit_vs_signal_eff.pdf")
