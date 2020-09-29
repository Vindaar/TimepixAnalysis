import strscans, os, strformat
import ggplotnim, nimhdf5
import ingrid / tos_helpers
import cligen

proc readCdlDsets(h5f: H5FileObj, fwKind: static FrameworkKind): DataFrame =
  var tfKindStr: string
  var year: string
  echo "Reading file: ", h5f.name
  for group in h5f:
    var mgrp = group
    if scanf(group.name, "/calibration-cdl-$+-$+kV", year, tfKindStr):
      let eccName = igEccentricity.toDset(fwKind)
      let frmsName = igFractionInTransverseRms.toDset(fwKind)
      let chName = igTotalCharge.toDset(fwKind)
      when fwKind == fkMarlin:
        let lenName = igLength.toDset(fwKind)
        let transRmsName = igRmsTransverse.toDset(fwKind)
        var df = seqsToDf({ "Eccentricity" : h5f[group.name / $eccName, float32],
                            "Length" : h5f[group.name / $lenName, float32],
                            "RmsTrans" : h5f[group.name / $transRmsName, float32],
                            "fracRmsTrans" : h5f[group.name / $frmsName, float32],
                            "totalCharge" : h5f[group.name / $chName, float32]})
        df = df.mutate(f{"L / RMS_trans" ~ `Length` / `RmsTrans`})
          .select("Eccentricity", "L / RMS_trans", "fracRmsTrans", "totalCharge")
      else:
        let ldivRmsName = igLengthDivRmsTrans.toDset(fwKind)
        var df = seqsToDf({ "Eccentricity" : h5f[group.name / $eccName, float64],
                            "L / RMS_trans" : h5f[group.name / $ldivRmsName, float64],
                            "fracRmsTrans" : h5f[group.name / $frmsName, float64],
                            "totalCharge" : h5f[group.name / $chName, float64]})
      df["Dataset"] = constantColumn(tfKindStr, df.len)
      result.add df
  result["Framework"] = constantColumn($fwKind, result.len)

proc main(marlin, tpa: string) =
  discard existsOrCreateDir("out")
  let marH5 = H5file(marlin, "r")
  let tpaH5 = H5file(tpa, "r")

  var df = marH5.readCdlDsets(fkMarlin)
  df.add tpaH5.readCdlDsets(fkTpa)
  # there are some events  with ridiculous charge values, filter them
  df = df.filter(f{`totalCharge` < 1e8},
                 f{`Eccentricity` < 1e2})
  df = df.group_by("Dataset")
  for (tup, subDf) in groups(df):
    let tfKind = tup[0][1]
    echo "Generating plots for ", tfKind
    # TODO: set the ranges based on percentiles?
    ggplot(subDf, aes("Eccentricity", fill = "Framework")) +
      geom_histogram(position = "identity", alpha = some(0.5)) +
      ggtitle(&"Eccentricity comparison of Marlin/TPA, dset: {tfKind}") +
      ggsave(&"out/eccentricity_marlin_tpa_comp_{tfKind}.pdf",
              width = 800, height = 480)
    ggplot(subDf, aes("L / RMS_trans", fill = "Framework")) +
      geom_histogram(position = "identity", alpha = some(0.5)) +
      ggtitle(&"L / RMS_trans comparison of Marlin/TPA, dset: {tfKind}") +
      ggsave(&"out/lengthDivRmsTrans_marlin_tpa_comp_{tfKind}.pdf",
              width = 800, height = 480)
    ggplot(subDf, aes("fracRmsTrans", fill = "Framework")) +
      geom_histogram(position = "identity", alpha = some(0.5)) +
      ggtitle(&"fracRmsTrans comparison of Marlin/TPA, dset: {tfKind}") +
      ggsave(&"out/fracRmsTrans_marlin_tpa_comp_{tfKind}.pdf",
              width = 800, height = 480)
    ggplot(subDf, aes("totalCharge", fill = "Framework")) +
      geom_histogram(position = "identity", alpha = some(0.5), bin_width = 1e4) +
      ggtitle(&"Total charge comparison of Marlin/TPA, dset: {tfKind}") +
      ggsave(&"out/total_charge_marlin_tpa_comp_{tfKind}.pdf",
              width = 800, height = 480)

when isMainModule:
  dispatch main
