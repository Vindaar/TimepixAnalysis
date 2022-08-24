import ingrid / [tos_helpers, ingrid_types]
import helpers / utils
import nimhdf5, seqmath, ggplotnim
import std / [unittest, strutils, os, strformat]

proc checkSeq(s1, s2: seq[float], isBins = true) =
  check s1.len == s2.len
  for i in 0 ..< s1.len:
    if isBins:
      check almostEqual(s1[i], s2[i])
    else:
      ## NOTE: the pre-generated CDL data and the on-the-fly computed data differ
      ## because the cuts are used slightly differently. In case of the pregenerated
      ## data the cuts were applied using `cutOnProperties` which effectively leaves
      ## all data *including* the bounds, whereas the new code uses the
      ## `withCdlData` template, which cuts *excluding* the bound values.
      let sMax = max(s1[i], s2[i])
      let sMin = min(s1[i], s2[i])
      if sMax > 0.0:
        check abs(1.0 - sMin / sMax) < 0.075 or sMax < 100.0
      else:
        check sMin == sMax

proc checkRefData(cdlFile, refFile: string, yearKind: YearKind) =
  ## reads the reference datasets from the `refFile` and returns them.
  var h5ref = H5open(refFile, "r")
  # create a table, which stores the reference datasets from the ref file
  const xray_ref = getXrayRefTable()

  # check if `FrameworkKind` defined in root of `h5ref` and if it is use the naming
  # scheme from there, otherwise use fkMarlin
  var frameworkKind: FrameworkKind = fkMarlin
  if "FrameworkKind" in h5ref.attrs:
    frameworkKind = parseEnum[FrameworkKind](h5ref.attrs["FrameworkKind", string])

  var
    ecc_ref = initTable[string, histTuple]()
    lengthDivRmsTrans_ref = initTable[string, histTuple]()
    fracRmsTrans_ref = initTable[string, histTuple]()

  var
    ## TODO: the following uses the `toDset` proc, which does not make sense for the original
    ## `XrayReferenceFile.h5` file, since that uses notation different from both normal Marlin
    ## and TPA. Also the `toDset` proc correctly fails for `igLengthDivRmsTrans` field, because
    ## it does not exist in Marlin. That of course is not the case for the reference file, where
    ## it is stored as ``lengthdivbyrmsy``.
    eccStr: string
    ldivrmsStr: string
    frmstStr: string
  case frameworkKind
  of fkMarlin:
    eccStr = "excentricity"
    ldivrmsStr = "lengthdivbyrmsy"
    frmstStr = "fractionwithinrmsy"
  of fkTpa:
    eccStr = igEccentricity.toDset(frameworkKind)
    ldivrmsStr = igLengthDivRmsTrans.toDset(frameworkKind)
    frmstStr = igFractionInTransverseRms.toDset(frameworkKind)

  var df: DataFrame
  for dset_name in values(xray_ref):
    # naming scheme does not depend on the actual data being processed, but only on what was used to
    # generate the `XrayReferenceFile.h5`
    var
      ecc = h5ref[(dset_name / eccStr).dset_str]
      ldivrms = h5ref[(dset_name / ldivrmsStr).dset_str]
      frmst = h5ref[(dset_name / frmstStr).dset_str]

    # to get the reference datasets, we read from the H5DataSet, reshape it to
    # a (N, 2) nested seq and finally split the two columns into two individual
    # sequences converted to float64
    ecc_ref[dset_name] = ecc.readAs(float64).reshape2D(ecc.shape).splitSeq(float64)
    lengthDivRmsTrans_ref[dset_name] = ldivrms.readAs(float64).reshape2D(ldivrms.shape).splitSeq(float64)
    fracRmsTrans_ref[dset_name] = frmst.readAs(float64).reshape2D(frmst.shape).splitSeq(float64)

    # now test refset stuff
    let (ecc_man, ldiv_man, frac_man) = genRefData(cdlFile, dsetName, yearKind, igEnergyFromCharge)

    template addDf(e, l, f: histTuple, key: string) =
      var dfDset = toDf({ "Eccentricity" : e.bins, "Ecc #" : e.hist,
                          "L / RMS_trans" : l.bins,
                          "L / RMS_trans #" : l.hist,
                          "fracRmsTrans" : f.bins,
                          "fracRmsTrans #" : f.hist })
      dfDset["Dset"] = dset_name
      dfDset["Algo"] = key
      df.add dfDset
    addDf(ecc_ref[dset_name], lengthDivRmsTrans_ref[dset_name], fracRmsTrans_ref[dset_name], "precomp")
    addDf((bins: ecc_man[0], hist: ecc_man[1]), (bins: ldiv_man[0], hist: ldiv_man[1]),
          (bins: frac_man[0], hist: frac_man[1]), "manual")


    ## The actual tests
    checkSeq ecc_ref[dset_name].bins, ecc_man[0]
    checkSeq lengthDivRmsTrans_ref[dset_name].bins, ldiv_man[0]
    checkSeq fracRmsTrans_ref[dset_name].bins, frac_man[0]

    checkSeq ecc_ref[dset_name].hist, ecc_man[1], false
    checkSeq lengthDivRmsTrans_ref[dset_name].hist, ldiv_man[1], false
    checkSeq fracRmsTrans_ref[dset_name].hist, frac_man[1], false


  #let xrayRef = getXrayRefTable()
  var labelOrder = initTable[Value, int]()
  for idx, el in xrayRef:
    labelOrder[%~ el] = idx
  ggplot(df, aes("Eccentricity", "Ecc #", fill = "Algo")) +
    ggridges("Dset", overlap = 1.75, labelOrder = labelOrder) +
    geom_histogram(stat = "identity", position = "identity", alpha = 0.5, hdKind = hdOutline) +
    ggtitle(&"Eccentricity of reference file, year: {yearKind}") +
    ggsave(&"out/eccentricity_ridgeline_{refFile.extractFilename}_{yearKind}.pdf",
            width = 800, height = 480)
  ggplot(df, aes("L / RMS_trans", "L / RMS_trans #", fill = "Algo")) +
    ggridges("Dset", overlap = 1.75, labelOrder = labelOrder) +
    geom_histogram(stat = "identity", position = "identity", alpha = 0.5, hdKind = hdOutline) +
    ggtitle(&"L / RMS_trans of reference file, year: {yearKind}") +
    ggsave(&"out/lengthDivRmsTrans_ridgeline_{refFile.extractFilename}_{yearKind}.pdf",
            width = 800, height = 480)
  ggplot(data = df.filter(f{Value: isNull(df["fracRmsTrans"][idx]) == (%~ false)}),
         aes("fracRmsTrans", "fracRmsTrans #", fill = "Algo")) +
    ggridges("Dset", overlap = 1.75, labelOrder = labelOrder) +
    geom_histogram(stat = "identity", position = "identity", alpha = 0.5, hdKind = hdOutline) +
    ggtitle(&"fracRmsTrans of reference file, year: {yearKind}") +
    ggsave(&"out/fracRmsTrans_ridgeline_{refFile.extractFilename}_{yearKind}.pdf",
            width = 800, height = 480)

let cdlFile = "/mnt/1TB/CAST/CDL_2019/calibration-cdl-2018.h5"
let refFile = "/mnt/1TB/CAST/CDL_2019/XrayReferenceFile2018.h5"
checkRefData(cdlFile, refFile, yr2018)
