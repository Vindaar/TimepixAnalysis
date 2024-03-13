import ingrid / [tos_helpers, ingrid_types]
import nimhdf5, ggplotnim, sequtils, cligen, os

proc main(fname: string) =
  let runNumber = 101
  let searchRadius = 50
  let epsilon = 65
  let clusterAlgo = caDefault


  let h5f = H5file(fname, "r")
  const grpName = "/reconstruction/run_101"
  let group = h5f[grpName.grp_str]
  let centerChip = group.attrs["centerChip", int]
  let numChips = group.attrs["numChips", int]
  let septemDf = h5f.getSeptemEventDF(runNumber)


  # now filter events for `centerChip` from and compare with `passedInds`
  let centerDf = septemDf.filter(f{int: `chipNumber` == centerChip})
  #echo "Center df ", centerDf
  var
    allDataX: seq[seq[seq[uint8]]]
    allDataY: seq[seq[seq[uint8]]]
    allDataToT: seq[seq[seq[uint16]]]
    allDataCh: seq[seq[seq[float]]]
  let vlenXY = special_type(uint8)
  let vlenCh = special_type(float64)
    #allData
  for i in 0 ..< numChips:
    allDataX.add h5f[group.name / "chip_" & $i / "x", vlenXY, uint8]
    allDataY.add h5f[group.name / "chip_" & $i / "y", vlenXY, uint8]
    allDataToT.add h5f[group.name / "chip_" & $i / "ToT", vlenCh, uint16]
    allDataCh.add h5f[group.name / "chip_" & $i / "charge", vlenCh, float]
  let lhoodCenter = h5f[group.name / "chip_3" / "likelihood", float]
  let energyCenter = h5f[group.name / "chip_3" / "energyFromCharge", float]
  let cXCenter = h5f[group.name / "chip_3" / "centerX", float]
  let cYCenter = h5f[group.name / "chip_3" / "centerY", float]
  let hitsCenter = h5f[group.name / "chip_3" / "hits", int]
  let rmsTCenter = h5f[group.name / "chip_3" / "rmsTransverse", float]
  let rmsLCenter = h5f[group.name / "chip_3" / "rmsLongitudinal", float]

  let chips = toSeq(0 .. 6)
  let gains = chips.mapIt(h5f[(group.name / "chip_" & $it / "gasGainSlices"), GasGainIntervalResult])
  let energies = h5f[group.name / "chip_3" / "energyFromCharge", float]
  let septemGrouped = septemDf.group_by("eventNumber")
  for (pair, evGroup) in groups(septemGrouped):
    let evNum = pair[0][1]
    if not (evNum.toInt == 5339 and runNumber == 101): continue
    echo evNum
    echo runNumber
    # then grab all chips for this event
    var septemFrame: PixelsInt
    for row in evGroup:
      let chip = row["chipNumber"].toInt
      let idx = row["eventIndex"].toInt
      let
        chX = allDataX[chip][idx]
        chY = allDataY[chip][idx]
        chToT = allDataToT[chip][idx]
        chCh = allDataCh[chip][idx]
      let numPix = chX.len
      var chpPix = newSeq[Pix](numPix)
      for i in 0 ..< numPix:
        chpPix[i] = (x: chX[i], y: chY[i], ch: chToT[i])
      # convert to septem coordinate and add to frame
      septemFrame.add chpPix.chpPixToSeptemPix(chip)

    # given the full frame run through the full reconstruction for this cluster
    # here we give chip number as -1, indicating "Septem"
    let recoEv = recoEvent((septemFrame, evNum.toInt.int), -1,
                           runNumber, searchRadius = searchRadius,
                           dbscanEpsilon = epsilon.float,
                           clusterAlgo = clusterAlgo)[]
    echo recoEv
    echo recoEv.cluster.len
    doAssert recoEv.cluster.len == 3, "It's actually only 3 clusters and the plotting is at fault here!"

when isMainModule:
  dispatch main
