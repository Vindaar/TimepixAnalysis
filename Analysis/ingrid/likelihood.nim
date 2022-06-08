import std / [os, tables, strutils, strformat, algorithm, sets,
              stats, sequtils, typetraits]
import docopt except Value
import nimhdf5, tos_helpers, seqmath, arraymancer
import helpers/utils
import ingrid / [ingrid_types, calibration]
import ingrid/calibration/fit_functions
from ingrid / private / geometry import recoEvent
import ingrid / private / likelihood_utils
import ingridDatabase / [databaseRead, databaseDefinitions]
import parsetoml except `{}`

let doc = """
InGrid likelihood calculator. This program is run after reconstruction is finished.
It calculates the likelihood values for each reconstructed cluster of a run and
writes them back to the H5 file

Usage:
  likelihood <HDF5file> [options]
  likelihood [options]
  likelihood <HDF5file> --h5out <outfile> [options]
  likelihood <HDF5file> --extract=FOLDER --to=OUTFOLDER

Options:
  --h5out <outfile>      The H5 file in which we store the events passing logL cut
  --extract=FOLDER       Given a H5 file created by this program under the
                         h5out argument, we will extract all raw event files,
                         which are contained in it from the given FOLDER.
                         Useful to e.g. plot passed events with the event display.
  --to=FOLDER            Output location of all extracted events. Events will just
                         be copied there.
  --region=REGION        The chip region to which we cut.
  --cdlYear=YEAR         The year from which to use the CDL data (2014, 2018).
                         Default for now is *2014*!
  --altCdlFile=CDLFILE   Alternative CDL file to use instead of `calibration-cdl`
                         in `resources`.
  --altRefFile=RefFILE   Alternative XrayRef file to use instead of `XrayReferenceFile`
                         in `resources`.
  --energyDset=NAME      Name of the energy dataset to use. By default `energyFromCharge`.
  --tracking             If flag is set, we only consider solar trackings (signal like)
  --scintiveto           If flag is set, we use the scintillators as a veto
  --fadcveto             If flag is set, we use the FADC as a veto
  --septemveto           If flag is set, we use the Septemboard as a veto
  --lineveto             If flag is set, we use an additional septem veto based on eccentric clusters
  --aggressive           If set, use aggressive veto. DO NOT USE (unless as a *reference*. Requires deep thought
                         about random coincidences & dead time of detector!)
  --plotSeptem           If flag is set, plots the SeptemEvents of all center clusters passing logL cut
  --createRocCurve       If flag is set, we create ROC curves for all energy bins. This
                         requires the input to already have a `likelihood` dataset!
  --plotLogL             If flag is set, we only plot the signal logL distributions.
  --computeLogL          If flag is set, we compute the logL dataset for each run in the
                         the input file. This is only required once or after changes to the
                         property datasets (e.g. energy calibration changed).
  -h --help              Show this help
  --version              Show version.
"""

const h5cdl_file = currentSourcePath() / "../../../resources/calibration-cdl.h5"
const XrayRefFile = currentSourcePath() / "../../../resources/XrayReferenceDataSet.h5"
const cdlExists = fileExists(h5cdl_file)
when not cdlExists:
  {.fatal: "CAST CDL reference file `calibration-cdl.h5` does not exist at: " &
    $h5cdl_file.}
const refExists = fileExists(XrayRefFile)
when not refExists:
  {.fatal: "X-ray reference file `XrayReferenceDataSet.h5` does not exist at: " &
    $XrayRefFile.}

# cut performed regardless of logL value on the data, since transverse
# rms > 1.5 cannot be a physical photon, due to diffusion in 3cm drift
# distance
const RmsCleaningCut = 1.5

type
  FlagKind = enum
    fkTracking, fkFadc, fkScinti, fkSeptem, fkLineVeto, fkAggressive, fkRocCurve, fkComputeLogL, fkPlotLogL, fkPlotSeptem

template withConfig(body: untyped): untyped =
  const sourceDir = currentSourcePath().parentDir
  let config {.inject.} = parseToml.parseFile(sourceDir / "config.toml")
  body

proc readSignalEff(): float =
  ## reads the `signalEfficiency` field from the TOML file
  withConfig:
    result = config["Likelihood"]["signalEfficiency"].getFloat

proc readMorphKind(): MorphingKind =
  ## reads the `morphingKind` field from the TOML file
  withConfig:
    result = parseEnum[MorphingKind](config["Likelihood"]["morphingKind"].getStr)

proc readSearchRadius(): int =
  ## reads the cluster finding `searchRadius` field from the TOML file
  withConfig:
    result = config["Likelihood"]["searchRadius"].getInt

proc readClusterAlgo(): ClusteringAlgorithm =
  ## Reads the clustering algorithm to use for the septem veto clustering
  withConfig:
    result = parseEnum[ClusteringAlgorithm](config["Likelihood"]["clusterAlgo"].getStr)

proc readDbscanEpsilon(): float =
  ## Reads the `ε` to use for the DBSCAN algorithm septem veto clustering
  withConfig:
    result = config["Likelihood"]["epsilon"].getFloat

proc writeLogLDsetAttributes[T: H5DataSet | H5Group](dset: var T,
                             cdlFile, refFile: string,
                             year: YearKind) =
  ## writes information about what datasets were used to calculate the likelihood
  ## dataset
  dset.attrs["year of CDL data"] = $year
  dset.attrs["calibration CDL file"] = cdlFile
  dset.attrs["X-ray reference file"] = refFile

proc calcLogLikelihood*(h5f: var H5File,
                        cdlFile, refFile: string,
                        year: YearKind,
                        energyDset: InGridDsetKind # pixel / charge energy
                       ) =
  ## - read all data of single run
  ## - get energy dataset
  ## - create energy bins
  ## - use digitize to get the correct energy bins for each cluster
  ## - create individual seq's for each energy bin (containing all
  ##   clusters and the corresponding properties)
  ## - create histogram for each of these energy binned properties

  ## in XrayLikelihoodProcessor we have
  ## - 1 histogram for each property and each energy bin
  ## - fill them with the log likelihood of all events falling into that histogram
  ##   - logLikelihood:
  ##     get number of elements in histogram at the bin for the element for which we get
  ##     the logL
  ##     # elements / # total in histogram = likelihood. Take log
  let morphKind = readMorphKind()
  for (num, group) in runs(h5f):
    echo &"Start logL calc of run {group}"
    # get number of chips from attributes
    var run_attrs = h5f[group.grp_str].attrs
    let nChips = run_attrs["numChips", int]

    var logL_chips = newSeq[seq[float64]](nChips)

    for (_, chipNumber, grp) in chipGroups(h5f, group):
      # iterate over all chips and perform logL calcs
      var attrs = h5f[grp.grp_str].attrs
      let logL = calcLikelihoodDataset(h5f, refFile, grp, year, morphKind, energyDset)
      # after walking over all events for this chip, add to correct
      # index for logL
      logL_chips[chipNumber] = logL
    # after we added all logL data to the seqs, write it to the file
    var logL_dsets = mapIt(toSeq(0..<nChips), h5f.create_dataset((group / &"chip_{it}/likelihood"),
                                                                 (logL_chips[it].len, 1),
                                                                 float64,
                                                                 overwrite = true))
    # write the data to the file
    echo &"Writing data of run {group} to file {h5f.name}"
    for tup in zip(logL_dsets, logL_chips):
      var (dset, logL) = tup
      dset[dset.all] = logL
      dset.writeLogLDsetAttributes(cdlFile, refFile, year)

proc writeLikelihoodData(h5f: var H5File,
                         h5fout: var H5File,
                         group: var H5Group,
                         cdlFile, refFile: string,
                         year: YearKind,
                         chipNumber: int,
                         cutTab: CutValueInterpolator,
                         passedInds: OrderedSet[int],
                         energyDset: InGridDsetKind, # pixel or charge
                         version: TimepixVersion # determines if we need to write ToA data
                        ) =
                         #durations: (float64, float64)) =
  ## writes all relevant cluster data of events corresponding to `passedInds` in
  ## the group given by `group` for chip `chipNumber` to the output file `h5f`

  when false:
    let (totalDuration, totalDurationPassed) = durations

  # TODO: add copy of attributes from energyFromPixel dataset, which contains
  # the calibration factor!
  # read all float datasets, which we want to write to the output file
  var dsetNames = @(getFloatDsetNames())
  # add the final two datasets, which we'd like to write
  dsetNames.add @["likelihood", energyDset.toDset, "x", "y", "ToT"]
  if version == Timepix3:
    dsetNames.add getFloatToANames()
    dsetNames.add @["ToA", "ToACombined", "toaMin"]

  let chpGrpName = group.name / &"chip_{chipNumber}"
  # get mutable group for this chip to copy attributes
  var chpGrpIn = h5f[chpGrpName.grp_str]

  # now get the event numbers to compare against the indices
  let evNumbers = h5f[(chpGrpName / "eventNumber"), int64]
  # then write to the new likelihood group for this chip
  # create the groups for the likelihood file, run group and chip group
  let
    runGrpName = group.name.replace("reconstruction", "likelihood")
    # group name of this chip's group
    logLgroup = runGrpName / &"chip_{chip_number}"
  var
    runGrp = h5fout.create_group(runGrpName)
    chpGrpOut = h5fout.create_group(logLgroup)
  ## walk all datasets and read data, filter to indices and write to output
  for dsetName in dsetNames: # all datasets including vlen
    let h5dset = h5f[(chpGrpName / dsetName).dset_str]
    withDset(h5dset): # reads full dataset into injected `dset` variable
      var data = newSeqOfCap[elementType(dset)](passedInds.card)
      for idx in passedInds:
        data.add dset[idx]
      var outDset = h5fout.create_dataset((logLgroup / dsetName),
                                          (passedInds.len, 1),
                                          elementType(dset))
      outDset[outDset.all] = data

  # get all event numbers from hash set by using elements as indices for event numbers
  let evNumsPassed = mapIt(passedInds, evNumbers[it]).sorted do (x, y: int64) -> int:
    result = cmp(x, y)
  # create dataset for allowed indices
  var evDset = h5fout.create_dataset((logLgroup / "eventNumber"),
                                     (evNumsPassed.len, 1),
                                     int)
  # write event numbers
  evDset[evDset.all] = evNumsPassed

  # finally write all interesting attributes
  chpGrpOut.attrs["MorphingKind"] = $cutTab.kind
  if cutTab.kind == mkNone:
    for key, val in pairs(cutTab.cutTab):
      chpGrpOut.attrs[&"logL cut value: {key}"] = val
      chpGrpOut.attrs["SpectrumType"] = "background"
  else:
    ## TODO: add dataset of energies & cut values for cut
    echo "WARN: Add information about used energies and cut values for mkLinear!"
    discard

  # write total number of clusters
  chpGrpOut.attrs["Total number of cluster"] = evNumbers.len
  when false:
    # still missing per chip information
    runGrp.attrs["totalDurationRun"] = totalDuration
    runGrp.attrs["totalDurationPassed"] = totalDurationPassed
  # copy attributes over from the input file
  runGrp.copy_attributes(group.attrs)
  chpGrpOut.copy_attributes(chpGrpIn.attrs)
  runGrp.writeLogLDsetAttributes(cdlFile, refFile, year)

func isVetoedByFadc(eventNumber: int, fadcTrigger, fadcEvNum: seq[int64],
                    fadcRise, fadcFall: seq[uint16]): bool =
  ## returns `true` if the event of `ind` is vetoed by the FADC based on cuts on
  ## the rise and fall time. Vetoed means the event must be thrown out
  ## because it does ``not`` conform to the X-ray hypothesis.
  ## ------------ NOTE --------------
  ## looking at ~/org/Figs/SPSC_Jan_2019/test/Run3/{rise,fall}Time_normalizedPDF_run3.pdf
  ## makes it seem like anything above ~130 is probably not an X-ray. Take that for
  ## now.
  ## TODO: CHOOSE THESE VALUE MORE WISELY!!!!!!
  const cutRiseHigh = 130'u16
  const cutRiseLow = 40'u16
  const cutFallLow = 400'u16
  const cutFallHigh = 600'u16
  result = false
  let fIdx = fadcEvNum.lowerBound(eventNumber)
  if fIdx < fadcEvNum.high and
     fadcEvNum[fIdx] == eventNumber and
     fadcTrigger[fIdx] == 1:
    # thus we know that `fIdx` points to an event with an FADC trigger
    # corresponding to `eventNumber`
    if fadcRise[fIdx] >= cutRiseLow and
       fadcRise[fIdx] >= cutRiseHigh and
       fadcFall[fIdx] >= cutFallLow and
       fadcFall[fIdx] >= cutFallHigh:
      result = false
    else:
      # outside either of the cuts, throw it out
      result = true

func isVetoedByScintis(eventNumber: int,
                       scintEvNum: seq[int64],
                       scinti1, scinti2: seq[int64]): bool =
  ## returns `true` if the event of `ind` is vetoed by the scintillators based
  ## on the fact that one of the two scintillators had a non trivial scintillator
  ## count value ( > 0 and < 4095; in practice ~< 400).
  ## Vetoed means the event must be thrown out, because the event was most
  ## likely induced by a muon
  result = false
  # throw out any event with a non trivial (> 0 and < 4095)
  # scintillator trigger
  const low = 0
  const high = 400 # be pessimistic about high
  let sIdx = scintEvNum.lowerBound(eventNumber)
  if sIdx < scintEvNum.high and
     scintEvNum[sIdx] == eventNumber and
     ((scinti1[sIdx] > low and scinti1[sIdx] < high) or
      (scinti2[sIdx] > low and scinti2[sIdx] < high)):
    # had a non trivial trigger, throw out
    result = true

proc writeVetoInfos(grp: H5Group, fadcVetoCount, scintiVetoCount: int,
                    flags: set[FlagKind]) =
  ## writes information about used vetoes and the number of events removed by
  ## the vetos
  var mgrp = grp
  mgrp.attrs["FADC Veto"] = $(fkFadc in flags)
  mgrp.attrs["Scinti Veto"] = $(fkScinti in flags)
  mgrp.attrs["# removed by FADC veto"] = fadcVetoCount
  mgrp.attrs["# removed by scinti veto"] = scintiVetoCount

# subtract from 768 as we do the same in `applyPitchConversion` (why?)
# normalize by full width (14 mm * 3 chips) and scale to all pixels
proc toXPix(x: float): int =
  clamp((768 - (x / (TimepixSize * 3.0)) * 768.0).int, 0, 767)

proc toYPix(y: float): int =
  clamp(((y / (TimepixSize * 3.0)) * 768.0).int, 0, 767)

type
  SeptemFrame = object
    pixels: PixelsInt     ## pure pixel data of all pixels in the septem frame (zero suppressed)
    charge: Tensor[float] ## charge values of all pixels in a [768,768] tensor
    centerEvIdx: int      ## index of the center event
    numRecoPixels: int    ## number of pixels ``after`` reconstruction. To truncate `pixels` after recoEvent

  AllChipData = object
    x: seq[seq[seq[uint8]]]
    y: seq[seq[seq[uint8]]]
    ToT: seq[seq[seq[uint16]]]
    charge: seq[seq[seq[float]]]

  CenterChipData = object
    lhoodCenter: seq[float]
    energies: seq[float]
    energyCenter: seq[float]
    cXCenter: seq[float]
    cYCenter: seq[float]
    hitsCenter: seq[int]
    rmsTCenter: seq[float]
    rmsLCenter: seq[float]

  ## The name might not be the most descriptive: This object stores information about different
  ## geometric properties of the separate clusters found in the septem event. This includes
  ## information about the original "center" event (the one that passed the logL cut on the center
  ## chip) (`center` in the field name) and sequences for the centers of each other cluster
  ## as well as lines that go through the centers of each cluster along the long axis.
  SeptemEventGeometry = object
    lines: seq[tuple[m, b: float]]   ## lines along the long axis through the center of each cluster
    centers: seq[tuple[x, y: float]] ## centers of each found cluster in the septem event
    xCenter: int ## x center pixel of the original "center event" (that passed logL)
    yCenter: int ## y center pixel of the original "center event" (that passed logL)
    centerRadius: float ## "radius" of the original "center event" (based on 3 * (RMS_T + RMS_L)/2 in pixel)

  DataView[T] = distinct ptr UncheckedArray[T]

func `[]`[T](dv: DataView[T], idx: int): T = cast[ptr UncheckedArray[T]](dv)[idx]
#proc `=copy`[T](dv1, dv2: DataView[T]) {.error: "Copying a data view is forbidden!".}
func toDataView[T](p: ptr seq[T]): DataView[T] = DataView[T](cast[ptr UncheckedArray[T]](p[0].addr))

proc evaluateCluster(clTup: (int, ClusterObject[PixInt]),
                     rs: var RunningStat, septemFrame: var SeptemFrame,
                     septemGeometry: var SeptemEventGeometry,
                     centerData: CenterChipData,
                     gainVals: seq[float],
                     calibTuple: tuple[b, m: float], ## contains the parameters required to perform energy calibration
                     region: ChipRegion,
                     refSetTuple: tuple[ecc, ldivRms, fracRms: Table[string, histTuple]],
                     flags: set[FlagKind]
                    ): tuple[logL, energy: float, lineVetoPassed: bool] =
  # total charge for this cluster
  let clusterId = clTup[0]
  let cl = clTup[1]
  let clData = cl.data
  # reset running stat for each cluster
  rs.clear()
  var pixIdx = septemFrame.numRecoPixels
  var totCharge = 0.0
  for pix in clData:
    # take each pixel tuple and reconvert it to chip based coordinates
    # first determine chip it corresponds to
    let pixChip = determineChip(pix)

    # for this pixel and chip get the correct gas gain and push it to the `RunningStat`
    rs.push gainVals[pixChip]
    # taken the chip of the pixel, reconvert that to a local coordinate system
    # given charge of this pixel, assign it to some intermediate storage
    totCharge += septemFrame.charge[pix.y, pix.x]
    # overwrite the `septemFrame` by `pix` and cluster id
    septemFrame.pixels[pixIdx] = (x: pix.x, y: pix.y, ch: clusterId)
    inc pixIdx
  septemFrame.numRecoPixels = pixIdx

  # determine parameters of the lines through the cluster centers
  # invert the slope
  let slope = tan(Pi - cl.geometry.rotationAngle)
  let cX = toXPix(cl.centerX)
  let cY = toYPix(cl.centerY)
  let intercept = cY.float - slope * cX.float
  septemGeometry.centers.add (x: cX.float, y: cY.float)
  septemGeometry.lines.add (m: slope, b: intercept)

  # using total charge and `RunningStat` calculate energy from charge
  ## TODO: need to look *only* at gains from the corresponding chips each
  ## and compute the energy of the part of each chip indidually and add
  let (b, m) = calibTuple
  let energy = totCharge * linearFunc(@[b, m], rs.mean) * 1e-6

  #let energy = cl.hits.float * 26.0
  let logL = calcLikelihoodForEvent(energy, # <- TODO: hack for now!
                                    cl.geometry.eccentricity,
                                    cl.geometry.lengthDivRmsTrans,
                                    cl.geometry.fractionInTransverseRms,
                                    refSetTuple)
  ## Check if the current cluster is in input chip region. If it is, either it is part of something
  ## super big that makes the center still fall into the gold region or it remains unchanged.
  ## In the unchanged case, let's compare the energy and cluster pixels
  let inRegionOfInterest = inRegion(cl.centerX - TimepixSize, cl.centerY - TimepixSize, region)

  var lineVetoPassed = true #
  if fkLineVeto in flags and not inRegionOfInterest: ## perform the line cut
    # if this is not in region of interest, check its eccentricity and compute if line points to original
    # center cluster
    let centerEvIdx = septemFrame.centerEvIdx
    # compute line orthogonal to this cluster's line
    let orthSlope = -1.0 / slope
    ## TODO: XXX: eccentricity check. Using a line of a cluster that is not very eccentric makes no sense!
    # determine y axis intersection
    septemGeometry.xCenter = clamp(256 - (centerData.cXCenter[centerEvIdx] / TimepixSize * 256.0).int, 0, 255) + 256
    septemGeometry.yCenter = clamp((centerData.cYCenter[centerEvIdx] / TimepixSize * 256).int, 0, 256) + 256

    let centerInter = septemGeometry.yCenter.float - orthSlope * septemGeometry.xCenter.float

    # use orthogonal line to determine crossing point between it and this cluster's line
    let xCut = (intercept - centerInter) / (orthSlope - slope)
    let yCut = orthSlope * xCut + centerInter
    septemGeometry.centers.add (x: xCut, y: yCut)
    # compute distance between two points
    let dist = sqrt( (septemGeometry.xCenter.float - xCut)^2 + (septemGeometry.yCenter.float - yCut)^2 )
    # if distance is smaller than radius
    septemGeometry.centerRadius = ((centerData.rmsTCenter[centerEvIdx] +
                                    centerData.rmsLCenter[centerEvIdx]) /
                                   2.0 * 3.0) /
                                   TimepixSize * 256.0
    if dist < septemGeometry.centerRadius:
      lineVetoPassed = false
    septemGeometry.lines.add (m: orthSlope, b: centerInter)
  result = (logL: logL,
            energy: energy,
            lineVetoPassed: lineVetoPassed) ## No line veto is equivalent to passing the line veto.
                                            ## Init to `true`, just use variable

proc readAllChipData(h5f: H5File, group: H5Group, numChips: int): AllChipData =
  ## Read all data for all chips of this run that we need for the septem veto
  let vlenXY = special_type(uint8)
  let vlenCh = special_type(float64)
  result = AllChipData(x: newSeq[seq[seq[uint8]]](numChips),
                       y: newSeq[seq[seq[uint8]]](numChips),
                       ToT: newSeq[seq[seq[uint16]]](numChips),
                       charge: newSeq[seq[seq[float]]](numChips))
  for i in 0 ..< numChips:
    result.x[i] =      h5f[group.name / "chip_" & $i / "x", vlenXY, uint8]
    result.y[i] =      h5f[group.name / "chip_" & $i / "y", vlenXY, uint8]
    result.ToT[i] =    h5f[group.name / "chip_" & $i / "ToT", vlenCh, uint16]
    result.charge[i] = h5f[group.name / "chip_" & $i / "charge", vlenCh, float]

proc readCenterChipData(h5f: H5File, group: H5Group, energyDset: InGridDsetKind): CenterChipData =
  ## Read all data of this run for the center chip that we need in the rest of
  ## the septem veto logic
  result.lhoodCenter = h5f[group.name / "chip_3" / "likelihood", float]
  result.energies = h5f[group.name / "chip_3" / energyDset.toDset, float]
  result.energyCenter = h5f[group.name / "chip_3" / energyDset.toDset, float]
  result.cXCenter = h5f[group.name / "chip_3" / "centerX", float]
  result.cYCenter = h5f[group.name / "chip_3" / "centerY", float]
  result.hitsCenter = h5f[group.name / "chip_3" / "hits", int]
  result.rmsTCenter = h5f[group.name / "chip_3" / "rmsTransverse", float]
  result.rmsLCenter = h5f[group.name / "chip_3" / "rmsLongitudinal", float]

proc buildSeptemEvent(evDf: DataFrame,
                      lhoodCenter, energies: seq[float],
                      cutTab: CutValueInterpolator,
                      allChipData: AllChipData,
                      centerChip: int): SeptemFrame =
  ## Given a
  result = SeptemFrame(centerEvIdx: -1, numRecoPixels: 0)
  # assign later to avoid indirection for each pixel
  var chargeTensor = zeros[float]([768, 768])
  var pixels: PixelsInt = newSeqOfCap[PixInt](5000) # 5000 = 5k * 26eV = 130keV. Plenty to avoid reallocations
  let chipData = evDf["chipNumber", int] # get data as tensors and access to avoid `Value` wrapping
  let idxData = evDf["eventIndex", int]
  for i in 0 ..< evDf.len:
    # get the chip number and event index, dump corresponding event pixel data
    # onto the "SeptemFrame"
    let chip = chipData[i]
    let idx = idxData[i]
    if chip == centerChip and lhoodCenter[idx.int] < cutTab[energies[idx.int]]:
      # assign center event index if this is the cluster that passes logL cut
      result.centerEvIdx = idx.int

    ## Use `idx` of this cluster & event to look up all x, y, ToT and charge data of the cluster
    let
      chX = allChipData.x[chip][idx].unsafeAddr.toDataView()
      chY = allChipData.y[chip][idx].unsafeAddr.toDataView()
      chToT = allChipData.ToT[chip][idx].unsafeAddr.toDataView()
      chCh = allChipData.charge[chip][idx].unsafeAddr.toDataView()
    let numPix = allChipData.x[chip][idx].len # data view has no length knowlegde!
    var chpPix = newSeq[PixInt](numPix)
    for j in 0 ..< numPix:
      let pix = (x: chX[j], y: chY[j], ch: chToT[j]).chpPixToSeptemPix(chip)
      # add current charge into full septem tensor
      chargeTensor[pix.y, pix.x] += chCh[j]
      chpPix[j] = pix
    # convert to septem coordinate and add to frame
    pixels.add chpPix
  result.charge = chargeTensor
  result.pixels = pixels

proc applySeptemVeto(h5f, h5fout: var H5File,
                     cdlFile, refFile: string,
                     runNumber: int,
                     year: YearKind,
                     energyDset: InGridDsetKind, # pixel or charge
                     passedInds: var OrderedSet[int],
                     cutTab: CutValueInterpolator,
                     version: TimepixVersion,
                     region: ChipRegion,
                     flags: set[FlagKind]) =
  ## Applies the septem board veto to the given `passedInds` in `runNumber` of `h5f`.
  ## Writes the resulting clusters, which pass to the `septem` subgroup (parallel to
  ## the `chip_*` groups into `h5fout`.
  ## If an event does not pass the septem veto cut, it is excluded from the `passedInds`
  ## set.
  const PlotCutEnergy = 5.0
  echo "Passed indices before septem veto ", passedInds.card
  let group = h5f[(recoBase() & $runNumber).grp_str]
  let centerChip = group.attrs["centerChip", int]
  let numChips = group.attrs["numChips", int]
  let septemDf = h5f.getSeptemEventDF(runNumber)

  # read clustering configuration
  let clusterAlgo = readClusterAlgo()
  # search radius for regular algo
  let searchRadius = readSearchRadius()
  # epsilon fro DBSCAN
  let epsilon = readDbscanEpsilon()

  # now filter events for `centerChip` from and compare with `passedInds`
  let centerDf = septemDf.filter(f{int: `chipNumber` == centerChip})
  #echo "Center df ", centerDf
  ## TODO: this is super slow!!
  let passedEvs = passedInds.mapIt(centerDf["eventNumber", it]).sorted.toOrderedSet

  ## Read all the pixel data for all chips
  let allChipData = readAllChipData(h5f, group, numChips)

  let centerData = readCenterChipData(h5f, group, energyDset)

  let refSetTuple = readRefDsets(refFile, year)
  let chips = toSeq(0 .. 6)
  let gains = chips.mapIt(h5f[(group.name / "chip_" & $it / "gasGainSlices"), GasGainIntervalResult])
  let septemHChips = chips.mapIt(getSeptemHChip(it))
  let toTCalibParams = septemHChips.mapIt(getTotCalibParameters(it, runNumber))
  let calibTuple = getCalibVsGasGainFactors(septemHChips[centerChip], runNumber)
  var rs: RunningStat

  # for the `passedEvs` we have to read all data from all chips
  let septemGrouped = septemDf.group_by("eventNumber")
  for (pair, evGroup) in groups(septemGrouped):
    let evNum = pair[0][1]
    if evNum in passedEvs:
      # then grab all chips for this event
      var septemFrame = buildSeptemEvent(evGroup, centerData.lhoodCenter, centerData.energies,
                                         cutTab, allChipData, centerChip)
      if septemFrame.centerEvIdx == -1:
        echo "Broken event! ", evGroup.pretty(-1)
        quit(1)
      # given the full frame run through the full reconstruction for this cluster
      # here we give chip number as -1, indicating "Septem"
      ## XXX: for full ToA support in Timepix3, need to also read the `toa` data and insert
      ## it here!
      let inp = (pixels: septemFrame.pixels, eventNumber: evNum.toInt.int,
                 toa: newSeq[uint16](), toaCombined: newSeq[uint64]())
      let recoEv = recoEvent(inp, -1,
                             runNumber, searchRadius = searchRadius,
                             dbscanEpsilon = epsilon,
                             clusterAlgo = clusterAlgo)

      # extract the correct gas gain slices for this event
      var gainVals: seq[float]
      for chp in 0 ..< numChips:
        let gainSlices = gains[chp]
        let gainEvs = gainSlices.mapIt(it.sliceStartEvNum)
        let sliceIdx = gainEvs.lowerBound(evNum.toInt)
        gainVals.add gainSlices[min(gainSlices.high, sliceIdx)].G # add the correct gas gain

      # calculate log likelihood of all reconstructed clusters
      var passed = false
      var lineVetoRejected = false
      var septemGeometry: SeptemEventGeometry # no need for constructor. `default` is fine
      if fkAggressive in flags:
        # if there's more than 1 cluster, remove
        if recoEv.cluster.len > 1:
          passedInds.excl septemFrame.centerEvIdx
          continue # skip to next iteration
      for clusterTup in pairs(recoEv.cluster):
        let (logL, energy, lineVetoPassed) = evaluateCluster(clusterTup, rs, septemFrame, septemGeometry,
                                                             centerData,
                                                             gainVals,
                                                             calibTuple,
                                                             region,
                                                             refSetTuple, flags)
        let cX = toXPix(clusterTup[1].centerX)
        let cY = toYPix(clusterTup[1].centerY)
        let chipClusterCenter = (x: cX, y: cY, ch: 0).determineChip(allowOutsideChip = true)
        if chipClusterCenter == centerChip and # this cluster's center is on center chip
           logL < cutTab[energy]:              # cluster passes logL cut
             passed = true
        if not lineVetoRejected and not lineVetoPassed:
          lineVetoRejected = true
      if not passed or lineVetoRejected:
        ## If `passed` is still false, it means *no* reconstructed cluster passed the logL now. Given that
        ## the original cluster which *did* pass logL is part of the septem event, the conclusion is that
        ## it was now part of a bigger cluster that did *not* pass anymore.
        passedInds.excl septemFrame.centerEvIdx

      if fkPlotSeptem in flags:
        if septemFrame.centerEvIdx < 0:
          doAssert false, "this cannot happen. it implies no cluster found in the given event"
        if centerData.energies[septemFrame.centerEvIdx] < PlotCutEnergy:
          # shorten to actual number of stored pixels. Otherwise elements with ToT / charge values will remain
          # in the `septemFrame`
          septemFrame.pixels.setLen(septemFrame.numRecoPixels)
          if fkLineVeto notin flags: # if no line veto, don't draw lines
            septemGeometry.lines = newSeq[tuple[m, b: float]]()
          plotSeptemEvent(septemFrame.pixels, runNumber, evNum.toInt,
                          lines = septemGeometry.lines,
                          centers = septemGeometry.centers,
                          passed = passed,
                          lineVetoRejected = lineVetoRejected,
                          xCenter = septemGeometry.xCenter,
                          yCenter = septemGeometry.yCenter,
                          radius = septemGeometry.centerRadius,
                          energyCenter = centerData.energies[septemFrame.centerEvIdx],
                          useTeX = true)
  echo "Passed indices after septem veto ", passedInds.card

proc filterClustersByLogL(h5f: var H5File, h5fout: var H5File,
                          cdlFile, refFile: string,
                          year: YearKind,
                          energyDset: InGridDsetKind, # pixel or charge
                          version: TimepixVersion,
                          flags: set[FlagKind],
                          region = crGold) =
  ## filters all clusters with a likelihood value in the given `h5f` by
  ## the logL cut values returned by `calcCutValueTab`
  ## clusters passing the cuts are stored in `h5fout`
  ## The `flags` argument decides what kind of filtering we perform aside from
  ## the logL cuts
  ## - fkTracking: only tracking data considered
  ## The remaining flags may not be available for all datasets!
  ## - fkFadc: FADC used as veto
  ## - fkScinti: Scintillators used as veto
  ## - fkSeptem: Septemboard used as veto
  # TODO: should the argument to calcCutValueTab not be crGold all the time?
  # We want to extract that data from the CDL data that most resembles the X-rays
  # we measured. This is guaranteed by using the gold region.
  let cutTab = calcCutValueTab(cdlFile, refFile, year, energyDset)
  # get the likelihood and energy datasets
  # get the group from file
  when false:
    # not yet supported, since no event duration for
    # each chip individually
    # Note: could only be calc'ed from correlating event numbers
    # with event durations, too much of a hassle RIGHT NOW
    # total duration in seconds of run time
    var
      totalDurations = initTable[int, float64]()
      totalDurationsPassed = initTable[int, float64]()
  else:
    # alternative, total duration of whole run
    var totalDuration: float64 = 0.0

  var useFadcVeto = fkFadc in flags
  var useScintiVeto = fkScinti in flags

  var
    totalScintiRemoveCount = 0
    totalScintiRemovedNotLogRemoved = 0
    totalEvCount = 0
    totalLogLCount = 0
  for num, group in runs(h5f):
    ## XXX: turn the following into a proper feature
    #if num != 261: continue
    echo &"Start logL cutting of run {group}"
    # get number of chips from attributes
    var mgrp = h5f[group.grp_str]
    var run_attrs = mgrp.attrs
    let nChips = run_attrs["numChips", int]
    let centerChip = run_attrs["centerChip", int]
    # get timestamp for run
    let tstamp = h5f[(group / "timestamp"), int64]
    let evDurations = h5f[group / "eventDuration", float64]
    let eventNumbers = h5f[group / "eventNumber", int64]
    # add sum of event durations
    totalDuration += evDurations.foldl(a + b, 0.0)

    var
      fadcTrigger: seq[int64]
      fadcRise: seq[uint16]
      fadcFall: seq[uint16]
      fadcEvNum: seq[int64]
      scinti1Trigger: seq[int64]
      scinti2Trigger: seq[int64]
    if fkFadc in flags:
      try:
        fadcTrigger = h5f[group / "fadcReadout", int64]
        fadcRise = h5f[group / "fadc/riseTime", uint16]
        fadcFall = h5f[group / "fadc/fallTime", uint16]
        fadcEvNum = h5f[group / "fadc/eventNumber", int64]
      except KeyError:
        echo "Run ", num, " has no FADC datasets!"
        useFadcVeto = false
    if fkScinti in flags:
      scinti1Trigger = h5f[group / "szint1ClockInt", int64]
      scinti2Trigger = h5f[group / "szint2ClockInt", int64]

    for (_, chipNumber, chipGroup) in chipGroups(h5f, group):
      if energyDset.toDset notin h5f[chipGroup.grp_str]:
        raise newException(IOError, "The input file " & $h5f.name & " does not contain the dataset " &
          energyDset.toDset() & " in the group: " & $chipGroup & ".")

      var fadcVetoCount = 0
      var scintiVetoCount = 0
      let chpGrp = h5f[chipGroup.grp_str]
      # iterate over all chips and perform logL calcs
      var attrs = chpGrp.attrs
      let
        # get chip specific dsets
        chipNumber = attrs["chipNumber", int]

      when false:
        # add duration for this chip to Duration table
        totalDurations[chipNumber] = 0.0
        totalDurationsPassed[chipNumber] = 0.0
        # vars for this chips durations
        # currently unsupported, since not part of reco files for
        # each chip group yet
        var
          totalDurationRun = 0.0
          totalDurationRunPassed = 0.0
      let
        # get the datasets needed for LogL
        energy = h5f[(chipGroup / energyDset.toDset), float64]
        logL = h5f[(chipGroup / "likelihood"), float64]
        centerX = h5f[(chipGroup / "centerX"), float64]
        centerY = h5f[(chipGroup / "centerY"), float64]
        rmsTrans = h5f[(chipGroup / "rmsTransverse"), float64]
        evNumbers = h5f[(chipGroup / "eventNumber"), int64].asType(int)

      # get event numbers corresponding to tracking (non tracking)
      var eventsInTracking: seq[int]
      if fkTracking in flags:
        eventsInTracking = h5f.getTrackingEvents(mgrp, tracking = true)
      else:
        eventsInTracking = h5f.getTrackingEvents(mgrp, tracking = false)

      # get all events part of tracking (non tracking)
      let indicesInTracking = filterTrackingEvents(evNumbers, eventsInTracking)
      if chipNumber == centerChip:
        totalEvCount += indicesInTracking.len

      # hash set containing all indices of clusters, which pass the cuts
      var passedInds = initOrderedSet[int]()
      # iterate through all clusters not part of tracking and apply logL cut
      for ind in indicesInTracking:
        when false:
          # see above, not yet implemented
          # add current ind to total duration; note before cut, since we want
          # total time detector was alive!
          totalDurationRun += evDurations[ind]

        var fadcVeto = false
        var scintiVeto = false
        if useFadcVeto:
          fadcVeto = isVetoedByFadc(evNumbers[ind], fadcTrigger, fadcEvNum,
                                    fadcRise, fadcFall)
        if fadcVeto:
          # increase if FADC vetoed this event
          inc fadcVetoCount
        if useScintiVeto:
          scintiVeto = isVetoedByScintis(evNumbers[ind], eventNumbers, scinti1Trigger,
                                         scinti2Trigger)
        if scintiVeto:
          # increase if Scintis vetoed this event
          inc scintiVetoCount

        # given datasest add element to dataset, iff it passes logL, region and
        # cleaning cut
        let inCutRegion = inRegion(centerX[ind], centerY[ind], region)
        if logL[ind] <= cutTab[energy[ind]] and
           inCutRegion and
           rmsTrans[ind] <= RmsCleaningCut and
           not fadcVeto and # if veto is true, means throw out!
           not scintiVeto:
          # include this index to the set of indices
          when false:
            totalDurationRunPassed += evDurations[ind]
          passedInds.incl ind

        if logL[ind] <= cutTab[energy[ind]] and inCutRegion and rmsTrans[ind] <= RmsCleaningCut and
           scintiVeto and chipNumber == centerChip:
          # only those events that otherwise wouldn't have made it by logL only
          inc totalScintiRemovedNotLogRemoved

      chpGrp.writeVetoInfos(fadcVetoCount, scintiVetoCount, flags)
      if chipNumber == centerChip:
        totalScintiRemoveCount += scintiVetoCount

      # create dataset to store it
      if passedInds.card > 0:
        # now in a second pass perform a septem veto if desired
        # If there's no events left, then we don't care about
        if fkSeptem in flags and chipNumber == centerChip:
          # read all data for other chips ``iff`` chip == 3 (centerChip):
          h5f.applySeptemVeto(h5fout, cdlFile, refFile, num, year,
                              energyDset,
                              passedInds,
                              cutTab = cutTab,
                              version = version,
                              region = region,
                              flags = flags)

        if passedInds.card > 0:
          # call function which handles writing the data
          h5f.writeLikelihoodData(h5fout,
                                  mgrp,
                                  cdlFile, refFile,
                                  year,
                                  chipNumber,
                                  cutTab,
                                  passedInds,
                                  energyDset,
                                  version)

        if chipNumber == centerChip:
          totalLogLCount += passedInds.card

        when false:
          (totalDurationRun, totalDurationRunPassed)
      else:
        var mchpGrp = chpGrp
        mchpGrp.attrs["LogLSpectrum"] = "No events passed cut"
        echo "No clusters found passing logL cut"

      when false:
        # finally add totalDuration to total duration vars
        totalDurations[chipNumber] += totalDurationRun
        totalDurationsPassed[chipNumber] += totalDurationRunPassed

  # once done write total duration as attributes to `likelihood` group
  var lhGrp = h5fout.getOrCreateGroup(likelihoodGroupGrpStr().string)
  when false:
    for key, val in totalDurations:
      lhGrp.attrs["totalDurationChip_" & $key] = val
    for key, val in totalDurationsPassed:
      lhGrp.attrs["totalDurationPassedChip_" & $key] = val
  else:
    # simply take full duration of all events
    lhGrp.attrs["totalDuration"] = totalDuration
    lhGrp.attrs["totalEvents"] = totalEvCount
    lhGrp.attrs["totalPassedEvents"] = totalLogLCount
    lhGrp.attrs["totalCutByScinti"] = totalScintiRemoveCount
    lhGrp.attrs["onlyCutByScinti"] = totalScintiRemovedNotLogRemoved
    lhGrp.attrs["MorphingKind"] = cutTab.kind
    ## TODO: add morphing kind to output!

  # write year and CDL and reference file used
  lhGrp.writeLogLDsetAttributes(cdlFile, refFile, year)

proc extractEvents(h5f: var H5File, extractFrom, outfolder: string) =
  ## extracts all events passing the likelihood cut from the folder
  ## ``extractFrom`` and copies them (plus potential FADC files) to
  ## the ``outfolder``
  # - iterate over groups
  # - get path of run folder
  # - iterate over chips
  # - get dataset
  # - for each dataset have set of event numbers
  for grp in items(h5f, start_path = "likelihood", depth = 1):
    var mgrp = grp
    let path = mgrp.attrs["pathName", string]
    let (head, tail) = path.splitPath
    let runFolder = joinPath(extractFrom, tail)
    if dirExists runFolder:
      # get all event numbers of all chips
      var eventNumbers = initSet[int64]()
      for chipGroup in items(h5f, start_path = mgrp.name, depth = 1):
        var mchip = chipGroup
        let dsetName = joinPath(mchip.name, "eventNumber").dset_str
        var dset = h5f[dsetName]
        let events = dset[int64]
        for ev in events:
          eventNumbers.incl ev
      # now copy over all files
      # iterate event numbers and copy files
      for ev in eventNumbers:
        let evFile = getFilenameFromEventNumber(ev)
        let infile = joinPath(runFolder, evFile)
        # check existence of corresponding FADC event
        let fadcExists = existsFile(infile & "-fadc")
        let outfile = joinPath(outfolder, evFile)
        # echo &"Copying {infile} to {outfile}"
        copyFile(infile, outfile)
        if fadcExists:
          copyFile(infile & "-fadc", outfile & "-fadc")

proc readLikelihoodDsets(h5f: H5File, energyDset: InGridDsetKind): DataFrame =
  ## reads all likelihood data in the given `h5f` file as well as the
  ## corresponding energies. Flattened to a 1D seq.
  ## This proc is for TPA generated H5 files! (i.e. containing run_* groups, ...)
  # iterate over all groups, read all likelihood and energy dsets
  var energies = newSeqOfCap[float](1_000_000)
  var logLs = newSeqOfCap[float](1_000_000)
  for run, grp in runs(h5f):
    let group = h5f[grp.grp_str]
    let centerChip = "chip_" & $group.attrs["centerChip", int]
    doAssert grp / centerChip / "likelihood" in h5f,
      "likelihood dataset must exist in input H5 file! Does not exist at path: " & $(grp / centerChip / "likelihood")
    let energy = h5f[grp / centerChip / energyDset.toDset, float64]
    let logL = h5f[grp / centerChip / "likelihood", float64]
    doAssert energy.len == logL.len
    energies.add energy
    logLs.add logL
  let bin_back = energies.mapIt(it.toRefDset)
  result = seqsToDf({ "Bin" : bin_back,
                      "Energy" : energies,
                      "Likelihood" : logLs })

proc readLikelihoodDsetsCdl(cdlFile, refFile: string,
                            yearKind: YearKind,
                            energyDset: InGridDsetKind,
                            region: ChipRegion): DataFrame =
  ## reads a CDL like H5 file and returns a DataFrame of the energies,
  ## likelihood values and categories (of the energy bin)
  # iterate over all groups, read all likelihood and energy dsets
  const xray_ref = getXrayRefTable()
  var
    energies: seq[float]
    logLs: seq[float]
    bins: seq[string]
  for bin in values(xray_ref):
    let (logL, energy) = buildLogLHist(cdlFile, refFile, bin, yearKind,
                                       energyDset, crGold)
    logLs.add logL
    energies.add energy
    bins.add sequtils.repeat(bin, energy.len)
  result = seqsToDf({ "Bin" : bins,
                      "Energy" : energies,
                      "Likelihood" : logLs })
  when false:
    # code to simply read all data w/o filtering.
    var energies = newSeqOfCap[float32](1_000_000)
    var logLs = newSeqOfCap[float32](1_000_000)
    var bins = newSeqOfCap[string](1_000_000)
    for grp in items(h5f):
      let energy = h5f[grp.name / "EnergyFromCharge", float32]
      let logL = h5f[grp.name / "LikelihoodMarlin", float32]
      doAssert energy.len == logL.len
      energies.add energy
      logLs.add logL
      proc removePref(s, prefix: string): string =
        result = s
        result.removePrefix(prefix)
      bins.add repeat(removePref(grp.name, "/" & cdlPrefix("2014")),
                      energy.len)
    let bin_back = energies.mapIt(it.toRefDset)
    result = seqsToDf({ "energy" : energies,
                        "logL" : logLs,
                        "bin" : bin_back })

proc determineEff(logLs: seq[float], cutVal: float,
                  isBackground = true): float =
  ## returns the efficiency given the sorted (!) `logLs`, a
  ## cut value `cutVal` and whether it's background or signal
  let cutIdx = logLs.lowerBound(cutVal)
  result = cutIdx.float / logLs.len.float
  if isBackground:
    result = 1.0 - result

proc calcSigEffBackRej(df: DataFrame, logLBins: seq[float],
                       isBackground = true): DataFrame =
  ## returns the signal eff and backround rej for all logLBins of the
  ## given data frame, split by the `bins` column (that is CDL classes)
  let dfG = df.group_by("Bin")
  for (pair, subDf) in groups(dfG):
    let logL = subDf.arrange("Likelihood")["Likelihood"].toTensor(float)
    var effs = newSeqOfCap[float](logLBins.len)
    for l in logLBins:
      let eff = determineEff(logL.toRawSeq, l, isBackground = isBackground)
      effs.add eff
    let binConst = toSeq(0 ..< effs.len).mapIt(pair[0][1].toStr)
    let effDf = seqsToDf({ "eff" : effs,
                           "cutVals" : logLBins,
                           "bin" : binConst })
    result.add effDf

proc calcRocCurve(dfSignal, dfBackground: DataFrame): DataFrame =
  # now use both to determine signal and background efficiencies
  # essentially have to generate some binning in `logL` we deem appropriate,
  const LogLBins = 500
  let logLBins = linspace(0.0, 40.0, LogLBins)
  let sigEffDf = calcSigEffBackRej(dfSignal, logLBins, isBackground = false)
    .rename(f{"sigEff" <- "eff"})
  let backRejDf = calcSigEffBackRej(dfBackground, logLBins, isBackground = true)
    .rename(f{"backRej" <- "eff"})
  result = innerJoin(sigEffDf, backRejDf, by = "cutVals")

proc createRocCurves(h5Back: H5File,
                     cdlFile, refFile: string,
                     yearKind: YearKind,
                     energyDset: InGridDsetKind, # pixel or charge
                     region: ChipRegion
                    ) =
  ## generates all ROC curves for the given two H5 files and the
  ## histograms of the likelihood distributions for the CDL data and
  ## the given background file.
  ## By default the file containing signal like events will be
  ## the X-ray reference file.
  let dfSignal = readLikelihoodDsetsCdl(cdlFile, refFile, yearKind, energyDset, region)
  let dfBack = readLikelihoodDsets(h5Back, energyDset)
    .filter(f{float: `Likelihood` != Inf})
  ggplot(dfBack, aes("Likelihood", fill = "Bin")) +
    geom_histogram(binWidth = 0.2) +
    ggtitle("-LnL distributions of non-tracking background, stacked",
            titlefont = font(11.0)) +
    ggsave("backgroundLogL.pdf")
  ggplot(dfSignal, aes("Likelihood", fill = "Bin")) +
    geom_histogram(binWidth = 0.2) +
    ggtitle("-LnL distributions of cdl calibration data, stacked",
            titlefont = font(11.0)) +
    ggsave("signalLogL.pdf")

  ggplot(dfBack, aes("Likelihood", fill = "Bin")) +
    geom_freqpoly(binWidth = 0.2,
                  position = "identity",
                  alpha = some(0.3)) +
    ggtitle("-LnL distributions of non tracking background as polygons, identity position",
            titlefont = font(11.0)) +
    ggsave("backgroundLogL_freqPoly.pdf")

  when false:
    # write the dfSignal data frame to file
    dfSignal.writeCsv("/tmp/dfSignal.csv")
  ggplot(dfSignal, aes("Likelihood", fill = "Bin")) +
    geom_freqpoly(binWidth = 0.2,
                  position = "identity",
                  alpha = some(0.3)) +
    ggtitle("-LnL distributions of cdl calibration data as polygons, identity position",
            titlefont = font(11.0)) +
    ggsave("signalLogL_freqPoly.pdf")

  let xrayRef = getXrayRefTable()
  var labelOrder = initTable[Value, int]()
  for idx, el in xrayRef:
    labelOrder[%~ el] = idx
  ggplot(dfBack, aes("Likelihood", fill = "Bin")) +
    geom_histogram(binWidth = 0.2,
                   position = "identity",
                   alpha = some(0.5)) +
    ggridges("Bin", overlap = 2.0,
             labelOrder = labelOrder) +
    ggtitle("-LnL distributions of non tracking background as ridgeline",
            titlefont = font(11.0)) +
    ggsave("backgroundLogL_ridgeline.pdf",
           height = 600.0)
  ggplot(dfSignal, aes("Likelihood", fill = "Bin")) +
    geom_histogram(binWidth = 0.2,
                   position = "identity",
                   alpha = some(0.5)) +
    ggridges("Bin", overlap = 2.0,
             labelOrder = labelOrder) +
    ggtitle("-LnL distributions of cdl calibration data as ridgeline",
            titlefont = font(11.0)) +
    ggsave("signalLogL_ridgeline.pdf",
           height = 600.0)


  ## TODO: IMPORTANT the results are still wrong!!! Especially the `0.9 Cu EPIC` line
  ## still is lower than it should be!
  # then determine efficiency in both signal and background
  # for computational efficiency reason only use raw, sorted `seq[float]`
  let res = calcRocCurve(dfSignal, readLikelihoodDsets(h5Back, energyDset))
  ggplot(res, aes("sigEff", "backRej", color = "bin")) +
    geom_line() +
    #ylim(0.8, 1.0) +
    ggtitle("ROC curves for likelihood method, 2014 data") +
    ggsave("roc_curves.pdf")
  #ggplot(effDf, aes("sigEff", "backRej")) +
  #  geom_line() +
  #  ggsave("roc_curve_full_range.pdf")

proc plotLogL(cdlFile, refFile: string,
              yearKind: YearKind,
              energyDset: InGridDsetKind,
              region: ChipRegion) =
  ## generates all ROC curves for the given two H5 files and the
  ## histograms of the likelihood distributions for the CDL data and
  ## the given background file.
  ## By default the file containing signal like events will be
  ## the X-ray reference file.
  let dfSignal = readLikelihoodDsetsCdl(cdlFile, refFile, yearKind, energyDset, region)
  ggplot(dfSignal, aes("Likelihood", fill = "Bin")) +
    geom_histogram(binWidth = 0.2) +
    ggtitle("-LnL distributions of cdl calibration data, stacked",
            titlefont = font(11.0)) +
    ggsave("signalLogL.pdf")

  when false:
    # write the dfSignal data frame to file
    dfSignal.writeCsv("/tmp/dfSignal.csv")
  ggplot(dfSignal, aes("Likelihood", fill = "Bin")) +
    geom_freqpoly(binWidth = 0.2,
                  position = "identity",
                  alpha = some(0.3)) +
    ggtitle("-LnL distributions of cdl calibration data as polygons, identity position",
            titlefont = font(11.0)) +
    ggsave("signalLogL_freqPoly.pdf")

  let xrayRef = getXrayRefTable()
  var labelOrder = initTable[Value, int]()
  for idx, el in xrayRef:
    labelOrder[%~ el] = idx
  ggplot(dfSignal, aes("Likelihood", fill = "Bin")) +
    geom_histogram(binWidth = 0.2,
                   position = "identity",
                   alpha = some(0.5)) +
    ggridges("Bin", overlap = 2.0,
             labelOrder = labelOrder) +
    ggtitle("-LnL distributions of cdl calibration data as ridgeline",
            titlefont = font(11.0)) +
    ggsave("signalLogL_ridgeline.pdf",
           height = 600.0)

# switch to cligen, then do:
# runs: seq[int] = @[]) = # `runs` allows to overwrite whihc run is logL cut'd
## Also do same for a `--useTeX` argument & flag! for Septem veto plots
proc main() =
  # create command line arguments
  let args = docopt(doc)
  echo args
  let
    h5f_file = $args["<HDF5file>"]
    extractFrom = $args["--extract"]

  var flags: set[FlagKind]
  if $args["--tracking"] == "true": flags.incl fkTracking
  if $args["--scintiveto"] == "true": flags.incl fkScinti
  if $args["--fadcveto"] == "true": flags.incl fkFadc
  if $args["--septemveto"] == "true": flags.incl fkSeptem
  if $args["--lineveto"] == "true": flags.incl fkLineVeto
  if $args["--aggressive"] == "true": flags.incl fkAggressive
  if $args["--createRocCurve"] == "true": flags.incl fkRocCurve
  if $args["--computeLogL"] == "true": flags.incl fkComputeLogL
  if $args["--plotLogL"] == "true": flags.incl fkPlotLogL
  if $args["--plotSeptem"] == "true": flags.incl fkPlotSeptem

  let cdlFile = if $args["--altCdlFile"] != "nil":
                  $args["--altCdlFile"]
                else:
                  h5cdl_file
  let refFile = if $args["--altRefFile"] != "nil":
                  $args["--altRefFile"]
                else:
                  XrayRefFile

  let region = if $args["--region"] != "nil":
                 parseEnum[ChipRegion]($args["--region"])
               else:
                 crGold

  let year = if $args["--cdlYear"] != "nil":
                   parseEnum[YearKind]($args["--cdlYear"])
                 else:
                   # default to 2014
                   yr2014
  let energyDset = if $args["--energyDset"] != "nil":
                     toIngridDset($args["--energyDset"])# , igEnergyFromCharge)
                   else:
                     igEnergyFromCharge
  doAssert energyDset != igInvalid, "Please enter a valid energy dataset. " &
    "Choices: {energyFromCharge, energyFromPixel}"

  var h5foutfile: string = ""
  if $args["--h5out"] != "nil":
    h5foutfile = $args["--h5out"]

  var h5f = H5open(h5f_file, "rw")
  h5f.visitFile
  let timepix = h5f.timepixVersion()

  if fkRocCurve in flags:
    ## create the ROC curves and likelihood distributios. This requires to
    ## previously run this tool with the default parameters
    createRocCurves(h5f, cdlFile, refFile, year, energyDset, region)
  if fkPlotLogL in flags:
    plotLogL(cdlFile, refFile, year, energyDset, region)
  if fkComputeLogL in flags:
    # perform likelihood calculation
    h5f.calcLogLikelihood(cdlFile, refFile, year, energyDset)
    h5f.flush() # flush so the data is written already
  if extractFrom == "nil" and h5foutfile.len > 0:
    var h5fout = H5open(h5foutfile, "rw", {akTruncate})

    if fkFadc in flags:
      echo "Using FADC as veto"
    if fkScinti in flags:
      echo "Using scintillators as veto"

    # now perform the cut on the logL values stored in `h5f` and write
    # the results to h5fout
    h5f.filterClustersByLogL(h5fout,
                             cdlFile,
                             refFile,
                             year, energyDset, timepix, flags, region)
    # given the cut values and the likelihood values for all events
    # based on the X-ray reference distributions, we can now cut away
    # all events not passing the cuts :)
  else:
    # extract all event numbers of the runs from the H5 file, check
    # existing in FOLDER/Run_???_* and copy to `outfolder`
    let outfolder = $args["--to"]
    h5f.extractEvents(extractFrom, outfolder)

  echo "Closing H5 file: ", h5f.name
  let err = h5f.close()
  if err != 0:
    echo &"Could not close h5 file properly! Return value was {err}"

  echo "Writing of all chips done"


when isMainModule:
  main()
