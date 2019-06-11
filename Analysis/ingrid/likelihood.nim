import os
import docopt
import tables
import strutils, strformat, ospaths
import algorithm
import sets
import nimhdf5
import tos_helpers
import helpers/utils
import sequtils
import seqmath
import loopfusion
import plotly
import ingrid_types

let doc = """
InGrid likelihood calculator. This program is run after reconstruction is finished.
It calculates the likelihood values for each reconstructed cluster of a run and
writes them back to the H5 file

Usage:
  likelihood <HDF5file> [options]
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
  --tracking             If flag is set, we only consider solar trackings (signal like)
  --scintiveto           If flag is set, we use the scintillators as a veto
  --fadcveto             If flag is set, we use the FADC as a veto
  --septemveto           If flag is set, we use the Septemboard as a veto
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
  histTuple = tuple[bins: seq[float64], hist: seq[float64]]

  FlagKind = enum
    fkTracking, fkFadc, fkScinti, fkSeptem

proc splitSeq[T, U](s: seq[seq[T]], dtype: typedesc[U]): (seq[U], seq[U]) =
  ## splits a (N, 2) nested seq into two seqs
  result[0] = newSeq[dtype](s.len)
  result[1] = newSeq[dtype](s.len)
  for i in 0..s.high:
    result[0][i] = s[i][0].U
    result[1][i] = s[i][1].U

proc buildLogLHist(h5file, dset: string, region: ChipRegion = crGold): seq[float] =
  ## given a file `h5file` containing a CDL calibration dataset
  ## `dset` apply the cuts on all events and build the logL distribution
  ## for the energy range
  ## `dset` needs to be of the elements contained in the returned map
  ## of tos_helpers.`getXrayRefTable`
  ## Default `region` is the gold region
  result = @[]
  var grp_name = cdlPrefix() & dset
  # create global vars for xray and normal cuts table to avoid having
  # to recreate them each time
  let xrayCutsTab {.global.} = getXraySpectrumCutVals()
  let cutsTab {.global.} = getEnergyBinMinMaxVals()

  withH5(h5file, "r"):
    # open h5 file using template
    let
      energy = h5f[(grp_name / "EnergyFromCharge"), float32]
      logL = h5f[(grp_name / "LikelihoodMarlin"), float32]
      centerX = h5f[(grp_name / "PositionX"), float32]
      centerY = h5f[(grp_name / "PositionY"), float32]
      ecc = h5f[(grp_name / "Excentricity"), float32]
      length = h5f[(grp_name / "Length"), float32]
      charge = h5f[(grp_name / "TotalCharge"), float32]
      rmsTrans = h5f[(grp_name / "RmsTransverse"), float32]
      npix = h5f[(grp_name / "NumberOfPixels"), float32]
      # get the cut values for this dataset
      cuts = cutsTab[dset]
      xrayCuts = xrayCutsTab[dset]
    for i in 0 .. energy.high:
      let
        # first apply Xray cuts (see C. Krieger PhD Appendix B & C)
        regionCut = inRegion(centerX[i], centerY[i], crSilver)
        xRmsCut = if rmsTrans[i] >= xrayCuts.minRms and
                     rmsTrans[i] <= xrayCuts.maxRms:
                    true
                  else:
                    false
        xLengthCut = if length[i] <= xrayCuts.maxLength: true else: false
        xEccCut = if ecc[i] <= xrayCuts.maxEccentricity: true else: false
        # then apply reference cuts
        chargeCut = if charge[i]   > cuts.minCharge and charge[i]   < cuts.maxCharge: true else: false
        rmsCut    = if rmsTrans[i] > cuts.minRms    and rmsTrans[i] < cuts.maxRms:    true else: false
        lengthCut = if length[i] < cuts.maxLength: true else: false
        pixelCut  = if npix[i]   > cuts.minPix:    true else: false
      # add event to likelihood if all cuts passed
      if allIt([regionCut, xRmsCut, xLengthCut, xEccCut, chargeCut, rmsCut, lengthCut, pixelCut], it):
        result.add logL[i]

  # now create plots of all ref likelihood distributions
  echo min(result), " ", max(result)
  if "4kV" in dset or "0.6kV" in dset:
    let binSize = (30.0 + 0.0) / 200.0
    histPlot(result)
      .binRange(0.0, 30.0)
      .binSize(binSize)
      .title(dset & " for region: " & $region)
      .show()

proc determineCutValue[T](hist: seq[T], eff: float): int =
  ## given a histogram `hist`, determine the correct bin to cut at to achieve
  ## a software efficiency of `eff`
  var
    cur_eff = 0.0
    last_eff = 0.0
  let hist_sum = hist.sum.float
  while cur_eff < eff:
    inc result
    last_eff = cur_eff
    cur_eff = hist[0..result].sum.float / hist_sum
  echo "Efficiency is at ", cur_eff, " and last Eff was ", last_eff

proc calcCutValueTab(region: ChipRegion = crGold): Table[string, float] =
  ## returns a table mapping the different CDL datasets to the correct cut values
  ## based on a chip `region`
  const
    xray_ref = getXrayRefTable()
    # software eff of 80%
    efficiency = 0.8
    # logL binning range
    nbins = 200 # NO CHANGE IF SET TO 200
    # range of histogram in logL
    logLrange = (0.0, 30.0)

  when true:
    let
      # get the raw log likelihood histograms (seq of individual values). Takes into account
      # cuts on properties and chip region
      rawLogHists = mapIt(toSeq(values(xray_ref)), buildLogLHist(h5cdl_file, it, region))
      # given raw log histograms, create the correctly binned histograms from it
      logHists = mapIt(rawLogHists, histogram(it, nbins, logLrange)[0])
      # get the cut value for a software efficiency of 80%
      cutVals = mapIt(logHists, determineCutValue(it, efficiency))
      # get the correct binning for the histograms
      bins = linspace(logLrange[0], logLrange[1], nbins + 1, endpoint = true)

    result = initTable[string, float]()
    for key, dset in pairs(xray_ref):
      # incl the correct values for the logL cut values
      result[dset] = bins[cutVals[key]]
    # some debugging output
  else:
    result = getChristophCutVals()
  when not defined(release) or defined(DEBUG):
    #echo logHists
    echo "Bins are ", bins
    echo "Cut values are ", cutVals
    echo mapIt(logHists, it.sum)
    echo "Corresponding to logL values of ", result

proc calcLogLikelihood*(h5f: var H5FileObj, ref_file: string) =
  ##
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

  var h5ref = H5file(ref_file, "r")
  # create a table, which stores the reference datasets from the ref file
  const xray_ref = getXrayRefTable()

  var
    ecc_ref = initTable[string, histTuple]()
    lengthDivRmsTrans_ref = initTable[string, histTuple]()
    fracRmsTrans_ref = initTable[string, histTuple]()
  for dset_name in values(xray_ref):
    var
      ecc = h5ref[(dset_name / "excentricity").dset_str]
      ldivrms = h5ref[(dset_name / "lengthdivbyrmsy").dset_str]
      frmst = h5ref[(dset_name / "fractionwithinrmsy").dset_str]

    # to get the reference datasets, we read from the H5DataSet, reshape it to
    # a (N, 2) nested seq and finally split the two columns into two individual
    # sequences converted to float64
    ecc_ref[dset_name] = ecc[float32].reshape2D(ecc.shape).splitSeq(float64)
    lengthDivRmsTrans_ref[dset_name] = ldivrms[float32].reshape2D(ldivrms.shape).splitSeq(float64)
    fracRmsTrans_ref[dset_name] = frmst[float32].reshape2D(frmst.shape).splitSeq(float64)

  # get the group from file
  for num, group in runs(h5f):
    echo &"Start logL calc of run {group}"
    # get number of chips from attributes
    var run_attrs = h5f[group.grp_str].attrs
    let nChips = run_attrs["numChips", int]

    var logL_chips = newSeq[seq[float64]](nChips)

    for grp in items(h5f, group):
      # iterate over all chips and perform logL calcs
      if "fadc" in grp.name:
        continue
      var attrs = grp.attrs
      let
        # get chip specific dsets
        chip_number = attrs["chipNumber", int]
        # get the datasets needed for LogL
        ecc = h5f[(grp.name / "eccentricity"), float64]
        lengthDivRmsTrans = h5f[(grp.name / "lengthDivRmsTrans"), float64]
        fracRmsTrans = h5f[(grp.name / "fractionInTransverseRms"), float64]
        # energy to choose correct bin
        energies = h5f[(grp.name / "energyFromCharge"), float64]

      # create seq to store data logL data for this chip
      var logL_s = newSeq[float64](ecc.len)
      for i in 0 .. ecc.high:
        var logL = 0.0'f64
        # try simple logL calc
        let refset = toRefDset(energies[i]) # / 1000.0) division needed for E from Pix,
                                            # since that is in eV inst of keV
        logL += logLikelihood(ecc_ref[refset][1],
                              ecc[i],
                              ecc_ref[refset][0])
        logL += logLikelihood(lengthDivRmsTrans_ref[refset][1],
                              lengthDivRmsTrans[i],
                              lengthDivRmsTrans_ref[refset][0])
        logL += logLikelihood(fracRmsTrans_ref[refset][1],
                              fracRmsTrans[i],
                              fracRmsTrans_ref[refset][0])
        logL *= -1.0
        # add logL to the sequence. May be Inf though
        logL_s[i] = logL
        # if logL != Inf:
        #   discard
      # after walking over all events for this chip, add to correct
      # index for logL
      logL_chips[chip_number] = logL_s
    # after we added all logL data to the seqs, write it to the file
    var logL_dsets = mapIt(toSeq(0..<nChips), h5f.create_dataset((group / &"chip_{it}/likelihood"),
                                                                 (logL_chips[it].len, 1),
                                                                 float64))
    # write the data to the file
    echo &"Writing data of run {group} to file {h5f.name}"
    #echo logL_dsets
    # forEach dset in var logL_dsets, logL in var logL_chips:
    #   dset[dset.all] = logL
    for tup in zip(logL_dsets, logL_chips):
      var (dset, logL) = tup
      dset[dset.all] = logL

proc writeLikelihoodData(h5f: var H5FileObj,
                         h5fout: var H5FileObj,
                         group: var H5Group,
                         chipNumber: int,
                         cutTab: Table[string, float],
                         passedInds: HashSet[int]) =
                         #durations: (float64, float64)) =
  ## writes all relevant cluster data of events corresponding to `passedInds` in
  ## the group given by `group` for chip `chipNumber` to the output file `h5f`

  when false:
    let (totalDuration, totalDurationPassed) = durations

  # TODO: add copy of attributes from energyFromPixel dataset, which contains
  # the calibration factor!
  # read all float datasets, which we want to write to the output file
  var float_dset_names = @(getFloatDsetNames())
  # add the final two datasets, which we'd like to write
  float_dset_names.add "likelihood"
  float_dset_names.add "energyFromCharge"
  var float_data_tab = initTable[string, seq[float]]()
  let chpGrpName = group.name / &"chip_{chipNumber}"
  # get mutable group for this chip to copy attributes
  var chpGrpIn = h5f[chpGrpName.grp_str]
  # fill table of float data sets
  for dset in float_dset_names:
    float_data_tab[dset] = h5f[(chpGrpName / dset), float64]

  # now get the event numbers to compare against the indices
  let evNumbers = h5f[(chpGrpName / "eventNumber"), int64]
  var float_data_passed = initTable[string, seq[float]]()
  # create new seqs of correct size in float_data_passed
  for dset in keys(float_data_tab):
    float_data_passed[dset] = newSeqOfCap[float](passedInds.card)
    for i in 0 .. float_data_tab[dset].high:
      if i in passedInds:
        # add element of index `i` to "passed data"
        float_data_passed[dset].add float_data_tab[dset][i]
  # then write to the new likelihood group for this chip
  # create the groups for the likelihood file, run group and chip group
  let
    runGrpName = group.name.replace("reconstruction", "likelihood")
    # group name of this chip's group
    logLgroup = runGrpName / &"chip_{chip_number}"
  var
    runGrp = h5fout.create_group(runGrpName)
    chpGrpOut = h5fout.create_group(logLgroup)

  # got all datasets ready for write
  for dset_name in keys(float_data_passed):
    var dset = h5fout.create_dataset((logLgroup / dset_name),
                                     (float_data_passed[dset_name].len, 1),
                                     float64)
    # write the data to the file
    dset[dset.all] = float_data_passed[dset_name]

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
  for key, val in pairs(cutTab):
    chpGrpOut.attrs[&"logL cut value: {key}"] = val
    chpGrpOut.attrs["SpectrumType"] = "background"
  # write total number of clusters
  chpGrpOut.attrs["Total number of cluster"] = evNumbers.len
  when false:
    # still missing per chip information
    runGrp.attrs["totalDurationRun"] = totalDuration
    runGrp.attrs["totalDurationPassed"] = totalDurationPassed
  # copy attributes over from the input file
  runGrp.copy_attributes(group.attrs)
  chpGrpOut.copy_attributes(chpGrpIn.attrs)

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

proc filterClustersByLogL(h5f: var H5FileObj, h5fout: var H5FileObj,
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
  let cutTab = calcCutValueTab(region)
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

  var logLSeq = newSeq[float64]()
  var logLSeqLow = newSeq[float64]()
  for num, group in runs(h5f):
    echo &"Start logL cutting of run {group}"
    # get number of chips from attributes
    var mgrp = h5f[group.grp_str]
    var run_attrs = mgrp.attrs
    let nChips = run_attrs["numChips", int]
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
    if fkSeptem in flags:
      raise newException(Exception, "Septemboard cuts not implemented at the " &
        "moment. To use it use the --septemVeto switch on the background rate " &
        "plot creation tool!")

    for chpGrp in items(h5f, group):
      if "fadc" in chpGrp.name:
        continue

      var fadcVetoCount = 0
      var scintiVetoCount = 0

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
        energy = h5f[(chpGrp.name / "energyFromCharge"), float64]
        logL = h5f[(chpGrp.name / "likelihood"), float64]
        centerX = h5f[(chpGrp.name / "centerX"), float64]
        centerY = h5f[(chpGrp.name / "centerY"), float64]
        rmsTrans = h5f[(chpGrp.name / "rmsTransverse"), float64]
        evNumbers = h5f[(chpGrp.name / "eventNumber"), int64].asType(int)

      # get event numbers corresponding to tracking (non tracking)
      var eventsInTracking: seq[int]
      if fkTracking in flags:
        eventsInTracking = h5f.getTrackingEvents(mgrp, tracking = true)
      else:
        eventsInTracking = h5f.getTrackingEvents(mgrp, tracking = false)

      # get all events part of tracking (non tracking)
      let indicesInTracking = filterTrackingEvents(evNumbers, eventsInTracking)

      if chipNumber == 0:
        for i, e in energy:
          let energyDset = e.toRefDset
          if "4kV" in energyDset:
            logLSeq.add logL[i]
          if "0.6kV" in energyDset:
            logLSeqLow.add logL[i]
      # hash set containing all indices of clusters, which pass the cuts
      var passedInds = initSet[int]()
      # iterate through all clusters not part of tracking and apply logL cut
      for ind in indicesInTracking:
        let dset = energy[ind].toRefDset
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
        if logL[ind] <= cutTab[dset] and
           inCutRegion and
           rmsTrans[ind] <= RmsCleaningCut and
           not fadcVeto and # if veto is true, means throw out!
           not scintiVeto:
          # include this index to the set of indices
          when false:
            totalDurationRunPassed += evDurations[ind]
          passedInds.incl ind

      chpGrp.writeVetoInfos(fadcVetoCount, scintiVetoCount, flags)
      # create dataset to store it
      if passedInds.card > 0:
        # call function which handles writing the data
        h5f.writeLikelihoodData(h5fout,
                                mgrp,
                                chipNumber,
                                cutTab,
                                passedInds)
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

  let binSize = (30.0 + 0.0) / 200.0
  let bins = (0.0, 30.0)
  histPlot(logLSeq)
    .binRange(bins[0], bins[1])
    .binSize(binSize)
    .title("Background logL for all runs 4kV Al")
    .show()
  histPlot(logLSeqLow)
    .binRange(bins[0], bins[1])
    .binSize(binSize)
    .title("Background logL for all runs 0.6kV C")
    .show()

  # once done write total duration as attributes to `likelihood` group
  var lhGrp = h5fout[likelihoodGroupGrpStr()]
  when false:
    for key, val in totalDurations:
      lhGrp.attrs["totalDurationChip_" & $key] = val
    for key, val in totalDurationsPassed:
      lhGrp.attrs["totalDurationPassedChip_" & $key] = val
  else:
    # simply take full duration of all events
    lhGrp.attrs["totalDuration"] = totalDuration

proc extractEvents(h5f: var H5FileObj, extractFrom, outfolder: string) =
  ## extracts all events passing the likelihood cut from the folder
  ## ``extractFrom`` and copies them (plus potential FADC files) to
  ## the ``outfolder``
  # - iterate over groups
  # - get path of run folder
  # - iterate over chips
  # - get dataset
  # - for each dataset have set of event numbers
  for grp in items(h5f, start_path = "likelihood", depth = 1):
    echo grp.name
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

proc main() =
  # create command line arguments
  let args = docopt(doc)
  let
    h5f_file = $args["<HDF5file>"]
    extractFrom = $args["--extract"]

  var flags: set[FlagKind]
  if $args["--tracking"] == "true": flags.incl fkTracking
  if $args["--scintiveto"] == "true": flags.incl fkScinti
  if $args["--fadcveto"] == "true": flags.incl fkFadc
  if $args["--septemveto"] == "true": flags.incl fkSeptem

  let region = if $args["--region"] != "nil":
                 parseEnum[ChipRegion]($args["--region"])
               else:
                 crGold

  var h5foutfile: string = ""
  if $args["--h5out"] != "nil":
    h5foutfile = $args["--h5out"]

  var h5f = H5file(h5f_file, "rw")
  h5f.visitFile

  if extractFrom == "nil":
    var h5fout: H5FileObj
    if h5foutfile != "":
      h5fout = H5file(h5foutfile, "rw")
    else:
      # in case no outfile given, we write to the same file
      h5fout = h5f

    if fkFadc in flags:
      echo "Using FADC as veto"
    if fkScinti in flags:
      echo "Using scintillators as veto"

    # perform likelihood calculation
    h5f.calcLogLikelihood(XrayRefFile)
    # now perform the cut on the logL values stored in `h5f` and write
    # the results to h5fout
    h5f.filterClustersByLogL(h5fout, flags, region)
    # given the cut values and the likelihood values for all events
    # based on the X-ray reference distributions, we can now cut away
    # all events not passing the cuts :)
  else:
    # extract all event numbers of the runs from the H5 file, check
    # existing in FOLDER/Run_???_* and copy to `outfolder`
    let outfolder = $args["--to"]
    h5f.extractEvents(extractFrom, outfolder)

  let err = h5f.close()
  if err != 0:
    echo &"Could not close h5 file properly! Return value was {err}"

  echo "Writing of all chips done"


when isMainModule:
  main()
