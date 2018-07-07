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
import ingrid_types

let doc = """
InGrid likelihood calculator. This program is run after reconstruction is finished.
It calculates the likelihood values for each reconstructed cluster of a run and
writes them back to the H5 file

Usage:
  likelihood <HDF5file> --reference <ref_file> [options]
  likelihood <HDF5file> --reference <ref_file> --h5out <outfile> [options]

Options:
  --reference <ref_file> The H5 file containing the X-ray reference spectra
  --h5out <outfile>      The H5 file in which we store the events passing logL cut
  --tracking             If flag is set, we only consider solar trackings (signal like)
  -h --help              Show this help
  --version              Show version.
"""

const h5cdl_file = "/mnt/Daten/Uni/CAST/data/CDL-reference/calibration-cdl.h5"

type
  histTuple = tuple[bins: seq[float64], hist: seq[float64]]

proc splitSeq[T, U](s: seq[seq[T]], dtype: typedesc[U]): (seq[U], seq[U]) =
  ## splits a (N, 2) nested seq into two seqs
  result[0] = newSeq[dtype](s.len)
  result[1] = newSeq[dtype](s.len)
  for i in 0..s.high:
    result[0][i] = s[i][0].U
    result[1][i] = s[i][1].U

proc cutPosition(centerX, centerY: float, region: ChipRegion): bool =
  ## returns the result of a cut on a certain chip `region`. Inputs the
  ## `centerX` and `centerY` position of a cluster and returns true if
  ## the cluster is within the region
  const centerChip = 7.0
  # make sure this is only initialized once somehow...
  let regCut = getRegionCut(region)
  case region
  of crGold:
    result = if centerX >= regCut.xMin and
                centerX <= regCut.xMax and
                centerY >= regCut.yMin and
                centerY <= regCut.yMax:
               true
             else:
               false
  of crAll:
    # simply always return good
    result = true
  else:
    # silver and bronze region only different by radius
    let
      xdiff = (centerX - centerChip)
      ydiff = (centerY - centerChip)
      radius = distance(xdiff, ydiff)
    # TODO: should gold cut be allowed? i.e. is gold region part of silver region?
    result = if radius <= regCut.radius: true else : false


proc buildLogLHist(h5file, dset: string, region: ChipRegion = crGold): seq[float] =
  ## given a file `h5file` containing a CDL calibration dataset
  ## `dset` apply the cuts on all events and build the logL distribution
  ## for the energy range
  ## `dset` needs to be of the elements contained in the returned map
  ## of tos_helpers.`getXrayRefTable`
  ## Default `region` is the gold region
  result = @[]
  var grp_name = cdlPrefix() & dset
  withH5(h5file, "r"):
    # open h5 file using template
    let
      energy = h5f[(grp_name / "EnergyFromCharge").dset_str][float32]
      logL = h5f[(grp_name / "LikelihoodMarlin").dset_str][float32]
      centerX = h5f[(grp_name / "PositionX").dset_str][float32]
      centerY = h5f[(grp_name / "PositionY").dset_str][float32]
      length = h5f[(grp_name / "Length").dset_str][float32]    
      charge = h5f[(grp_name / "TotalCharge").dset_str][float32]
      rmsTrans = h5f[(grp_name / "RmsTransverse").dset_str][float32]
      npix = h5f[(grp_name / "NumberOfPixels").dset_str][float32]
      # get the cut values for this dataset
      cuts = getEnergyBinMinMaxVals()[dset]

    for i in 0 .. energy.high:
      let
        regionCut = cutPosition(centerX[i], centerY[i], region)
        chargeCut = if charge[i]   > cuts.minCharge and charge[i]   < cuts.maxCharge: true else: false
        rmsCut    = if rmsTrans[i] > cuts.minRms    and rmsTrans[i] < cuts.maxRms:    true else: false
        lengthCut = if length[i] < cuts.maxLength: true else: false
        pixelCut  = if npix[i]   > cuts.minPix:    true else: false
      # add event to likelihood if all cuts passed
      if allIt([regionCut, chargeCut, rmsCut, lengthCut, pixelCut], it):
        result.add logL[i]

proc determineCutValue[T](hist: seq[T], eff: float): int =
  ## given a histogram `hist`, determine the correct bin to cut at to achieve
  ## a software efficiency of `eff`
  var
    cur_eff = 0.0
  let hist_sum = hist.sum.float
  while cur_eff < eff:
    inc result
    cur_eff = hist[0..result].sum.float / hist_sum

proc calcCutValueTab(region: ChipRegion = crGold): Table[string, float] =
  ## returns a table mapping the different CDL datasets to the correct cut values
  ## based on a chip `region`
  const
    xray_ref = getXrayRefTable()
    # software eff of 80%
    efficiency = 0.8
    # logL binning range
    nbins = 500
    # range of histogram in logL
    logLrange = (5.0, 20.0)

  let
    # get the raw log likelihood histograms (seq of individual values). Takes into account
    # cuts on properties and chip region
    rawLogHists = mapIt(toSeq(values(xray_ref)), buildLogLHist(h5cdl_file, it, region))
    # given raw log histograms, create the correctly binned histograms from it
    logHists = mapIt(rawLogHists, histogram(it, nbins, logLrange))
    # get the cut value for a software efficiency of 80%
    cutVals = mapIt(logHists, determineCutValue(it, efficiency))
    # get the correct binning for the histograms
    bins = linspace(logLrange[0], logLrange[1], nbins + 1, endpoint = true)

  result = initTable[string, float]()
  for key, dset in pairs(xray_ref):
    # incl the correct values for the logL cut values
    result[dset] = bins[cutVals[key]]
  # some debugging output
  when not defined(release):
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
        ecc = h5f[(grp.name / "eccentricity").dset_str][float64]
        lengthDivRmsTrans = h5f[(grp.name / "lengthDivRmsTrans").dset_str][float64]
        fracRmsTrans = h5f[(grp.name / "fractionInTransverseRms").dset_str][float64]
        # energy to choose correct bin
        energies = h5f[(grp.name / "energyFromPixel").dset_str][float64]

      # create seq to store data logL data for this chip
      var logL_s = newSeq[float64](ecc.len)
      for i in 0 .. ecc.high:
        var logL = 0.0'f64
        # try simple logL calc
        let refset = toRefDset(energies[i] / 1000.0)
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
  ## writes all relevant cluster data of events corresponding to `passedInds` in
  ## the group given by `group` for chip `chipNumber` to the output file `h5f`

  # TODO: add copy of attributes from energyFromPixel dataset, which contains
  # the calibration factor!
  
  # read all float datasets, which we want to write to the output file
  var float_dset_names = @(getFloatDsetNames())
  # add the final two datasets, which we'd like to write
  float_dset_names.add "likelihood"
  float_dset_names.add "energyFromPixel"
  var float_data_tab = initTable[string, seq[float]]()
  let chpGrpName = group.name / &"chip_{chipNumber}"
  # get mutable group for this chip to copy attributes
  var chpGrpIn = h5f[chpGrpName.grp_str]
  # fill table of float data sets
  for dset in float_dset_names:
    float_data_tab[dset] = h5f[(chpGrpName / dset).dset_str][float64]

  # now get the event numbers to compare against the indices
  let evNumbers = h5f[(chpGrpName / "eventNumber").dset_str][int64]
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
  # copy attributes over from the input file
  runGrp.copy_attributes(group.attrs)
  chpGrpOut.copy_attributes(chpGrpIn.attrs)


proc filterClustersByLogL(h5f: var H5FileObj, h5fout: var H5FileObj, tracking = false) =
  ## filters all clusters with a likelihood value in the given `h5f` by 
  ## the logL cut values returned by `calcCutValueTab`
  ## clusters passing the cuts are stored in `h5fout`
  ## if `tracking == false` we cut on non tracking data
  const region = crGold

  let cutTab = calcCutValueTab(region)
  # get the likelihood and energy datasets
  # get the group from file
  for num, group in runs(h5f):
    echo &"Start logL cutting of run {group}"
    # get number of chips from attributes
    var mgrp = h5f[group.grp_str]
    var run_attrs = mgrp.attrs
    let nChips = run_attrs["numChips", int]
    # get timestamp for run
    let tstamp = h5f[(group / "timestamp").dset_str][int64]
    for chpGrp in items(h5f, group):
      if "fadc" in chpGrp.name:
        continue
      # iterate over all chips and perform logL calcs
      var attrs = chpGrp.attrs
      let
        # get chip specific dsets
        chipNumber = attrs["chipNumber", int]
        # get the datasets needed for LogL
        energy = h5f[(chpGrp.name / "energyFromPixel").dset_str][float64]
        logL = h5f[(chpGrp.name / "likelihood").dset_str][float64]
        centerX = h5f[(chpGrp.name / "centerX").dset_str][float64]
        centerY = h5f[(chpGrp.name / "centerY").dset_str][float64]        
        evNumbers = h5f[(chpGrp.name / "eventNumber").dset_str][int64].asType(int)
        # get indices (= event numbers) corresponding to no tracking
        tracking_inds = h5f.getTrackingEvents(mgrp, tracking = tracking)
        # get all events part of tracking
        tracking_events = filterTrackingEvents(evNumbers, tracking_inds)

      # hash set containing all indices of clusters, which pass the cuts
      var passedInds = initSet[int]()
      # iterate through all clusters not part of tracking and apply logL cut
      for ind in tracking_events:
        let dset = energy[ind].toRefDset
        # given datasest add element to dataset, iff it passes cut
        let regionCut = cutPosition(centerX[ind], centerY[ind], region)
        if logL[ind] <= cutTab[dset] and regionCut == true:
          # include this index to the set of indices
          passedInds.incl ind

      # create dataset to store it
      if passedInds.card > 0:
        # call function which handles writing the data
        h5f.writeLikelihoodData(h5fout, mgrp, chipNumber, cutTab, passedInds)
      else:
        var mchpGrp = chpGrp
        mchpGrp.attrs["LogLSpectrum"] = "No events passed cut"
        echo "No clusters found passing logL cut"

proc main() =
  # create command line arguments
  let args = docopt(doc)
  let
    h5f_file = $args["<HDF5file>"]
    ref_file = $args["--reference"]

  let tracking_flag = if $args["--tracking"] == "true": true else: false

  var h5foutfile: string = ""
  if $args["--h5out"] != "nil":
    h5foutfile = $args["--h5out"]

  var h5f = H5file(h5f_file, "rw")
  h5f.visitFile

  var h5fout: H5FileObj
  if h5foutfile != "":
    h5fout = H5file(h5foutfile, "rw")
  else:
    # in case no outfile given, we write to the same file
    h5fout = h5f
  
  # perform likelihood calculation
  h5f.calcLogLikelihood(ref_file)
  # now perform the cut on the logL values stored in `h5f` and write
  # the results to h5fout
  h5f.filterClustersByLogL(h5fout, tracking = tracking_flag)
  # given the cut values and the likelihood values for all events
  # based on the X-ray reference distributions, we can now cut away
  # all events not passing the cuts :)
  let err = h5f.close()
  if err != 0:
    echo &"Could not close h5 file properly! Return value was {err}"
  
  echo "Writing of all chips done"


when isMainModule:
  main()
