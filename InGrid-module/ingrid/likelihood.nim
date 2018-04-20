import docopt
import tables
import strutils, strformat, ospaths
import nimhdf5
import tos_helper_functions
import sequtils
import seqmath
import loopfusion

let doc = """
InGrid likelihood calculator. This program is run after reconstruction is finished.
It calculates the likelihood values for each reconstructed cluster of a run and
writes them back to the H5 file

Usage:
  likelihood <HDF5file> --reference <ref_file>

"""

type histTuple = tuple[bins: seq[float64], hist: seq[float64]]

proc splitSeq[T, U](s: seq[seq[T]], dtype: typedesc[U]): (seq[U], seq[U]) =
  ## splits a (N, 2) nested seq into two seqs
  result[0] = newSeq[dtype](s.len)
  result[1] = newSeq[dtype](s.len)
  for i in 0..s.high:
    result[0][i] = s[i][0].U
    result[1][i] = s[i][1].U

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

  echo ecc_ref[xray_ref[0]]
  
  # get the group from file
  for num, group in runs(h5f):
    echo &"Start logL calc of run {group}"
    # get number of chips from attributes
    var run_attrs = h5f[group.grp_str].attrs
    let nChips = run_attrs["numChips", int]
    # get timestamp for run
    let tstamp = h5f[(group / "timestamp").dset_str][int64]

    var logL_chips = newSeq[seq[float64]](nChips)
    
    for grp in items(h5f, group):
      # iterate over all chips and perform logL calcs
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
                              ecc[i] + 1.0,
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
    echo logL_dsets
    # forEach dset in var logL_dsets, logL in var logL_chips:
    #   dset[dset.all] = logL
    for tup in zip(logL_dsets, logL_chips):
      var (dset, logL) = tup
      dset[dset.all] = logL
    
                                                                 
                                                                  
      

proc main() =
  # create command line arguments
  let args = docopt(doc)
  let
    h5f_file = $args["<HDF5file>"]
    ref_file = $args["<ref_file>"]

  var h5f = H5file(h5f_file, "rw")
  h5f.visitFile
  # perform likelihood calculation
  h5f.calcLogLikelihood(ref_file)
  echo "Writing of all chips done"
  

when isMainModule:
  main()
