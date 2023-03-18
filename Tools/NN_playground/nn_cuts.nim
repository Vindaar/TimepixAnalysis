#from std / algorithm import lowerBound
#import pkg / datamancer except linspace
#import seqmath
#
#import flambeau/[flambeau_nn, tensors]
## local
#import ./io_helpers, ./nn_types

## Include this file!! After `io_helpers` imported!

proc determineEff*(pred: seq[float], cutVal: float,
                  isBackground = true): float =
  ## returns the efficiency given the sorted (!) predictions, a
  ## cut value `cutVal` and whether it's background or signal
  let cutIdx = pred.lowerBound(cutVal)
  result = cutIdx.float / pred.len.float
  if isBackground:
    result = (1.0 - result)

proc calcSigEffBackRej*(df: DataFrame, bins: seq[float],
                       isBackground = true): DataFrame =
  ## returns the signal eff and backround rej for all bins of the
  ## given data frame, split by the `bins` column (that is CDL classes)
  let vals = df.arrange("predictions")["predictions", float]
  var effs = newSeqOfCap[float](bins.len)
  for l in bins:
    let eff = determineEff(vals.toSeq1D, l, isBackground = isBackground)
    effs.add eff
  result = toDf({ "eff" : effs,
                  "cutVals" : bins })

proc calcRocCurve*(predictions: seq[float], targets: seq[int]): DataFrame =
  # now use both to determine signal and background efficiencies
  # essentially have to generate some binning in `logL` we deem appropriate,
  let dfPlt = toDf(predictions, targets)
    .mutate(f{"isSignal" ~ `targets` == 1})
  const nBins = 1000
  let bins = linspace(predictions.min, predictions.max, nBins)
  let dfSignal = dfPlt.filter(f{`isSignal` == true})
  let dfBackground = dfPlt.filter(f{`isSignal` == false})

  var
    sigEffDf = newDataFrame()
    backRejDf = newDataFrame()
  if dfSignal.len > 0:
    sigEffDf = calcSigEffBackRej(dfSignal, bins, isBackground = false)
      .rename(f{"sigEff" <- "eff"})
  if dfBackground.len > 0:
    backRejDf = calcSigEffBackRej(dfBackground, bins, isBackground = true)
      .rename(f{"backRej" <- "eff"})
  if sigEffDf.len > 0 and backRejDf.len > 0:
    result = innerJoin(sigEffDf, backRejDf, by = "cutVals")
  elif sigEffDf.len > 0:
    result = sigEffDf
  elif backRejDf.len > 0:
    result = backRejDf
  else:
    doAssert false, "Both signal and background dataframes are empty!"

proc determineCutValue(model: AnyModel, device: Device, df: DataFrame, ε: float): float =
  let (cdlInput, cdlTarget) = df.toInputTensor()
  let cdlPredict = model.predict(cdlInput, device)

  # Compute alternative cut value by percentile
  # E.g. `ε = 0.8` means we need `1.0 - ε = 0.2` because the network predicts the
  # signal data on the positive side of the background (if we were to cut on the RHS
  # we'd keep all background!). Times 100 as int because `percentile` takes an integer
  let εLocal = 1.0 - ε
  result = cdlPredict.percentile((εLocal * 100.0).round.int)

  when false:
    ## XXX: this is not possible here anymore!
    # and now based on 'ROC curve'
    # Note: `predTargets` & `cdlTarget` contain same data. May not have same order though!
    let dfMLPRoc = calcRocCurve(cdlPredict, predTargets)
    dfMlpRoc.showBrowser()

    let dfme = dfMLPRoc.filter(f{`sigEff` >= εLocal}).head(20) # assumes sorted by sig eff (which it is)
    let cutVal = dfme["cutVals", float][0]
    dfMLPRoc.filter(f{`sigEff` >= ε}).showBrowser()
  echo "Cut value: ", result

proc determineCutValue*(model: AnyModel, device: Device, ε: float, readRaw: bool): float =
  ## Returns the required cut value for a desired software efficiency `ε`.
  ## This reads the (cleaned) CDL data, computes their output and returns the
  ## cut value at the percentile corresponding to `ε`.
  ##
  ## XXX: `readRaw` currently ignored!
  # get the CDL data
  let dfCdl = prepareCdl(readRaw)
  result = model.determineCutValue(device, dfCdl, ε)

import ingrid / ingrid_types
proc determineLocalCutValue*(model: AnyModel, device: Device, ε: float, readRaw: bool): OrderedTable[string, float] =
  ## Returns the required cut value for a desired software efficiency `ε`, which
  ## will be calculated per target/filter combination (instead of globally)
  ##
  ## This reads the (cleaned) CDL data, computes their output and returns the
  ## cut value at the percentile corresponding to `ε`.
  # get the CDL data
  result = initOrderedTable[string, float]()
  let dfCdl = prepareCdl(readRaw)
  for (tup, subDf) in groups(dfCdl.group_by("Target")):
    let target = tup[0][1].toStr
    result[target] = model.determineCutValue(device, subDf, ε)
  echo "Local cut values: ", result

proc calcCutValues*(model: AnyModel, device: Device, ε: float, readRaw: bool,
                    morphKind: MorphingKind): CutValueInterpolator =
  ## Important note: the `MorphingKind` is not actually doing any morphing in the NN
  ## case! (we only have a single distribution, equivalent to the likelihood distribution
  ## in the lnL cut methode case).
  ## As such all we can do is interpolate between the cut values from energy to energy.
  ## For the time being
  discard
