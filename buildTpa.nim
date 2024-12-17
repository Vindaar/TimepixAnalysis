# Simple script to build all necessary tools of the TimepixAnalysis as well as
# the TPA tools itself.
import shell
import std / [strutils, terminal, os, sequtils]
# The script will make no attempt to install the HDF5 library (only the Nim lib)!


var toContinue = true
template shellCheck(actions: untyped): untyped =
  if toContinue:
    let res = shellVerbose:
      actions
    toContinue = res[1] == 0
  if not toContinue:
    stdout.styledWrite(fgRed, "[ERROR]: Building TimepixAnalysis failed.\n")
    quit(1)

proc verifyTool(tool: string) =
  var res = ""
  shellAssign:
    res = which ($tool)
  if "not found" in res:
    stdout.styledWrite(fgRed, "[ERROR]: " & tool & " not found. Please install it!\n")
    quit(1)

proc findLib(lib: string, tool: string): bool =
  ## Checks if the given shared library `lib` is found
  ## using `locate`. Feel free to overwrite using `-d:LocatTool=foo`.
  var res = ""
  shellAssign:
    res = ($tool) ($lib)
  result = res.len > 0
  if result:
    doAssert lib in res, "`" & $tool & "` returned: " & $res & "\nbut did not find the library: " & $lib

proc checkNLopt(allowClone: bool, clonePath: string,
                locateTool: string): bool =
  ## Returns a bool indicating whether user action is required, i.e.
  ## setting up the path to the shared lib.
  echo "Checking for NLopt"
  let nloptCheck = findLib("libnlopt.so", locateTool)
  if nloptCheck:
    echo "[INFO] NLopt shared library found."
  elif allowClone:
    echo "[INFO] Clone and build NLopt"
    verifyTool("git")
    verifyTool("cmake")
    verifyTool("make")
    let path = clonePath.expandTilde()
    shellCheck:
      one:
        mkdir -p ($path)
        cd ($path)
        git clone "https://github.com/vindaar/nimnlopt"
        cd nimnlopt/c_header
        git clone "https://github.com/stevengj/nlopt"
        cd nlopt
        mkdir build
        cd build
        cmake "-DNLOPT_PYTHON=OFF -DNLOPT_GUILE=OFF -DNLOPT_SWIG=OFF -DNLOPT_MATLAB=OFF -DNLOPT_OCTAVE=OFF .."
        make
    stdout.styledWrite(fgYellow, "\n\nNOTE: Please call `sudo make install` in " &
      clonePath & "/nimnlopt/c_header/nlopt/build! Alternatively, copy the " &
      "`libnlopt.so` file into an appropriate location scanned by `ld.so`.\n\n")
    result = true # return true, because user action is required
  else:
    stdout.styledWrite(fgRed, "\n\nNOTE: NLopt shared library is not found using " & $locateTool &
      " and cloning is forbidden. Please install it yourself or launch this program with " &
      "--allowClone=true\n\n")

proc checkMPFIT(allowClone: bool, clonePath: string,
                locateTool: string): bool =
  ## Returns a bool indicating whether user action is required, i.e.
  ## setting up the path to the shared lib.
  echo "Checking for mpfit"
  let mpfitCheck = findLib("libmpfit.so", locateTool)
  if mpfitCheck:
    echo "[INFO] MPFIT shared library found."
  elif allowClone:
    echo "[INFO] Clone and build MPFIT"
    let path = clonePath.expandTilde()
    verifyTool("git")
    shellCheck:
      one:
        mkdir -p ($path)
        cd ($path)
        git clone "https://github.com/vindaar/nim-mpfit"
        cd "nim-mpfit"
        cd "c_src"
        gcc "-c -O3 -Wall -Werror -fpic" mpfit.c mpfit.h
        gcc "-shared -o" libmpfit.so mpfit.o
    stdout.styledWrite("\n\nNOTE: Please copy `libmpfit.so` located in " &
      clonePath & "/nim-mpfit/c_src into an appropriate location, " &
      "e.g. '/usr/local/lib' scanned by `ld.so`. \n\n")
    result = true # return true, because user action is required
  else:
    stdout.styledWrite(fgRed, "\n\nNOTE: MPFIT shared library is not found using " & $locateTool &
      " and cloning is forbidden. Please install it yourself or launch this program with " &
      "--allowClone=true\n\n")

proc generatePathsFile() =
  ## This generates a `nimble.paths` file from the `nimble.lock` file in
  ## `Analysis`. This file is crucial, because it sets up the Nim compiler
  ## such that it uses exactly the dependencies defined in the lockfile.
  ##
  ## This is done via the `--path:/foo/bar/` mechanism of the compiler.
  shellCheck:
    one:
      cd Analysis
      nimble setup

proc compileBinaries(compileLikelihood: bool) =
  proc compile(bin: string, flags: seq[string]) =
    let f = @flags.mapIt("-d:" & it).join(" ")
    let cmd = &"nim c {f} {bin}"
    shellCheck:
      ($cmd)

  var bins = @[
    # Base programs
    ("Analysis/ingrid/parse_raw_tpx3", @["danger", "blosc"]),
    ("Analysis/ingrid/raw_data_manipulation", @["danger", "blosc"]),
    ("Analysis/ingrid/reconstruction", @["danger"]),
    # Similarly the limit calculation can be compiled manually.
    # ("Analysis/ingrid/mcmc_limit_calculation", @["danger"]),
    ("Analysis/ingrid/fake_event_generator", @["danger"]),
    # Basic limit calcultion tool
    ("Analysis/ingrid/basic_limit_calc", @["danger"]),
    # Processing pipeline helpers
    ("Analysis/ingrid/runAnalysisChain", @["release"]),
    ("Analysis/createAllLikelihoodCombinations", @[]),
    # Plotting helpers
    ("Plotting/karaPlot/plotData", @["danger"]),
    ("Plotting/plotBackgroundRate/plotBackgroundRate", @["release"]),
    ("Plotting/plotBackgroundClusters/plotBackgroundClusters", @["release"]),
  ]

  let InCI = getEnv("IN_CI", "false").parseBool
  let ignoreCdl = "IgnoreCdlFile=" & $InCI
  if compileLikelihood:
    # NOTE: For the moment we don't build `likelihood`, because it requires some data files that are
    # too big to be hosted in the TPA repository. It is usually not required for most users anyway.
    bins.add ("Analysis/ingrid/likelihood", @["danger", "useMalloc", ignoreCdl])

  for (b, f) in bins:
    compile(b, f)

proc main(locateTool = "locate",
          allowClone = false,
          clonePath  = "~/src",
          compileLikelihood = false,
          args = "") =
  let dir = getCurrentDir()
  doAssert dir.endsWith("TimepixAnalysis"), "`buildTpa` must be run from the TimepixAnalysis root directory!"

  # 0. verify `locateTool` can be found
  verifyTool(locateTool)

  # 1. Set up NLopt and MPFIT dependencies
  var action = false
  action = action or checkNLopt(allowClone, clonePath, locateTool)
  action = action or checkMPFIT(allowClone, clonePath, locateTool)
  if action:
    stdout.styledWrite(fgYellow, "[WARNING]: Continuing with further setup after compilation of " &
      "NLopt and/or MPFIT. You will need to perform manual steps to make the shared libraries available " &
      "before you can run any of the compiled binaries.\n\n")

  # 2. Generate a `nimble.paths` file from the lockfile
  generatePathsFile()

  # 3. Compile all relevant binaries
  compileBinaries(compileLikelihood)

when isMainModule:
  import cligen
  dispatch main, help = {
    "locateTool" : "Program to use to detect installed shared libraries on the system.",
    "allowClone" : "If true will automatically clone a git repository and build shared library dependencies.",
    "clonePath"  : "Base path in which cloned directories will be installed.",
    "args"       : "An additional command line argument string passed to all programs being compiled."
    }
