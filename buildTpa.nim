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
    quit("Building TimepixAnalysis failed!")

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
                locateTool: string) =
  echo "Checking for NLopt"
  let nloptCheck = findLib("libnlopt.so", locateTool)
  if nloptCheck:
    echo "[INFO] NLopt shared library found."
  elif allowClone:
    echo "[INFO] Clone and build NLopt"
    let path = clonePath.expandFilename()
    shellCheck:
      one:
        mkdir -p ($path)
        cd ($path)
        git clone "https://github.com/vindaar/nimnlopt"
        cd nimnlopt/c_header
        git clone "git://github.com/stevengj/nlopt"
        cd nlopt
        mkdir build
        cd build
        cmake ".."
        make
    stdout.styledWrite(fgYellow, "\n\nNOTE: Please call `sudo make install` in " &
      "$HOME/src/nimnlopt/c_header/nlopt/build! Alternatively, copy the " &
      "`libnlopt.so` file into an appropriate location scanned by `ld.so`.\n\n")
  else:
    stdout.styledWrite(fgRed, "\n\nNOTE: NLopt shared library is not found using " & $locateTool &
      " and cloning is forbidden. Please install it yourself or launch this program with " &
      "--allowClone=true")

proc checkMPFIT(allowClone: bool, clonePath: string,
                locateTool: string) =
  echo "Checking for mpfit"
  let mpfitCheck = findLib("libmpfit.so", locateTool)
  if mpfitCheck:
    echo "[INFO] MPFIT shared library found."
  elif allowClone:
    echo "[INFO] Clone and build MPFIT"
    shellCheck:
      one:
        mkdir -p "$HOME/src"
        cd "$HOME/src"
        git clone "https://github.com/vindaar/nim-mpfit"
        cd "nim-mpfit"
        cd "c_src"
        gcc "-c -O3 -Wall -Werror -fpic" mpfit.c mpfit.h
        gcc "-shared -o" libmpfit.so mpfit.o
    stdout.styledWrite("\n\nNOTE: Please copy `libmpfit.so` located in " &
      "$HOME/src/nim-mpfit/c_src into an appropriate location, " &
      "e.g. '/usr/local/lib' scanned by `ld.so`. \n\n")
  else:
    stdout.styledWrite(fgRed, "\n\nNOTE: MPFIT shared library is not found using " & $locateTool &
      " and cloning is forbidden. Please install it yourself or launch this program with " &
      "--allowClone=true")

proc generatePathsFile() =
  shellCheck:
    one:
      cd Analysis
      nimble setup


proc compileBinaries() =
  proc compile(bin: string, flags: seq[string]) =
    let f = @flags.mapIt("-d:" & it).join(" ")
    let cmd = &"nim c {f} {bin}"
    shellCheck:
      ($cmd)

  let bins = @[
    ("Analysis/ingrid/parse_raw_tpx3", @["danger", "blosc"]),
    ("Analysis/ingrid/raw_data_manipulation", @["danger", "blosc"]),
    ("Analysis/ingrid/reconstruction", @["danger"]),
    ("Analysis/ingrid/likelihood", @["danger"]),
    ("Analysis/ingrid/runAnalysisChain", @["release"]),
    ("Analysis/ingrid/fake_event_generator", @["danger"]),
    ("Analysis/createAllLikelihoodCombinations", @[""]),
    ("Plotting/karaPlot/plotData", @["danger"]),
    ("Plotting/plotBackgroundRate/plotBackgroundRate", @["release"]),
    ("Plotting/plotBackgroundClusters/plotBackgroundClusters", @["release"]),
  ]
  for (b, f) in bins:
    compile(b, f)

echo "Finished building TimepixAnalysis!"

proc main(locateTool = "locate",
          allowClone = true,
          clonePath  = "~/src",
          args = "") =
  let dir = getCurrentDir()
  doAssert dir.endsWith("TimepixAnalysis"), "`buildTpa` must be run from the TimepixAnalysis root directory!"

  # 1. Set up NLopt and MPFIT dependencies
  checkNLopt(allowClone, clonePath, locateTool)
  checkMPFIT(allowClone, clonePath, locateTool)

  # 2. Generate a `nimble.paths` file from the lockfile
  generatePathsFile()

  # 3. Compile all relevant binaries
  compileBinaries()

when isMainModule:
  import cligen
  dispatch main, help = {
    "locateTool" : "Program to use to detect installed shared libraries on the system.",
    "allowClone" : "If true will automatically clone a git repository and build shared library dependencies.",
    "clonePath"  : "Base path in which cloned directories will be installed.",
    "args"       : "An additional command line argument string passed to all programs being compiled."
    }
