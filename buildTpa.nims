# Simple script to build all necessary tools of the TimepixAnalysis as well as
# the TPA tools itself.
import shell, strutils
# The script will make no attempt to install the HDF5 library (only the Nim lib)!

echo "Adding helpers to nimble"
let dir = getCurrentDir()
shell:
  one:
    cd `$dir`/NimUtil
    nimble develop "-y"

echo "Checking for NLopt"
var nloptCheck = ""
shellAssign:
  nloptCheck = nimble path nlopt
if "Error" in nloptCheck:
  echo "Clone and build NLopt"
  shell:
    one:
      cd "$HOME/src"
      git clone "https://github.com/vindaar/nimnlopt"
      cd nimnlopt/c_header
      git clone "git://github.com/stevengj/nlopt"
      cd nlopt
      mkdir build
      cd build
      cmake ".."
      make
      cd "../../../"
      nimble develop "-y"
  echo "\n\nNOTE: Please call `sudo make install` in " &
    "$HOME/src/nimnlopt/c_header/nlopt/build! \n\n"

echo "Checking for mnpfit"
var mpfitCheck = ""
shellAssign:
  mpfitCheck = nimble path mpfit
if "Error" in mpfitCheck:
  echo "Clone and build mpfit"
  shell:
    one:
      cd "$HOME/src"
      git clone "https://github.com/vindaar/nim-mpfit"
      cd "nim-mpfit"
      cd "c_src"
      gcc "-c -Wall -Werror -fpic" mpfit.c mpfit.h
      gcc "-shared -o" libmpfit.so mpfit.o
      cd ".."
      nimble develop "-y"
  echo "\n\nNOTE: Please copy `libmpfit.so` located in " &
    "$HOME/src/nim-mpfit/c_src into an appropriate location, " &
    "e.g. '/usr/local/lib'! \n\n"

echo "Adding ingridDatabase to nimble"
shell:
  one:
    cd `$dir`/InGridDatabase
    nimble develop "-y"
    nimble koch

echo "Adding ingrid to nimble"
shell:
  one:
    cd `$dir`/Analysis
    nimble develop "-y"

echo "Now add the InGrid-Python module via `setup.py develop` to your " &
  "correct Python installation"

echo "Compiling Nim procs for Python"
shell:
  one:
    cd `$dir`/Analysis/ingrid
    nim c "-d:release --app:lib --out:procsForPython.so" procsForPython.nim
    mv procsForPython.so "../../InGrid-Python/ingrid"
