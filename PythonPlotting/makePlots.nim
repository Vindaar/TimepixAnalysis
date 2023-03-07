import shell
import strutils, os, strformat, sequtils

const kinds = ["All", "Bronze", "Silver", "Gold"]
const
  cmd = "./Plotting/PyS_createBackgroundRate.py"
  options = "--ck_binning --fancy --dontShow --title \"\""

proc getFiles(k, startPath: string): seq[string] =
  for f in walkFiles(&"{startPath}/*{k}*.h5"):
    if "septem" notin f:
      result.add f
  echo "Paths ", result


func noneFile(f: seq[string]): string =
  result = f.filterIt("fadc" notin it and "scinti" notin it)[0]
func fadcFile(f: seq[string]): string =
  result = f.filterIt("fadc" in it and "scinti" notin it)[0]
func scintiFile(f: seq[string]): string =
  result = f.filterIt("fadc" notin it and "scinti" in it)[0]
func bothFile(f: seq[string]): string =
  result = f.filterIt("fadc" in it and "scinti" in it)[0]

proc makeVetoComparisons =

  let startPath = paramStr(1)

  template onlyFadc(files: seq[string]): untyped {.dirty.} =
    let fadc = files.fadcFile
    shell:
      `$cmd` `$none` "--files" `$fadc` `$options`
    shell:
      `$cmd` `$none` "--files" `$fadc` `$options` "--log"
  template onlyScinti(files: seq[string]): untyped {.dirty.} =
    let scinti = files.scintiFile
    shell:
      `$cmd` `$none` "--files" `$scinti` `$options`
    shell:
      `$cmd` `$none` "--files" `$scinti` `$options` "--log"
  template onlyBoth(files: seq[string]): untyped {.dirty.} =
    let both = files.bothFile()
    shell:
      `$cmd` `$none` "--files" `$both` `$options`
    shell:
      `$cmd` `$none` "--files" `$both` `$options` "--log"
  template all(files: seq[string]): untyped {.dirty.} =
    shell:
      `$cmd` `$none` "--files" `$fadc` `$scinti` `$both` `$options`
    shell:
      `$cmd` `$none` "--files" `$fadc` `$scinti` `$both` `$options` "--log"

  for k in kinds:
    echo "-----"
    let files = getFiles(k, startPath)

    let none = files.filterIt("fadc" notin it and "scinti" notin it)[0]
    onlyFadc(files)
    onlyScinti(files)
    onlyBoth(files)
    all(files)

proc makeCdlComparisons =
  template compareNone(f2014, f2019: seq[string]): untyped {.dirty.} =
    let none2014 = f2014.noneFile
    let none2019 = f2019.noneFile
    shell:
      `$cmd` `$none2014` "--files" `$none2019` `$options`
    shell:
      `$cmd` `$none2014` "--files" `$none2019` `$options` "--log"
  template compareFadc(f2014, f2019: seq[string]): untyped {.dirty.} =
    let fadc2014 = f2014.fadcFile
    let fadc2019 = f2019.fadcFile
    shell:
      `$cmd` `$fadc2014` "--files" `$fadc2019` `$options`
    shell:
      `$cmd` `$fadc2014` "--files" `$fadc2019` `$options` "--log"
  template compareScinti(f2014, f2019: seq[string]): untyped {.dirty.} =
    let scinti2014 = f2014.scintiFile
    let scinti2019 = f2019.scintiFile
    shell:
      `$cmd` `$scinti2014` "--files" `$scinti2019` `$options`
    shell:
      `$cmd` `$scinti2014` "--files" `$scinti2019` `$options` "--log"
  template compareBoth(f2014, f2019: seq[string]): untyped {.dirty.} =
    let both2014 = f2014.bothFile
    let both2019 = f2019.bothFile
    echo "Both 2014 ", both2014
    echo "Both 2019 ", both2019
    shell:
      `$cmd` `$both2014` "--files" `$both2019` `$options`
    shell:
      `$cmd` `$both2014` "--files" `$both2019` `$options` "--log"

  for k in kinds:
    echo "-----"
    let files2014 = getFiles(k.toLowerAscii, "likelihoods_cdl2014")
    let files2019 = getFiles(k, "likelihoods_cdl2019")

    compareNone(files2014, files2019)
    compareFadc(files2014, files2019)
    compareScinti(files2014, files2019)
    compareBoth(files2014, files2019)
    #all(files)

makeCdlComparisons()
