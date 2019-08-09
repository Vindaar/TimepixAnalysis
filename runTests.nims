# executes all tests
import shell

shell:
  nim c "-r --threads:on" Tests/io/tReadFadc.nim
  nim c "-r --threads:on" Tests/tMatchAnyRunFolder.nim
  nim c "-r --threads:on" Tests/tOldTosFolder.nim
