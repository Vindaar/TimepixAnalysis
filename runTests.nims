# executes all tests
import shell

let dir = getCurrentDir()

shell:
  nim c "-r --threads:on" `$dir`/Tests/io/tReadFadc.nim
  nim c "-r --threads:on" `$dir`/Tests/tMatchAnyRunFolder.nim
  nim c "-r --threads:on" `$dir`/Tests/tOldTosFolder.nim
