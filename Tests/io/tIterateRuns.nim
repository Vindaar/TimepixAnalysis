import unittest, sets
import ingrid / tos_helpers
import nimhdf5

#[
In the past (and as a regression now) there was a bug which caused us to
skip groups when iterating over all runs in an input H5 file.
In the past that was because we modified the group table while iterating
over it (we added groups).
That was fixed (commit: ?), but it seems to be back...
Let's write a test, to avoid such a regression again in the future.
]#

suite "Run iteration test":
  test "Iterate all runs in file":
    # first create a H5 file with a bunch of groups
    var h5f = H5file("test.h5", "rw")
    var groupNames = initHashSet[string]()
    for i in 0 .. 500:
      let name = recoBase() & $i
      groupNames.incl name
      discard h5f.create_group(name)

    var err = h5f.close()
    check err == 0

    ## now reopen the file and first check if all groups exist
    h5f = H5file("test.h5", "r")
    # first check if all groups exist
    for grp in groupNames:
      check grp in h5f

    err = h5f.close()
    check err == 0

    # finally reopen for the last time and iterate groups again.
    # Create second set and compare sets
    var readGroups = initHashSet[string]()
    h5f = H5file("test.h5", "r")

    for (num, grp) in runs(h5f):
      readGroups.incl grp

    check readGroups.card == 501
    check readGroups == groupNames

    err = h5f.close()
    check err == 0

    removeFile("test.h5")
