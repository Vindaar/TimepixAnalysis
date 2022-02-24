import seqmath ## `zipEm` has since been moved to `seqmath`
import unittest

suite "Helper tests":
  test "zipEm macro":
    let a = @[1,2,3]
    let b = @[1'f64, 1.5, 3.0]
    let c = @[0'u16, 4, 8]
    let f = zipEm(a, b, c)
    let exp = @[(1, 1'f64, 0'u16), (2, 1.5, 4'u16), (3, 3.0, 8'u16)]
    check f == exp
