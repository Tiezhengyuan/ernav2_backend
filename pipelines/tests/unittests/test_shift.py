from unittest import TestCase
from ddt import ddt, unpack, data
from pipelines.match.shift import Shift

@ddt
class TestShift(TestCase):
  @data(
    ['cdd', 'acdd', [1,]],
    ['nina', 'ninjaninaaninann', [5, 10]],
    ['cdc', 'acddcd', []],
    ['cdc', 'cdc', [0]],
    ['cdc', 'acdcdc', [1,3]],
    ['NNN', 'NNNNN', [0,1,2]]
  )
  @unpack
  def test_shift_and(self, pattern, text, expect):
    res = Shift(pattern).shift_and(text)
    assert res == expect
    

