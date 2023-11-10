from unittest import TestCase
from ddt import ddt, unpack, data
from pipelines.sequence.sequence import Sequence

@ddt
class TestSequence(TestCase):
  @data(
    ['ATCG', {'A': '0001', 'T': '0010', 'C': '0100', 'G': '1000'} ],
    ['nina', {'n': '0101', 'i': '0010', 'a': '1000'}],
  )
  @unpack
  def test_to_mask(self, input, expect):
    _, res = Sequence(input).nt_bitvector()
    # print(res)
    assert res == expect
    

