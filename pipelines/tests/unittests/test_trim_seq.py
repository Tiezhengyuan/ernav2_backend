from unittest import TestCase
from ddt import ddt, unpack, data
from pipelines.sequence.trim_seq import TrimSeq

@ddt
class TestTrimSeq(TestCase):
  @data(
    #trim: exact 3'-end match
    ['ATTTGCGGGACTACGTTACAAG', 'ACGTTACAAG', 6, 0, 'ATTTGCGGGACT'],
    #trim: 3-end longer than adapter, and exact match
    ['ATTTGCGGGACTACGTTACAAGTTT', 'ACGTTACAAG', 6, 0, 'ATTTGCGGGACT'],
    #trim: 3-end longer than adapter, with 1 mismatch
    ['ATTTGCGGGACTACGTTACAAGTTT', 'ACGTTACAAT', 6, 1, 'ATTTGCGGGACT'],
    # no trimming: not exact match
    ['ATTTGCGGGACTACGTTACACG', 'ACGTTACAAG', 6, 0, 'ATTTGCGGGACTACGTTACACG'],
    # no trimming: seed is not exact match
    ['ATTTGCGGGACTACGTTACACG', 'TCGTTACACG', 6, 1, 'ATTTGCGGGACTACGTTACACG'],
  )
  @unpack
  def test_trim_3end(self, read, adapter, min_match, max_err, expect):
    t = TrimSeq(adapter, min_match, None, max_err)
    res, pos = t.trim_3end(read)
    assert res == expect
    