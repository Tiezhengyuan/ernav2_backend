from unittest import TestCase
from ddt import ddt, unpack, data
from pipelines.match.kmer import Kmer

@ddt
class TestKmer(TestCase):
  @data(
    ['ATCGG', 2, None, ['AT', 'TC', 'CG', 'GG']],
    ['ATCGG', 3, None, ['ATC', 'TCC', 'CGG']],
    ['ATCGG', 6, None, ['ATCGG']],
    ['ATCGG', 3, 2, ['CGG']],
    ['ATCGG', 3, -4, ['TCG', 'CGG']],
  )
  @unpack
  def test_kmer(self, input, k, start, expect):
    iter = Kmer(input).kmer(k, start)
    res = [i for i in iter]
    print(res)
    # assert res == expect
    