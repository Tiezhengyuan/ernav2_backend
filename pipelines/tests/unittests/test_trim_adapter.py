import os
from unittest import TestCase
from ddt import ddt, data, unpack

from pipelines.tests import DIR_DATA, DIR_TEMP
from pipelines.process.trim_adapter import TrimAdapter

@ddt
class TestTrimAdapter(TestCase):

    @data(
        ['reads_1.fq','reads_1.fq',],
        ['reads_1.fq.gz','reads_1.fq',],
    )
    @unpack
    def test_3end(self, infile, outfile):
        params = {
            'input': os.path.join(DIR_DATA, infile),
            'output': os.path.join(DIR_TEMP, outfile),
            'adapter': 'ATGCTG',
            'trim_end': '3end',
            'min_match': 6,
            'max_err': 0,
        }
        res = TrimAdapter(params)()
        assert res.get('trim_percentage') > 0

