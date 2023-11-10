import os
from unittest import TestCase
from ddt import ddt, data, unpack

from pipelines.tests import DIR_DATA
from pipelines.biofile.fastq import FASTQ

@ddt
class TestFASTQ(TestCase):

    @data(
        ['reads_1.fq', 'r1']
    )
    @unpack
    def test_parse_records(self, input, expect):
        infile = os.path.join(DIR_DATA, input)
        fastq_iter = FASTQ(infile).parse_records()
        record = next(fastq_iter)
        assert record.id == expect
        print(record)
