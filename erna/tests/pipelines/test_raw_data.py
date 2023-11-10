from unittest import TestCase

from pipelines.process.process_raw_data import ProcessRawData
from pipelines.tests import DIR_RAW_DATA

class TestRawData(TestCase):

  def test_scan_raw_data(self):
    ProcessRawData(DIR_RAW_DATA).scan_raw_data()
