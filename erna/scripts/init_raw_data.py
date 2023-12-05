'''
initialize models
example:
    python3 erna/manage.py shell < erna/scripts/init_raw_data.py

- Delete all data in RawData, Sample, SampleFile, SampleProject
- Uncompress *gz files
- Reload raw data into RawData
'''
# import sys
# print(sys.path)

from pipelines.process.process_raw_data import ProcessRawData

client = ProcessRawData()
print('Search and uncompress gz files...')
client.uncompress_raw_data()

print('Refresh RawData...')
res = client.reset_sample()


