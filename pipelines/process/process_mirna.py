
from connector.connect_mirbase import ConnectMirbase

class ProcessMiRNA:

    def __init__(self):
        pass

    def load_mirbase(self, overwrite:bool=None):
        '''
        download data from miRBase
        load data into eRNA
        '''
        local_path, local_files = ConnectMirbase().download_mirbase(overwrite)
        print(local_path, local_files)
