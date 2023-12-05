"""
connect FTP
"""
import os
import re
from sh import gunzip
from ftplib import FTP
from pipelines.utils import Dir

class ConnectFTP:
    def __init__(self, url:str, username:str=None, password:str=None):
        self.url = url

    def is_dir(self, ftp, name=str):
        origin_dir = ftp.pwd()
        try:
            ftp.cwd(name)
            ftp.cwd(origin_dir)
            return True
        except Exception as e:
            # print(e)
            pass
        # print('origin dir: ', ftp.pwd())
        return False

    def list_files(self, endpoint:str=None, match:str=None):
        '''
        list ftp path of files
        not recursive
        '''
        ftp_files = []
        # connect FTP
        ftp = FTP(self.url)
        ftp.login()
        if endpoint:
            ftp.cwd(endpoint)
        for file_name in ftp.nlst():
            if not self.is_dir(ftp, file_name):
                if match is None or re.search(match, file_name):
                    ftp_files.append((endpoint, file_name))
        return ftp_files

    def download_file(
            self,
            endpoint:str,
            file_name:str,
            local_path:str,
            run_gunzip:bool=None,
            overwrite:bool=None
        ):
        '''
        download one file from FTP
        one download one connection avoiding timeout
        '''
        if run_gunzip is None: run_gunzip = True
        if overwrite is None: overwrite = False
        Dir(local_path).init_dir()
        local_file = os.path.join(local_path, file_name)
        unzip_file = local_file.replace('.gz', '')
        # doesn't download if file exists and overwrite is False
        if os.path.isfile(unzip_file) and overwrite is False:
            return unzip_file
        
        # connect FTP
        ftp = FTP(self.url)
        ftp.login()
        if endpoint:
            ftp.cwd(endpoint)
        ftp_file = f"{self.url}/{endpoint}/{file_name}"
        
        # download
        try:
            with open(local_file, 'wb') as f:
                ftp.retrbinary(f"RETR {file_name}", f.write)
                print(f"Download {ftp_file}")
            # unzip .gz file
            if run_gunzip and local_file.endswith('gz'):
                print(f"decompress {local_file} to {unzip_file}")
                gunzip(local_file, '-f')
                return unzip_file
            return local_file
        except Exception as e:
            print('Failure: download data from FTP, error=', e)
            os.remove(local_file)
        return None
    
    def download_files(
            self,
            endpoint:str=None,
            match:str=None,
            local_path:str=None,
            overwrite:bool=None
        ):
        '''
        download files from FTP path
        That isnot recursive
        '''
        # list all files
        ftp_files = self.list_files(endpoint, match)

        # download files
        local_files = []
        for current_endpoint, file_name in ftp_files:
            local_file = self.download_file(
                current_endpoint,
                file_name,
                local_path,
                overwrite
            )
            if local_file and local_file not in local_files:
                local_files.append(local_file)
        return local_files

    def download_tree(self, local_path:str, endpoint:str=None,\
            match:str=None):
        '''
        arg: local_name is determined by os.path.join()
        Download FTP directory recursively
        '''
        # initialize endpoint
        ftp = FTP(self.url)
        ftp.login()
        origin_endpoint = ftp.pwd()

        # scan ftp path
        local_files = []
        pool = [(endpoint, local_path)]
        while pool:
            ftp.cwd(origin_endpoint)
            _endpoint, _local_path = pool.pop(0)
            if _endpoint:
                ftp.cwd(_endpoint)
            print(ftp.pwd(), ftp.nlst())
            # check if name is directory
            for name in ftp.nlst():
                if self.is_dir(ftp, name):
                    sub_endpoint = f"{_endpoint}/{name}"
                    sub_local_path = os.path.join(_local_path, name)
                    Dir(sub_local_path).init_dir()
                    pool.append((sub_endpoint, sub_local_path))
            #download files
            local_files += self.download_files(
                None, match, _local_path)
        return local_files

