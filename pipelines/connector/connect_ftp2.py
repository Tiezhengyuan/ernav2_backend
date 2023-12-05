"""
connect FTP using static methods
1. tremendous files have to be download
2. FTP login can't be alive long time.
3. internet connection is low-speed
"""
import os
import time
from typing import Callable
import ftplib
from django.conf import settings

from pipelines.utils.dir import Dir

class ConnectFTP2:
    ref_dir = getattr(settings, 'REFERENCES_DIR')
    dir_download =  ref_dir if os.path.isdir(ref_dir) \
        else os.environ.get('DIR_DOWNLOAD', '')

    @staticmethod
    def connect_ftp(ftp_endpoint:str, ftp_path:str, func:Callable=None,
                    try_times:str=None, **kwargs):
        '''
        connection handler
        '''
        error = ''
        if try_times is None: try_times = 3
        for i in range(1, try_times+1):
            try:
                res = None
                ftp = ftplib.FTP(ftp_endpoint)
                ftp.login()
                ftp.cwd(ftp_path)
                if func:
                    res = func(ftp, **kwargs)
                ftp.close()
                return { 'return': res }
            except ftplib.error_perm as err1:
                print(f"Failure: {ftp_endpoint} is ok. but the path {ftp_path} is wrong.")
                return {
                        'errno': 550,
                        'error': err1,
                    }
            except ftplib.all_errors as err2:
                if i == try_times:
                    print(f"Failure: can't target {ftp_endpoint}/{ftp_path}. " + 
                            f"Try {try_times} times. error={error}")
                    return {
                        'error': err2,
                    }
                time.sleep(2)
            except Exception as e:
                print(f"Error: error={e}")
                return {
                        'error': e,
                    }

    @staticmethod
    def _download(ftp, **kwargs):
        ftp_file_name = kwargs.get('ftp_file_name')
        local_file = kwargs.get('local_file')
        if ftp_file_name and local_file:
            try:
                with open(local_file, 'wb') as f:
                    ftp.retrbinary(f"RETR {ftp_file_name}", f.write)
                    print(f"Download FTP data {ftp_file_name} to {local_file}.")
                return True
            except Exception as e:
                print(f"Failure: can't download FTP data: {ftp_file_name}. error={e}")
                if os.path.isfile(local_file):
                    os.remove(local_file)
        return False
    
    @staticmethod
    def list_contents(ftp_endpoint:str, ftp_path:str):
        def _func(ftp):
            contents = []
            ftp.dir(contents.append)
            return contents
        res = ConnectFTP2.connect_ftp(ftp_endpoint, ftp_path, _func)
        if 'return' in res:
            return res['return']


    @staticmethod
    def parse_contents(ftp_endpoint:str, ftp_path:str, pattern:str=None):
        parser = {
            'dirs': [],
            'files': [],
            'current_path': ftp_path,
        }
        contents = ConnectFTP2.list_contents(ftp_endpoint, ftp_path)
        for i in contents:
            items = i.split(' ')
            if items[0].startswith('d'):
                parser['dirs'].append(f"{ftp_path}/{items[-1]}")
            elif items[0].startswith('-'):
                if pattern is None or (pattern is not None \
                        and items[-1].endswith(pattern)):
                    parser['files'].append(items[-1])
        return parser

    @staticmethod
    def scan_tree(ftp_endpoint:str, ftp_path:str, pattern:str=None)->list:
        '''
        anonymous user login FTP site
        walk through the directory recrusively
        '''
        path = []
        pool = [ftp_path,]
        while pool:
            current_ftp_path = pool.pop(0)
            list_contents = ConnectFTP2.parse_contents(
                ftp_endpoint,  current_ftp_path, pattern
            )
            pool += list_contents['dirs']
            for name in list_contents['files']:
                path.append((list_contents['current_path'], name))
        return path


    @staticmethod
    def retrieve_file_names(ftp_endpoint:str, ftp_path:str, name_pattern:str=None)->list:
        '''
        Suppose that large number of files exists in one directory.
        Don't consider sub-directory
        '''
        # scan ftp path
        res = ConnectFTP2.connect_ftp(
            ftp_endpoint,
            ftp_path,
            lambda ftp: ftp.nlst()
        )
        if 'return' in res:
            file_names = []
            for name in res['return']:
                if name_pattern is None or (name_pattern is not None\
                            and name.endswith(name_pattern)):
                    # print(ftp_path, name)
                    file_names.append(name)
            return file_names
        return []


    @staticmethod
    def download_file(ftp_endpoint:str, ftp_path:str, \
            ftp_file_name:str, local_path:str)->bool:
        '''
        anonymous user login FTP site
        Download one file at a time. Try 3 times in default
        local file would be covered if that exists.
        '''
        local_file = os.path.join(local_path, ftp_file_name)
        res = ConnectFTP2.connect_ftp(
            ftp_endpoint=ftp_endpoint,
            ftp_path=ftp_path,
            func=ConnectFTP2._download,
            ftp_file_name=ftp_file_name,
            local_file=local_file,
        )
        if 'return' in res:
            return res['return']
        return False
    

    @staticmethod
    def download_tree(ftp_endpoint:str, ftp_path:str, pattern:str,
                    local_path:str)->bool:
        '''
        anonymous user login FTP site
        Download a large number of files, but Download one file at a time. 
        local file would be covered if that exists.
        '''
        tree_path = ConnectFTP2.scan_tree(ftp_endpoint, ftp_path, pattern)
        for current_ftp_path, file_name in tree_path:
            if current_ftp_path.startswith('/'):
                current_ftp_path = current_ftp_path[1:]
            dir_names = [local_path,] + current_ftp_path.split('/')
            current_local_path = os.path.join(*dir_names)
            Dir(current_local_path).init_dir()
            # print(current_ftp_path, current_local_path, file_name)
            # download file
            ConnectFTP2.download_file(
                ftp_endpoint,
                current_ftp_path,
                file_name,
                current_local_path
            )


