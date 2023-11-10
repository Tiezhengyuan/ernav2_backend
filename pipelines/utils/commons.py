
import os
import xml.dom.minidom

class Commons:
    cascade_num = 2
    
    def __init__(self):
        # default directory
        self.dir_download = os.environ.get('DIR_DOWNLOAD', '')
        self.dir_cache = os.environ.get('DIR_CACHE', '')
        self.dir_map = os.path.join(self.dir_cache, 'map')
        self.dir_chunk_map = os.path.join(self.dir_cache, 'chunk_map')
        self.dir_bin = os.environ.get('DIR_BIN', '')

        # default file
        # format: {<class name>:{<method name>:<local path>}}
        self.json_cache = os.path.join(self.dir_cache, 'cache_local_path.json')
        self.json_download = os.path.join(self.dir_cache, 'download_local_path.json')
    
    def print_xml(self, xml_str:str):
        temp = xml.dom.minidom.parseString(xml_str)
        new_xml = temp.toprettyxml()
        print(new_xml)
    
    def print_dict(self, indict):
        '''
        print dictionary to stdout for debugging
        '''
        n = 1
        for key in sorted(indict.keys()):
            print('{:5}: {:10}\t{}'.format(n, key, indict[key]))
            n += 1
