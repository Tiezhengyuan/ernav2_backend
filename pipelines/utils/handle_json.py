"""
Note: file size of JSON might be <200MB
"""
from typing import Iterable
import json
import os
import re
import pandas as pd
from .commons import Commons
from .utils import Utils

class HandleJson(Commons):
    def __init__(self, infile:str=None):
        super(HandleJson, self).__init__()
        self.infile = infile

    def from_text(self):
        outfile = re.sub('txt$', 'json', self.infile)
        df = pd.read_csv(self.infile, sep="\t", skiprows=0, header=1)
        df.rename(columns={list(df)[0]: 'assembly_accession'}, inplace=True)
        df.to_json(outfile, orient='index', indent=2)
        return outfile
    
    def to_dict(self)->dict:
        try:
            with open(self.infile, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(e)
            return {}

    def read_json(self)->Iterable:
        try:
            with open(self.infile, 'rb') as f:
                data = json.load(f)
                for k,v in data.items():
                    yield (k, v)
        except Exception as e:
            print('###error: ', e)
            yield

    def search_value(self, keys:list)->list:
        if not keys:  return []
        try:
            for k,v in self.read_json():
                if k == keys[0]:
                    return Utils.get_deep_value({k:v}, keys)
        except Exception as e:
            pass


    def update_json(self, input_dict:dict)->None:
        new = {}
        if os.path.isfile(self.infile):
            try:
                with open(self.infile, 'r') as f:
                    origin = json.load(f)
                    new.update(origin)
            except Exception as e:
                print(e)
        print(new)
        # input_dict override origin if there are overlapped
        new.update(input_dict)
        print(new)
        return self.save_json(new)

    def save_json(self, input:dict):
        try:
            with open(self.infile, 'w') as f:
                json.dump(input, f, indent=4, sort_keys=True)
                print(f"{self.infile} is saved.")
        except Exception as e:
            print(e)
            return False
        return True


