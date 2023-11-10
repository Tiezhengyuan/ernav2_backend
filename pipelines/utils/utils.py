"""

"""
from copy import deepcopy
import gzip
import os
import re
import pandas as pd 


class Utils:

    @staticmethod
    def sort_array(arr):
        '''
        element: characters+numerics
        for example: sort chromosome name
        '''
        arr_len = len(arr)
        for i in range(0, arr_len-1):
            for j in range(1, arr_len):
                a=arr[j-1][3:]
                a_char= re.findall(r"[^0-9]", a)
                b=arr[j][3:]
                b_char= re.findall(r"[^0-9]", b)
                #a is char
                if a_char:
                    if b_char==[] or (b_char and a>b):
                        arr[j], arr[j-1] = arr[j-1], arr[j]
                #a_char==[]
                else:
                    if b_char==[] and int(a)>int(b):
                        arr[j], arr[j-1] = arr[j-1], arr[j]
        return arr
    
    @staticmethod
    def init_dict(input:dict, keys:list, default_val=None):
        '''
        arg: default_val = '', [], {}
        '''
        curr = input
        if isinstance(input, dict):
            for k in keys[:-1]:
                if k not in curr:
                    curr[k] = {}
                curr = curr[k]
            if keys[-1] not in curr:
                curr[keys[-1]] = default_val if \
                    default_val is not None else ''


    @staticmethod
    def update_dict(input:dict, key, val):
        if key not in ('', '-', None):
            if key not in input:
                input[key] = []
            else:
                if not isinstance(input[key], list):
                    input[key] = [input[key],]
            tmp = val if isinstance(val, list) else [val,]
            for t in tmp:
                if t not in input[key]:
                    input[key].append(t)
   
    @staticmethod
    def merge_dict(d1:dict, d2:dict)->dict:
        '''
        values of corresponding keys between d1 and d2
        should match to each other
        '''
        merged = deepcopy(d1)
        for k in d1:
            if k in d2:
                Utils.update_dict(merged, k, d2[k])
                del d2[k]
        if d2:
            merged.update(d2)
        return merged

    @staticmethod
    def get_deep_value(input:dict, keys:list):
        if not keys:
            return []
        val = []
        pool = [(keys, input),]
        while pool:
            curr_keys, curr_input = pool.pop(0)
            # print(curr_keys, curr_input)
            if isinstance(curr_input, dict):
                key = curr_keys[0]
                if key in curr_input:
                    if len(curr_keys) == 1:
                        tmp = []
                        if isinstance(curr_input[key], list):
                            tmp += curr_input[key]
                        else:
                            tmp = [curr_input[key]]
                        for t in tmp:
                            if t not in val and t not \
                                in (None, '', [], {}, '-'):
                                val.append(t)
                    else:
                        pool.append((curr_keys[1:], curr_input[key]))
            elif isinstance(curr_input, list):
                for item in curr_input:
                    pool.append((curr_keys, item))
        return val
    
    @staticmethod
    def switch_key_value(input:dict)->dict:
        '''
        switch key ~ value of a dictionary
        '''
        if not isinstance(input, dict):
            return None
        new = {}
        for key,val in input.items():
            for v in val if isinstance(val, list) else [val,]:
                if type(v) in (str, int, tuple):
                    Utils.update_dict(new, v, key)
        return new

    @staticmethod
    def search_series(series:pd.Series, index)->list:
        '''
        Note: index name of series must be unique
        '''
        val = series.get(index)
        if val is not None:
            series[index], series.iloc[-1] = series.iloc[-1], series[index]
            # print('##', series)
            return list(val) if type(val) == pd.Series else [val,] 
        return []

    @staticmethod
    def parse_ncbi_acc(infile)->dict:
        '''
        value: NCBI_protein_accession, index: UniProtKB_protein_accession
        source file: *_gene_refseq_uniprotkb_collab.gz
        '''
        accessions = {}
        # infile = os.path.join(self.dir_source, "gene_refseq_uniprotkb_collab.gz")
        df = pd.read_csv(infile, sep='\t', header=0)
        index_name = 'UniProtKB_protein_accession'
        df['acc_group'] = df[index_name].str[:2]
        for name, sub_df in df.groupby(['acc_group']):
            series = None
            if len(sub_df) > 1:
                series = pd.Series(sub_df.iloc['#NCBI_protein_accession'].squeeze())
                series.index = sub_df[index_name]
            elif len(sub_df) == 1:
                series = pd.Series(sub_df.iat[0,0], index=[sub_df.iat[0,1],])
            if series is not None:
                series.index.name = index_name
                series.name = name
                accessions[name] = series
            print(series.shape, series)
            print('\n\n')
        return accessions