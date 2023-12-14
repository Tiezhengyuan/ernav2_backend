'''
GTF/GFF
'''
import json
import re
from typing import Iterable
import numpy as np
import pandas as pd
import anndata as ad
import os

from pipelines.utils.dir import Dir

class Annot:
    fields = ['seqname','source','feature','start','end','score','strand','frame',]

    def __init__(self, annot_file):
        self.annot_file = annot_file
        self.file_type = annot_file[-3:]
        self.outdir =  os.path.join(os.path.dirname(annot_file), \
            f'{self.file_type}_features')

    def __call__(self):
        # skip read gff if features have been decomposed
        if os.path.isdir(self.outdir):
            res = {}
            for file_name in os.listdir(self.outdir):
                if file_name.endswith('json'):
                    feature, _ = os.path.splitext(file_name)
                    res[feature] = os.path.join(self.outdir, file_name)
            return res
        else:
            Dir(self.outdir).init_dir()
            return self.decompose_by_feature()

    def decompose_by_feature(self):
        annot = {}
        for rec in self.parse_annot_file():
            feature = rec.get('feature') or 'unknown'
            if feature not in annot:
                annot[feature] = []
            annot[feature].append(rec)

        res = {}
        for feature, records in annot.items():
            outfile = os.path.join(self.outdir, f"{feature}.json")
            with open(outfile, 'w') as f:
                json.dump(records, f, indent=4)
            res[feature] = outfile
        return res
    
    def parse_annot_file(self) -> Iterable:
        with open(self.annot_file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    line = line.rstrip()
                    items = line.split('\t')
                    rec = dict([(k, v) for k, v in zip(self.fields, items[:-1])])

                    attr = self.gff_attribute(items[-1]) if self.file_type == 'gff' \
                        else self.gtf_attribute(items[-1])
                    rec.update(attr)
                    yield rec

    def gff_attribute(self, attribute:str):
        attr = {}
        for item in attribute.split(';'):
            name1, val1 = item.split('=')
            if ',' in val1:
                attr[name1] = val1.split(',')
            else:
                attr[name1] = val1
        return attr

    def gtf_attribute(self, attribute:str):
        attr = {}
        attribute = re.sub('\";$', '', attribute)
        for item in attribute.split("\"; "):
            name1, val1 = item.split(" \"")
            if name1 not in attr:
                attr[name1] = []
            attr[name1].append(val1)
        return attr

    def get_feature(self, feature:str):
        '''
        return dataframe object
        '''
        infile = os.path.join(self.outdir, f'{feature}.json')
        if os.path.isfile(infile):
            with open(infile, 'r') as f:
                data = pd.read_json(f)
                data.index = data['ID']
                data = data.transpose()
                return data
        print(f"{infile} does not exit.")
        return None

