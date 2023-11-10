import os, sys
import re
import gzip
import json

class File:
    def __init__(self, infile):
        self.infile = infile

    def readonly_handle(self):
        '''
        ouput is read-only file handle
        '''
        if self.infile.endswith('.gz'):
            in_obj = gzip.open(self.infile, 'rt')
        else:
            in_obj = open(self.infile, 'rt')
        return in_obj

    def read_slice(self, rows:int=None, skip:int=None):
        '''
        read file and return some rows in a block
        '''
        if rows is None: rows = 100
        if skip is None: skip = 0
        in_obj = gzip.open(self.infile, 'rt') if self.infile.endswith('.gz') \
            else open(self.infile, 'rt')
        with in_obj as f:
            block = []
            # skip some top rows
            for i in range(skip):
                next(f)
            for line in f:
                line = line.rstrip()
                block.append(line)
                if len(block) == rows:
                    yield block
                    block = []
            if block:
                yield block            

    def read_top_lines(self, rows:int=None):
        '''
        read file and return some rows in a block
        '''
        if rows is None: rows = 1
        in_obj = gzip.open(self.infile, 'rt') if self.infile.endswith('.gz') \
            else open(self.infile, 'rt')
        with in_obj as f:
            lines = []
            for i in range(rows):
                line = next(f)
                lines.append(line.rstrip())
            return lines


    def list_to_file(self, inlist, out_file):
        '''
        export list to a text file seperated by return
        '''
        out_obj = open(out_file, 'wt')
        for key in inlist:
            out_obj.write(str(key)+'\n')
        out_obj.close()
        print('write a list to ', out_file)

    def file_to_list(self):
        '''
        read certain column in a file
        '''
        outlist=[]
        try: 
            in_obj = self.readonly_handle(self.infile)
            for line in in_obj:
                line = line.strip()
                outlist.append(line)
            in_obj.close()
        except FileNotFoundError:
            print(self.infile, 'didnot exit!')
            pass
        return outlist

    def file_to_dict(self, pattern=","):
        '''
        read text file into dictionary
        '''
        outdict={}
        with open(self.infile, 'r') as f:
            for line in f:
                line = line.rstrip()
                if not line.startswith('#'):
                    items=line.split(pattern)
                    outdict[items[0]]=items[1]
        return outdict

  
    def read_dump_file(self):
        '''
        *.dmp is exported from Oracle database
        '''
        with open(self.infile, 'r') as f:
            for line in f:
                items = re.split(r'\t*\|\t*', line.rstrip())
                yield items

    def unzip_gz(self):
        '''
        unzip .gz file using gunzip in linux
        '''
        pass