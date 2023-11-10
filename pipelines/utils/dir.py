
import os, sys
import re

class Dir:
    def __init__(self, indir):
        self.indir = indir

    def init_dir(self):
        '''
        create the directory if that doesn't exists
        '''
        pool = [self.indir,]
        while pool:
            curr_dir = pool.pop(0)
            if not os.path.isdir(curr_dir):
                parent_dir = os.path.dirname(curr_dir)
                if os.path.isdir(parent_dir):
                    try:
                        os.mkdir(curr_dir, 0o777)
                    except Exception as e:
                        print(e)
                        return False
                else:
                    pool = [parent_dir, curr_dir] + pool
        return True

    def format_dir(self):
        '''
        format directory
        '''
        #judge if absolute dir
        if self.indir.find('/')==0:
            pass
        elif self.indir.find('~')==0:
            self.indir = re.sub('~', os.getenv('HOME'), self.indir)
        elif indir is None:
            self.indir = os.getcwd()
        else: # ./test or test
            self.indir = os.path.abspath(self.indir)
        #alway followed by '/'
        if not self.indir[-1]=='/':
            self.indir += '/'
            
        #create directory if not exists
        if not os.path.isdir(self.indir):
            os.mkdir(self.indir, 0o777)
        return self.indir


    def recrusive_files(self): 
        '''
        list all files with a given directory and sub directories
        get all files
        '''
        all_files=[]
        for root, dirs, files in os.walk(self.indir):
            #print('########', root, dirs, len(all_files) )
            for filename in list(files):
                out_file = os.path.join(root, filename)
                if os.path.isfile(out_file) and out_file.find('/.') == -1:
                    all_files.append(out_file)
                    #print(len(all_files), out_file)
        return all_files

    def seek_fq(self):
        """
        get all fastq files given a certain directory
        """
        print('Retrieve all *.fastq files under', self.indir)
        raw_files = []
        #get all files
        all_files = self.recrusive_files() 
        #print(all_files)
        #find file with .fastq or .fq
        for af in all_files:
            m = re.search(r'fastq$|fq$|fastq.gz$|fq.gz$', af)
            if m:
                #print('raw data:',af)
                raw_files.append(af)
        #print(raw_files)
        R2_files=[x for x in raw_files if '_R2' in x]
        R1_files=list(set(raw_files)-set(R2_files))
        #print(R2_files)
        return R1_files, R2_files, raw_files

    @staticmethod
    def cascade_dir(parent_dir:str, id_str:str, num:int):
        '''
        restrict numbers of files or subdirectory for tuning up
        '''
        if len(id_str) < num:
            return parent_dir
        next_dir = os.path.join(parent_dir, id_str[:num])
        return Dir.cascade_dir(next_dir, id_str[num:], num)
        
