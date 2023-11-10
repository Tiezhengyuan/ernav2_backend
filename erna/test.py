from ftplib import FTP

class Test:
  def download_file(self):
    ftp = FTP('ftp.ncbi.nlm.nih.gov')
    ftp.login()
    ftp.cwd('pubmed/')

class Test2(Test):
  def __init__(self):
    super(Test, self).__init__()

Test2().download_file()