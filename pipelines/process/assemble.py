
class Assemble:
  def __init__(self, params:dict):
    self.params = params
  
  def align_transcriptome(self):
    cmd = self.cmd_stringtie()
    print(cmd)


  def cmd_stringtie(self):
    cmd = [
      self.params['tool'].exe_path,
    ]
    return cmd