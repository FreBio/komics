'''
TO DO:
  * include version check for trimmomatic
  * include java check
'''

import os
import re
import sys
import subprocess


class Error (Exception): pass


class Trimmomatic:
  def __init__(self, 
    out,
    reads1,
    reads2,
    threads=1,
    dir=None,
    maxlength=125,
    minlength=100,
    minqual=30,
    windowbp=10,
    ):
      self.out = out
      self.reads1 = os.path.abspath(reads1)
      self.reads2 = os.path.abspath(reads2)
      self.threads = threads
      self.jarfile = self._build_jar(dir)
      self.adapters = self._build_adapters(dir)
      self.maxlength = maxlength
      self.minlength = minlength
      self.minqual = minqual
      self.windowbp = windowbp
      
      if not os.path.exists(self.reads1):
        sys.stderr.write('\nERROR: reads1 file not found: "' + self.reads1 + '"\n')
        sys.exit(0)
      if not os.path.exists(self.reads2):
        sys.stderr.write('\nERROR: reads2 file not found: "' + self.reads2 + '"\n')
        sys.exit(0)
      if not os.path.exists(self.jarfile):
        sys.stderr.write('\nERROR: TRIMMOMATIC jarfile not found: "' + self.jarfile + '"\n')
        sys.exit(0)
      if not os.path.exists(self.adapters):
        sys.stderr.write('\nERROR: TRIMMOMATIC adapters not found: "' + self.adapters + '"\n')
        sys.exit(0)


  def _build_jar(self, dir):
    if dir is None:
      try:
        path_trim=os.environ['DIR_TRIMMOMATIC']
        path_exe=re.sub("/$","",re.sub(".*/Trimmo", "trimmo", path_trim))
        #directory = os.path.abspath(os.environ['DIR_TRIMMOMATIC'] + '/trimmomatic-0.36.jar')
        directory = os.path.abspath(os.environ['DIR_TRIMMOMATIC'] + '/' + path_exe + '.jar')
      except KeyError:
        sys.stderr.write('\nPlease set the environment variable "DIR_TRIMMOMATIC" or set the directory name with the option "--dir"\n')
        sys.exit(0)   
    else:
      directory = os.path.abspath(dir + '/' + re.sub("/$","",re.sub(".*/Trimmo", "trimmo", dir)) + '.jar')
    return directory


  def _build_adapters(self, dir):
    if dir is None:
      try:
       os.environ['DIR_TRIMMOMATIC']
       adapters = os.path.abspath(os.environ['DIR_TRIMMOMATIC'] + '/adapters/TruSeq3-PE.fa')
      except KeyError:
        sys.stderr.write('\nPlease set the environment variable "DIR_TRIMMOMATIC""\n')
        sys.exit(0)   
    else:
      adapters = os.path.abspath(dir + '/adapters/TruSeq3-PE.fa')
    return adapters


  def trim_reads(self):
    sys.stderr.write('\nTrimming reads\n')
    sys.stderr.write('==============\n')

    trim_command = [
    "java -jar", self.jarfile,
    "PE",
    "-threads", str(self.threads),
    self.reads1,
    self.reads2,
    "-baseout", self.out + "_trimmed.fq.gz",
    "LEADING:" + str(self.minqual),
    "TRAILING:" + str(self.minqual),
    "ILLUMINACLIP:" + self.adapters + ":" + "2" + ":" + str(self.minqual) + ":" + "15" + ":" + "1" + ":true",
    "SLIDINGWINDOW:" + str(self.windowbp) + ":" + str(self.minqual),
    "MINLEN:" + str(self.minlength),
    "CROP:" + str(self.maxlength)
    ]
    
    subprocess.call(' '.join(trim_command), shell = True)
    os.remove(self.out + "_trimmed_1U.fq.gz")
    os.remove(self.out + "_trimmed_2U.fq.gz")
    sys.stderr.write('\nkomics trimfq successfully completed.\n\n')
