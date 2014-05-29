#! /usr/bin/python

# hacky way to add in job retries for ER5
# we didn't want to alter release code so this
# scripts adds job retries

# also add +DailyIhopeJob = True to submit files

from tempfile import mkstemp
from shutil import move
from os import remove, close
from sys import argv

import glob

flist = glob.glob('*.sub')
for f in flist:
  #Create temp file
  fh, abs_path = mkstemp()
  new_file = open(abs_path,'w')
  old_file = open(f)

  for line in old_file:
    if line.startswith('queue'):
      line = '+DailyIhopeJob = True\nqueue 1'
      new_file.write(line)
    else:
      new_file.write(line)

  #close temp file
  new_file.close()
  close(fh)
  old_file.close()

  #Remove original file
  remove(f)

  #Move new file
  move(abs_path, f)

