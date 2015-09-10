import os
from subprocess import call

def copy_file_to_output_directory(output_directory, filename, source_timestamp):
  if os.path.exists(output_directory + "/" + filename):
    stat = os.stat(output_directory + "/" + filename)
    if int(source_timestamp) <= int(stat.st_mtime):
      return False  
  return True

def get_stuff_to_sync(input_directory, output_directory):
  
  files_to_sync = []
  directories_to_sync = []
  
  for f in os.listdir(input_directory):
    stat = os.stat(input_directory + "/" + f)
    if copy_file_to_output_directory(output_directory, f, stat.st_mtime):
      
      if os.path.isfile(input_directory + "/" + f):
        files_to_sync.append(input_directory + "/" + f)
      else:
        directories_to_sync.append(input_directory + "/" + f)

  return {'files': files_to_sync, 'directories': directories_to_sync}


def sync_files(input_directory, output_directory):
  things_to_sync = get_stuff_to_sync(input_directory, output_directory)

  files = things_to_sync['files']
  dirs = things_to_sync['directories']

  print files
  print
  print dirs

  for f in files:
    call(["cp", "-pr", f, output_directory + "/"])

  for d in dirs:
    call(["cp", "-pr", "-R", d, output_directory + "/"])



input_directory = "/home/aashish/root/macros/MODAnalyzer/plots/Version 3/"
output_directory = "/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Plots/Version 3/"

sync_files(input_directory, output_directory)