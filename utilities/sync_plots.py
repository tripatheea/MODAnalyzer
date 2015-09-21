import os
from subprocess import call

def copy_file_to_output_directory(output_directory, filename, source_timestamp):
  if os.path.exists(output_directory + "/" + filename):
    stat = os.stat(output_directory + "/" + filename)
    if int(source_timestamp) <= int(stat.st_mtime):
      return False  
  return True

def sync_files(input_directory, output_directory):

  files_to_sync = []
  directories_to_sync = []

  if not os.path.isdir(output_directory):
    call(["mkdir", output_directory])
  
  for f in os.listdir(input_directory):
    stat = os.stat(input_directory + "/" + f)
    
      
    if os.path.isfile(input_directory + "/" + f):
      if copy_file_to_output_directory(output_directory, f, stat.st_mtime):
        files_to_sync.append(input_directory + "/" + f)
    else:
      directories_to_sync.append((input_directory + "/" + f, output_directory + "/" + f))

  for f in files_to_sync:
    call(["cp", "-pr", f, output_directory + "/"])

  for d in directories_to_sync:
    input_dir = d[0]
    output_dir = d[1]
    sync_files(input_dir, output_dir)
    

input_directory = "/home/aashish/root/macros/MODAnalyzer/plots/"
output_directory = "/home/aashish/Dropbox (MIT)/Research/CMSOpenData/Plots/"

sync_files(input_directory, output_directory)