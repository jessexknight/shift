import os

# full path to parent of /code/utils/fio.py
therootpath = os.path.abspath(__file__).replace(os.path.join('code','utils','fio.py'),'')

def rootpath(*args):
  # e.g. rootpath('a','b') returns therootpath/a/b
  return os.path.join(therootpath,*args)

def genpath(fname):
  # create path to fname if needed
  path = os.path.dirname(fname)
  if path: os.makedirs(path,exist_ok=True)
  return fname
