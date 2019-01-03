import os
import re
import sys
import argparse

#Regex to pulls routines from text file.
def gen_subroutines():
  """
  Regexes to pull subroutines and functions from .txt file and write to .f files.
  """
  subroutines_re = re.compile(r"./ ADD NAME=([A-Z1-9]+)\n\s+(SUBROUTINE .*?\s+END)", re.S)
  functions_re = re.compile(r"./ ADD NAME=([A-Z]+)\n\n(\s+FUNCTION .*?\s+END)", re.S)
  int_functions_re = re.compile(r"./ ADD NAME=([A-Z]+)\n\n\s+(INTEGER FUNCTION.*?\s+END)", re.S)
  logical_functions_re = re.compile(r"./ ADD NAME=([A-Z]+)\n\n\s+(LOGICAL FUNCTION.*?\s+END)", re.S)
  complex_functions_re = re.compile(r"./ ADD NAME=([A-Z]+)\n\n\s+(COMPLEX FUNCTION.*?\s+END)", re.S)
  f_extension = '.f'
  
  with open("reclib_subroutines.txt", 'r') as f:
    subroutines_file = f.read()
  
  #grab all subroutines
  subroutines = subroutines_re.findall(subroutines_file)
  
  #grab all the different function types (INTEGER, LOGICAL, COMPLEX, FUNCTION)
  functions = int_functions_re.findall(subroutines_file)
  functions += logical_functions_re.findall(subroutines_file)
  functions += functions_re.findall(subroutines_file)
  functions += complex_functions_re.findall(subroutines_file)
  
  for i, x in enumerate(subroutines):
    f_name = x[0].lower() + f_extension
    print i, f_name
    #write with correct initial spacing.
    with open(f_name,'w') as f:
      print >> f, 6*' '+x[1]
  
  for i, y in enumerate(functions):
    f_name = y[0].lower() + f_extension
    print i, f_name
    #write with correct initial spacing.
    with open(f_name,'w') as f:
      print >> f, 6*' '+y[1]
  print 'Found', len(subroutines), 'subroutines'
  print 'Found', len(functions), 'functions'

def gen_programs():
  f_extension = '.f'
  program_re = re.compile(r"./ ADD NAME=([A-Z1-9]+)\n\s+(RUNFILE .*?\s+)ENDPROG", re.S)
  with open('reclib_examples.txt','r') as f:
    program_file = f.read()
  programs = program_re.findall(program_file)  
  print len(programs)
  for i, program in enumerate(programs):
    f_name = program[0].lower() + f_extension
    with open(f_name,'w') as f:
      prog = program[1].split('\n')
      print >> f, 6*' '+'PROGRAM '+program[0] + '\n' + '\n'.join(prog[5:])

def build_lib():
  pass

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--gen_code', '-gc', help='Parse text file with subroutines and functions.', action="store_true")
  parser.add_argument('--gen_prog', '-gp', help='Parse text file with main programs.', action="store_true")
  parser.add_argument('--build_lib', '-b', help='Build shared library.', action="store_true")
  args = parser.parse_args()

  if args.gen_code:
    gen_subroutines()    

  if args.build_lib:
    build_lib()

  if args.gen_prog:
    gen_programs()

