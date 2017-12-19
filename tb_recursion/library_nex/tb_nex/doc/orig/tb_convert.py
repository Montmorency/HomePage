import os
import sys
import string

def swap_rst(fname='preface.rst'):
  with open(fname, 'r') as f:
    file_str = f.read()
  new_str = '. '.join(map(string.capitalize, [x.strip().replace('\n\n','\n') for x in file_str.split('.')]))

  cap_names = ['heine', 'nex', 'bullett', 'haydock', 'hamiltonian', 'slater', 'koster']
  for cp in cap_names:
    new_str = new_str.replace(cp, string.capitalize(cp))

#removes spurious white space around punctuation marks
  punc_spaces = [['( ', '('], [' )',')'],[' , ',', ']]
  for ps in punc_spaces:
    new_str = new_str.replace(ps[0], ps[1])
    

  with open('new_'+fname, 'w') as f:
    print >> f, new_str

if __name__ == '__main__':
  swap_rst(fname=sys.argv[1])
