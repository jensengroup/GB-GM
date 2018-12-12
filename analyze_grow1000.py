'''
Written by Jan H. Jensen 2018
'''

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
#import pickle
import operator
import collections
import re

def read_file(file_name):
  smiles_list = []
  with open(file_name,'r') as file:
    for smiles in file:
      smiles_list.append(smiles)

  return smiles_list

def get_probs(smarts_list,smiles_list,ring=False):
  bonds = []
  probs = collections.OrderedDict()

  for smarts in smarts_list:
    probs[smarts] = 0

  number_of_molecules = 0
  tot = 0
  for smiles in smiles_list:
    #print smiles
    number_of_molecules += 1
    mol = Chem.MolFromSmiles(smiles)
    Chem.Kekulize(mol)
    for smarts in smarts_list:
      matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts),uniquify=ring)
      num_bonds = len(matches) 
      probs[smarts] += num_bonds
      tot += num_bonds
   
  tot = 0
  probs2 = collections.OrderedDict()
  for key in probs:
    if probs[key] > 0:
      #print key, probs[key]
      tot += probs[key]
      probs2[key] = probs[key]
            
  return tot, probs2

def clean_probs(probs):
  exceptions = ['[#7]#','[#8]=','[#9]','[#17]','[#35]','[#53]']
  probs2 = collections.OrderedDict()
  # for key in probs:
    # skip = False
    # for exception in exceptions:
    #   if exception in key:
    #     tokens = re.split('\[|\]|;',key)
    #     alt_key = '['+tokens[3]+']'+tokens[2]+'['+tokens[1]+';!R]'
    #     probs[alt_key] += probs[key]

  for key in probs:
    skip = False
    for exception in exceptions:
      if exception in key:
        skip = True
    if not skip:
      probs2[key] = probs[key]
  
  tot = 0
  for key in probs2:
      tot += probs2[key]
      
  return tot, probs2

  return probs

def get_p(probs):
  p = []
  for key in probs:
    p.append(float(probs[key])/tot)
    
  return p

def get_rxn_smarts_make_rings(probs):
  X = {'[#6R': 'X4', '[#7R': 'X3'}
  rxn_smarts = []
  for key in probs:
    tokens = key.split(']')

    smarts = ''
    if '=' in key:
      smarts += tokens[0][:-1] + X[tokens[0]] + ';!R:1]'
    else:
      smarts += tokens[0][:-1] + ';!R:1]=,'
    smarts += tokens[2][:-1] + ';!R:2]>>'
    smarts += '[*:1]1' + tokens[1] + '][*:2]1'

    rxn_smarts.append(smarts)
    
  return rxn_smarts

def get_rxn_smarts_rings(probs):
  X = {'[#6R': 'X4', '[#7R': 'X3'}
  rxn_smarts = []
  for key in probs:
    tokens = key.split(']')

    smarts = ''
    if '=' in key:
      smarts += tokens[0] + X[tokens[0]] + ';!r6;!r7;!R2:1]'
    else:
      smarts += tokens[0] + ';!r6;!r7;!R2:1]'

    smarts += tokens[2] + ';!r6;!r7:2]>>'
    smarts += '[*:1]' + tokens[1] + '][*:2]'

    rxn_smarts.append(smarts)
    
  return rxn_smarts

def get_rxn_smarts(probs):
  rxn_smarts = []
  for key in probs:
    smarts = ''
    tokens = key.split(']')
    smarts = tokens[0]
    if '-' in key and '#16' not in smarts:  # key <-> smarts
      smarts += ';!H0:1]>>[*:1]'
    if '=' in key and '#16' not in smarts:  # key <-> smarts
      smarts += ';!H1;!H0:1]>>[*:1]'
    if ']#[' in key:
      smarts += ';H3:1]>>[*:1]'
    if '#16' in smarts:  # key <-> smarts
      smarts += ':1]>>[*:1]'
      
    smarts += tokens[-2] + ']'
    rxn_smarts.append(smarts)
    
  return rxn_smarts

def get_mean_size(smiles_list):
  size = []
  for smiles in smiles_list:
    mol = Chem.MolFromSmiles(smiles)
    num_atoms = mol.GetNumAtoms()
    size.append(num_atoms)
 
  return np.mean(size), np.std(size)

def count_macro_cycles(smiles_list,smarts_list,tot,probs):  
  #probs = collections.OrderedDict()
  for smarts in smarts_list:
    probs[smarts] = 0
    
  for smiles in smiles_list:
    for smarts in smarts_list:
      mol = Chem.MolFromSmiles(smiles)
      Chem.Kekulize(mol)
      matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts),uniquify=True)
      if len(matches) > 0:
        probs[smarts] += 1
        tot += 1
      
  return tot,probs


#################################


file_name = 'grow1000.smi'

elements = ['#5','#6','#7','#8','#9','#14','#15','#16','#17','#35','#53']
bonds = ['-','=','#']

smiles_list = read_file(file_name)


mean_size, size_stdv = get_mean_size(smiles_list)
print('mean number of non-H atoms',mean_size,'+/-', size_stdv)
print('')

smarts = ['[*]','[R]','[!R]','[R2]']

tot,probs = get_probs(smarts,smiles_list)
print('Probability of ring atoms',float(probs['[R]'])/probs['[*]'])
print('Probability of non-ring atoms',float(probs['[!R]'])/probs['[*]'])
print('Probability of fused-ring atoms',float(probs['[R2]'])/probs['[*]'])
print('')

smarts = ['[R]~[R]~[R]','[R]-[R]-[R]','[R]=[R]-[R]']

tot,probs = get_probs(smarts,smiles_list,ring=True)

#print(tot,probs)

print('Probability of [R]-[R]-[R]',float(probs['[R]-[R]-[R]'])/probs['[R]~[R]~[R]'])
print('Probability of [R]=[R]-[R]',float(probs['[R]=[R]-[R]'])/probs['[R]~[R]~[R]'])
print('')

smarts = []
for element in elements:
  smarts.append('['+element+']')

#print(get_probs(smarts,smiles_list))

smarts = []
for element in elements:
  smarts.append('['+element+'R]')

tot_Ratoms,probs_Ratoms = get_probs(smarts,smiles_list)
#print (tot_Ratoms,probs_Ratoms)

R_elements = []
for key in probs_Ratoms:
  R_elements.append(key)
  
#print (R_elements)

smarts = []

for i,e1 in enumerate(R_elements):
  for e2 in R_elements:
    for j,e3 in enumerate(R_elements):
      if j >= i:
        sm_s = e1 + '-' + e2 + '-' + e3
        if sm_s not in smarts:
          smarts.append(sm_s)
      sm_d = e1 + '=' + e2 + '-' + e3
      if sm_d not in smarts:
        smarts.append(sm_d)
        
#print (len(smarts),smarts)

tot,probs = get_probs(smarts,smiles_list,ring=True)

#print (tot,probs)

sorted_x = sorted(probs.items(), key=operator.itemgetter(1), reverse=True)
count = 0
for i in range(len(sorted_x)):
  print (sorted_x[i][0],sorted_x[i][1]/tot)

print('')
#print (count)

rxn_smarts_rings = get_rxn_smarts_rings(probs)
#print (rxn_smarts_rings)

rxn_smarts_make_rings = get_rxn_smarts_make_rings(probs)
#print (rxn_smarts_make_rings)

p_rings = get_p(probs)
#print (p_rings)

#pickle.dump(p_rings,open('p_ring.p','wb')) 
#pickle.dump(rxn_smarts_rings,open('rs_ring.p','wb'))
#pickle.dump(rxn_smarts_make_rings,open('rs_make_ring.p','wb'))

smarts = []

for bond in bonds:
  for element1 in elements:
    for element2 in elements:
      smarts.append('['+element1+']'+bond+'['+element2+';!R]')

#print (len(smarts))
tot,probs = get_probs(smarts,smiles_list)
tot,probs = clean_probs(probs)
#print (tot, probs)
p = get_p(probs)
#print (p)

sorted_x = sorted(probs.items(), key=operator.itemgetter(1), reverse=True)
count = 0
for i in range(len(sorted_x)):
  print (sorted_x[i][0],sorted_x[i][1]/tot)

rxn_smarts = get_rxn_smarts(probs)
#print (rxn_smarts)

#pickle.dump(p,open('p1.p','wb')) 
#pickle.dump(rxn_smarts,open('r_s1.p','wb')) 


smarts_list = ['[*]1-[*]-[*]-1','[*]1-[*]=[*]-1','[*]1-[*]-[*]-[*]-1','[*]1=[*]-[*]-[*]-1','[*]1=[*]-[*]=[*]-1',
          '[*]1-[*]-[*]-[*]-[*]-1','[*]1=[*]-[*]-[*]-[*]-1','[*]1=[*]-[*]=[*]-[*]-1',
          '[*]1-[*]-[*]-[*]-[*]-[*]-1','[*]1=[*]-[*]-[*]-[*]-[*]-1','[*]1=[*]-[*]=[*]-[*]-[*]-1',
          '[*]1=[*]-[*]-[*]=[*]-[*]-1','[*]1=[*]-[*]=[*]-[*]=[*]-1']

smarts_macro = ['[r;!r3;!r4;!r5;!r6;!r8;!r9;!r10;!r11;!r12]','[r;!r3;!r4;!r5;!r6;!r7;!r9;!r10;!r11;!r12]',
         '[r;!r3;!r4;!r5;!r6;!r7;!r8;!r10;!r11;!r12]','[r;!r3;!r4;!r5;!r6;!r7;!r8;!r9;!r11;!r12]',
         '[r;!r3;!r4;!r5;!r6;!r7;!r8;!r9;!r10;!r12]','[r;!r3;!r4;!r5;!r6;!r7;!r8;!r9;!r10;!r11]']


tot, probs = get_probs(smarts_list,smiles_list,ring=True)

tot, probs = count_macro_cycles(smiles_list,smarts_macro,tot,probs)

num_rings = 0
print('')
for key in probs:
  print(key,probs[key])
  num_rings += probs[key]

print('')
print('number of rings',num_rings)
