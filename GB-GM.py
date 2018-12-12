'''
code written by Jan H. Jensen 2018
'''
from rdkit import Chem
from rdkit.Chem import AllChem
import random
import numpy as np
import pickle
from io import StringIO
import sys
#sio = sys.stderr = StringIO()
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

def valences_not_too_large(mol):
    valence_dict = {5:3, 6:4, 7:3, 8:2, 9:1, 16:6, 17:1, 35:1, 53:1}
    atomicNumList = [a.GetAtomicNum() for a in mol.GetAtoms()]
    valences = [valence_dict[atomic_num] for atomic_num in atomicNumList]
    BO = Chem.GetAdjacencyMatrix(mol,useBO=True)
    number_of_bonds_list = BO.sum(axis=1)
    for valence, number_of_bonds in zip(valences,number_of_bonds_list):
        if number_of_bonds > valence:
            return False

    return True

def run_rxn(rxn_smarts,mol):
  new_mol_list = []
  patt = rxn_smarts.split('>>')[0]
  # work on a copy so an un-kekulized version is returned
  # if the molecule is not changed
  mol_copy = Chem.Mol(mol)
  try:
    Chem.Kekulize(mol_copy)
  except:
    pass
  if mol_copy.HasSubstructMatch(Chem.MolFromSmarts(patt)):
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    new_mols = rxn.RunReactants((mol_copy,))
    for new_mol in new_mols:
      try:
        Chem.SanitizeMol(new_mol[0])
        #new_mol_list.append(Chem.MolFromSmiles(Chem.MolToSmiles(new_mol[0])))
        new_mol_list.append(new_mol[0])
      except:
        pass
    if len(new_mol_list) > 0:
      new_mol = random.choice(new_mol_list) 
      return new_mol
    else:
      return mol
  else:
      return mol

def add_atom(mol):
  if np.random.random() < 0.63: # probability of adding ring atom
    rxn_smarts = np.random.choice(rxn_smarts_ring_list, p=p_ring)
    if not mol.HasSubstructMatch(Chem.MolFromSmarts('[r3,r4,r5]'))\
       or AllChem.CalcNumAliphaticRings(mol) == 0:
      rxn_smarts = np.random.choice(rxn_smarts_make_ring, p=p_make_ring)
      if np.random.random() < 0.036: # probability of starting a fused ring
        rxn_smarts = rxn_smarts.replace("!", "")
  else:
    if mol.HasSubstructMatch(Chem.MolFromSmarts('[*]1=[*]-[*]=[*]-1')):
      rxn_smarts = '[r4:1][r4:2]>>[*:1]C[*:2]'
    else:
      rxn_smarts = np.random.choice(rxn_smarts_list, p=p)
    
  mol = run_rxn(rxn_smarts,mol)
  smiles = Chem.MolToSmiles(mol)

  return mol, smiles

def expand_small_rings(mol):  
  Chem.Kekulize(mol,clearAromaticFlags=True)
  rxn_smarts = '[*;r3,r4;!R2:1][*;r3,r4:2]>>[*:1]C[*:2]'
  count = 0
  while mol.HasSubstructMatch(Chem.MolFromSmarts('[r3,r4]=[r3,r4]')):
    mol = run_rxn(rxn_smarts,mol)
    #print('expanding ring',Chem.MolToSmiles(mol))
#    count += 1
#    if count > 5:
#    	print('breaking in expand_small_rings')
#    	sys.exit()
    
  return mol

def generate_mol(smiles,max_atoms,average_size,size_stdev):
  #np.random.seed(0)
  mol = Chem.MolFromSmiles(smiles)
  num_atoms = mol.GetNumAtoms()
 
  target_size = size_stdev*np.random.randn() + average_size

  count = 0
  while num_atoms < max_atoms and count < max_atoms:
    count += 1
    mol, smiles = add_atom(mol)
    num_atoms = mol.GetNumAtoms()
    if num_atoms > target_size:
      break
  
  if valences_not_too_large(mol):   
    mol = expand_small_rings(mol)
    return mol
  else:
    return None

def scale_p_ring(rxn_smarts_ring_list,p_ring,new_prob_double):
  p_single = []
  p_double = []
  for smarts,p in zip(rxn_smarts_ring_list,p_ring):
    if '=' in smarts:
      p_double.append(p)
    else:
      p_single.append(p)
    
  prob_double, prob_single = sum(p_double), sum(p_single)
  scale_double = new_prob_double/prob_double
  scale_single = (1.0 - new_prob_double)/(1-prob_double)
  for i, smarts in enumerate(rxn_smarts_ring_list):
    if '=' in smarts:
      p_ring[i] *= scale_double
    else:
      p_ring[i] *= scale_single
      
  print(scale_double,scale_single*prob_single,sum(p_ring))
  
  return p_ring

#########################

import time

average_size, size_stdev = 23.2, 4.4
p_ring = pickle.load(open('p_ring.p','rb'))
p_make_ring = p_ring
rxn_smarts_make_ring = pickle.load(open('rs_make_ring.p','rb'))
rxn_smarts_ring_list = pickle.load(open('rs_ring.p','rb'))

rxn_smarts_list = pickle.load(open('r_s1.p','rb'))
p = pickle.load(open('p1.p','rb'))

prob_double = 0.8
p_ring = scale_p_ring(rxn_smarts_ring_list,p_ring,prob_double)
p_make_ring = p_ring

t0 = time.time()

smiles = "CC"

mol_list = []

with open('grow1000.smi','w') as file:
  count = 1
  while count <= 1000:
    mol = generate_mol(smiles,50,average_size,size_stdev)
    if mol:
      new_smiles = Chem.MolToSmiles(mol)
      file.write(new_smiles+'\n')
      count += 1

 
t1 = time.time()

print (t1-t0)