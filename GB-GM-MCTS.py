'''
code adapted from https://github.com/haroldsultan/MCTS/blob/master/mcts.py by Jan H. Jensen 2018

'''
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdmolops

import time

import numpy as np
import random
import pickle

#import networkx as nx
import sascorer

from io import StringIO
import sys
#sio = sys.stderr = StringIO()

from rdkit import rdBase
rdBase.DisableLog('rdApp.error')

def logP_score(m):
  logp = Descriptors.MolLogP(m)
  SA_score = -sascorer.calculateScore(m)
  #cycle_list = nx.cycle_basis(nx.Graph(rdmolops.GetAdjacencyMatrix(m)))
  cycle_list = m.GetRingInfo().AtomRings() #remove networkx dependence
  if len(cycle_list) == 0:
      cycle_length = 0
  else:
      cycle_length = max([ len(j) for j in cycle_list ])
  if cycle_length <= 6:
      cycle_length = 0
  else:
      cycle_length = cycle_length - 6
  cycle_score = -cycle_length
  #print cycle_score
  #print SA_score
  #print logp
  SA_score_norm=(SA_score-SA_mean)/SA_std
  logp_norm=(logp-logP_mean)/logP_std
  cycle_score_norm=(cycle_score-cycle_mean)/cycle_std
  score_one = SA_score_norm + logp_norm + cycle_score_norm
  
  return score_one

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
  old_mol = Chem.Mol(mol)
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
  if valences_not_too_large(mol):
    return mol
  else:
    return old_mol


def expand_small_rings(mol): 
  Chem.Kekulize(mol,clearAromaticFlags=True)     
  rxn_smarts = '[*;r3,r4;!R2:1][*;r3,r4:2]>>[*:1]C[*:2]'
  while mol.HasSubstructMatch(Chem.MolFromSmarts('[r3,r4]=[r3,r4]')):
    mol = run_rxn(rxn_smarts,mol)

  return mol

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

# code modified from https://github.com/haroldsultan/MCTS/blob/master/mcts.py
import random
import math
import hashlib
#import logging
#import argparse
import numpy as np
 
SCALAR=1/math.sqrt(2.0)


class State():
  NUM_TURNS = 60	# max number of atoms. Later: average_size + 3*size_stdev?
  max_children = 25
  def __init__(self, mol=None, smiles='', turn=NUM_TURNS):
    self.mol = mol
    self.turn = turn
    self.smiles = smiles
    
  def next_state(self):
    #nextmove=random.choice(self.MOVES)
    smiles = self.smiles
    for i in range(100):
      mol = add_atom(self.mol)
      smiles = Chem.MolToSmiles(mol)
      if smiles != self.smiles:
        break

    next = State(mol, smiles, self.turn-1)
    return next

  def terminal(self):
    target_size = size_stdev*np.random.randn() + average_size
    if self.mol == None:
      num_atoms = 0
    else:
      num_atoms = self.mol.GetNumAtoms()

    if self.turn == 0 or num_atoms > target_size:
      self.mol = expand_small_rings(self.mol)
      self.smiles = Chem.MolToSmiles(self.mol)
      return True
    
    return False

  def reward(self):
    global max_score
    global count
    count += 1
    #logP = Descriptors.MolLogP(self.mol)
    logP = logP_score(self.mol)
 
    if logP > max_score[0]:
      max_score = [logP,self.smiles]

      return 1.0
    else:
      return 0.0
    #return logP/(1+abs(logP))
 
  def __hash__(self):
    return int(hashlib.md5(str(self.smiles).encode('utf-8')).hexdigest(),16)
  def __eq__(self,other):
    if hash(self)==hash(other):
      return True
    return False
  def __repr__(self):
    s="Value: %d; Moves: %s; Turn %s"%(self.value,self.moves,self.turn)
    return s
	

class Node():
  def __init__(self, state, parent=None):
    self.visits=1
    self.reward=0.0	
    self.state=state
    self.children=[]
    self.parent=parent	
  def add_child(self,child_state):
    child=Node(child_state,self)
    self.children.append(child)
  def update(self,reward):
    self.reward+=reward
    self.visits+=1
  def fully_expanded(self):
    if len(self.children)==self.state.max_children:
      return True
    return False
  def __repr__(self):
    #s="Node; children: %d; visits: %d; reward: %f; value %d"%(len(self.children),self.visits,self.reward,self.state.value)
    s=str(self.state.smiles)
    return s
		

def UCTSEARCH(budget,root):
  for iter in range(int(budget)):
    #print iter
    front=TREEPOLICY(root)
    #print 'len(node.children)',len(front.children)
    for child in front.children:
      #print "computing reward"
      reward=DEFAULTPOLICY(child.state)
      BACKUP(child,reward)
  return BESTCHILD(root,0)

def TREEPOLICY(node):
  #a hack to force 'exploitation' in a game where there are many options, and you may never/not want to fully expand first
    #print 'len(node.children)',len(node.children)
    while node.fully_expanded():
      node = BESTCHILD(node,SCALAR)
    
    #print 'best child',node.state.moves
    if node.state.terminal():
      return node
    else:
      #print 'expanding',node.state.moves
      node = EXPAND_ALL(node)
      return node
    
def EXPAND_ALL(node):
  lcount = 0
  while not node.fully_expanded() and lcount < node.state.max_children:
    lcount += 1
    node = EXPAND(node)
    #print 'not enough children, expanding', node.state.smiles, len(node.children)
  return node

def EXPAND(node):
  tried_children=[c.state for c in node.children]
  new_state=node.state.next_state()
  lcount = 0
  while new_state in tried_children and lcount < new_state.max_children:
    lcount += 1
    new_state=node.state.next_state()
  node.add_child(new_state)
  return node

#current this uses the most vanilla MCTS formula it is worth experimenting with THRESHOLD ASCENT (TAGS)
def BESTCHILD(node,scalar):
  bestscore=0.0
  bestscore=-99.
  bestchildren=[]
  for c in node.children:
    exploit=c.reward/c.visits
    explore=math.sqrt(2.0*math.log(node.visits)/float(c.visits))	
    score=exploit+scalar*explore
    #print score, bestscore
    if score==bestscore:
      bestchildren.append(c)
    if score>bestscore:
      bestchildren=[c]
      bestscore=score
  if len(bestchildren)==0:
    print ("OOPS: no best child found, probably fatal")
  return random.choice(bestchildren)

def DEFAULTPOLICY(state):
	while state.terminal()==False:
		state=state.next_state()
	return state.reward()

def BACKUP(node,reward):
	while node != None:
		node.visits+=1
		node.reward+=reward
		node=node.parent
	return

global max_logP, count
global p_ring, p_make_ring, rxn_smarts_make_ring, rxn_smarts_ring_list, rxn_smarts_list, p
global average_size, size_stdev

logP_values = np.loadtxt('logP_values.txt')
SA_scores = np.loadtxt('SA_scores.txt')
cycle_scores = np.loadtxt('cycle_scores.txt')
SA_mean =  np.mean(SA_scores)
SA_std=np.std(SA_scores)
logP_mean = np.mean(logP_values)
logP_std= np.std(logP_values)
cycle_mean = np.mean(cycle_scores)
cycle_std=np.std(cycle_scores)

average_size, size_stdev = 39.15, 3.50

p_ring = pickle.load(open('p_ring.p','rb'))
p_make_ring = p_ring
rxn_smarts_make_ring = pickle.load(open('rs_make_ring.p','rb'))
rxn_smarts_ring_list = pickle.load(open('rs_ring.p','rb'))

rxn_smarts_list = pickle.load(open('r_s1.p','rb'))
p = pickle.load(open('p1.p','rb'))

prob_double = 0.8
p_ring = scale_p_ring(rxn_smarts_ring_list,p_ring,prob_double)
p_make_ring = p_ring

num_sims = 40 # 40 = 1000 logP evaluations

print('num_sims',num_sims)
print('average_size, size_stdev', average_size, size_stdev)
print('max_children = 25')
print('')

results = []
size = []
t0 = time.time()
for i in range(10):
  max_score = [-99999.,'']
  count = 0

  smiles = 'CC'
  mol = Chem.MolFromSmiles(smiles)

  current_node = Node(State(mol,smiles))
  current_node = UCTSEARCH(num_sims,current_node)
  print(i, max_score[0], max_score[1], Chem.MolFromSmiles(max_score[1]).GetNumAtoms())
  results.append(max_score[0])
  size.append(Chem.MolFromSmiles(max_score[1]).GetNumAtoms())

t1 = time.time()
print('')
print('time, count ',t1-t0, count)
print(max(results),np.array(results).mean(),np.array(results).std())
print(max(size),np.array(size).mean(),np.array(size).std())


