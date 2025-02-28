# Selection methods.
# Copyright (c)
# Authors: Jonatan Gomez and Elizabeth León  
# E-mails: jgomezpe@unal.edu.co  and eleonguz@unal.edu.co
# All rights reserved.
# Licence
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
# in the documentation and/or other materials provided with the distribution.
# Neither the name of the copyright owners, their employers, nor the names of its contributors may be used to endorse or 
# promote products derived from sgoal.this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import random as rand

############### UTILITY FUNCTIONS ################
# A permutation of n elements
def permutation(n):
  x = [i for i in range(0,n)]
  rand.shuffle(x)
  return x

# Normalizes a vector of weights
def normalize(weight):
  c = 0
  for w in weight: 
    c += w
  return [ w/c for w in weight ]

# Gest the best index (position with the best quality value)
# Minimization version
def max_best(quality):
  k = 0
  for i in range(1,len(quality)):
    if(quality[i] >= quality[k]): 
      k = i
  return k

# Maximization version
def min_best(quality):
  k = 0
  for i in range(1,len(quality)):
    if(quality[i] <= quality[k]): 
      k = i
  return k
 
############### SELECTION MECHANISMS ################
#### Uniform selection. Picks N elements (indices) with the same probability ####
def uniform(quality, N):
  return [rand.randint(0,len(quality)-1) for i in range(N)]

#### Tournament selection. 
# Picks 4 individuals at random and returns the best one

# Selects 1 individual. ####
# Minimization version
def min_tournament1(quality):
  m = 4 #Tournament's size
  candidate = uniform(quality, m)
  x = 0
  for k in range(1,m):
    if(quality[candidate[x]] >= quality[candidate[k]]):
      x = k
  return candidate[x]

# Tournament selection. Selects N individuals.
# Uses tournament1 for each individual to be selected
# Minimization version
def min_tournament(quality, N):
  return [min_tournament1(quality) for i in range(N)]

# Tournament 1. Maximization version
def max_tournament1(quality):
  m = 4 #Tournament's size
  candidate = uniform(quality, m)
  x = 0
  for k in range(1,m):
    if(quality[candidate[x]] <= quality[candidate[k]]):
      x = k
  return candidate[x]

# Tournament selection. Selects N individuals.
# Uses tournament1 for each individual to be selected
# Maximization version
def max_tournament(quality, N):
  return [max_tournament1(quality) for i in range(N)]

#### Weighted selection: p_i is the probability of selecting element i ####
def weighted(p):
  y = rand.random()
  k=0
  while k<len(p) and y>=p[k]:
    y -= p[k]
    k+=1
  return k

# Adjust function values (quality arrays) to real quality measures (q_i > 0)
def adjustquality(quality):
  m = min(quality)
  i=0
  n = len(quality)
  while(i<n and m==quality[i]):
    i+=1
  if(i==n):
    return [1 for i in range(n)]
  m2 = quality[i]
  i+=1
  while(i<n):
    if(m < quality[i] and m2>quality[i]):
      m2 = quality[i]
    i+=1
  d = m2 - m
  return normalize([q - m + d for q in quality])

# Adjust function values (quality arrays) minimization version
def min_adjustquality(quality): return adjustquality([-q for q in quality])

# Adjust function values (quality arrays) maximization version
def max_adjustquality(quality): return adjustquality(quality.copy())

#### Roulette wheel selection. Selects N individuals. ####
# Minimization version
def min_roulette(quality, N):
  p = min_adjustquality(quality)
  return [weighted(p) for i in range(N)]

# Maximization version
def max_roulette(quality, N):
  p = max_adjustquality(quality)
  return [weighted(p) for i in range(N)]