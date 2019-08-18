import numpy as np
import math
import argparse
import itertools
import copy
import math

def_common_FN = 0.01 #false negative rate

class myMultipleDetectorData:

   def setObs( self, obs, testOnly = False ):
      if len(obs) != self.N:
         raise ValueError("The numbers of elements for false positives and observations are not the same.")
      if not testOnly:
         self.obs = obs

   def __init__( self, FPs, FN = def_common_FN, obs = None ):
      self.FPs = FPs
      self.FN = FN
      self.N = len(self.FPs)
      if min(self.FPs)<0 or max(self.FPs)>1:
         raise ValueError("Invalid false positive value is found. Must be within [0,1]")
      if self.FN<0 or self.FN>1:
         raise ValueError("Invalid inefficiency value is found. Must be within [0,1]")
      if obs is not None:
         self.setObs( obs )
      self.PartialMooN_FP = None
      self.PartialMooN_FN = None
      
   def getPartialMooNProb( self ):
      self.PartialMooN_FP = []
      self.PartialMooN_FN = []
      indexList = [ i for i in range(self.N) ]
      for M in range(self.N+1):
         partial_FP = 0
         partial_FN = 0
         for hitList in itertools.combinations(indexList,M):
            singlePat_FP = 1 
            singlePat_FN = 1
            for i in range(self.N):
               FP = self.FPs[i]
               FN = self.FN
               singlePat_FP *= FP if i in hitList else (1-FP)
               singlePat_FN *= (1-FN) if i in hitList else FN
            partial_FP += singlePat_FP
            partial_FN += singlePat_FN
         self.PartialMooN_FP.append( partial_FP )
         self.PartialMooN_FN.append( partial_FN )

   def getCombinedMooNProb( self, M ):
      if M<0 or M>self.N:
         raise IndexError( f"M={M} is not valid, must be within [0,{self.N}]" )
      if self.PartialMooN_FP is None:
         self.getPartialMooNProb()

      combined_FP = 0
      combined_FN = 0
      for m in range(0,self.N+1):
         if m<M:
            combined_FN += self.PartialMooN_FN[m]
         else:
            combined_FP += self.PartialMooN_FP[m]

      return combined_FP, combined_FN

   def getLikelihood( self, obs=None ):
      if obs is None:
         if self.obs is None:
            print( "No hit infomation is given." )
            return None
         else:
            obs = self.obs
      else: 
         self.setObs( obs, testOnly=True )
         
      positive_NLL = 0 #Negative Log Likelihood for anomaly assumption
      negative_NLL = 0 #Negative Log Likelihood for normal assumption
      for i in range(self.N):
         positive_NLL -= math.log(1-self.FN) if obs[i]>0.5 else math.log(self.FN)
         negative_NLL -= math.log(self.FPs[i]) if obs[i]>0.5 else math.log(1-self.FPs[i])

      return positive_NLL, negative_NLL

   def getLikelihoodProb( self ):
      combined_FP = 0
      combined_FN = 0
      for i in range(2**self.N):
         obs = []
         tmp = i  
         positive_prob = 1
         negative_prob = 1
         for j in range(self.N):
            obs.append( tmp%2 )
            positive_prob *= (1-self.FN) if obs[-1]>0.5 else self.FN
            negative_prob *= self.FPs[j] if obs[-1]>0.5 else (1-self.FPs[j])
            tmp = int(tmp/2)
         positive_NLL, negative_NLL = self.getLikelihood( obs )
         if positive_NLL<negative_NLL:
            combined_FP += negative_prob
         else:
            combined_FN += positive_prob
      return combined_FP, combined_FN

def findBest( FPs, FN = def_common_FN, method = "MooN", verbose = False ):

   all_results = {}
   reduced_FP = FPs.copy()
   best_FP = min( FPs )
   for i in range(len(FPs)-1): #exclude M=1 case 
      N = len(reduced_FP)
      detectorData = myMultipleDetectorData( reduced_FP, FN = FN )
      for m in range(N+1): 
         combined_FP, combined_FN = ( detectorData.getCombinedMooNProb( m )
                                      if method=="MooN" else detectorData.getLikelihoodProb() )
         if combined_FP>best_FP or combined_FN>FN:
            continue
         Fvalue = 2*(1-combined_FP)*(1-combined_FN)/(2-combined_FP-combined_FN)
         all_results[(m,N)] = (combined_FP, combined_FN, Fvalue)
         if method!="MooN":
            break
      combined_FP, combined_FN = detectorData.getCombinedMooNProb( m )
      iWorstFP = reduced_FP.index(max(reduced_FP))
      reduced_FP.pop(iWorstFP)

   if len(all_results) == 0:
      return None

   best_Fvalue = None
   best_key = None
   for key in all_results.keys():
      fvalue = all_results[key][2]
      if verbose:
         print( f"{key} : {all_results[key]}" )
      if best_Fvalue is None or best_Fvalue<fvalue:
         best_Fvalue = fvalue
         best_key = key

   return best_key[0], best_key[1], all_results[best_key][0], all_results[best_key][1] #bestN, bestM, combined_FP, FN
