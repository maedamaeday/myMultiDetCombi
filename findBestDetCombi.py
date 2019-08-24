import numpy as np
import math
import argparse
import itertools
import copy
import math

import myMultiDetCombi as mdc

def findBest( FPs, FN = mdc.def_common_FN, method = "MooN", verbose = False ):

   all_results = {}
   reduced_FP = FPs.copy()
   best_FP = min( FPs )
   for i in range(len(FPs)-1): #exclude M=1 case 
      N = len(reduced_FP)
      detectorData = mdc.myMultipleDetectorData( reduced_FP, FN = FN )
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

                                  
def main( FPs, FN = mdc.def_common_FN, obs = None ):

   if obs is not None:
      detectorData = mdc.myMultipleDetectorData( FPs, FN = FN, obs = obs )
      print( detectorData.getLikelihood() )
      return

   print( f"input false positive rates : {FPs}" )
   print( f"input false negative rate : {FN} (common for all the detectors)" )
   print( "****M-out-of-N logic****" )
   try:
      bestM, bestN, best_FP, best_FN = findBest( FPs, FN, verbose = True )
   except TypeError:
      print( "No improvement is expected." )
   else:
      print( f"best condition : N = {bestN}, M = {bestM}, FP = {best_FP}, FN = {best_FN}" )

   print( "****likelihood methods****" )
   try:
      bestM, bestN, best_FP, best_FN = findBest( FPs, FN, method = 'likelihood',
                                                 verbose = True )
   except TypeError:
      print( "No improvement is expected." )
   else:
      print( f"best condition : N = {bestN}, FP = {best_FP}, FN = {best_FN}" )

   return

if __name__ == '__main__':

   parser = argparse.ArgumentParser(description='find the best combination to maximize detector performance with M-out-of-N logic and likelihood method',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   parser.add_argument('-p', '--FPs', '--false_positives', nargs='+', type=float, required=True,
                       help='list of false positive rates for each detector with fixed false negative rate given by the "--FN" option.')
   parser.add_argument('-n', '--FN', '--false_negative', type=float, default=mdc.def_common_FN,
                       help='false negative rate, which is assumed to be common among detectors')
   parser.add_argument('-o', '--obs', nargs='+', type=int, default=None,
                       help='list of observations for each detector. 1 means positive detection and 0 means negative detection. When this option is specified, only a result of negative log likelihood for positive/negative assumption is returned.')
   args = parser.parse_args()
   FPs = args.FPs
   FN = args.FN
   obs = args.obs
   main( FPs = FPs, FN = FN, obs = obs )
                       
