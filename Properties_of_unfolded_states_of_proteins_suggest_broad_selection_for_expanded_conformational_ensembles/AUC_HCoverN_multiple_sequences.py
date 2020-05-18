#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 14:25:28 2018

@author: anabel
"""

'This script calculates global z-scores and HC/N (area under curve)'
'Input file is a csv table with sequences of interest'
'This script will read through each line (sequence) and perform the calculation'


import pandas as pd
import numpy as np

filename = input('Name your file:')

f1 = input("Drag and drop the csv file into command line:")
file1 = f1.strip('\'')   #removes quotes from string input
df = pd.read_csv(file1)

MJ_HW = {"A":-0.0645, "C":0.502, "D":-1.338, "E":-1.348, "F":1.483, "G":-0.382, "H":-0.0737, "I":1.144, "K":-1.472, "L":1.329,
      "M":0.767, "N":-0.551, "P":-0.424, "Q":-0.518, "R":-1.209, "S":-0.540, "T":-0.216, "V":0.811, "W":1.289, "Y":0.812
      }

######## removes symbol * used to designate stop codon  #########
valid_amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M',
                     'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

raw_AA_seq = str(input("Write column name containing sequences to analyze:"))
contains_X = df[raw_AA_seq] #pulls out sequences as a series to loop over it

all_only_aa = []
for sequence in contains_X:
    only_real_aa = []
    for aa in sequence:
        if aa in valid_amino_acids:     #if aa is not in the valid aa list, it will be removed. If an entire row contains X, the row is kept, but will be empty
            only_real_aa.append(aa)
    all_only_aa.append(only_real_aa)

#format the list to be a string of characters, which is needed for script below
aalist_join = []
for aalist in all_only_aa:
    join = ''.join(aalist)   #removes quotes and commas, creates one big list, each protein sequence is now separated by a comma
    aalist_join.append(join)

df['AA_seq_noSTOPcodon'] = aalist_join    #I don't need to convert this list to pd.Series because of the format I created in for loop
        



##########################################################################################
######################## assign z-scores to all amino acids #############################
##########################################################################################

#retreive amino acid sequence columns from DataFrame as a series to loop over it
np_AAseq = df['AA_seq_noSTOPcodon']

#add a column with protein length
protein_length = []
for AAseq in np_AAseq:
    length = len(AAseq)
    protein_length.append(length)
df['protein length'] = pd.Series(protein_length)    

#go through every line in PDB file, and assign z scores to each amino acid for every protein
Zscores_proteins = []
for protein in np_AAseq:
    ListofValues = []
    #for each protein, assign z score to each amino acid
    for aa in protein :
        hydrophobicity = MJ_HW[aa]
        ListofValues.append(hydrophobicity)
    Zscores_proteins.append(ListofValues)

#add a column with z scores to existing dataframe for every amino acid for all proteins
df['aa z scores'] = pd.Series(Zscores_proteins) 



##########################################################################################
############ calculate global (average) z-scores and add as a new column #################
##########################################################################################

np_Zscores = df['aa z scores'] #pulls out z scores as a series to do math

globalZscores = []
for Zlist in np_Zscores:
    avgZ = np.average(Zlist)
    globalZscores.append(avgZ)

df['Global z score'] = pd.Series(globalZscores)



##########################################################################################
###### calculate average Z score in a defined window and add as a new column #############
##########################################################################################


# set up parameters: window size and step size needed for user-defined function: calculateSlidingWindow    
windowSize = int(input('Enter sliding window size: '))
while windowSize%2 != 1:
    windowSize = int(input('**ERROR**: Please enter an odd numbered window size: '))
if windowSize >= 20:
    windowSize = int(input('**ERROR**: Window size can not be larger than 20 amino acids. Please re-enter window size: '))
#the nr dataset does not have proteins smaller than 20 amino acids, so this is upper limit of window size

stepSize = int(input('Enter step size: '))
if stepSize > windowSize:
    stepSize = int(input('**ERROR**: Step size can not not be larger than window size. Please re-enter step size: '))    

print("Computing...")
print("This might take a while...")

#define the function that will average z scores across a user-defined window
def calculateAvgZscore(zscores_list, windowSize, stepSize):
    "This function calculates the average z-score for a defined amino acid window and step size" 
    windowZscoreValues = []
    avgZscore = []
    
    numOfchunks = int(((len(zscores_list)-windowSize)/stepSize)+1)

    #creates the sliding window
    for values in range(0,numOfchunks*stepSize,stepSize):
        window = zscores_list[values:values+windowSize]
        windowZscoreValues.append(window)
    
    #calculates the average of the defined sliding window at window center, and adds the value NaN to amino acids that are not the center of window
    for aa in windowZscoreValues :
        np_window = np.array(aa)
        avgaaWindow = np.average(np_window)
        avgZscore.append(avgaaWindow)
    
    return avgZscore

avgWindowZscore = []
for zscores_list in np_Zscores:        
    windowZscore = calculateAvgZscore(zscores_list, windowSize, stepSize)
    avgWindowZscore.append(windowZscore)
df['avg window z score'] = pd.Series(avgWindowZscore)



##########################################################################################
###### calculate the area under the curve, generated from avg Z scores #############
##########################################################################################

np_windowZscores = df['avg window z score']  #pull out column as a series to do math

def abstrapz(y, x=None, dx=1.0):
    "This function re-writes the source code for np.trapz() to remove contributions of negative values" 
    y = np.asanyarray(y)
    if x is None: #when x is None, points are assumed to be equally spaced dx apart
        d = dx
    else:
        x = np.asanyarray(x)
        d = np.diff(x)
    ret = (d * (y[1:] +y[:-1]) / 2.0)
    return ret[ret>0].sum()  #The important line, only positive areas are summed


pos_AUCs = []
HCoverN_list = []
for y in np_windowZscores :   
    pos_AUC = abstrapz(y, x=None, dx = 1)
    pos_AUC = round(float(pos_AUC), 4)
    HCoverN = pos_AUC/len(y)
    HCoverN = round(float(HCoverN),4)
    pos_AUCs.append(pos_AUC)
    HCoverN_list.append(HCoverN)


df['pos_AUCs'] = pd.Series(pos_AUCs)

df['HC/N'] = pd.Series(HCoverN_list)

#new_col = raw_AA_seq_noSTOPcodon + ' HC/N'
#df.rename(columns={'HC/N':new_col}, inplace=True) #renames column  to specify HC/N

#before file export, remove rows with NaN values since this will return error when plotting
#and does not provide any information
df = df.dropna(axis=0, inplace = False)  #this drops the entire row if NaN is present, inplace false means it's dropped instead of replacing NaN with 'None')
df = df.drop('avg window z score', axis = 1, inplace = False)
#export csv file
df.set_index('EcoGene_ID')  #sets EcoGene ID as index
df.to_csv(str(filename) + '.csv', index=False)


