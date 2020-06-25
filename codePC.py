#!/usr/bin/env python3

import os
import sys
import shutil
import time
import math
import csv
import statistics
import threading
import functools
import itertools
import scipy
import numpy as np
import pandas as pd
import argparse as ap

from numpy import *
from scipy import stats
from scipy.stats import norm





def outputFile(position, alphaV, filterF, countSizeTile, ind):
    filter_pos = list()
    filter_pos = set(position)
    file_name = " "
    
    if alphaV == alpha and len(filter_pos) != 0:
        file_name = "pc_alpha_{}.txt".format(alphaV)
        print("Generating output file >> {}".format(file_name))

    with open(file_name, "w+") as file:
        for i in filter_pos:
            file.write("{} ".format(i))
    

    filter_func(filterF, filter_pos, countSizeTile, ind)


##-- Information about the time execution
def info(s, init_new_line=False):
    if s:
        nfo = '\n' if init_new_line else ''
        nfo += '[i] '
        sys.stdout.write(nfo + str(s) + '\n')
        sys.stdout.flush()

def getCorrelation(data, rows, columns, alpha, genesFrom, filter_file):
    refRow =  0 
    countCor = 0

    genes = []
    while(i <= rows):
        genes.append(csv_file[ref].iloc[refRow])
        refRow += 1
        i += 1

    indP = 0
    if genesFrom in genes and genesFrom != "":
        indP = genes.index(genesFrom)
    
    X = []
    Y = []
    Final_pos = list()
    IDx = ''
    IDy = ''

    txtGenX = "Name of Gene {}: {}"
    txtGenY = "Name of Gene {}: {}"
    txtPear = "The correlation value of these genes is: {}"
    txtPV = "The p-value of these genes is: {}"

    for i in csv_file.iloc[indP]:
        if i in genes:
            IDx = i
            continue
        try:
            X.append(float(i))
        except ValueError:
            v = i.replace(',', '.')
            X.append(float(v))

    ps = 0
    while ps < refRow -1:
        ps += 1

        for i in csv_file.iloc[ps-1]:
            if i in genes:
                IDy = i
                continue
            try:
                Y.append(float(i))
            except ValueError:
                v = i.replace(',', '.')
                Y.append(float(v))
        
        #pearson correlation 
        correlation, p_value = stats.pearsonr(X, Y)

        fisher_z = (1/2) * np.log((1+correlation)/(1-correlation))

        z_decision = sqrt(columns - 3) * abs(fisher_z)

        # inverse cumulative distribution functions
        inv_cdf = norm.ppf(1-alpha/2)

        if z_decision > inv_cdf :
            countCor += 1
            Final_pos.append(ps-1)
        
            print("\nPositions of the genes: [{},{}]".format(indP, ps-1), '\n')

            print(txtGenX.format(indP, IDx))
            print(txtGenY.format(ps-1, IDy),'\n')
            print(txtPear.format(correlation))
            print(txtPV.format(p_value))
            print("Fisher's z-transform: {}".format(fisher_z))
            print("Equation to reject (kalisch): {}".format(z_decision))
            print("inv_cdf: ", inv_cdf)
            print('_______________________________________________________')
        
        Y.clear()
        
    outputFile(Final_pos, alpha, filter_file, countCor, indP)


   
def read_params():
    par = ap.ArgumentParser()
    arg = par.add_argument

    arg('-m', '--mode', type=int, required=False, help='>> mode = 2 (filter file) mode = 1 (complete analysis)')
    arg('-csv', '--csvFile', type=str,required=False, help='')
    arg('-lgn', '--lgnFile', type=str,required=False, help='')
    arg('-a', '--alpha', type=float, required=False, help='')
    arg('-t', '--tileFile', type=str, required=True, help='tile file to filter')
    arg('-s', '--sizeTile', type=int, required=False, help='')
    arg('-df', '--dataFile', type=str,required=False, help='>> list of genes (not filtered)')
    
       

    args = par.parse_args()

    return vars(par.parse_args())


def filter_func(filterF, data_file, countSize_tile, ind):
        
    f2 = open(filterF, 'r')

    if(size_tile != None):
        countSize_tile = size_tile

    data = [[int(v) for v in line.split()] for line in f2]
    
    outfile = "filtered_file.txt"      
    f3 = open(outfile, "w")

    print("Generating output file >> {}".format(outfile))

    f3.write("{}\n".format(csvFile))
    for j in data:
        for i in j:
            if i in data_file:
                f3.write("{} ".format(i))
        f3.write("\n")

 
def filter_mode(filterF, dataF, size,csv):
       
    dataFile = open(dataF, 'r')

    inputData = [[int(v) for v in line.split()] for line in dataFile]

    filterFile = open(filterF, 'r')

    filter_Data = [[int(v) for v in line.split()] for line in filterFile]

    list_data = []
    for i in inputData:
        for j in i:
            list_data.append(j)
    
    out_name = "provafiltered_file.txt"
    out = open(out_name, "w+")
    print("Generating output file >> {}".format(out_name))
    
    out.write("{}\n".format(csv))
    for j in filter_Data: 
        for i in j:
            if i in list_data:
                out.write("{} ".format(i))
        out.write("\n")
        
    

if __name__ == '__main__':
    t0 = time.time()
    pars = read_params()

    exec_mode = 1
    alpha = pars['alpha']
    filter_file = pars['tileFile']
    size_tile = pars['sizeTile']
    exec_mode = pars['mode']
    dataFile = pars['dataFile']
    csvFile = pars['csvFile']
    lgnFile = pars['lgnFile']

    if(exec_mode != 2):

        #------------------------------------ IMPORT csv file

        csv_file = pd.read_csv(csvFile, delimiter= ';', header=None, skiprows=1 , low_memory=False)
        
        if(lgnFile != None):
            lgn_input = pd.read_csv(lgnFile, delimiter= ',')
            lgn_from = lgn_input.loc[0, 'from']
            lgn_to = lgn_input.loc[0, 'to']
        else:
            lgn_from = ""
            lgn_to = ""
        

        #-----------------------------------------------------
        Ncolumns = len(csv_file.columns)-1
        Nrows = csv_file.shape[0]-1
        #-----------------------------------------------------
        
        print(">> Computation for \u03B1 = {}\n".format(alpha))
        getCorrelation(csv_file, Nrows, Ncolumns, alpha, lgn_from, filter_file)

    else:
        if(csvFile == 0):
            print("\nMissing parameters: -csv <name csv>\n")
            sys.exit(0)
        print("Filtering analysis")
        filter_mode(filter_file, dataFile, size_tile, csvFile)

    t1 = time.time()   
    print('\n')
    info('Computation end. Execution time: {}s'.format(int(t1-t0)))


sys.exit(0)
