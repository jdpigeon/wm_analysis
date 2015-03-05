# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 13:58:44 2015

@author: Dano

Purpose: Load Water Maze data into python arrays

Input: a csv file containing either training or probe data
"""
import pandas as pd
import scipy as sp
import scipy.stats
import pickle #To unpickle the tank coordinates out of the .ann file
import numpy as np
import glob #For use of wild cards in getting .ann files
import matplotlib.pylab as plt
import matplotlib.mlab as mlab

pd.set_option('display.precision',5)

def load_data(filename):
    """
    loads data from csv file, cleaning up rows and columns that include only nan entries
    
    Input: filename    
    Output: dataframe (df) object of data
    """
    probe = pd.read_csv(filename, index_col="Animal")
    probe = probe.dropna(axis=1, how='all')
    probe = probe.dropna(axis=0, how='any')
    return probe
    
    
def label_groups(frame):
    
    """
    labels every other 2 mice as GFP or APP (ie. m1256 = GFP, m3478= APP)
    """
    frame['Group'] = 'GFP'
    i = 0
    j = 0
    for index, row in frame.iterrows():
        if j % 2 == 0:
            frame['Group'][index] = 'GFP'
            i = i + 1
            if i % 2 == 0:
                j = j + 1
            continue
        frame['Group'][index] = 'APP'
        i = i + 1
        if i % 2 == 0:
            j = j + 1
    
    return frame
    
##########Graphing Functions##############
    
def single(frame):
    """
    produces a graph of time in zone for an individual mouse's probe trial
    In the future, may be able to take string input and find appropriate data from master frame
    as of now, requires individual mouse's zone data as input
   
    """
    plt.figure(figsize=(12,12))
    frame.plot(kind='bar')
    plt.title('Time spent in Zone', fontsize='xx-large')
    plt.ylabel('Time (sec)', fontsize='xx-large')
    plt.xlabel('')
    plt.yticks(fontsize='xx-large')
    locs, labels = plt.xticks()
    plt.setp(labels,rotation=0,fontsize='xx-large')
    plt.legend(['Zone','1','T','2','4'], fontsize='xx-large')
    fig = matplotlib.pyplot.gcf()
    fig.set_size_inches(8,5)
    fig.savefig('test2png.png',dpi=100)

def zoneplot(frame, err):
    """
    produces a graph for zone plots. Works for both all zone and zone vs other
    """
    frame.plot(kind='bar',yerr=err)
    plt.title('Time spent in Zone')
    plt.ylabel('Time (sec)')
    plt.xlabel('Zone')
    locs, labels = plt.xticks()
    plt.setp(labels,rotation=0)
    
    
def quadplot(frame, err):
    """
    produces a graph for zone plots. Works for both all zone and zone vs other
    """
    frame.plot(kind='bar',yerr=err)
    plt.title('Time spent in Quadrant')
    plt.ylabel('Time (%)')
    plt.xlabel('Quadrant')
    locs, labels = plt.xticks()
    plt.setp(labels,rotation=0)
    
def confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m-h


    """
    

To do:
1. write exclude mouse function
2. Change position of bars in zoneplot so that manipulation groups are seperate, with their distinct T O bars
3. x
"""
##########################  

if __name__ == "__main__":
    
    probe = load_data('C1C2probe.csv')
    
    grouped = label_groups(probe)
    
    gfpmice = probe[probe['Group']=='GFP']
    appmice = probe[probe['Group']=='APP']
    
    ####ZONe
    """
    gfpmice = gfpmice[['Zone (% time)   Radius:   15.0','Unnamed: 41','Unnamed: 42','Unnamed: 43']]
    appmice = appmice[['Zone (% time)   Radius:   15.0','Unnamed: 41','Unnamed: 42','Unnamed: 43']]
    #gfpmice = gfpmice.drop('m2')
    
    
    #error
    gfperr = gfpmice.apply(scipy.stats.sem, axis=0)
    gfperr.name = "GFP"
    apperr = appmice.apply(scipy.stats.sem, axis=0)
    apperr.name = "APP"
    err = pd.concat([gfperr,apperr], axis=1)
    err.index = ['1','T','3','4']
    #mean
    avggfpzone = gfpmice.mean(axis=0)
    avgappzone = appmice.mean(axis=0)
    avggfpzone.name = "GFP"
    avgappzone.name = "APP"
    zonestats = pd.concat([avggfpzone, avgappzone], axis=1)
    zonestats.index = ['1','T','3','4']
    
    zoneplot(zonestats,err)
    
    #mean
    avgothergfp = avggfpzone[['Zone (% time)   Radius:   15.0','Unnamed: 42','Unnamed: 43']].mean(axis=0)
    avgotherapp = avgappzone[['Zone (% time)   Radius:   15.0','Unnamed: 42','Unnamed: 43']].mean(axis=0)
    otherstats = zonestats[1:3]
    otherstats['GFP'][1] = avgothergfp
    otherstats['APP'][1] = avgotherapp
    otherstats.index = ['T','O']
    #error
    othererr = err[1:3]
    othererr[1:2] = [err.drop('T').mean(axis=0)]
    othererr.index = ['T','O']
    
    zoneplot(otherstats,othererr)
    """
   
    
    
    
    gfpmice = gfpmice[['Quadrant time (%)','Unnamed: 29','Unnamed: 30','Unnamed: 31']].astype(float)
    appmice = appmice[['Quadrant time (%)','Unnamed: 29','Unnamed: 30','Unnamed: 31']].astype(float)
    

        
    
    #error
    gfperr = gfpmice.apply(scipy.stats.sem, axis=0)
    gfperr.name = "GFP"
    apperr = appmice.apply(scipy.stats.sem, axis=0)
    apperr.name = "APP"
    err = pd.concat([gfperr,apperr], axis=1)
    err.index = ['R','T','L','O']
    #mean
    avggfpzone = gfpmice.mean(axis=0)
    avgappzone = appmice.mean(axis=0)
    avggfpzone.name = "GFP"
    avgappzone.name = "APP"
    zonestats = pd.concat([avggfpzone, avgappzone], axis=1)
    zonestats.index = ['R','T','L','O']
    
    quadplot(zonestats,err)



    
    

4