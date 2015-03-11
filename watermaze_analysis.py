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
import sys
import matplotlib.patches as mpatches

pd.set_option('display.precision',5)

def load_data(filename):
    """
    loads data from csv file, cleaning up rows and columns that include only nan entries
    
    Input: filename    
    Output: dataframe (df) object of data
    """
    if filename == 'trainingdata.csv':
        probe = pd.read_csv(filename)
        
    else:
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
    
def zone(frame):
    """
    Creates 2 plots showing 1.) average time in all zones and 2.) target zone vs average time spent in all other zones
    
    Input: probe data
    Output: zone_all and zone_other plot 
    """
    gfpmice = frame[frame['Group']=='GFP']
    appmice = frame[frame['Group']=='APP']
    gfpmice = gfpmice[['Zone (% time)   Radius:   20.0','Unnamed: 41','Unnamed: 42','Unnamed: 43']]
    appmice = appmice[['Zone (% time)   Radius:   20.0','Unnamed: 41','Unnamed: 42','Unnamed: 43']]
    
    #mean
    avggfpzone = gfpmice.mean(axis=0)
    avgappzone = appmice.mean(axis=0)
    avggfpzone.name = "GFP"
    avgappzone.name = "APP"
    zone_all = pd.concat([avggfpzone, avgappzone], axis=1)
    zone_all.index = ['1','T','3','4']
    
    #error
    gfperr = gfpmice.apply(scipy.stats.sem, axis=0)
    gfperr.name = "GFP"
    apperr = appmice.apply(scipy.stats.sem, axis=0)
    apperr.name = "APP"
    err = pd.concat([gfperr,apperr], axis=1)
    err.index = ['1','T','3','4']

    zoneplot(zone_all,err)
    
    #mean
    avgothergfp = avggfpzone[['Zone (% time)   Radius:   20.0','Unnamed: 42','Unnamed: 43']].mean(axis=0)
    avgotherapp = avgappzone[['Zone (% time)   Radius:   20.0','Unnamed: 42','Unnamed: 43']].mean(axis=0)
    zone_other = pd.concat([avggfpzone[1:3],avgappzone[1:3]],axis=0)
    zone_other[1] = avgothergfp
    zone_other[3] = avgotherapp
    zone_other.index = ['T','O','T','O']
    
    #error
    avggfperr = gfperr[['Zone (% time)   Radius:   20.0','Unnamed: 42','Unnamed: 43']].mean(axis=0)
    avgapperr = apperr[['Zone (% time)   Radius:   20.0','Unnamed: 42','Unnamed: 43']].mean(axis=0)
    othererr = pd.concat([gfperr[1:3],apperr[1:3]], axis=0)
    othererr[1] = avggfperr
    othererr[3] = avgapperr
    othererr.index = ['T','O','T','O']
    
    otherplot(zone_other,othererr)

def quadrant(frame):
    """
    Creates 1 plot showing average time in all quadrants
    
    Input: probe data
    Output: quadrant plot
    """
    gfpmice = frame[frame['Group']=='GFP']
    appmice = frame[frame['Group']=='APP']
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
    plt.show()
    
def otherplot(frame, err):
    """
    produces a bar graph where T and O indices have distinct columns for all groups (ie. GFP, APP)
    """
    frame.plot(kind='bar', yerr=err, color = ['blue','blue','green','green'])
    blue_patch = mpatches.Patch(color='blue', label='GFP')
    green_patch = mpatches.Patch(color='green', label='APP')
    plt.legend(handles=[blue_patch,green_patch])
    plt.title('Time spent in Zone')
    plt.ylabel('Time (sec)')
    plt.xlabel('Zone')
    locs, labels = plt.xticks()
    plt.setp(labels,rotation=0)
    plt.show()
    
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
4. write advanced exclude mouse function
2. Finish up training data graphing function
3. Alter zone_other graphing so that bars are grouped together
"""
##########################  

if __name__ == "__main__":
    
    probe = load_data('probe_1-4_20cm.csv')
    """
    trainingdata = load_data('trainingdata.csv')
    trainingdata = trainingdata[['Animal','Date','Trial duration', 'Distance travelled (cm)', 'Average speed','% time near walls']]
    grouped = trainingdata.groupby(['Animal','Date']).mean()
    grouped = grouped.reset_index()
    appmice = ['m3','m4','m7','m8']
    gfpmice = ['m1','m2','m5','m6']
    app = grouped[grouped['Animal'].isin(appmice)]
    app.name=('APP')
    gfp = grouped[grouped['Animal'].isin(gfpmice)]
    gfp.name=('GFP')
    app = app.groupby('Date').mean()
    app.index = ['1','2','3']
    gfp = gfp.groupby('Date').mean()
    gfp.index = ['1','2','3']
    
    grouperr = trainingdata.groupby(['Animal','Date']).sem()
    grouperr = grouperr.reset_index()
    apperr = grouperr[grouperr['Animal'].isin(appmice)]
    apperr.name=('APP')
    gfperr = grouperr[grouperr['Animal'].isin(gfpmice)]
    gfperr.name=('GFP')
    apperr = apperr.groupby('Date').mean()
    gfperr = gfperr.groupby('Date').mean()
    apperr.index = ['1','2','3']
    gfperr.index = ['1','2','3']
    
    duration = pd.concat([gfp['Trial duration'],app['Trial duration']],axis=1)
    duration.columns = ['GFP','APP']
    durationerr = pd.concat([gfperr['Trial duration'],apperr['Trial duration']],axis=1)
    durationerr.columns = ['GFP','APP']
    duration.plot(yerr=durationerr)
    plt.title('Trial Duration')
    plt.ylabel('Time (sec)')
    plt.xlabel('Day')
    plt.show()    
    
    distance = pd.concat([gfp['Distance travelled (cm)'],app['Distance travelled (cm)']],axis=1)
    distance.columns = ['GFP','APP']
    distanceerr = pd.concat([gfperr['Distance travelled (cm)'],apperr['Distance travelled (cm)']],axis=1)
    distanceerr.columns = ['GFP','APP']
    distance.plot(yerr=distanceerr)
    plt.title('Distance travelled')
    plt.ylabel('Distance (cm)')
    plt.xlabel('Day')
    plt.show()   
    
    speed = pd.concat([gfp['Average speed'],app['Average speed']],axis=1)
    speed.columns = ['GFP','APP']
    speederr = pd.concat([gfperr['Average speed'],apperr['Average speed']],axis=1)
    speederr.columns = ['GFP','APP']
    speed.plot(yerr=speederr)
    plt.title('Average Swim Speed')
    plt.ylabel('Speed (cm/sec)')
    plt.xlabel('Day')
    plt.show()   
    
    thigmo = pd.concat([gfp['% time near walls'],app['% time near walls']],axis=1)
    thigmo.columns = ['GFP','APP']
    thigmoerr = pd.concat([gfperr['% time near walls'],apperr['% time near walls']],axis=1)
    thigmoerr.columns = ['GFP','APP']
    thigmo.plot(yerr=thigmoerr)
    plt.title('Thigmotaxis')
    plt.ylabel('Time (sec)')
    plt.xlabel('Day')
    plt.show()   
    """
    
    label_groups(probe)
    probe = probe.drop('m3')
    probe = probe.drop('m7')
    probe = probe.drop('m12')
    probe = probe.drop('m15')
    probe = probe.drop('m16')

    zone(probe)
    t = scipy.stats.ttest_ind(gfpmice['Unnamed: 41'],appmice['Unnamed: 41'])
    a = scipy.stats.f_oneway(gfpmice['Unnamed: 41'],appmice['Unnamed: 41'])
    """
    gfpmice = probe[probe['Group']=='GFP']
    appmice = probe[probe['Group']=='APP']
   
    prox = pd.DataFrame()
    prox['GFP'] = [gfpmice['Average Proximity'].astype(float).mean()]
    prox['APP'] = [appmice['Average Proximity'].astype(float).mean()]
    proxerr = pd.DataFrame()
    proxerr['GFP'] = [gfpmice['Average Proximity'].astype(float).sem()]
    proxerr['APP'] = [appmice['Average Proximity'].astype(float).sem()]
    prox.plot(kind='bar', yerr=proxerr)
    plt.title('Average Proximity')
    plt.ylabel('Proximity (cm)')
    plt.xlabel('')
    locs, labels = plt.xticks()
    plt.setp(labels,rotation=0)
    plt.show()
    """
    
    #quadrant(probe)
 
  
