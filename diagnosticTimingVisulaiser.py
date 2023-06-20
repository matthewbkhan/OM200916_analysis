import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
from scipy.io import netcdf_file
import matplotlib as mpl
from pyhdf.SD import *
import pandas as pd
import numpy as np
import sys,csv,math

saving = False

shotAvailable = [98246,98248,98249,98250,98251,98252,98253,98254,98257,98262,98263,98264]

shotNums      = shotAvailable
sydorDict     = {98246:[0.0,0.77,2.06,2.28],98248:[0.0,1.09,2.12,2.34],98249:[0.0,1.42,2.17,2.39],98250:[0.0,1.73,2.33,2.55],98251:[0.0,0.75,2.06,2.28],98252:[0.0,1.39,2.27,2.49],98253:[0.0,0.75,2.16,2.38],98254:[0.0,0.75,2.16,2.38],98257:[0.0,1.07,2.22,2.44],98262:[0.0,1.39,2.27,2.49],98263:[0.0,1.07,2.22,2.44],98264:[0.0,1.07,2.22,2.44]}
mmiT0Dict     = {98246:2.600, 98248:2.610, 98249:2.620, 98250:2.580, 98251:2.600, 98252:2.630, 98253:2.600, 98254:2.600, 98257:2.620, 98262:2.630, 98263:2.620, 98264:2.620}
kbfT1Dict     = {98246:2.650, 98248:2.660, 98249:2.760, 98250:2.630, 98251:2.650, 98252:2.680, 98253:2.650, 98254:2.650, 98257:2.740, 98262:2.750, 98263:2.810, 98264:2.670}
trxiTimesDict = {98246:2.773, 98248:2.806, 98249:2.822, 98250:2.783, 98251:2.801, 98252:2.833, 98253:2.801, 98254:2.801, 98257:2.821, 98262:2.833, 98263:2.810, 98264:2.821}
laserNameDict = {98246:[1,79370],98248:[2,79371],98249:[3,79372],98250:[4,79373],98251:[16,79385],98252:[18,79387],98253:[1,79370],98254:[1,79370],98257:[2,79371],98262:[3,79372],98263:[2,79371],98264:[2,79371]}

# shotNums = [98246,98248,98249]

for shotNum in shotNums:
    shotNumStr = str(shotNum)
    print("\n##########")
    print(shotNumStr)
    print("##########")

    """ Begin the plotting """
    # plt.figure(figNum,figsize=(16,9))
    mpl.rcParams['font.size']=22
    fig, ax1 = plt.subplots(figsize=(14,9))
    # plt.title(shotNumStr)
    ax2 = ax1.twinx()
    ax2.set_ylabel('Normalised Neutron Rate')
    ax1.set_zorder(ax2.get_zorder()+1)
    ax1.set_frame_on(False)
    ax1.set_xlim([-100,4100])
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Power (TW/beam)')
    lns = []

    """ Load in Laser Profile """
    laserProfiles  = pd.read_csv("./shotData/OM200916_shotsheet-pulse_shape_data_200914b_allArDD_200914.csv")
    laserName      = laserNameDict[shotNum]
    requestedPower = laserProfiles["power_"+str(laserName[0])+"_SRF"+str(laserName[1])]
    requestedTime  = laserProfiles["times_"+str(laserName[0])+"_SRF"+str(laserName[1])]

    """ Interpolate the requested time, up to index 16 (not sure why 16?) """
    times       = np.linspace(0.0,0.15,num=150)
    interp_func = interpolate.interp1d(requestedTime[:16], requestedPower[:16], kind='linear')
    interp      = interp_func(times)
    p0_index    = (np.abs(0.02*np.max(requestedPower)-interp)).argmin()
    t0          = times[p0_index]

    """ Load in the actual laser powers """
    hdfFile     = "./shotData/P510/P510-p510_data_"+shotNumStr+".hdf"
    file        = SD(hdfFile, SDC.READ)
    time_obj    = file.select(10)
    masterTime  = time_obj.get()
    pow_obj     = file.select(2)
    actualPower = pow_obj.get()

    """ Plot the lasers """
    lns1 = ax1.plot([(i-t0)*1e3 for i in requestedTime],[i*(1./60.) for i in requestedPower],label='Requested Power',linewidth=3,color='tab:orange')
    lns2 = ax1.plot(masterTime,                         actualPower,                         label='Actual Power',   linewidth=3,color='tab:blue')
    lns += lns1
    lns += lns2

    """ KBF """
    kbfT1              = kbfT1Dict[shotNum]
    kbfTimes           = [0,100,200,300]#kbfTimesDict[shotNum]
    kbfAcquisitionTime = np.array([(kbfT1*1000.)+i for i in kbfTimes])
    kbfHeight          = np.max(actualPower)*np.ones(len(kbfAcquisitionTime))
    lns4 = ax1.plot(kbfAcquisitionTime,kbfHeight,marker='o',linestyle=' ',markersize=10,label='KBF',color='tab:red',zorder=5)
    lns += lns4

    """ TRXI """
    trxiT0             = trxiTimesDict[shotNum]*1000.
    lns5 = ax1.plot([trxiT0,trxiT0+120],0.75*np.max(actualPower)*np.ones(2),marker=' ',linewidth=7,markersize=10,label='TRXI',color='tab:olive',zorder=5)
    lns += lns5

    """ SYDOR """
    sydorAquisitionTime = sydorDict[shotNum]
    sydorAquisitionTime = 1000. * np.array(sydorAquisitionTime)
    # lns6 = ax1.plot(sydorAquisitionTime,0.5*np.max(actualPower)*np.ones(len(sydorAquisitionTime)),marker='o',linestyle=' ',markersize=10,label='SYDOR',color='tab:cyan')
    lns6 = ax1.plot([sydorAquisitionTime[0],sydorAquisitionTime[0]+215],0.5*np.max(actualPower)*np.ones(2),marker=' ',linewidth=7,label='SYDOR',color='tab:cyan')
    for t,h in zip(sydorAquisitionTime[1:],[0.515,0.53,0.545]):
        ax1.plot([t,t+215],h*np.max(actualPower)*np.ones(2),marker=' ',linewidth=7,color='tab:cyan',zorder=5)
    lns += lns6

    """ MMI """
    mmiT0              = mmiT0Dict[shotNum]
    mmiTimes           = [0.0,0.1,0.2,0.3]#mmi[1]
    mmiAcquisitionTime = np.array([(mmiT0+i)*1000. for i in mmiTimes])
    mmiHeight          = 0.25*np.ones(len(mmiAcquisitionTime))
    lns7 = ax1.plot(mmiAcquisitionTime,mmiHeight,marker='o',linestyle=' ',markersize=10,label='MMI',color='tab:orange',zorder=5)
    lns += lns7

    """ Experimental Neutrons """
    """ TODO: Remove the normalisation once the simulation neutrons are correct.
    """
    neutrons  = np.loadtxt('./shotData/cryontd/'+shotNumStr+'_cryontd.txt')
    neutTime  = neutrons[:,0]
    neutYield = neutrons[:,1]
    lns3 = ax2.plot(neutTime,neutYield/np.max(neutYield),label='Exp. Yield',linewidth=3,color='tab:green')
    lns += lns3

    """ Simulated Neutrons """
    """ TODO: I am not certain that I have got the timing correct for the simulations (the definition of t0 differs from sim to exp)
        It would also be good to compare the actual neutron values and not the normalised numbers.
    """
    fileName = 'hyChDD'
    hyPath   = "../hyadesSims/230607a_m3b,75m4,50dt1,2FlxB/230607a_m3b,75m4,50dt1,2FlxB/"
    nfile    = netcdf_file(hyPath+str(shotNum)+'/'+fileName+'.cdf', mode='r', mmap=False) # Load in the output from the simulation
    simTime  = nfile.variables['DumpTimes'] # The times that the simulation generated an output
    tn       = nfile.variables['Bpeprdr'] # Fusion energy emission rate
    tnEmR    = np.sum(tn[:],axis=1)
    lnsX     = ax2.plot((simTime[:]*1e12)-t0,tnEmR[:]/np.max(tnEmR[:]),label='Sim. Yield',linewidth=3,color='tab:purple',linestyle="--")
    lns     += lnsX

    labs = [l.get_label() for l in lns]
    ax1.legend(lns, labs, loc='upper left', prop={"size":14})
    if(saving):
        plt.savefig("./pulseShape_"+shotNumStr+".png")
    plt.show()
