#!/usr/bin/env python3.5

import h5py
import os
import time
import epics
import argparse
import pickle
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from datetime import datetime, timedelta

ngs2ip = '172.17.65.31'
ngs2user = 'aouser'
ngs2pass = 'ngs2@CerroP'
host = '"tcp://localhost:45000"'
utlPath = '/opt/ao/bin/aocmd'
# filePath = '/home/aouser/focusSeq/'
filePath = '/home/iarriagada/work/projects/ngs2/autofocus/'
# filePath = './'
fNameRoot = 'autoFcsSeq-'
numSamp = 100
maxWait = 60
arcppix = 0.229
os.environ["EPICS_CA_ADDR_LIST"] = ngs2ip
# EPICS Records definitions
w1Rec = 'ngs2:wfsCentroids.VALB'
w2Rec = 'ngs2:wfsCentroids.VALH'
w3Rec = 'ngs2:wfsCentroids.VALN'
focusDmdRec = 'ngs2:ngsFocus.TGTX'
focusPosRec = 'ngs2:ngsFocus.POSX'
focusDirRec = 'ngs2:ngsFocus.DIR'
focusModRec = 'ngs2:ngsFocus.MODE'
axisNames = ['w1FWHMave', 'w1FWHMave', 'w1FWHMave', 'Focus Pos']

def parse_args():
    '''
    This routines parses the arguments used when running this script
    '''
    parser = argparse.ArgumentParser(
        description='Use this script for NGS2 autofocus')

    parser.add_argument('-sp',
                        '--startpos',
                        dest='spos',
                        default=0.0,
                        help='Focus start position e.g.: -sf 0.5')

    parser.add_argument('-ep',
                        '--endpos',
                        dest='epos',
                        default=6.0,
                        help='Focus end position e.g.: -ef 4.2')

    parser.add_argument('-stp',
                        '--stepsize',
                        dest='step',
                        default=0.5,
                        help='Step size between each focus position\
                        e.g.: -stp 0.5')

    parser.add_argument('-pkl',
                        '--pickle',
                        dest='pklFile',
                        default='',
                        help='Name of the pickle file to be plotted and used\
                        e.g.: -pkl autofocus-20190919T185622.pkl')

    parser.add_argument('-naf',
                        '--noAF',
                        dest='noAF',
                        action='store_true',
                        help='Use this option to disable Auto Focus')

    args = parser.parse_args()
    return args

def onPosChange(pvname=None, value=None, char_value=None, host=None, **kws):
    global posFlag # Global in position variable
    global focus # Focus stage position demand
    # Check if current pos is within tolerance from target position only if
    # position hasn't been reached
    if (abs(focus-value) < 0.001) and not(posFlag):
        print("Focus position {0} reached".format(str(value)))
        posFlag = True

def captureFWHM(focusPos, focusDmd, focusMod, focusDir,
                spos, epos, step, fname):
    global posFlag
    global focus
    wFWHMAvgList = []
    # Start connection with epics records
    w1fwhm = epics.PV(w1Rec)
    w2fwhm = epics.PV(w2Rec)
    w3fwhm = epics.PV(w3Rec)
    # Handles the case in which the starting position is greater than end
    # position
    if spos > epos:
        posList = list(np.arange(float(epos),
                                 float(spos)+float(step),
                                 float(step)))
    else:
        posList = list(np.arange(float(spos),
                                 float(epos)+float(step),
                                 float(step)))
    print('This is the range of positions to move:')
    print(posList)
    with open(fname, 'ab') as f:
        # write channels names at the top of pickle file
        pickle.dump(axisNames, f)
    focusMod.put('MOVE')
    for focus in posList:
        wFWHM = []
        # Sets focus demand, then executes focus command, then wait until
        # focus demand is reached
        posFlag = False
        focusDmd.put(focus)
        print('Position demand is {0}'.format(focus))
        focusDir.put('START')
        while not(posFlag):
            continue
        # Captures the FWHM value for each probe and stores it into array
        sample = 0
        while sample < (numSamp + 1):
            wFWHM.append([w1fwhm.value/arcppix,
                          w2fwhm.value/arcppix,
                          w3fwhm.value/arcppix,
                          focus])
            sample += 1
            time.sleep(0.2)
        # Transpose matrix, calculate the average of each column then append
        # results to final list
        wFWHMMat = np.array(wFWHM).T
        wFWHMAvg = list(np.average(wFWHMMat, axis=1))
        with open(fname, 'ab') as f:
            # write channels names at the top of pickle file
            pickle.dump(wFWHMAvg, f)
        wFWHMAvgList.append(wFWHMAvg)
    fwhmAvgData = np.array(wFWHMAvgList).T
    print(np.array(wFWHMAvgList).T)
    print(fwhmAvgData.min(axis=1))
    return fwhmAvgData

def analyzeFile(fName):
    print('Reading data from ' + fName)
    fwhmData = []
    with open(fName, 'rb') as f:
        while(True):
            try:
                fwhmData.append(pickle.load(f))
            except (EOFError):
                break
    # print("sample of FWHM array:")
    # print(fwhmData[0])
    # print(fwhmData[1])
    # print("size of pickle: " + str(len(fwhmData)))
    return fwhmData

if __name__ == '__main__':
    # Capture command line arguments
    args = parse_args()
    # Define lists and global variables
    focus = 0
    posFlag = True
    startTime = datetime.now() # starting time of the capture
    startDateStr = datetime.strftime(startTime, '%Y%m%dT%H%M%S')
    fileName = 'autofocus-'+startDateStr+'.pkl' # define file name
    fDmd = epics.PV(focusDmdRec)
    fDir = epics.PV(focusDirRec)
    fMod = epics.PV(focusModRec)
    fPos = epics.PV(focusPosRec, callback=onPosChange)
    if args.pklFile == '':
        fwhmAverage = captureFWHM(fPos, fDmd, fMod, fDir, args.spos,
                                  args.epos, args.step, fileName)
    else:
        fwhmFromFile = analyzeFile(args.pklFile)
        fwhmAverageList = fwhmFromFile[1:]
        fwhmAverage = np.array(fwhmAverageList).T
    pOrder = 10
    cW1 = np.polyfit(fwhmAverage[3], fwhmAverage[0], pOrder)
    cW2 = np.polyfit(fwhmAverage[3], fwhmAverage[1], pOrder)
    cW3 = np.polyfit(fwhmAverage[3], fwhmAverage[2], pOrder)
    # print(cW1)
    # print(cW1[0] + cW1[1]*(fwhmAverage[3][0]**1) + cW1[2]*(fwhmAverage[3][0]**2) + cW1[3]*(fwhmAverage[3][0]**3) + cW1[4]*(fwhmAverage[3][0]**4))
    # pnW1 = np.poly1d(cW1)
    # print (pnW1(fwhmAverage[3][0]))
    # print(fwhmAverage[3][0])
    # print(fwhmAverage[0][0])



    pfW1 = [np.sum([cW1[pOrder-j] * x**j for j in list(range(pOrder+1))])
                                                for x in fwhmAverage[3]]
    pfW2 = [np.sum([cW2[pOrder-j] * x**j for j in list(range(pOrder+1))])
                                                for x in fwhmAverage[3]]
    pfW3 = [np.sum([cW3[pOrder-j] * x**j for j in list(range(pOrder+1))])
                                                for x in fwhmAverage[3]]
    minIndW1 = list(pfW1).index(np.amin(pfW1))
    minFW1 = fwhmAverage[3][minIndW1]
    labelW1 = 'OFP = ' + str(minFW1)
    minIndW2 = list(pfW2).index(np.amin(pfW2))
    minFW2 = fwhmAverage[3][minIndW2]
    labelW2 = 'OFP = ' + str(minFW2)
    minIndW3 = list(pfW3).index(np.amin(pfW3))
    minFW3 = fwhmAverage[3][minIndW3]
    labelW3 = 'OFP = ' + str(minFW3)

    optFocus = np.average([minFW1, minFW2, minFW3])
    print('Optimal Focus Position (OFP) is {0} [mm]'.format(optFocus))
    print('Moving focus to position')
    fMod.put('MOVE')
    posFlag = False
    focus = optFocus
    fDmd.put(optFocus)
    fDir.put('START')
    while not(posFlag):
        continue
    plt.plot(fwhmAverage[3], fwhmAverage[0], label='W1',
                 color='#0000FF', linestyle = '', marker='x')
    plt.plot(fwhmAverage[3], pfW1, label='Polyfit W1',
                 color='#00AAFF', linestyle = '--', marker='')
    plt.axvline(minFW1, label=labelW1, color='#0000FF')
    plt.plot(fwhmAverage[3], fwhmAverage[1], label='W2',
                 color='#00FF00', linestyle = '', marker='x')
    plt.plot(fwhmAverage[3], pfW2, label='Polyfit W2',
                 color='#00FFAA', linestyle = '--', marker='')
    plt.axvline(minFW2, label=labelW2, color='#00FF00')
    plt.plot(fwhmAverage[3], fwhmAverage[2], label='W3',
                 color='#FF0000', linestyle = '', marker='x')
    plt.plot(fwhmAverage[3], pfW3, label='Polyfit W3',
                 color='#FF00AA', linestyle = '--', marker='')
    plt.axvline(minFW3, label=labelW3, color='#FF0000')
    plt.ylabel('FWHM [pix]')
    plt.xlabel('Focus Pos [mm]')
    plt.title('Focus Optimization')
    plt.legend()
    plt.show()

