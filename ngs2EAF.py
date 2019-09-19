#!/usr/bin/env python3.5

import h5py
import os
import time
import epics
import argparse
import pickle
import numpy as np
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

    parser.add_argument('-naf',
                        '--noAF',
                        dest='noAF',
                        action='store_true',
                        help='Use this option to disable Auto Focus')

    args = parser.parse_args()
    return args

def onPosChange(pvname=None, value=None, char_value=None, host=None, **kws):
    global posFlag
    global focus
    if (abs(focus-value) < 0.001) and not(posFlag):
        print("Focus position {0} reached".format(str(focusPos.value)))
        posFlag = True

if __name__ == '__main__':
    # Capture command line arguments
    args = parse_args()
    # Define lists and global variables
    wFWHMAvgList = []
    focus = 0
    posFlag = True
    # Start connection with epics records
    w1fwhm = epics.PV(w1Rec)
    w2fwhm = epics.PV(w2Rec)
    w3fwhm = epics.PV(w3Rec)
    focusDmd = epics.PV(focusDmdRec)
    focusPos = epics.PV(focusPosRec, callback=onPosChange)
    focusDir = epics.PV(focusDirRec)
    focusMod = epics.PV(focusModRec)
    # Handles the case in which the starting position is greater than end
    # position
    if args.spos > args.epos:
        posList = list(np.arange(float(args.epos),
                                 float(args.spos)+float(args.step),
                                 float(args.step)))
    else:
        posList = list(np.arange(float(args.spos),
                                 float(args.epos)+float(args.step),
                                 float(args.step)))
    print('This is the range of positions to move:')
    print(posList)
    focusMod.put('MOVE')
    startTime = datetime.now() # starting time of the capture
    startDateStr = datetime.strftime(startTime, '%Y%m%dT%H%M%S')
    fileName = 'autofocus-'+startDateStr+'.pkl' # define file name
    with open(fileName, 'ab') as f:
        # write channels names at the top of pickle file
        pickle.dump(axisNames, f)
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
        # waitTime = 0
        # startTime = datetime.now()
        # # Wait until file becomes available
        # while not(os.path.exists(filePath + fileName)):
            # time.sleep(0.1)
            # currTime = datetime.now()
            # waitTime = (currTime - startTime).total_seconds()
            # if waitTime > maxWait:
                # break
        # Transpose matrix, calculate the average of each column then append
        # results to final list
        wFWHMMat = np.array(wFWHM).T
        wFWHMAvg = list(np.average(wFWHMMat, axis=1))
        with open(fileName, 'ab') as f:
            # write channels names at the top of pickle file
            pickle.dump(wFWHMAvg, f)
        wFWHMAvgList.append(wFWHMAvg)
    fwhmAvgData = np.array(wFWHMAvgList).T
    print(np.array(wFWHMAvgList).T)
    print(fwhmAvgData.min(axis=1))
    fwhmData = []
    with open(fileName, 'rb') as f:
        while(True):
            try:
                fwhmData.append(pickle.load(f))
            except (EOFError):
                break
    print("sample of FWHM array:")
    print(fwhmData[0])
    print(fwhmData[1])
    print("size of pickle: " + str(len(fwhmData)))
    plt.plot(fwhmAvgData[3], fwhmAvgData[0], label='W1')
    plt.plot(fwhmAvgData[3], fwhmAvgData[1], label='W2')
    plt.plot(fwhmAvgData[3], fwhmAvgData[2], label='W3')
    plt.ylabel('FWHM [pix]')
    plt.xlabel('Focus Pos [mm]')
    plt.title('Focus Optimization')
    plt.legend()
    plt.show()

