#!/usr/bin/env python3.5

import numpy as np
import h5py
import os
from datetime import datetime, timedelta
import time
import paramiko
import matplotlib.pyplot as plt

host = '"tcp://localhost:45000"'
utlPath = '/opt/ao/bin/aocmd'
# filePath = '/home/aouser/focusSeq/'
filePath = '/home/iarriagada/work/projects/ngs2/autofocus/'
# filePath = './'
fNameRoot = 'autoFcsSeq-'
numSamp = 1000
maxWait = 60

def parse_args():
    '''
    This routines parses the arguments used when running this script
    '''
    parser = argparse.ArgumentParser(
        description='Use this script for NGS2 autofocus')

    parser.add_argument('-sp',
                        '--startpos',
                        dest='sPos',
                        default=0.0,
                        help='Focus start position e.g.: -sf 0.5')

    parser.add_argument('-ef',
                        '--endpos',
                        dest='ePos',
                        default=6.0,
                        help='Focus end position e.g.: -ef 4.2')

    parser.add_argument('-stp',
                        '--steps',
                        dest='steps',
                        default=10,
                        help='Number of steps between start and end position\
                        e.g.: -stp 5')

    parser.add_argument('-naf',
                        '--noAF',
                        dest='noAF',
                        action='store_true',
                        help='Use this option to disable Auto Focus')

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    focus = 0.0
    wFWHMAvgList = []
    while (focus < 6.5):
        startDate = datetime.now() # starting time of the capture
        startDateStr = datetime.strftime(startDate, '%Y%m%d-')
        # File name definition
        fileName = \
        fNameRoot + startDateStr + str(focus).replace('.', 'p') + '.h5' # define file name
        fileName = fileName.replace('17', '18')
        # Focus command string definition
        focusCmd = '"FOCUS SETPOS ' + str(focus) + '"'
        # Sequence command string definition
        seqCmd = '"SEQUENCE ' + str(numSamp) + ' ' + filePath + fileName + '"'
        # Execute focus command, then wait until focus demand is reached
        # os.system(utlPath + ' ' + host + ' ' + focusCmd)
        # time.sleep(5) # needs to be replaced with EPICS Focus position value
        print("Focus position {0} reached".format(str(focus)))
        # Execute sequence command
        # os.system(utlPath + ' ' + host + ' ' + seqCmd)
        print("Writing " + filePath + fileName)
        waitTime = 0
        startTime = datetime.now()
        # Wait until file becomes available
        while not(os.path.exists(filePath + fileName)):
            time.sleep(0.1)
            currTime = datetime.now()
            waitTime = (currTime - startTime).total_seconds()
            if waitTime > maxWait:
                break
        # Load created file and store FWHM values on variable as a numpy array
        seqFile = h5py.File(filePath + fileName, 'r')
        wFWHM = np.array(seqFile['fwhm_focus'])
        # Transpose matrix, calculate the average of each column then append
        # results to final list
        wFWHMMat = wFWHM.T
        wFWHMAvg = list(np.average(wFWHMMat, axis=1))
        wFWHMAvgList.append(wFWHMAvg)
        focus += 0.5
    fwhmAvgData = np.array(wFWHMAvgList).T
    print(np.array(wFWHMAvgList).T)
    print(fwhmAvgData.min(axis=1))
    plt.plot(fwhmAvgData[3], fwhmAvgData[0], label='W1')
    plt.plot(fwhmAvgData[3], fwhmAvgData[1], label='W2')
    plt.plot(fwhmAvgData[3], fwhmAvgData[2], label='W3')
    plt.ylabel('FWHM [pix]')
    plt.xlabel('Focus Pos [mm]')
    plt.title('Focus Optimization')
    plt.legend()
    plt.show()

