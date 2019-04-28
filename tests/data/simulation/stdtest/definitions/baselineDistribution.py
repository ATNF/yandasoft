#!/usr/bin/env python
import numpy as np
import pylab as plt
from configs import *
import math
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--config', '-c', default='A27CR3P6B', type=str, help='Configuration to plot')
parser.add_argument('--nbins', '-n', default=30, type=int, help='Number of bins (default : %(default)s)')
parser.add_argument('--maxB', '-m', default=-1, type=int, help='Max value of baseline in plot (-1, the default, to set automatically)')
options = parser.parse_args()

config=getConfigs()

confName=options.config

antList = AntennaList(config[string.upper(confName)], 50, "south")

baselines=antList.baselineDistribution()

print baselines.size
if options.maxB > 0.:
    plt.hist(baselines,bins=options.nbins,range=(0,options.maxB))
    plt.xlim(0,options.maxB)
else:
    plt.hist(baselines,bins=options.nbins)
plt.title(confName)
plt.savefig('baselineDist.png')

