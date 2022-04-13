import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib
import sys
import os
import numba as nb
from matplotlib.ticker import FormatStrFormatter

'''
# russian letters
from matplotlib import rc
rc('font',**{'family':'sansserif'})
rc('text', usetex=True)
rc('text.latex',unicode=True)
rc('text.latex',preamble='\usepackage[utf8]{inputenc}')
rc('text.latex',preamble='\usepackage[russian]{babel}')
'''
from matplotlib import rcParams
rcParams.update({'figure.autolayout': False})

font = {'family': 'normal', 'weight': 'normal', 'size': 30}
matplotlib.rc('font', **font)
