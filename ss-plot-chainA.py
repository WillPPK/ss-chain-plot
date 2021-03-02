import matplotlib as mpl
mpl.use('Agg')
import pylab
import pandas as pd
import numpy as np
import gromacs.formats
import sys
import numpy as np
from matplotlib import rc, rcParams
from matplotlib.pyplot import *
from gromacs.fileformats import *
from matplotlib import gridspec
import MDAnalysis as mda

mpl.rcParams.update({'font.size': 18})
rc('axes', linewidth=2)

U = mda.Universe(sys.argv[1])
A = U.select_atoms("(resname MCYS or protein) and segid A and name CA")
P = U.select_atoms("(resname MCYS or protein) and segid A")
All = U.atoms

#residue1=1 # edit here

pylab.figure(figsize=(10,10))
gs = gridspec.GridSpec(4, 1)
pylab.title('')
subplots_adjust(wspace=0,hspace=0)

files = []
         
for x in range(1, 4):
                files.append('MD_%d/ss_chain_a.xpm' % (x))

data = ([(f,gromacs.formats.XPM(f, reverse=True, autoconvert=False)) for f in files])

values = (np.concatenate([xvg.array for xvg in list(zip(*data))[1]]))

values = pd.DataFrame(values)

values = values.T

AHelix = (values == "A-Helix").sum(axis=1)
BHelix = (values == "3-Helix").sum(axis=1)
CHelix = (values == "5-Helix").sum(axis=1)

BBridge=(values == "B-Bridge").sum(axis=1)
BSheet=(values == "B-Sheet").sum(axis=1)

helix = pd.concat([AHelix,BHelix,CHelix],axis=1)
strand = pd.concat([BSheet,BBridge],axis=1)

alphanorm = np.true_divide(helix.sum(axis=1),len(values.T))
betanorm = np.true_divide(strand.sum(axis=1),len(values.T))

ax = subplot(gs[0])
v=[A.resnums[0], A.resnums[-1], 0, 1.01]
pylab.axis(v)

alphanorm.plot.bar(color='red',linewidth=0)
betanorm.plot.bar(color='navy',linewidth=0)

pylab.ylabel( "SS Propensity", fontsize='18' )

ax.xaxis.set_major_formatter( NullFormatter() )
ax.yaxis.set_ticks(np.arange(0, 1.1, 0.5))

ax = subplot(gs[1:,0])
v=[A.resnums[0], A.resnums[-1], 0, 8.95]
pylab.axis(v)

files = []

for x in range(1, 4):
                files.append('MD_%d/rmsf.xvg' % (x))

data = ([(f,gromacs.formats.XVG(f)) for f in files])

ABC = (np.array([xvg.array[1] for xvg in list(zip(*data))[1]]).T)*10

ABCmean = np.mean(ABC,axis=1)
ABCstd = np.std(ABC,axis=1)
ABCmax = np.amax(ABC,axis=1)
ABCmin = np.amin(ABC,axis=1)

pylab.fill_between(A.resnums,ABCmean-ABCstd,ABCstd+ABCmean, color="grey", alpha=0.5)
pylab.plot(A.resnums,ABCmean, linestyle="-",lw=2,color='black')

pylab.ylabel( "RMSF ($\AA$)", fontsize='18' )
pylab.xlabel( "Residue Number", fontsize='18' )

pylab.savefig("ss-rmsf-A.pdf", format='pdf')
pylab.clf()

