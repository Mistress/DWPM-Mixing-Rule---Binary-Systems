#!/usr/bin/env python

from scipy import *
from pylab import *
from os import mkdir, path, remove, listdir
from mpl_toolkits.mplot3d import Axes3D
import scipy.optimize
import scipy.linalg
import scipy.interpolate 

import glob
import tables
import vanDerWaals
import ErrorClasses
import GibbsClasses
import AnalyticSystemsClasses
import AnalyticPhaseCompSystemsClasses

class Mixture:
        
    def __init__(self, Compounds, MixtureDataDir, PureDataDir):
        
        MixtureName = Compounds[0]
        self.Compounds = Compounds
        self.Data =dict((Compound, {}) for Compound in Compounds)
        for Compound in Compounds:
            h5file = tables.openFile(PureDataDir+'/'+Compound+'.h5', 'r')
            properties = h5file.root.Properties
            self.Data[Compound] = dict(((field, row[field]) for row in properties.iterrows() for field in properties.colnames))        
            h5file.close()
        self.Name = '-'.join(Compounds)
        h5file = tables.openFile(MixtureDataDir+'/'+ self.Name +'.h5', 'r')
        self.M = dict(Compounds = h5file.root.Compounds.read(), ExpComp = h5file.root.ExperimentalData.TielineData.ExpComp.read(), T =h5file.root.ExperimentalData.TielineData.T.read())
        for row in h5file.root.UNIQUACParams.iterrows():
            for field in h5file.root.UNIQUACParams.colnames:
                self.M[field] = row[field]
        self.vdWaalsInstance = vanDerWaals.Slope(self.Data, Compounds)
        [Slope, VPData] = self.vdWaalsInstance.BestFit()
        h5file.close()
        self.Binaries = zeros((size(self.M['T']), 2))

    def OverallSOpt(self, Cell_s, Model, ModelInstance, R):
        
        h5file = tables.openFile('Results/'+self.Name+'/'+ Model+'/Cells('+str(Cell_s[0])+','+str(Cell_s[1])+')/'+ self.Name +'.h5', 'r')
        DWPMParams = array([row['ModelParams'] for row in h5file.root.Outputs.iterrows()])
        Convergence = [row['Converged'] for row in h5file.root.Outputs.iterrows()]
        PureCompParams = array([row['PureCompParams'] for row in h5file.root.Outputs.iterrows()])
        Temps = array([row['T'] for row in h5file.root.Outputs.iterrows()])
        Test = [x==1 for x in Convergence]
        h5file.close()
        
        Tc = reshape(array([self.Data[Compound]['Tc'] for Compound in self.Compounds]), -1)
       
        TrTerms = array([Temps/Tc[0], Temps/Tc[1]])
               
        Constants = array([log(-1*DWPMParams[:,0])/(1-TrTerms[0,:]),log(-1*DWPMParams[:,1])/(1-TrTerms[1,:])])
                
        if all(Test):
            ScaledSV = [std(Constants[0,:])/average(Constants[0,:]), std(Constants[1,:])/average(Constants[1,:])]
            return (ScaledSV)
        else:
            return (NaN,NaN)

##========================================================================================================================================

Model = 'DWPM'
ModelInstance = GibbsClasses.DWPM
MixtureDataDir = 'Data/Mixtures'
PureDataDir = 'Data/PureComps'
R = 8.314

if not(path.exists('Results/')):
    mkdir('Results/')   

for file in listdir(MixtureDataDir):

    s = arange(0.1, 1.0, 0.1)

    h5file = tables.openFile(MixtureDataDir+'/'+file, 'r')
    Compounds = h5file.root.Compounds.read()
    h5file.close()

    Name = '-'.join(Compounds)

    Optimization = Mixture(Compounds, MixtureDataDir, PureDataDir)
    
    Z = abs(array([Optimization.OverallSOpt(array([s1, s2]), Model, ModelInstance, R) for s1 in s for s2 in s]))
    plotZ1 = transpose(reshape(Z[:,0], (size(s), -1)))
    plotZ2 = transpose(reshape(Z[:,1], (size(s), -1)))
    s_1, s_2 = meshgrid(s, s)
                           
    fig = matplotlib.pyplot.figure()
    subfigures = fig.add_subplot(1, 2, 1).get_position()
    ax = Axes3D(fig, subfigures)
    ax.plot_wireframe(s_1, s_2, plotZ1, rstride=1, cstride=1) 
    matplotlib.pyplot.title(r"Relative Variation of $m_{1}$", fontsize=14) 
    subfigures = fig.add_subplot(1, 2, 2).get_position()
    ax = Axes3D(fig, subfigures)
    ax.plot_wireframe(s_1, s_2, plotZ2, rstride=1, cstride=1) 
    matplotlib.pyplot.title(r"Relative Variation of $m_{2}$", fontsize=14) 
    fig.savefig('Results/'+Name+'/'+Model+'/StdVarOfParametersWireframe.png')
    matplotlib.pyplot.close()

    fig2 = matplotlib.pyplot.figure()
    ax = Axes3D(fig2)
    plotTotalZ = plotZ1 + plotZ2
    ax.plot_wireframe(s_1, s_2, plotTotalZ, rstride=1, cstride=1) 
    matplotlib.pyplot.title(r"Sum of Variations of $m_{1}$ and $m_{2}$ ", fontsize=14) 
    fig2.savefig('Results/'+Name+'/'+Model+'/SumStdVarOfParametersWireframe.png')
    matplotlib.pyplot.close()
    

