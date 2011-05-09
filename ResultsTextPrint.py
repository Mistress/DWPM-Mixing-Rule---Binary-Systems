#!/usr/bin/env python

from scipy import *
from tables import *
from pickle import load
from os import mkdir, path, remove, listdir
from pylab import *

ResultsDir = 'Results'
MixtureDataDir = 'Data/Mixtures'
PlotPredX =dict()

for model in ('DWPM', 'NRTL', 'UNIQUAC'):
        
    if path.exists('Results/TextFiles/'):
        fileList = listdir('Results/TextFiles')
        for fn in fileList: 
            remove(path.join('Results/TextFiles/', fn))
    else:
        mkdir('Results/TextFiles/')
      
for Mixture in listdir(MixtureDataDir):
    OutTextFile = open('Results/TextFiles/'+Mixture[0:-3]+'.tex', 'w')
  
    
    AllModelParams = dict()
    for model in ('DWPM', 'NRTL', 'UNIQUAC'):
        h5file = openFile(ResultsDir+'/'+ Mixture[0:-3]+'/'+model+'/IndividualT/'+ Mixture, 'r')
        table = h5file.root.Outputs
        PlotPredX[model] = array([row['Predicted'] for row in h5file.root.Outputs.iterrows()])
        AllModelParams['T'] = [row['T'] for row in table.iterrows()]
        AllModelParams[model]= [row['ModelParams'] for row in table.iterrows()]
        if model == 'DWPM':
            AllModelParams['DWPMPureCompParams'] = [row['PureCompParams'] for row in table.iterrows()]
        h5file.close()

    for line in arange(size(AllModelParams['T'])):
        OutTextFile.write(r'%.2f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f  ' %(AllModelParams['T'][line], AllModelParams['DWPM'][line][0], AllModelParams['DWPM'][line][1],AllModelParams['DWPM'][line][2],AllModelParams['DWPM'][line][3],AllModelParams['DWPMPureCompParams'][line][0],AllModelParams['DWPMPureCompParams'][line][1]))
        OutTextFile.write('\\')
        OutTextFile.write('\\')
        OutTextFile.write('\n')
        
  
    OutTextFile.flush
    OutTextFile.close()
            
           

    h5file = openFile(MixtureDataDir+'/'+ Mixture, 'r')
    Compounds = h5file.root.Compounds.read()
    PlotExpX = h5file.root.ExperimentalData.TielineData.ExpComp.read()
    PlotT = reshape(h5file.root.ExperimentalData.TielineData.T.read(),(-1, 1))
    print PlotT
    print PlotExpX
    h5file.close()
    matplotlib.rc('text', usetex = True)
    fig = matplotlib.pyplot.figure()
    l1, l1a = matplotlib.pyplot.plot(PlotPredX['DWPM'],PlotT,'k-')
    l2, l2a = matplotlib.pyplot.plot(PlotPredX['NRTL'],PlotT , 'k:')
    l3, l3a = matplotlib.pyplot.plot(PlotPredX['UNIQUAC'],PlotT, 'k--')
    l4  = matplotlib.pyplot.plot(PlotExpX[0], PlotT, 'ko')
    matplotlib.pyplot.plot(PlotExpX[1], PlotT, 'ko')
    matplotlib.pyplot.xlabel(r'Mole Fraction of '+Compounds[0].capitalize(),fontsize = 14)
    matplotlib.pyplot.ylabel(r'Temperature', fontsize = 14)
    leg = matplotlib.pyplot.legend((l1, l2,l3, l4,), ('DWPM','NRTL','UNIQUAC', 'Experimental'), loc=0)
    for t in leg.get_texts():
        t.set_fontsize('small')
    savefig('Results/TextFiles/'+Mixture[0:-3]+'PhaseDiagram.ps')
    matplotlib.pyplot.close()
