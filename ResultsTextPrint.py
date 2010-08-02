#!/usr/bin/env python

from scipy import *
from tables import *
from pickle import load
from os import mkdir, path, remove, listdir
from pylab import *
import GibbsClasses

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
            AllModelParams['PureCompParams'] = [row['PureCompParams'] for row in table.iterrows()]
        h5file.close()

    for line in arange(size(AllModelParams['T'])):
        DWPM0 = (AllModelParams['PureCompParams'][line][0]/AllModelParams['DWPM'][line][0])**0.5
        DWPM1 = (AllModelParams['PureCompParams'][line][1]/AllModelParams['DWPM'][line][1])**0.5
        OutTextFile.write(r'%.2f & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f ' %(AllModelParams['T'][line],DWPM0 ,DWPM1, AllModelParams['NRTL'][line][0], AllModelParams['NRTL'][line][1], AllModelParams['UNIQUAC'][line][0], AllModelParams['UNIQUAC'][line][1]))
        OutTextFile.write('\\')
        OutTextFile.write('\\')
        OutTextFile.write('\n')
        
  
    OutTextFile.flush
    OutTextFile.close()
            
           

    h5file = openFile(MixtureDataDir+'/'+ Mixture, 'r')
    Compounds = h5file.root.Compounds.read()
    PlotExpX = h5file.root.ExperimentalData.TielineData.ExpComp.read()
    PlotT = reshape(h5file.root.ExperimentalData.TielineData.T.read(),(-1, 1))
    h5file.close()
    
    matplotlib.rc('text', usetex = True)
    
    fig = matplotlib.pyplot.figure()
    l1, l1a = matplotlib.pyplot.plot(PlotPredX['DWPM'],PlotT,'k-')
    l2, l2a = matplotlib.pyplot.plot(PlotPredX['NRTL'],PlotT , 'k:')
    l3, l3a = matplotlib.pyplot.plot(PlotPredX['UNIQUAC'],PlotT, 'k--')
    l4  = matplotlib.pyplot.plot(PlotExpX[0], PlotT, 'ko')
    matplotlib.pyplot.plot(PlotExpX[1], PlotT, 'ko')
    matplotlib.pyplot.xlabel(r'Mole Fraction of '+Compounds[0],fontsize = 14)
    matplotlib.pyplot.ylabel(r'Temperature', fontsize = 14)
    leg = matplotlib.pyplot.legend((l1, l2,l3, l4,), ('Wilson3','NRTL','UNIQUAC', 'Experimental'), loc=0)
    for t in leg.get_texts():
        t.set_fontsize('small')
    savefig('Results/TextFiles/'+Mixture[0:-3]+'PhaseDiagram.ps')
    matplotlib.pyplot.close()

    
    PlotT = reshape(PlotT, -1)
    if size(PlotT)<= 4:
        fig = matplotlib.pyplot.figure(frameon=False,figsize=(5, 4 ))
        ax = fig.add_subplot(111)
       
        matplotlib.pyplot.xlabel(r'Mole Fraction of '+Compounds[0].capitalize())
        matplotlib.pyplot.ylabel(r'$\displaystyle\frac{\Delta G_{mix}}{RT}$')
        matplotlib.pyplot.title(r'\textbf{'+Compounds[0].capitalize()+'-'+Compounds[1].capitalize()+'}')
        matplotlib.pyplot
        for T in PlotT:
            Actual = reshape( (PlotExpX[0][where(PlotT==T)], PlotExpX[1][where(PlotT==T)]), -1) 
            
            TangentComps = reshape(PlotPredX['DWPM'][where(PlotT==T)], -1)

            ModelInstance = GibbsClasses.DWPM((AllModelParams['DWPM'][ int(reshape(where(PlotT==T),-1))][0], AllModelParams['DWPM'][ int(reshape(where(PlotT==T),-1))][1],0.5))
            c = (AllModelParams['PureCompParams'][ int(reshape(where(PlotT==T),-1))][0], AllModelParams['PureCompParams'][ int(reshape(where(PlotT==T),-1))][1])

            print c
            M = dict()
   
            Gmix = array([ModelInstance.deltaGmix(x, T, c, M)for x in arange(0.001, 1, 0.001)])
            TangentGibbs = [ModelInstance.deltaGmix(x, T, c, M)for x in TangentComps]
            ActualPlot = array([ModelInstance.deltaGmix(x, T, c, M) for x in Actual])
            ax.plot(arange(0.001, 1, 0.001), Gmix, 'b-',TangentComps, TangentGibbs, 'r-', Actual, ActualPlot, 'k.',linewidth=2)          
           
        ax.yaxis.set_ticks_position('right')
        savefig('Results/TextFiles/'+Mixture[0:-3]+'Equilibrium.pdf')
        

