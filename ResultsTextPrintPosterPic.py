#!/usr/bin/env python

from scipy import *
from tables import *
from pickle import load
from os import mkdir, path, remove, listdir
from pylab import *

ResultsDir = 'Results'
MixtureDataDir = 'Data/Mixtures'
PlotPredX =dict()
Mixtures = ('1-butanol-water.h5','aniline-water.h5', 'methanenitro-cyclohexane.h5','methanol-hexane.h5')

matplotlib.rc('text', usetex = True)
fig = matplotlib.pyplot.figure(frameon=False)

for model in ('DWPM', 'NRTL', 'UNIQUAC'):
        
    if path.exists('Results/TextFiles/'):
        fileList = listdir('Results/TextFiles')
        for fn in fileList: 
            remove(path.join('Results/TextFiles/', fn))
    else:
        mkdir('Results/TextFiles/')
      
for Mixture in Mixtures :
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
   
    index = '2'.join(('2', str(int(reshape(where(array(Mixtures)==Mixture), -1))+1)))
    ax = subplot(int(index))
    ax.yaxis.set_major_locator(MaxNLocator(6))
    l1, l1a = matplotlib.pyplot.plot(PlotPredX['DWPM'],PlotT,'k-',linewidth=2)
    l2, l2a = matplotlib.pyplot.plot(PlotPredX['NRTL'],PlotT , 'b--', linewidth=2)
    l3, l3a = matplotlib.pyplot.plot(PlotPredX['UNIQUAC'],PlotT, 'r--', linewidth=2)
    l4  = matplotlib.pyplot.plot(PlotExpX[0], PlotT, 'k.')
    matplotlib.pyplot.plot(PlotExpX[1], PlotT, 'k.')
    if Mixture == '1-butanol-water.h5':
        matplotlib.pyplot.title(r'\textbf{1-Butanol'+'-'+Compounds[1].capitalize()+'}')
    else:
        matplotlib.pyplot.title(r'\textbf{'+Compounds[0].capitalize()+'-'+Compounds[1].capitalize()+'}')
    matplotlib.pyplot.xlabel(r'Mole Fraction of '+Compounds[0].capitalize())
    matplotlib.pyplot.ylabel(r'Temperature (K)',fontsize = 14)
    #matplotlib.pyplot.xlim([0, 1])
    matplotlib.pyplot.subplots_adjust(hspace=0.35,wspace=0.3)

    
   

matplotlib.pyplot.subplots_adjust(bottom=0.15)
leg = fig.legend((l1, l2,l3, l4,), (r'Wilson3',r'NRTL',r'UNIQUAC', r'Experimental'), loc=8, mode="expand", ncol =4)
leg.draw_frame(False)
for t in leg.get_texts():
        t.set_fontsize('small')

savefig('Results/TextFiles/'+'Figure.pdf')
matplotlib.pyplot.close()
