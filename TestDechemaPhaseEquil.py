 
#!/usr/bin/env python

from MainBinaryLLE import *
import tables
from os import mkdir, path, remove, listdir
from scipy import *


Models = ('NRTL', 'UNIQUAC')
ModelInstances = (GibbsClasses.NRTL, GibbsClasses.UNIQUAC)
MixtureDataDir = 'Data/Mixtures'
PureDataDir = 'Data/PureComps'
#Bounds = [((-1000, 0), (-1000, 0), (0.5, 0.5)), ((-800, 3000), (-800, 3000)), ((-800, 3000), (-800, 3000))]
R = 8.314

for file in listdir(MixtureDataDir):
    Fit = 'Test'
    h5file = tables.openFile(MixtureDataDir+'/'+file, 'r')
    Compounds = h5file.root.Compounds.read()
    params = (h5file.root.DechemaParams.NRTL.read(), h5file.root.DechemaParams.UNIQUAC.read())
    System = Mixture(Compounds, MixtureDataDir, PureDataDir)  
    if not(path.exists('Results/'+System.Name)):
        mkdir('Results/'+System.Name)
    h5file.close()
    
    for i in arange(size(Models)):
        Model = Models[i]
        if not(path.exists('Results/'+System.Name+'/'+Model)):
            mkdir('Results/'+System.Name+'/'+Model)        
        if path.exists('Results/'+System.Name+'/'+Model+'/'+Fit):
            fileList = listdir('Results/'+System.Name+'/'+Model+'/'+Fit)
            for fn in fileList: 
                remove(path.join('Results/'+System.Name+'/'+Model+'/'+Fit, fn))
        else:
            mkdir('Results/'+System.Name+'/'+Model+'/'+Fit)
        
        for T in System.M['T']:
           
            Actual  = array([interp(T,cast['f'](System.M['T']), cast['f'](System.M['ExpComp'][0])),interp(T,cast['f'](System.M['T']), cast['f'](System.M['ExpComp'][1]))])
            CompC = System.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in Compounds]
            System.Plotter(reshape(params[i][:, where(System.M['T']==T)], 2), Models[i], Fit, ModelInstances[i](reshape(params[i][:, where(System.M['T']==T)], 2)),Actual, c, T)
 
 
