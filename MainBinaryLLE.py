#!/usr/bin/env python

from scipy import *
from pylab import *
from os import mkdir, path, remove, listdir
import scipy.optimize
import scipy.linalg
import glob
import tables
import vanDerWaals
import ErrorClasses
import GibbsClasses
import AnalyticSystemsClasses


class ResultsFile1(tables.IsDescription):
    T = tables.Float64Col()
    ModelParams = tables.Float64Col(shape = (4)) 
    Actual = tables.Float64Col(shape = (2))
    PureCompParams = tables.Float64Col(shape = (2))
    Converged = tables.Float64Col(shape = (1))

class ResultsFile2(tables.IsDescription):
    T = tables.Float64Col()
    ModelParams = tables.Float64Col(shape = (2)) 
    Actual = tables.Float64Col(shape = (2))
    Converged = tables.Float64Col(shape = (1))
    
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
    
    def Plotter(self, Params, m, b, Cell_s, Model, ModelInstance, Actual, c, T, Converged):
       
        deltaGibbsmix = array([ModelInstance.deltaGmix(x, T, c, self.M) for x in arange(0.001, 1, 0.001)])
        ActualPlot = array([ModelInstance.deltaGmix(x, T, c, self.M) for x in Actual])
        Tieline = [m*x + b for x in arange(0.001, 1, 0.001)]
        
        matplotlib.rc('text', usetex = True)
        fig = matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(arange(0.001, 1, 0.001), deltaGibbsmix, 'b-', arange(0.001, 1, 0.001), Tieline, 'r-', Actual, ActualPlot, 'ko')
        matplotlib.pyplot.xlabel(r'Mole Fraction of '+self.Compounds[0].capitalize(), fontsize = 14)
        matplotlib.pyplot.ylabel(r'$\displaystyle\frac{\Delta G_{mix}}{RT}$', fontsize = 14)
        matplotlib.pyplot.title(r'\textbf{Predicted $\displaystyle\frac{\Delta G_{mix}}{RT}$ %3.2f}'%T, fontsize = 14)
        ax = fig.add_axes([0,0,1,1])
        if Model == 'DWPM':
            ax.text(0,0, 'Model: %s, Params: %.3f, %.3f, %.3f, %.3f, %.1f'%(Model, Params[0], Params[1], Cell_s[0], Cell_s[1], Converged), fontsize=12, transform=ax.transAxes)
        else:
            ax.text(0,0, 'Model: %s, Params: %.3f, %.3f, %.1f '%(Model, Params[0], Params[1], Converged), fontsize=12, transform=ax.transAxes)
        ax.set_axis_off()
        #show()
        savefig('Results/'+self.Name+'/'+Model+'/T_'+str(T) +'.pdf')
        matplotlib.pyplot.close()
        
        if not(path.exists('Results/'+self.Name+'/'+ Model+'/'+self.Name +'.h5')):
            h5file = tables.openFile('Results/'+self.Name+'/'+Model+'/'+self.Name +'.h5', 'w', "Optimization Outputs")
            if Model=='DWPM':
                table = h5file.createTable("/", "Outputs", ResultsFile1, "Optimal model parameters, predicted phase equilibrium, errors etc")
            else:
                table = h5file.createTable("/", "Outputs", ResultsFile2, "Optimal model parameters, predicted phase equilibrium, errors etc")
        else:
            h5file = tables.openFile('Results/'+self.Name+'/'+Model+'/'+self.Name +'.h5', 'r+')
            table = h5file.root.Outputs

        table.row['T'] = T
        table.row['ModelParams'] = reshape(append(Params, Cell_s),-1)
        table.row['Actual'] = Actual
        if Model=='DWPM':
            table.row['PureCompParams']  = reshape(array(c),-1)
        table.row['Converged'] = Converged
        table.row.append()
        table.flush()
        h5file.close()

    def ParamCalc(self, Model, ModelInstance, InitParamsLimit, R):

        if path.exists('Results/'+self.Name+'/'+Model):
            fileList = listdir('Results/'+self.Name+'/'+Model)
            for fn in fileList: 
                remove(path.join('Results/'+self.Name+'/'+Model, fn))
        else:
            mkdir('Results/'+self.Name+'/'+Model)
            
        if Model == "DWPM":
            Cell_s = (0.3, 0.8)
            System = AnalyticSystemsClasses.SystemDWPM(Cell_s)
        elif Model == "NRTL":
            Cell_s = ()
            alpha = 0.2
            System = AnalyticSystemsClasses.SystemNRTL(alpha)
        elif Model == "UNIQUAC":
            Cell_s = ()
            z = 10
            System = AnalyticSystemsClasses.SystemUNIQUAC(self.M, z)

            
        InitParams = (InitParamsLimit[1]-InitParamsLimit[0])*random(2) + append(InitParamsLimit[0], InitParamsLimit[0])
        InitParamj = append(InitParams, array([-1., -3.]))      
                        
        for T in self.M['T']:
            Actual =  array([interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][0])),interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][1]))])
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in Compounds]
           
            MaxFEval = 1000
            TolParam = 1e-6
            Converge = 3
            NStarts = 1
            MaxNStarts = 10000

            while (Converge!=1 and NStarts<MaxNStarts):
                                  
                      [Paramj, infodict, Converge, Mesg] = scipy.optimize.fsolve(System.SystemEquations, InitParamj, (Actual, R, T), System.SystemEquationsJac, 1, 0, TolParam, MaxFEval, None, 0.0, 100, None)
                      print Mesg
                      NFCalls = infodict['nfev']
                      NJCalls = infodict['njev']
                      NStarts = NStarts + 1
                      
                      InitParams = (InitParamsLimit[1]-InitParamsLimit[0])*random(2) + append(InitParamsLimit[0], InitParamsLimit[0])
                      InitParamj = append(InitParams, array([1., -5.]))  
                      
            print NStarts
            print Paramj
            InitParamj = Paramj
                                        
            m = Paramj[2]
            b = Paramj[3]
            
            if Model == "DWPM":
                c12 = c[0]*Paramj[0]
                c21 = c[1]*Paramj[1]
            else:
                c12 = Paramj[0]
                c21 = Paramj[1]
                            
            Params = array([c12, c21])
                        
            ParamInstance = ModelInstance(append(Params, Cell_s))
            self.Plotter(Params, m, b, Cell_s, Model, ParamInstance, Actual, c, T, Converge)

        return Params

##=============================================================##
Models = ('DWPM', 'NRTL', 'UNIQUAC')
ModelInstances = (GibbsClasses.DWPM, GibbsClasses.NRTL, GibbsClasses.UNIQUAC)
MixtureDataDir = 'Data/Mixtures'
PureDataDir = 'Data/PureComps'
Bounds = [((-15, 0.001), (-15, 0.001)), ((-1000, 3000), (-1000, 3000)), ((-800, 3000), (-800, 3000))]
Scale = ((15, 15), (4000, 4000), (3800, 3800))
InitParamsLimit =((0.00, 0.15), (-100., 1200.), (-800., 3000.))
R = 8.314

if not(path.exists('Results/')):
    mkdir('Results/')   

for file in listdir(MixtureDataDir):

    h5file = tables.openFile(MixtureDataDir+'/'+file, 'r')
    Compounds = h5file.root.Compounds.read()
    h5file.close()

    Optimization = Mixture(Compounds, MixtureDataDir, PureDataDir)

    if not(path.exists('Results/'+Optimization.Name)):
        mkdir('Results/'+Optimization.Name)

    for i in arange(size(Models)): 
##    for i in array([1]): 
        Name = '-'.join(Compounds)
        print Models[i]

        Optimization.ParamCalc(Models[i], ModelInstances[i], InitParamsLimit[i], R)

        h5file = tables.openFile('Results/'+Name+'/'+ Models[i]+'/'+ Name +'.h5', 'r')
        PlotT = array([row['T'] for row in h5file.root.Outputs.iterrows()])
        PlotParams = array([row['ModelParams'] for row in h5file.root.Outputs.iterrows()])
        if Models[i] == "DWPM":
            PlotPureParams = array([row['PureCompParams'] for row in h5file.root.Outputs.iterrows()])
            PlotLambda = PlotParams[:,:2]/PlotPureParams
        h5file.close()

        matplotlib.rc('text', usetex = True)
        fig = matplotlib.pyplot.figure()
        if Models[i] == "DWPM":
            matplotlib.pyplot.plot(PlotT, PlotParams[:,:1],'ro', PlotT, PlotParams[:,1:2],'ro', PlotT, PlotPureParams[:,:1], 'ko', PlotT, PlotPureParams[:,1:2],\
                                       'ko',PlotT, PlotLambda[:,:1], 'r*', PlotT, PlotLambda[:,1:2], 'r*')
        else:
            matplotlib.pyplot.plot(PlotT, PlotParams[:,:1],'ro', PlotT, PlotParams[:,1:2],'ro')
##        matplotlib.pyplot.ylabel(r'Parameters for '+self.Name.capitalize(), fontsize = 14)
        matplotlib.pyplot.xlabel(r'Temperature', fontsize = 14)
        matplotlib.pyplot.title(r'\textbf{Parameters for '+Optimization.Name.capitalize()+'}', fontsize = 14)
        savefig('Results/'+Name+'/'+Models[i]+'/VariationOfParameters.pdf')
        matplotlib.pyplot.close()
