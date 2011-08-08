#!/usr/bin/env python

from scipy import *
from pylab import *
from os import mkdir, path, remove, listdir
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

    def ModelsPlotter(self, Models, ModelInstances, Cell_s, alpha, z):

        if path.exists('Results/'+self.Name+'/AllModelsGibbsPlots' ):
            fileList = listdir('Results/'+self.Name+'/AllModelsGibbsPlots')
            for fn in fileList: 
                remove(path.join('Results/'+self.Name+'/AllModelsGibbsPlots', fn))
        else:
            mkdir('Results/'+self.Name+'/AllModelsGibbsPlots')

            
        for i in arange(size(Models)):

            h5file = tables.openFile('Results/'+self.Name+'/'+ Models[i]+'/'+ self.Name +'.h5', 'r')
            if Models[i]=="DWPM":
                DWPMParams = array([row['ModelParams'] for row in h5file.root.Outputs.iterrows()])
            elif Models[i]=="NRTL":
                NRTLParams = array([row['ModelParams'] for row in h5file.root.Outputs.iterrows()])
            elif Models[i]=="UNIQUAC":
                UNIQUACParams = array([row['ModelParams'] for row in h5file.root.Outputs.iterrows()])
            h5file.close()

        for T in self.M['T']:
            
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in Compounds]
            Actaul = array([self.M['ExpComp'][0], self.M['ExpComp'][1]])
            ##ActualPlot = array([ModelInstance.deltaGmix(x, T, c, self.M) for x in Actual])

            print append(array([interp(T, self.M['T'], DWPMParams[:, 0]), interp(T, self.M['T'], DWPMParams[:, 1]) ]), Cell_s)
            DWPMInstance = ModelInstances[0](append(array([interp(T, self.M['T'], DWPMParams[:, 0]), interp(T, self.M['T'], DWPMParams[:, 1]) ]) , Cell_s))
            DWPMPlot = array([DWPMInstance.deltaGmix(x, T, c, self.M) for x in arange(0.001, 1, 0.001)])

            NRTLInstance = ModelInstances[1](array([interp(T, self.M['T'], NRTLParams[:, 0]), interp(T, self.M['T'], NRTLParams[:, 1]) ]))
            NRTLPlot = array([NRTLInstance.deltaGmix(x, T, c, self.M) for x in arange(0.001, 1, 0.001)])

            UNIQUACInstance = ModelInstances[2](array([interp(T, self.M['T'], UNIQUACParams[:, 0]), interp(T, self.M['T'], UNIQUACParams[:, 1]) ]))
            UNIQUACPlot = array([UNIQUACInstance.deltaGmix(x, T, c, self.M) for x in arange(0.001, 1, 0.001)])

            fig = matplotlib.pyplot.figure()
            matplotlib.pyplot.plot(arange(0.001, 1, 0.001), DWPMPlot, 'b-', arange(0.001, 1, 0.001), NRTLPlot, 'r-', arange(0.001, 1, 0.001), UNIQUACPlot, 'k-')
            matplotlib.pyplot.xlabel(r'Mole Fraction of '+self.Compounds[0].capitalize(), fontsize = 14)
            matplotlib.pyplot.ylabel(r'$\displaystyle\frac{\Delta G_{mix}}{RT}$', fontsize = 14)
            matplotlib.pyplot.title(r'\textbf{Predicted $\displaystyle\frac{\Delta G_{mix}}{RT}$ %3.2f}'%T, fontsize = 14)
            #ax = fig.add_axes([0,0,1,1])
            #ax.set_axis_off()

            savefig('Results/'+self.Name+'/AllModelsGibbsPlots/T_'+str(T) +'.pdf')
            matplotlib.pyplot.close()          
            
    def ParamCalc(self, Model, ModelInstance, InitParamsLimit, Bounds, R, System, Cell_s):

        if path.exists('Results/'+self.Name+'/'+Model):
            fileList = listdir('Results/'+self.Name+'/'+Model)
            for fn in fileList: 
                remove(path.join('Results/'+self.Name+'/'+Model, fn))
        else:
            mkdir('Results/'+self.Name+'/'+Model)
                        
        InitParams = (InitParamsLimit[1]-InitParamsLimit[0])*random(2) + append(InitParamsLimit[0], InitParamsLimit[0])
        InitParamj = append(InitParams, array([-1., -3.]))      
                        
        for T in self.M['T']:
            Actual =  array([interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][0])),interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][1]))])
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in Compounds]
           
            MaxFEval = 1000
            TolParam = 1e-4
            Converge = 3
            NStarts = 1
            MaxNStarts = 10000

            while (Converge!=1 and NStarts<MaxNStarts):
                                  
                      [Paramj, infodict, Converge, Mesg] = scipy.optimize.fsolve(System.SystemEquations, InitParamj, (Actual, R, T), System.SystemEquationsJac, 1, 0, TolParam, MaxFEval, None, 0.0, 100, None)
                     ## print Mesg
                      NFCalls = infodict['nfev']
                      NJCalls = infodict['njev']
                      NStarts = NStarts + 1
                      
                      if (Paramj[3]>0.01) or not(Bounds[0]<= Paramj[0]<=Bounds[1]) or not(Bounds[0]<= Paramj[1]<=Bounds[1]):
                          Converge = 3
                      
                      InitParams = (InitParamsLimit[1]-InitParamsLimit[0])*random(2) + append(InitParamsLimit[0], InitParamsLimit[0])
                      InitParamj = append(InitParams, array([-1., -5.]))  
                      
            ##print NStarts
            ##print Paramj
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

            if Model == "DWPM":
                ParamInstance = ModelInstance(append(Params, Cell_s))
                self.Plotter(Params, m, b, Cell_s, Model, ParamInstance, Actual, c, T, Converge)
            else:
                ParamInstance = ModelInstance(Params)
                self.Plotter(Params, m, b, (), Model, ParamInstance, Actual, c, T, Converge)            

             
    def PhaseDiagramCalc(self, Model, System, Params, R, T):

        CompC = self.vdWaalsInstance.CompC(T)
        c = [CompC[Compound] for Compound in Compounds]

        if Model == 'DWPM':
            Params[0] = Params[0]/c[0]
            Params[1] = Params[1]/c[1]
            
        MaxFEval = 1000
        TolPhase = 1e-4
        Converge = 3
        NStarts = 1
        MaxNStarts = 10000
        
        while (Converge!=1 and NStarts<MaxNStarts):
            
            InitPhase = append(array([0.1, -0.1])*random(2)+array([0.0, 1.0]), array([-2., -3]))
                        
            [Phase, infodict, Converge, Mesg] = scipy.optimize.fsolve(System.SystemEquations, InitPhase, (Params, R, T), System.SystemEquationsJac, 1, 0, TolPhase, MaxFEval, None, 0.0, 100, None)
            
            NFCalls = infodict['nfev']
            NJCalls = infodict['njev']
            NStarts = NStarts + 1
            
            if (Phase[3]>0.01) or not(0.0 <= Phase[0]<=1.0) or not(0.0 <= Phase[1]<= 1.0) or (abs(Phase[0]-Phase[1])<0.1):
                Converge = 3
        
        ##print NStarts
        ##print Phase
        InitPhase = Phase

        Predicted1 = Phase[0]
        Predicted2 = Phase[1]
        m = Phase[2]
        b = Phase[3]

        return array([Predicted1, Predicted2])

    def OverallSOpt(self, Cell_s, Model, ModelInstance, InitParamsLimit, Bounds, R):

        System = AnalyticSystemsClasses.SystemDWPM(Cell_s)
        self.ParamCalc(Model, ModelInstance, InitParamsLimit, Bounds, R, System, Cell_s)
        
        h5file = tables.openFile('Results/'+self.Name+'/'+ Model+'/'+ self.Name +'.h5', 'r')
        DWPMParams = array([row['ModelParams'] for row in h5file.root.Outputs.iterrows()])
        h5file.close()
        
        return sum(abs(std(DWPMParams[:,0])),abs( std(DWPMParams[:,1])))
                                   
##=============================================================##
Models = ('DWPM', 'NRTL', 'UNIQUAC')
ModelInstances = (GibbsClasses.DWPM, GibbsClasses.NRTL, GibbsClasses.UNIQUAC)
MixtureDataDir = 'Data/Mixtures'
PureDataDir = 'Data/PureComps'
Bounds = ((0.00, 100), (-1000, 3000), (-800, 3000))
InitParamsLimit =((0.00, 0.15), (-100., 1000.), (-800., 3000.))
alpha = 0.2
z = 10
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
        if Models[i] == "DWPM":
            InitCell_s = array([0.3, 0.8])
            Cell_s = scipy.optimize.fmin(Optimization.OverallSOpt, InitCell_s,(Models[i], ModelInstances[i], InitParamsLimit[i], Bounds[i], R), xtol=0.01, ftol=0.001, maxiter=None, maxfun=None, full_output=0, disp=1, retall=0, callback=None)
            System = AnalyticSystemsClasses.SystemDWPM(Cell_s)
            PhaseSystem = AnalyticPhaseCompSystemsClasses.SystemDWPM(Cell_s)
        elif Models[i] == "NRTL":
            System = AnalyticSystemsClasses.SystemNRTL(alpha)
            PhaseSystem = AnalyticPhaseCompSystemsClasses.SystemNRTL(alpha)
        elif Models[i] == "UNIQUAC":
            System = AnalyticSystemsClasses.SystemUNIQUAC(Optimization.M, z)
            PhaseSystem = AnalyticPhaseCompSystemsClasses.SystemUNIQUAC(Optimization.M, z)

        Name = '-'.join(Compounds)
        print Models[i]

        Optimization.ParamCalc(Models[i], ModelInstances[i], InitParamsLimit[i], Bounds[i], R, System, Cell_s)
        
        h5file = tables.openFile('Results/'+Name+'/'+ Models[i]+'/'+ Name +'.h5', 'r')
        PlotT = array([row['T'] for row in h5file.root.Outputs.iterrows()])
        PlotParams = array([row['ModelParams'] for row in h5file.root.Outputs.iterrows()])
        PlotActual = array([row['Actual'] for row in h5file.root.Outputs.iterrows()])
        
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

        LinearParams1 = scipy.interpolate.interp1d(PlotT, PlotParams[:,0])
        LinearParams2 = scipy.interpolate.interp1d(PlotT, PlotParams[:,1])
        
        InterpT = linspace(PlotT[0], PlotT[-1], 10)
        PredictedLinear = array([Optimization.PhaseDiagramCalc(Models[i], PhaseSystem, array([LinearParams1(T), LinearParams2(T)]), R, T) for T in InterpT])
      
        fig = matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(PredictedLinear[:,:1], InterpT, 'r--', PredictedLinear[:,1:2], InterpT, 'r--', PlotActual[:,:1], PlotT, 'ko', PlotActual[:,1:2], PlotT, 'ko')
        matplotlib.pyplot.ylabel(r'Temperature', fontsize = 14)
        matplotlib.pyplot.xlabel(r'Composition', fontsize = 14)
        matplotlib.pyplot.title(r'\textbf{Phase diagram for '+Optimization.Name.capitalize()+'}', fontsize = 14)
        savefig('Results/'+Name+'/'+Models[i]+'/PhaseDiagram.pdf')
        matplotlib.pyplot.close()

 ##   Optimization.ModelsPlotter(Models, ModelInstances, Cell_s, alpha, z)
