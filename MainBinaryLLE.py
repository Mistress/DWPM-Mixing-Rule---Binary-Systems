#!/usr/bin/env python

from scipy import *
from pylab import *
from os import mkdir, path, remove, listdir
import scipy.optimize
import glob
import tables
import vanDerWaals
import ErrorClasses
import GibbsClasses
import PhaseStability


class ResultsFile1(tables.IsDescription):
    T = tables.Float64Col()
    ModelParams = tables.Float64Col(shape = (3)) 
    Predicted = tables.Float64Col(shape = (2))
    Actual = tables.Float64Col(shape = (2))
    SumSqrError = tables.Float64Col()
    AbsError = tables.Float64Col(shape = (2))

class ResultsFile2(tables.IsDescription):
    T = tables.Float64Col()
    ModelParams = tables.Float64Col(shape = (2)) 
    Predicted = tables.Float64Col(shape = (2))
    Actual = tables.Float64Col(shape = (2))
    SumSqrError = tables.Float64Col()
    AbsError = tables.Float64Col(shape = (2))


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
        self.AdachiLuParam = dict((Compound,  Slope[Compound]['BestFit']) for Compound in Compounds)
        h5file.close()
    
    def Plotter(self, BestParams, Model, Fit, ModelInstance, Actual, c, T):
       
        deltaGibbsmix = array([ModelInstance.deltaGmix(x, T, c, self.M) for x in arange(0.001, 1, 0.001)])
        ActualPlot = array([ModelInstance.deltaGmix(x, T, c, self.M) for x in Actual])
        TangentComps = PhaseStability.CalcPhaseStability(ModelInstance, T, c, self.M)
        TangentGibbs = [ModelInstance.deltaGmix(x, T, c, self.M)for x in TangentComps]
        AbsError = ErrorClasses.SumAbs(TangentComps, Actual).Error()
        matplotlib.rc('text', usetex = True)
        fig = matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(arange(0.001, 1, 0.001), deltaGibbsmix, 'b-',TangentComps, TangentGibbs, 'r-', Actual, ActualPlot, 'ko')
        matplotlib.pyplot.xlabel(r'Mole Fraction of '+self.Compounds[0].capitalize(), fontsize = 14)
        matplotlib.pyplot.ylabel(r'$\displaystyle\frac{\Delta G_{mix}}{RT}$', fontsize = 14)
        matplotlib.pyplot.title(r'\textbf{Predicted Phase Equilibrium at %3.2f}'%T, fontsize = 14)
        ax = fig.add_axes([0,0,1,1])
        if Model == 'DWPM':
            ax.text(0,0, 'Model: %s, Params: %.3f, %.3f, %.3f, Abslote Error(sum): %4.4e'%(Model,BestParams[0],BestParams[1], BestParams[2],AbsError), fontsize=12, transform=ax.transAxes)
        else:
            ax.text(0,0, 'Model: %s, Params: %.3f, %.3f, Abslote Error(sum): %4.4e'%(Model,BestParams[0],BestParams[1] ,AbsError), fontsize=12, transform=ax.transAxes)
        ax.set_axis_off()
        #show()
        savefig('Results/'+self.Name+'/'+Model+'/'+ Fit+'/T_'+str(T) +'.pdf')
        matplotlib.pyplot.close()
        
        if not(path.exists('Results/'+self.Name+'/'+ Model+'/'+Fit+'/'+self.Name +'.h5')):
            h5file = tables.openFile('Results/'+self.Name+'/'+Model+'/'+Fit+'/'+self.Name +'.h5', 'w', "Optimization Outputs")
            if Model=='DWPM':                                             
                table = h5file.createTable("/", "Outputs", ResultsFile1, "Optimal model parameters, predicted phase equilibrium, errors etc")
            else:
                table = h5file.createTable("/", "Outputs", ResultsFile2, "Optimal model parameters, predicted phase equilibrium, errors etc")

            table.row['T'] = T
            table.row['ModelParams'] = array(BestParams)
            table.row['Predicted'] = TangentComps
            table.row['Actual'] = Actual
            table.row['SumSqrError'] = ErrorClasses.SumSquare(TangentComps, Actual).Error()
            table.row['AbsError'] = ErrorClasses.AbsError(TangentComps, Actual).Error()

            table.row.append()
            table.flush()
            h5file.close()
        else:
            h5file = tables.openFile('Results/'+self.Name+'/'+Model+'/'+Fit+'/'+self.Name +'.h5', 'r+')
            table = h5file.root.Outputs
            table.row['T'] = T
            table.row['ModelParams'] = array(BestParams)
            table.row['Predicted'] = TangentComps
            table.row['Actual'] = Actual
            table.row['SumSqrError'] = ErrorClasses.SumSquare(TangentComps, Actual).Error()
            table.row['AbsError'] = ErrorClasses.AbsError(TangentComps, Actual).Error()
            
            table.row.append()
            table.flush()
            h5file.close()

           
    def OptFunctionIndvT(self, params, ModelInstance, Actual, T, c, Scale):
        
        UnScaledParams = params*Scale
        Predicted = PhaseStability.CalcPhaseStability(ModelInstance(UnScaledParams) , T, c, self.M)
        Error = ErrorClasses.SumSquare(Predicted ,Actual).Error()

        return Error
  
    def BestFitParamsIndvT(self,Model, ModelInstance, InitParams, Bounds, Scale):
        Fit = 'IndividualT'
        if not(path.exists('Results/'+self.Name+'/'+Model)):
            mkdir('Results/'+self.Name+'/'+Model)        
        if path.exists('Results/'+self.Name+'/'+Model+'/'+Fit):
            fileList = listdir('Results/'+self.Name+'/'+Model+'/'+Fit)
            for fn in fileList: 
                remove(path.join('Results/'+self.Name+'/'+Model+'/'+Fit, fn))
        else:
            mkdir('Results/'+self.Name+'/'+Model+'/'+Fit)
            

        ScaledInitParams = array(InitParams)/array(Scale)
        BoundsList = [array(item) for item in Bounds]
        ScaledBoundsArrays = [BoundsList[i]/array(Scale)[i] for i in arange(size(Scale))]
        ScaledBounds = [tuple(item) for item in ScaledBoundsArrays]
        
        for T in self.M['T']:
                    
            Actual =  array([interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][0])),interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][1]))])
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in self.Compounds]
            
            [params, fx, its, imode, smode] = scipy.optimize.fmin_slsqp(self.OptFunctionIndvT, ScaledInitParams,[], None, [], None, ScaledBounds, None, None, None, (ModelInstance, Actual, T, c, array(Scale)), 1000, 1e-6, 1, 1, 1e-6)
            #params, fx, dict = scipy.optimize.fmin_l_bfgs_b(self.OptFunctionIndvT, InitParams, None, (ModelInstance, Actual, T, c),1, Bounds, 10, 1000, 1e-6, 1e-6, 1,1500)
            ScaledInitParams = params    
            self.Plotter(params*array(Scale), Model, Fit, ModelInstance(params*array(Scale)), Actual, c, T)
            
        return params*array(Scale)
    
    def OptFunctionOvrlT(self, params, ModelInstance, Scale):

        Errors = zeros(size(self.M['T']))
        UnScaledParams = params*Scale

        for T in self.M['T']:
            Actual =  array([interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][0])),interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][1]))])
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in Compounds]
            Predicted = PhaseStability.CalcPhaseStability(ModelInstance(UnScaledParams), T, c, self.M)
            Errors[where(self.M['T']== T)] = ErrorClasses.SumAbs(Predicted ,Actual).Error()
        
        OvrlError = sum(Errors**2)

        return OvrlError
        
    def BestFitParamsOvrlT(self, Model, ModelInstance, InitParams, Bounds, Scale):
        
        Fit = 'OverallT'
        if not(path.exists('Results/'+self.Name+'/'+Model)):
            mkdir('Results/'+self.Name+'/'+Model)        
        if path.exists('Results/'+self.Name+'/'+Model+'/'+Fit):
            fileList = listdir('Results/'+self.Name+'/'+Model+'/'+Fit)
            for fn in fileList: 
                remove(path.join('Results/'+self.Name+'/'+Model+'/'+Fit, fn))
        else:
            mkdir('Results/'+self.Name+'/'+Model+'/'+Fit)

        ScaledInitParams = array(InitParams)/array(Scale)
        BoundsList = [array(item) for item in Bounds]
        ScaledBoundsArrays = [BoundsList[i]/array(Scale)[i] for i in arange(size(Scale))]
        ScaledBounds = [tuple(item) for item in ScaledBoundsArrays]

        [params, fx, its, imode, smode] = scipy.optimize.fmin_slsqp(self.OptFunctionOvrlT, ScaledInitParams,[], None, [], None, ScaledBounds, None, None, None, (ModelInstance, Scale), 1000, 10e-6, 1, 1, 10e-6)
         
        for T in self.M['T']:
            Actual =  array([interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][0])),interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][1]))])
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in self.Compounds]
            self.Plotter(params*array(Scale), Model, Fit, ModelInstance(params*array(Scale)), Actual, c, T)

        return params*array(Scale)

    
       

##=============================================================##
Models = ('DWPM', 'NRTL', 'UNIQUAC')
ModelInstances = (GibbsClasses.DWPM, GibbsClasses.NRTL, GibbsClasses.UNIQUAC)
MixtureDataDir = 'Data/Mixtures'
PureDataDir = 'Data/PureComps'
Bounds = [((-1500, 0), (-1500, 0), (0.5, 0.5)), ((-1000, 3000), (-1000, 3000)), ((-800, 3000), (-800, 3000))]
Scale = ((1500, 1500, 1), (4000, 4000), (4000, 4000))

R = 8.314

if not(path.exists('Results/')):
    mkdir('Results/')   


 

for file in listdir(MixtureDataDir):

    h5file = tables.openFile(MixtureDataDir+'/'+file, 'r')
    Compounds = h5file.root.Compounds.read()
    InitUNIQUAC = tuple(h5file.root.DechemaParams.UNIQUAC.read()[:,0])
    InitNRTL = tuple(h5file.root.DechemaParams.NRTL.read()[:,0])
    h5file.close()
    Optimization = Mixture(Compounds, MixtureDataDir, PureDataDir) 
    InitParams =[(-1000.0,-100.0, 0.5),InitNRTL, InitUNIQUAC]

    if not(path.exists('Results/'+Optimization.Name)):
        mkdir('Results/'+Optimization.Name)

    
    for i in arange(size(Models)):  
        Name = '-'.join(Compounds)
        if path.exists('Results/'+Name+'/'+ Models[i]+'/IndividualT/'+Name +'.h5'):
            remove('Results/'+Name+'/'+ Models[i]+'/IndividualT/'+Name +'.h5') 

        Optimization.BestFitParamsIndvT(Models[i], ModelInstances[i], InitParams[i], Bounds[i], Scale[i])

        h5file = tables.openFile('Results/'+Name+'/'+ Models[i]+'/IndividualT/'+Name +'.h5', 'r')
        PlotT = array([row['T'] for row in h5file.root.Outputs.iterrows()])
        PlotExpX = array([row['Actual'] for row in h5file.root.Outputs.iterrows()])
        PlotPredX = array([row['Predicted'] for row in h5file.root.Outputs.iterrows()])
        h5file.close()
        matplotlib.rc('text', usetex = True)
        fig = matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(PlotPredX, PlotT, 'r-', PlotExpX, PlotT, 'ko')
        matplotlib.pyplot.xlabel(r'Mole Fraction of '+Compounds[0].capitalize(), fontsize = 14)
        matplotlib.pyplot.ylabel(r'Temperature', fontsize = 14)
        matplotlib.pyplot.title(r'\textbf{Predicted Phase Diagram}', fontsize = 14)
        savefig('Results/'+Name+'/'+Models[i]+'/IndividualT/PhaseDiagram.pdf')
        matplotlib.pyplot.close()

        if path.exists('Results/'+Name+'/'+ Models[i]+'/OverallT/'+Name +'.h5'):
            remove('Results/'+Name+'/'+ Models[i]+'/OverallT/'+Name +'.h5') 

        Optimization.BestFitParamsOvrlT(Models[i], ModelInstances[i], InitParams[i], Bounds[i],Scale[i])

        h5file = tables.openFile('Results/'+Name+'/'+ Models[i]+'/OverallT/'+Name +'.h5', 'r')
        PlotT = array([row['T'] for row in h5file.root.Outputs.iterrows()])
        PlotExpX = array([row['Actual'] for row in h5file.root.Outputs.iterrows()])
        PlotPredX = array([row['Predicted'] for row in h5file.root.Outputs.iterrows()])
        h5file.close()
        matplotlib.rc('text', usetex = True)
        fig = matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(PlotPredX, PlotT, 'r-', PlotExpX, PlotT, 'ko')
        matplotlib.pyplot.xlabel(r'Mole Fraction of '+Compounds[0].capitalize(), fontsize = 14)
        matplotlib.pyplot.ylabel(r'Temperature', fontsize = 14)
        matplotlib.pyplot.title(r'\textbf{Predicted Phase Diagram}', fontsize = 14)
        savefig('Results/'+Name+'/'+Models[i]+'/OverallT/PhaseDiagram.pdf')
        matplotlib.pyplot.close()


   
