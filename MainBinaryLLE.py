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
        self.M = dict(Compounds = h5file.root.Compounds.read(), ExpComp = h5file.root.ExperimentalData.T_CompData.ExpComp.read(), T =h5file.root.ExperimentalData.T_CompData.T.read())
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
        (TangentComps, Fvalue) = PhaseStability.CalcPhaseStability(ModelInstance, T, c, self.M)
        TangentGibbs = [ModelInstance.deltaGmix(x, T, c, self.M)for x in TangentComps]
        AbsError = ErrorClasses.SumAbs(TangentComps, Actual).Error()
        matplotlib.rc('text', usetex = True)
        fig = matplotlib.pyplot.figure()
        matplotlib.pyplot.plot(arange(0.001, 1, 0.001), deltaGibbsmix, 'b-',TangentComps, TangentGibbs, 'r-', Actual, ActualPlot, 'ko')
        matplotlib.pyplot.xlabel(r'Mole Fraction of '+self.Compounds[0].capitalize(), fontsize = 14)
        matplotlib.pyplot.ylabel(r'$\displaystyle\frac{\Delta G_{mix}}{RT}$', fontsize = 14)
        matplotlib.pyplot.title(r'\textbf{Predicted Phase Equilibrium at %3.2f}'%T, fontsize = 14)
        ax = fig.add_axes([0,0,1,1])
        ax.text(0,0, 'Model: %s, Params: %s, Abslote Error(sum): %4.4e, Termination: %s'%(Model, str(around(BestParams, 4)),AbsError, Fvalue), fontsize=12, transform=ax.transAxes)
        ax.set_axis_off()
        #show()
        savefig('Results/'+self.Name+'/'+Model+'/'+ Fit+'/T_'+str(T) +'.pdf')
        matplotlib.pyplot.close()
        print Actual, TangentComps
        print self.Name,Fvalue
                   
              
            
           
    def OptFunctionIndvT(self, params, ModelInstance, Actual, T, c):
        
        Predicted = PhaseStability.CalcPhaseStability(ModelInstance(params) , T, c, self.M)
        Error = ErrorClasses.SumSquare(Predicted ,Actual).Error()

        return Error
  
    def NonEqConstrIndvT(self, params, ModelInstance, Actual, T, c):
        
        deltaGibbsMixTest = array([ModelInstance(params).deltaGmix(x, T, c, self.M) for x in arange(0.001, 1, 0.001)])
        TangentTestComps = PhaseStability.CalcPhaseStability(ModelInstance(params), T, c, self.M)
        TangentTest = array([((ModelInstance(params).FirstDerivative(TangentTestComps[0], T, c, self.M))*(x - TangentTestComps[0])+ ModelInstance(params).deltaGmix(TangentTestComps[0], T, c, self.M))-ModelInstance(params).deltaGmix(x, T, c, self.M) for x in arange(0.001, 1, 0.001)])
        DistinctTest = -1*array([abs(TangentTestComps[0]-TangentTestComps[1])+0.01])
       
        return -1000*append(append(deltaGibbsMixTest, TangentTest), DistinctTest)

    def BestFitParamsIndvT(self,Model, ModelInstance, InitParams, Bounds):
        Fit = 'IndividualT'
        if not(path.exists('Results/'+self.Name+'/'+Model)):
            mkdir('Results/'+self.Name+'/'+Model)        
        if path.exists('Results/'+self.Name+'/'+Model+'/'+Fit):
            fileList = listdir('Results/'+self.Name+'/'+Model+'/'+Fit)
            for fn in fileList: 
                remove(path.join('Results/'+self.Name+'/'+Model+'/'+Fit, fn))
        else:
            mkdir('Results/'+self.Name+'/'+Model+'/'+Fit)
        
        for T in self.M['T']:
                    
            Actual =  array([interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][0])),interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][1]))])
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in self.Compounds]
            
            [params, fx, its, imode, smode] = scipy.optimize.fmin_slsqp(self.OptFunctionIndvT, InitParams,[], None, [], self.NonEqConstrIndvT, Bounds, None, None, None, (ModelInstance, Actual, T, c), 1000, 10e-6, 1, 1, 5)
            #InitParams = params    
            self.Plotter(params, Model, Fit, ModelInstance(params), Actual, c, T)
            
        return params
    
    def OptFunctionOvrlT(self, params, ModelInstance, R):

        Errors = zeros(size(self.M['T']))
        
        for T in self.M['T']:
            Actual =  array([interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][0])),interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][1]))])
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in Compounds]
            Predicted = PhaseStability.CalcPhaseStability(ModelInstance(params), T, c, self.M)
            Errors[where(self.M['T']== T)] = ErrorClasses.SumAbs(Predicted ,Actual).Error()
        
        OvrlError = sum(Errors**2)

        return OvrlError

    def NonEqConstrOvrlT(self, params, ModelInstance,R):
        
        Test = array([])
        
        for T in self.M['T']:
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in Compounds]
            deltaGibbsMixTest = array([ModelInstance(params).deltaGmix(x, T, c, self.M) for x in arange(0.001, 1, 0.001)])
            TangentTestComps = PhaseStability.CalcPhaseStability(ModelInstance(params), T, c, self.M)
            TangentTest = array([((ModelInstance(params).FirstDerivative(TangentTestComps[0], T, c, self.M))*(x - TangentTestComps[0])+ ModelInstance(params).deltaGmix(TangentTestComps[0], T, c, self.M))-ModelInstance(params).deltaGmix(x, T, c, self.M) for x in arange(0.001, 1, 0.001)])
            DistinctTest = -1*array([abs(TangentTestComps[0]-TangentTestComps[1])+0.01])
            Test = append(Test, append(append(deltaGibbsMixTest, TangentTest), DistinctTest))
       
        return -1*Test
        
    def BestFitParamsOvrlT(self, Model, ModelInstance, InitParams, Bounds):
        
        Fit = 'OverallT'
        if not(path.exists('Results/'+self.Name+'/'+Model)):
            mkdir('Results/'+self.Name+'/'+Model)        
        if path.exists('Results/'+self.Name+'/'+Model+'/'+Fit):
            fileList = listdir('Results/'+self.Name+'/'+Model+'/'+Fit)
            for fn in fileList: 
                remove(path.join('Results/'+self.Name+'/'+Model+'/'+Fit, fn))
        else:
            mkdir('Results/'+self.Name+'/'+Model+'/'+Fit)

        [params, fx, its, imode, smode] = scipy.optimize.fmin_slsqp(self.OptFunctionOvrlT, InitParams,[], None, [], self.NonEqConstrOvrlT, Bounds, None, None, None, (ModelInstance, R), 100, 10e-4, 1, 1, 10e-2)
         
        for T in self.M['T']:
            Actual =  array([interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][0])),interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][1]))])
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in self.Compounds]
            self.Plotter(params, Model, Fit, ModelInstance(params), Actual, c, T)

        return params, 

    
       

##=============================================================##
#Models = ('DWPM', 'NRTL', 'UNIQUAC')
#ModelInstances = (GibbsClasses.DWPM, GibbsClasses.NRTL, GibbsClasses.UNIQUAC)
#MixtureDataDir = 'Data/Mixtures'
#PureDataDir = 'Data/PureComps'
#Compounds = ('1-butanol', 'water')
#Bounds = [((-1000, 0), (-1000, 0), (0.5, 0.5)), ((-800, 3000), (-800, 3000)), ((-800, 3000), (-800, 3000))]
#InitParams =[(-200.0,-25.0, 0.5), (-250.0, 1500.0), (-20, 300.00)]
#R = 8.314

#if not(path.exists('Results/')):
#    mkdir('Results/')   


#Optimization = Mixture(Compounds, MixtureDataDir, PureDataDir)     
#if not(path.exists('Results/'+Optimization.Name)):
#        mkdir('Results/'+Optimization.Name)

#for i in arange(size(Models)):  
#    Name = '-'.join(Compounds)
#    if path.exists('Results/'+Name+'/'+ Models[i]+'/IndividualT/'+Name +'.h5'):
#        remove('Results/'+Name+'/'+ Models[i]+'/IndividualT/'+Name +'.h5') 
#    
#    Optimization.BestFitParamsIndvT(Models[i], ModelInstances[i], InitParams[i], Bounds[i])
#    
#    h5file = tables.openFile('Results/'+Name+'/'+ Models[i]+'/IndividualT/'+Name +'.h5', 'r')
#    PlotT = array([row['T'] for row in h5file.root.Outputs.iterrows()])
#    PlotExpX = array([row['Actual'] for row in h5file.root.Outputs.iterrows()])
#    PlotPredX = array([row['Predicted'] for row in h5file.root.Outputs.iterrows()])
#    h5file.close()
#    matplotlib.rc('text', usetex = True)
#    fig = matplotlib.pyplot.figure()
#    matplotlib.pyplot.plot(PlotPredX, PlotT, 'r-', PlotExpX, PlotT, 'ko')
#    matplotlib.pyplot.xlabel(r'Mole Fraction of '+Compounds[0].capitalize(), fontsize = 14)
#    matplotlib.pyplot.ylabel(r'Temperature', fontsize = 14)
#    matplotlib.pyplot.title(r'\textbf{Predicted Phase Diagram}', fontsize = 14)
#    savefig('Results/'+Name+'/'+Models[i]+'/IndividualT/PhaseDiagram.pdf')
#    matplotlib.pyplot.close()
    
   # if path.exists('Results/'+Name+'/'+ Models[i]+'/OverallT/'+Name +'.h5'):
   #     remove('Results/'+Name+'/'+ Models[i]+'/OverallT/'+Name +'.h5') 
    
   # Optimization.BestFitParamsOvrlT(Models[i], ModelInstances[i], InitParams[i], Bounds[i])

   # h5file = tables.openFile('Results/'+Name+'/'+ Models[i]+'/OverallT/'+Name +'.h5', 'r')
   # PlotT = array([row['T'] for row in h5file.root.Outputs.iterrows()])
   # PlotExpX = array([row['Actual'] for row in h5file.root.Outputs.iterrows()])
   # PlotPredX = array([row['Predicted'] for row in h5file.root.Outputs.iterrows()])
   # h5file.close()
   # matplotlib.rc('text', usetex = True)
   # fig = matplotlib.pyplot.figure()
   # matplotlib.pyplot.plot(PlotPredX, PlotT, 'r-', PlotExpX, PlotT, 'ko')
   # matplotlib.pyplot.xlabel(r'Mole Fraction of '+Compounds[0].capitalize(), fontsize = 14)
   # matplotlib.pyplot.ylabel(r'Temperature', fontsize = 14)
   # matplotlib.pyplot.title(r'\textbf{Predicted Phase Diagram}', fontsize = 14)
   # savefig('Results/'+Name+'/'+Models[i]+'/OverallT/PhaseDiagram.pdf')
   # matplotlib.pyplot.close()


#    for Method in Methods:
#        print('Fitting data using ', Method)
#        [Best_Params{m}, x_eq{m}] = CalcBestFitParam(R, T, M, Best_Params, c, delGmixfun, ddelGmix_dx, m )
        


#delGmix_plot(m,:) = arrayfun(@(x) delGmixfun{m}(R, T, M, Best_Params{m}, c, x), xplot);
#Gmix_eq(m,:) = arrayfun(@(x) delGmixfun{m}(R, T, M, Best_Params{m}, c, x), x_eq{m});

