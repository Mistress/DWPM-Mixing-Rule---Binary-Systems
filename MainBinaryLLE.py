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


class ResultsFile1(tables.IsDescription):
    T = tables.Float64Col()
    ModelParams = tables.Float64Col(shape = (4)) 
    Actual = tables.Float64Col(shape = (2))
    PureCompParams = tables.Float64Col(shape = (2))

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
        h5file.close()
        self.Binaries = zeros((size(self.M['T']), 2))
    
    def Plotter(self, Params, m, b, Cell_s, Model, ModelInstance, Actual, c, T):
       
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
            ax.text(0,0, 'Model: %s, Params: %.3f, %.3f, %.3f, %.3f'%(Model, Params[0], Params[1], Cell_s[0], Cell_s[1]), fontsize=12, transform=ax.transAxes)
        else:
            ax.text(0,0, 'Model: %s, Params: %.3f, %.3f, Abslote Error(sum): %4.4e'%(Model,BestParams[0],BestParams[1] ,AbsError), fontsize=12, transform=ax.transAxes)
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
        else:
            table.row['Predicted'] = TangentComps
            table.row['SumSqrError'] = SumSqrError
            table.row['AbsError'] = AbsError            
        table.row.append()
        table.flush()
        h5file.close()

    def SystemJac(self, Params, Actual, Cell_s):
        
        lambda_12 = Params[0]
        lambda_21 = Params[1]
        m = Params[2]
        b = Params[3]

        s_1 = Cell_s[0]
        s_2 = Cell_s[1]

        x = Actual[0]
                
        Row1 = [-lambda_12**(s_1-1)*(1-x)*x/(x+lambda_12**s_1*(1-x)), -lambda_21**(s_2-1)*(1-x)*x/(lambda_21**s_2*x-x+1), -x, -1.]
        Row3 = [lambda_12**(s_1-1)*x/(x+lambda_12**s_1*(1-x))-lambda_12**(s_1-1)*(1-x)/(x+lambda_12**s_1*(1-x))+\
        lambda_12**(s_1-1)*(1-lambda_12**s_1)*(1-x)*x/(x+lambda_12**s_1*(1-x))**2,\
        lambda_21**(s_2-1)*x/(lambda_21**s_2*x-x+1)-lambda_21**(s_2-1)*(1-x)/(lambda_21**s_2*x-x+1)+\
        lambda_21**(s_2-1)*(lambda_21**s_2-1)*(1-x)*x/(lambda_21**s_2*x-x+1)**2,-1., 0.]

        x = Actual[1]
        
        Row2 = [-lambda_12**(s_1-1)*(1-x)*x/(x+lambda_12**s_1*(1-x)), -lambda_21**(s_2-1)*(1-x)*x/(lambda_21**s_2*x-x+1), -x, -1.]
        Row4 = [lambda_12**(s_1-1)*x/(x+lambda_12**s_1*(1-x))-lambda_12**(s_1-1)*(1-x)/(x+lambda_12**s_1*(1-x))+\
        lambda_12**(s_1-1)*(1-lambda_12**s_1)*(1-x)*x/(x+lambda_12**s_1*(1-x))**2,\
        lambda_21**(s_2-1)*x/(lambda_21**s_2*x-x+1)-lambda_21**(s_2-1)*(1-x)/(lambda_21**s_2*x-x+1)+\
        lambda_21**(s_2-1)*(lambda_21**s_2-1)*(1-x)*x/(lambda_21**s_2*x-x+1)**2,-1., 0.]

        return array([Row1, Row2, Row3, Row4])

    def System(self, Params, Actual, Cell_s):
        
        lambda_12 = Params[0]
        lambda_21 = Params[1]
        m = Params[2]
        b = Params[3]
        
        s_1 = Cell_s[0]
        s_2 = Cell_s[1]

        x = Actual[0]
        
        Eq1 = -(1-x)*log(lambda_21**s_2*x-x+1)/s_2-x*log(x+lambda_12**s_1*(1-x))/s_1+x*log(x)-m*x+log(1-x)*(1-x)-b
        Eq3 = log(lambda_21**s_2*x-x+1)/s_2-log(x+lambda_12**s_1*(1-x))/s_1+log(x)-(lambda_21**s_2-1)*(1-x)/(s_2*(lambda_21**s_2*x-x+1))\
            -(1-lambda_12**s_1)*x/(s_1*(x+lambda_12**s_1*(1-x)))-log(1-x)-m

        x = Actual[1]
        
        Eq2 =  -(1-x)*log(lambda_21**s_2*x-x+1)/s_2-x*log(x+lambda_12**s_1*(1-x))/s_1+x*log(x)-m*x+log(1-x)*(1-x)-b
        Eq4 = log(lambda_21**s_2*x-x+1)/s_2-log(x+lambda_12**s_1*(1-x))/s_1+log(x)-(lambda_21**s_2-1)*(1-x)/(s_2*(lambda_21**s_2*x-x+1))\
            -(1-lambda_12**s_1)*x/(s_1*(x+lambda_12**s_1*(1-x)))-log(1-x)-m

        return array([Eq1, Eq2, Eq3, Eq4])

    def ParamCalc(self, Model, ModelInstance, InitParams):

        if path.exists('Results/'+self.Name+'/'+Model):
            fileList = listdir('Results/'+self.Name+'/'+Model)
            for fn in fileList: 
                remove(path.join('Results/'+self.Name+'/'+Model, fn))
        else:
            mkdir('Results/'+self.Name+'/'+Model)
            
        if Model == "DWPM":
            Cell_s = (0.5, 0.5)
        else:
            Cell_s = ()
        
        for T in self.M['T']:
            Actual =  array([interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][0])),interp(T,cast['f'](self.M['T']), cast['f'](self.M['ExpComp'][1]))])
            CompC = self.vdWaalsInstance.CompC(T)
            c = [CompC[Compound] for Compound in Compounds]

            InitParamj =  append(InitParams[0]/c[0],append(InitParams[1]/c[1],[1., -5.]))                                

            MaxFEval = 1000
            TolParam = 1e-6

            [Paramj, infodict, Nit, Mesg] = scipy.optimize.fsolve(self.System, InitParamj, (Actual, Cell_s), self.SystemJac, 1, 0, TolParam, MaxFEval, None, 0.0, 100, None)
            print Mesg
            print infodict['nfev']
            print infodict['njev']
            print infodict['fvec']
            print infodict['fjac']
        
##            j = 1
##            Norm = 10
##            Paramj = array([InitParams[0]/c[0], InitParams[1]/c[1], 1., -2.])
       
##            while (j<MaxIter)and(Norm>TolParam):
##                f = self.System(Paramj, Actual, Cell_s)
##                J = self.SystemJac(Paramj, Actual, Cell_s)
##                Parami = Paramj
##                Delta = scipy.linalg.solve(J,-1*f)
##                Paramj = Parami + Delta
##                Norm = scipy.linalg.norm(Delta, inf)
##                j = j+1

##            print(j)
##            print(Paramj)

            lambda_12 = Paramj[0]
            lambda_21 = Paramj[1]
            m = Paramj[2]
            b = Paramj[3]

##          Do a check using Actual[1]-> Should give the same parameter values
  
            c12 = c[0]*lambda_12
            c21 = c[1]*lambda_21

            Params = array([c12, c21])
            print Params
            
            ParamInstance = ModelInstance(append(Params, Cell_s))
            self.Plotter(Params, m, b, Cell_s, Model, ParamInstance, Actual, c, T)

        return Params

##=============================================================##
Models = ('DWPM', 'NRTL', 'UNIQUAC')
ModelInstances = (GibbsClasses.DWPM, GibbsClasses.NRTL, GibbsClasses.UNIQUAC)
MixtureDataDir = 'Data/Mixtures'
PureDataDir = 'Data/PureComps'
Bounds = [((-15, 0.001), (-15, 0.001)), ((-1000, 3000), (-1000, 3000)), ((-800, 3000), (-800, 3000))]
Scale = ((15, 15), (4000, 4000), (3800, 3800))
InitParams =((-9, -1), (300, 300), (300, 300))
#InitParams =((-1/22.5, -1/618.4), (300, 300), (300, 300))

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

##  for i in arange(size(Models)): 
    for i in arange(1): 
        Name = '-'.join(Compounds)
        print Models[i]

        Optimization.ParamCalc(Models[i], ModelInstances[i], InitParams[i])

##        h5file = tables.openFile('Results/'+Name+'/'+ Models[i]+'/'+ Name +'.h5', 'r')
##        PlotT = array([row['T'] for row in h5file.root.Outputs.iterrows()])
##        PlotExpX = array([row['Actual'] for row in h5file.root.Outputs.iterrows()])
##        PlotPredX = array([row['Predicted'] for row in h5file.root.Outputs.iterrows()])
##        h5file.close()

##        matplotlib.rc('text', usetex = True)
##        fig = matplotlib.pyplot.figure()
##        matplotlib.pyplot.plot(PlotPredX, PlotT, 'r-', PlotExpX, PlotT, 'ko')
##        matplotlib.pyplot.xlabel(r'Mole Fraction of '+Compounds[0].capitalize(), fontsize = 14)
##        matplotlib.pyplot.ylabel(r'Temperature', fontsize = 14)
##        matplotlib.pyplot.title(r'\textbf{Predicted Phase Diagram}', fontsize = 14)
##        savefig('Results/'+Name+'/'+Models[i]+'/PhaseDiagram.pdf')
##        matplotlib.pyplot.close()

       


   
