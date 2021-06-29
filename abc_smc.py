import numpy as np
from scipy.stats import norm, uniform, multivariate_normal
from scipy.optimize import minimize
import sys,ast
from random import choices,seed,random
from tqdm import tqdm
from p_tqdm import p_umap # parallel code
from functools import partial
import os

#x_data = np.array([0.1, 1, 2, 10, 100])
#y_data = np.array([2, 1, 0.9, 0.8, 0.7])
#sg2 data time 100 from javier
x_data=np.array([0, 3.125e-06, 6.250e-06, 1.250e-05, 2.500e-05, 5.000e-05,
                 1.000e-04, 2.000e-04, 2.000e-01])
y_data = np.array([ 541388.61, 499462.04, 450763.51, 295709.72,
                   141880.68,  69208.92, 39025.37,  29463.05,  12465.16])
y_data = y_data/np.max(y_data) #fix the max to 1

if os.path.isdir('smc') is False: ## if 'smc' folder does not exist:
    os.mkdir('smc') ## create it, the output will go there

def pars_to_dict(pars):
### This function is not necessary, but it makes the code a bit easier to read,
### it transforms an array of pars e.g. p[0],p[1],p[2] into a
### named dictionary e.g. p['k0'],p['B'],p['n'],p['x0']
### so it is easier to follow the parameters in the code
    dict_pars = {}
    for ipar,par in enumerate(parlist):
        dict_pars[par['name']] = pars[ipar] 
    return dict_pars

parlist = [ # list containing information of each parameter
    {'name' : 'F0', 'lower_limit':0.1,'upper_limit':10.0},# fluorescence in absence of repressor
    {'name' : 'Finf','lower_limit':0.1,'upper_limit':10.0},# fluorescence in abundance of repressor
    {'name' : 'xc','lower_limit':0.1,'upper_limit':100000.0},# characeristic concentration of the repressor
    {'name' : 'n','lower_limit':0.5,'upper_limit':4.0}# Hill exponent
]

def sampleprior():
### Generate a random parameter inside the limits stablished. The shape of the distribution can be changed if required
    prior = []
    for ipar,par in enumerate(parlist):
        prior.append(uniform.rvs(loc = par['lower_limit'],
                                 scale = par['upper_limit']))
    return prior

def evaluateprior(pars, model = None):
### Given a set of parameters return how probable is given the prior knowledge.
### If the priors are uniform distributions, this means that it will return 0 when a parameter
### outside the prior range is tested.    
    prior = 1
    for ipar,par in enumerate(parlist):
            prior *= uniform.pdf(pars[ipar],loc = par['lower_limit'],
                scale = par['upper_limit'])
    return prior

def model(x,pars):
### This is the function you are trying to fit, so given a set of parameters (pars)
### and a set of concentrations x of signal/inducer, it returns the predicted response
### For instance for a repressive hill function:
    p = pars_to_dict(pars)
    return p['Finf']+(p['F0']-p['Finf'])/(1+np.power(x/p['xc'],p['n'])) 

def distance(pars): #### 
### This function defines the distance between the data and the model for a given set of parameters
### for instance it can be a least squares minimization e.g. for a set of experimental observations
### stored as xdata and ydata:
    return np.sum(np.power(y_data-model(x_data,pars),2))


def GeneratePar(processcall = 0,  previousparlist = None, previousweights = None,
                eps_dist = 10000, kernel = None):
# This function generates a valid parameter given a set of sampled pars previousparlist,previousweights
# that is under a threshold eps_dist

# processcall is a dummy variable that can be useful when tracking the function performance
# it also allows the use of p_tqdm mapping that forces the use of an iterator 

        seed() # setting random seeds for each thread/process to avoid having the same random sequence in each thread
        np.random.seed()

        evaluated_distances = []
        d = eps_dist + 1 # original distance is beyond eps_dist
        while d > eps_dist: # until a suitable parameter is found
            if previousparlist is None: # if is the first iteration
                proposed_pars = sampleprior()
            else:
                selected_pars = choices(previousparlist, weights = previousweights)[0] # select a point from the previous sample
                proposed_pars = selected_pars + kernel.rvs() # and perturb it

            if (evaluateprior(proposed_pars,model) > 0):

                d = distance(proposed_pars) 
                evaluated_distances.append(d)

        # Calculate weight
        if previousparlist is None:
            weight = 1
        else:
            sum_denom = 0
            for ipars,pars in enumerate(previousparlist):
                kernel_evaluation = kernel.pdf(proposed_pars-pars)
                sum_denom += kernel_evaluation*previousweights[ipars]

            weight = evaluateprior(proposed_pars)/sum_denom
        return proposed_pars, d, weight, evaluated_distances

def GeneratePars(previousparlist = None, previousweights = None, eps_dist = 10000,
    Npars = 1000, kernelfactor = 1.0):
## Call to function GeneratePar in parallel until Npars points are accepted

    if previousparlist is not None:
        previouscovar = 2.0*kernelfactor*np.cov(np.array(previousparlist).T)
        # print('covariance matrix previous parset:',previouscovar)
        kernel = multivariate_normal(cov = previouscovar)

    else:
        kernel = None

    trials = 0

    # for N in tqdm(range(Npars)): ## not parallel version
    #   GenerateParstePar(0,model = 'Berg',gamma = gamma, previousparlist = previousparlist,
    #   previousweights = previousweights, eps_dist = eps_dist, kernel = kernel)

    results = p_umap(
        partial(GeneratePar, previousparlist = previousparlist,
        previousweights = previousweights, eps_dist = eps_dist, kernel = kernel), range(Npars))
    accepted_distances = [result[1] for result in results]
    evaluated_distances = [res for result in results for res in result[3]] # flatten list 
    newparlist = [result[0] for result in results]
    newweights = [result[2] for result in results]

    newweights /= np.sum(newweights) #Normalizing
    print("acceptance rate:", Npars/len(evaluated_distances))
    print("min accepted distance: ",np.min(accepted_distances))
    print("median accepted distance: ",np.median(accepted_distances))
    print("median evaluated distance: ",np.median(evaluated_distances))

    return(newparlist,newweights, accepted_distances, Npars/len(evaluated_distances) )

#final_dist =100
def Sequential_ABC(initial_dist = 10000, final_dist =0.01, Npars = 1000, prior_label = None,
                   adaptative_kernel = True):

## Sequence of acceptance threshold start with initial_dis and keeps on reducing until
## a final threshold final_dist is reached.

## prior_label can be used to restart a sampling from a previous prior distribution in case
## further exploration with a lower epsilon is needed    

    distance = initial_dist
    not_converged = True
    last_round = False
    kernelfactor = 1.0

    if prior_label is None:
        pars = None
        weights = None
        idistance = 0
    else: # a file with the label is used to load the posterior, use always numerical (not 'final')
        pars = np.loadtxt('smc/pars_{}_{}_{}_{}.out'.format(model,sto,gamma,prior_label))
        weights = np.loadtxt('smc/weights_{}_{}_{}_{}.out'.format(model,sto,gamma,prior_label))
        accepted_distances = np.loadtxt('smc/distances_{}_{}_{}_{}.out'.format(model,sto,gamma,prior_label))
        distance = np.median(accepted_distances)
        idistance = prior_label

    while not_converged:

        idistance += 1
        print("SMC step with target distance: {}".format(distance))
        pars,weights,accepted_distances,acceptance = GeneratePars( previousparlist = pars,
            previousweights = weights, eps_dist = distance,
            Npars = Npars, kernelfactor = kernelfactor)
        proposed_dist = np.median(accepted_distances)
        if last_round is True:
            not_converged = False
            label = 'final'
        else:
            label = idistance
        if proposed_dist<final_dist:
            distance = final_dist
            last_round = True
        else:
            distance = proposed_dist
        np.savetxt('smc/pars_{}.out'.format(label), pars)
        np.savetxt('smc/weights_{}.out'.format(label), weights)
        np.savetxt('smc/distances_{}.out'.format(label), accepted_distances)

        if acceptance < 0.1 and kernelfactor>0.1 and adaptative_kernel:
            kernelfactor = kernelfactor * 0.7
            print('Reducing kernel width to : ',kernelfactor)
        elif acceptance > 0.5 and kernelfactor<1 and adaptative_kernel:
            kernelfactor = kernelfactor / 0.7
            print('Increasing kernel width to : ',kernelfactor)


    

Sequential_ABC()



