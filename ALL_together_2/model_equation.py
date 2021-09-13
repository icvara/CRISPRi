#CRISPRi equation + fitting

import matplotlib.pyplot  as plt
import numpy as np
import os
import pandas as pd


initdist=1000#40.
finaldist=0.5#2.5
priot_label=None
 

parlist = [ # list containing information of each parameter

  #  {'name' : 'leak_1','lower_limit':0.0,'upper_limit':0.5},# fluorescence in abundance of repressor
    {'name' : 'log_k_1','lower_limit':-4.5,'upper_limit':0},# characeristic concentration of the repressor
    {'name' : 'n_sg','lower_limit':1.0,'upper_limit':3.5},# Hill exponent

  #  {'name' : 'leak_1t4','lower_limit':0.0,'upper_limit':0.5},# fluorescence in abundance of repressor
    {'name' : 'log_k_1t4','lower_limit':-4.5,'upper_limit':0},# characeristic concentration of the repressor
   # {'name' : 'n_1t4','lower_limit':1.0,'upper_limit':3.5},# Hill exponent

  #  {'name' : 'leak_2','lower_limit':0.0,'upper_limit':0.5},# fluorescence in abundance of repressor
    {'name' : 'log_k_2','lower_limit':-4.5,'upper_limit':0},# characeristic concentration of the repressor
  #  {'name' : 'n_2','lower_limit':1.0,'upper_limit':3.5},# Hill exponent

  #  {'name' : 'leak_3','lower_limit':0.0,'upper_limit':0.5},# fluorescence in abundance of repressor
    {'name' : 'log_k_3','lower_limit':-4.5,'upper_limit':0},# characeristic concentration of the repressor
   # {'name' : 'n_3','lower_limit':1.0,'upper_limit':3.5},# Hill exponent

  #  {'name' : 'leak_4','lower_limit':0.0,'upper_limit':0.5},# fluorescence in abundance of repressor
    {'name' : 'log_k_4','lower_limit':-4.5,'upper_limit':0},# characeristic concentration of the repressor
   # {'name' : 'n_4','lower_limit':1.0,'upper_limit':3.5},# Hill exponent

  #  {'name' : 'leak_4t4','lower_limit':0.0,'upper_limit':0.5},# fluorescence in abundance of repressor
    {'name' : 'log_k_4t4','lower_limit':-4.5,'upper_limit':0},# characeristic concentration of the repressor
   # {'name' : 'n_4t4','lower_limit':1.0,'upper_limit':3.5},# Hill exponent

  #  {'name' : 'leak_5','lower_limit':0.0,'upper_limit':0.5},# fluorescence in abundance of repressor
    {'name' : 'log_k_5','lower_limit':-4.5,'upper_limit':0},# characeristic concentration of the repressor
  #  {'name' : 'n_sg','lower_limit':1.0,'upper_limit':3.5},# Hill exponent

  #  {'name' : 'leak_6','lower_limit':0.0,'upper_limit':0.5},# fluorescence in abundance of repressor
    {'name' : 'log_k_6','lower_limit':-4.5,'upper_limit':0},# characeristic concentration of the repressor
   # {'name' : 'n_sg','lower_limit':1.0,'upper_limit':3.5},# Hill exponent

    #first node which is activate by ara(x) and repress the GFP
    {'name' : 'leak_ARA', 'lower_limit':0.0,'upper_limit':0.1},# fluorescence in abundance of inducer
    {'name' : 'B_ARA','lower_limit':0.5,'upper_limit':1.5},# fluorescence in absence of inducer
    {'name' : 'k_ARA','lower_limit':-5,'upper_limit':-3.5},# characeristic concentration of the activation
    {'name' : 'n_ARA','lower_limit':1.0,'upper_limit':3.5}# Hill exponent
]


def Get_data(df,s):
    sub_df=df[df['sample']==s]
    max_input= sub_df['GFP'][df['arabinose']=='off-target'].tolist()
    x_data = np.array(sub_df['arabinose'][1:len(sub_df['arabinose'])].astype(float).tolist())
    y_data = np.array(sub_df['GFP'][1:len(sub_df['GFP'])].tolist())
    return x_data,y_data, max_input

#equation
def model(x,max_input, p):
### This is the function you are trying to fit, so given a set of parameters (pars)
### and a set of concentrations x of signal/inducer, it returns the predicted response
### For instance for a repressive hill function:

    #p=pars
    node1 = p['leak_ARA'] + np.power((p['B_ARA']-p['leak_ARA'])*x,p['n_ARA'])/(np.power(10**p['k_ARA'],p['n_ARA'])+np.power(x,p['n_ARA']))
    #return p['Finf']+(p['F0']-p['Finf'])/(1+np.power(node1/10**p['log_xc'],p['n']))
    node2_sg1 = (max_input['sg1'])/(1+np.power(node1/10**p['log_k_1'],p['n_sg'])) 
    node2_sg1t4 = (max_input['sg1t4'])/(1+np.power(node1/10**p['log_k_1t4'],p['n_sg'])) 
    node2_sg2 = (max_input['sg2'])/(1+np.power(node1/10**p['log_k_2'],p['n_sg'])) 
    node2_sg3 = (max_input['sg3'])/(1+np.power(node1/10**p['log_k_3'],p['n_sg'])) 
    node2_sg4 = (max_input['sg4'])/(1+np.power(node1/10**p['log_k_4'],p['n_sg'])) 
    node2_sg4t4 = (max_input['sg4t4'])/(1+np.power(node1/10**p['log_k_4t4'],p['n_sg'])) 
    node2_sg5 = (max_input['sg5']-)/(1+np.power(node1/10**p['log_k_5'],p['n_sg'])) 
    node2_sg6 = (max_input['sg6'])/(1+np.power(node1/10**p['log_k_6'],p['n_sg'])) 
    
    return node2_sg1,node2_sg1t4,node2_sg2,node2_sg3,node2_sg4,node2_sg4t4,node2_sg5,node2_sg6

def distance(pars, x_data, y_data, max_input): #### 
### This function defines the distance between the data and the model for a given set of parameters
### for instance it can be a least squares minimization e.g. for a set of experimental observations
### stored as xdata and ydata:
    d1,d1t4,d2,d3,d4,d4t4,d5,d6 = model(x_data,max_input,pars)


    
    dist = np.sum(np.power(y_data['sg1']-d1,2)) + np.sum(np.power(y_data['sg1t4']-d1t4,2)) + np.sum(np.power(y_data['sg2']-d2,2)) + np.sum(np.power(y_data['sg3']-d3,2)) + np.sum(np.power(y_data['sg4']-d4,2)) + np.sum(np.power(y_data['sg4t4']-d4t4,2))+np.sum(np.power(y_data['sg5']-d5,2))+np.sum(np.power(y_data['sg6']-d6,2))
                                                                                                                                                                                                                                                                                                

    return dist
    
    
path='data/input.txt'

sample = ['sg1','sg1t4','sg2' ,'sg3' ,'sg4' ,'sg4t4','sg5','sg6']
path='data/data.txt'
#path='C:/Users/Administrator/Desktop/Modeling/CRISPRi/data/data.txt'
dataframe = pd.read_csv(path,sep='\t')
dataframe['arabinose']=dataframe['arabinose'].astype(str)
  
x_data = Get_data(dataframe, 'sg1')[0]  
y_data={}
max_input={}
for s in sample:
    x,y, max_i = Get_data(dataframe, s)
    y_data[s]=y
    max_input[s]=max_i[0]
    
#print(max_input)
#print(y_data)  
    
    
    
    
    
    
'''
def Flow_crispri(par,GFP,ARA):
    #crispri - repression
    dGFP_dt = par.B_gfp / (1+(par.K_crispri*ARA)**par.n_crispri) - par.delta_gfp*GFP

    return dGFP_dt

def Flow_crispra(par,GFP,ARA):
    #crispra -activation
    dGFP_dt = par.B_gfp*(par.K_crispra*ARA)**par.n_crispra / (1+(par.K_crispra*ARA)**par.n_crispra) - par.delta_gfp*GFP

    return dGFP_dt


def Integration(maxtime,dt, ARA):

    totalstep= int(maxtime/dt)+1


    GFP=np.zeros((totalstep,len(ARA)))
    ti=0
    i=0
    GFPi=np.zeros(len(ARA))

    while (ti<maxtime):
        GFP[i]=GFPi
        dGFP_dt=Flow_crispri(par,GFPi,ARA)
        GFPi= GFPi + dGFP_dt*dt
        
        ti=ti+dt
        i=i+1
        
    return GFP
'''   





