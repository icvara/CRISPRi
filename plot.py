#ploting parameter

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statistics
import os
from collections import Counter
import sys
from scipy.signal import argrelextrema
from matplotlib.colors import LogNorm, Normalize
import multiprocessing
import time
from functools import partial


filename="ALL_together_1"
n=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29',
    '30','31','32','33','34','35','36','37']
n=['1','10','20','30','37']
#n=['final']
#
#sys.path.insert(0, '/users/ibarbier/AC-DC/'+filename)
sys.path.insert(0, 'C:/Users/Administrator/Desktop/Modeling/CRISPRi/'+filename)
import model_equation as meq
  
parlist=meq.parlist


######################################################################33
#########################################################################
###########################################################################

def load(number= n,filename=filename,parlist=parlist):
    namelist=[]
    for i,par in enumerate(parlist):
        namelist.append(parlist[i]['name'])
        
    path = filename+'/smc/pars_' + number + '.out'
    dist_path = filename+'/smc/distances_' + number + '.out'

    raw_output= np.loadtxt(path)
    dist_output= np.loadtxt(dist_path)
    df = pd.DataFrame(raw_output, columns = namelist)
    df['dist']=dist_output
    df=df.sort_values('dist',ascending=False)
    distlist= sorted(df['dist'])
    p=[]
    for dist in distlist:
        
        p_0=df[df['dist']==dist]
        p0=[]
        for n in namelist:
          p0.append(p_0[n].tolist()[0])
   
        p0=pars_to_dict(p0)
        p.append(p0)

    
    return p, df 

def pars_to_dict(pars):
### This function is not necessary, but it makes the code a bit easier to read,
### it transforms an array of pars e.g. p[0],p[1],p[2] into a
### named dictionary e.g. p['k0'],p['B'],p['n'],p['x0']
### so it is easier to follow the parameters in the code
    dict_pars = {}
    for ipar,par in enumerate(parlist):
        dict_pars[par['name']] = pars[ipar] 
    return dict_pars
##plotting part

     
def plot_alltime(n,filename,parlist):
    namelist=[]
    for i,par in enumerate(parlist):
        namelist.append(parlist[i]['name'])
    parl = np.append(namelist,'dist')
    index=1
    size=round(np.sqrt(len(parl)))
    for i,name in enumerate(parl):
        plt.subplot(size,size,index)
       # plt.tight_layout()
        for ni,nmbr in enumerate(n):
            p,df= load(nmbr,filename,parlist)
            sns.kdeplot(df[name],bw_adjust=.8,label=nmbr,linewidth=.1)
        #plt.ylim(0,1)
        #if i < (len(parl)-2):
        #    plt.xlim((parlist[i]['lower_limit'],parlist[i]['upper_limit']))
        if index==size:       
          plt.legend(bbox_to_anchor=(1.05, 1))
        index=index+1
    plt.savefig(filename+"/plot/"+'ALLround_plot.pdf', bbox_inches='tight')
    plt.close()




##############################################################################################################3   

if __name__ == "__main__":
   
    if os.path.isdir(filename+'/plot') is False: ## if 'smc' folder does not exist:
        os.mkdir(filename+'/plot') ## create it, the output will go there
    
    X=meq.x_data
    X[X==0]=10e-7
    Y=meq.y_data

    print(Y['sg1'])


    #print(ARA)
    #plot_alltime(n,filename,parlist)
    p,pdf= load('37',filename,parlist)

    p1=pdf.iloc[0]
    for sgi,sg in enumerate(['sg1','sg1t4','sg2','sg3','sg4','sg4t4','sg5','sg6']):
      plt.subplot(2,4,sgi+1)
      plt.xscale("log")
      for pi,par in enumerate(p):
          par=pdf.iloc[pi]
          d1,d1t4,d2,d3,d4,d4t4,d5,d6 = meq.model(X,meq.max_input,par)
          plt.plot(X,d1)      
      plt.plot(X,Y[sg],'o')


    plt.show()
    #print(pdf)



  
 

    
    