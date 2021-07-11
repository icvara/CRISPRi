#ploting parameter

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import abc_smc
import statistics

from collections import Counter

final_par = pd.DataFrame(columns=['Finf','xc','n','L1','B1','log_k1','n1','distance'], index= ['sg1','sg2','sg3','sg4','sg5','sg6','sg1t4','sg4t4'])

number="final"


sample='sg1'
for sample in ['sg1','sg2','sg3','sg4','sg5','sg6','sg1t4','sg4t4']:
    path = 'smc_'+sample+'\pars_' + number + '.out'
    dist_path = 'smc_'+sample+'\distances_' + number + '.out'

    raw_output= np.loadtxt(path)
    dist_output= np.loadtxt(dist_path)
    df = pd.DataFrame(raw_output, columns = ['Finf','xc','n','L1','B1','log_k1','n1'])
    df['dist']=dist_output

    # print(df[df['dist']==min(df['dist'])])

    p_min=df[df['dist']==min(df['dist'])]

    

    p= [ p_min['Finf'].tolist()[0], p_min['xc'].tolist()[0], p_min['n'].tolist()[0],
          p_min['L1'].tolist()[0], p_min['B1'].tolist()[0], p_min['log_k1'].tolist()[0],
          p_min['n1'].tolist()[0] ]
    p2=p.copy()
    p2.append(min(df['dist']))
    final_par.loc[sample]=p2
    

    path='data/data.txt'
    dataframe = pd.read_csv(path,sep='\t')
    x_zero=1e-8
    x_data,y_data,max_input= abc_smc.Get_data(dataframe,sample)
    x_model=np.logspace(np.log10(x_zero), np.log10(max(x_data)),100,base=10)
    x_model[0]=0
    y_model=abc_smc.model(x_model,max_input,p)



    x=np.copy(x_data)
    x2=np.copy(x_model)

    x[0] = x2[0] = x_zero

    plt.close() #
    plt.plot(x2, y_model,'-',label="models")
    plt.plot(x,y_data,'o',label=sample)
    plt.xscale('log')
    plt.xlabel("arabinose [%]")
    plt.ylabel("Normalized GFP")
    plt.legend()
    plt.xscale('log')
    labels=x_data
    plt.xticks(x, labels, rotation='vertical')
    #plt.show()
    plt.savefig(sample+'data+fit.pdf', bbox_inches='tight') #

    plt.close()

    g=sns.pairplot(df[['Finf','xc','n','L1','B1','log_k1','n1']], kind='kde')
    #plt.show()
    plt.savefig(sample+'par_plot.pdf', bbox_inches='tight')
    print(sample)

#plt.savefig('ALL_data+fit.pdf', bbox_inches='tight')
final_par.to_csv('used_par.csv')
