#ploting parameter

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import abc_smc

number="final" 

for sample in ['sg1','sg2','sg3','sg4','sg5','sg6','sg1t4','sg4t4']:
    path = 'smc_'+sample+'\pars_' + number + '.out'
    raw_output= np.loadtxt(path)
    df = pd.DataFrame(raw_output, columns = ['Finf','xc','n','L1','B1','log_k1','n1'])
    p=np.mean(df)
    print(p)
    
    path='data/data.txt'
    dataframe = pd.read_csv(path,sep='\t')
    x_data,y_data,max_input= abc_smc.Get_data(dataframe,sample)
    x=np.copy(x_data)
    x[0]=1e-8

   # for index, p in dfmean.iterrows():
    y_model=abc_smc.model(x,max_input,p)
   # plt.plot(x_data, y_model,'-r')
    plt.close()
    plt.plot(x, y_model,'-r',label="models")
    plt.plot(x,y_data,'bo',label=sample)
    plt.xscale('log')
    plt.xlabel("arabinose [%]")
    plt.ylabel("Normalized GFP")
    plt.legend()
    plt.xscale('log')
    labels=x_data
    plt.xticks(x, labels, rotation='vertical')
    #plt.show()
    plt.savefig(sample+'data+fit.pdf', bbox_inches='tight')

    plt.close()
    sns.pairplot(df)
    plt.savefig(sample+'par_plot.pdf', bbox_inches='tight')


'''
path='data/data.txt'
dataframe = pd.read_csv(path,sep='\t')
dataframe['arabinose']=dataframe['arabinose'].astype(str)

x_data,y_data,max_input= abc_smc.Get_data(dataframe,sample)



for index, p in df1.iterrows():
    y_model=p['Finf']+(p['F0']-p['Finf'])/(1+np.power(x_data/p['xc'],p['n']))
    plt.plot(x_data, y_model,'--g')

dtp=df

plt.title("")
x=np.copy(x_data)
x[0]=1e-8

for index, p in dtp.iterrows():
    y_model=abc_smc.model(x,max_input,p)
    plt.plot(x, y_model,'-r')

plt.plot(x, y_model,'-r',label="models")

plt.plot(x,y_data,'-bo',label="sg2")
plt.xscale('log')
plt.xlabel("arabinose [%]")
plt.ylabel("Normalized GFP")
plt.legend()
plt.xscale('log')
labels=x_data
plt.xticks(x, labels, rotation='vertical')
#plt.show()
plt.savefig(sample+'model_par_abc_smc.pdf', bbox_inches='tight')


sns.pairplot(dtp)
plt.savefig(sample+'par_plot.pdf', bbox_inches='tight')



plt.title("")
plt.xlabel("arabinose [%]")
plt.ylabel("Normalized GFP")
plt.legend()

plt.xscale('log')
labels=x_data
plt.xticks(x, labels, rotation='vertical')



#plt.show()
plt.savefig('model_par_manuel.pdf', bbox_inches='tight')


plt.show()

'''
print('youpi')
