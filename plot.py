#ploting parameter

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

number="17" 

#raw_output= np.loadtxt('smc\pars_final.out')
path = 'smc\pars_' + number + '.out'
raw_output= np.loadtxt(path)

df = pd.DataFrame(raw_output, columns = ['F0','Finf','xc','n'])


#distance= np.loadtxt('smc\distances_final.out')
#distance= np.loadtxt('smc\distances_17.out')
#df['distance']=distance
#df_2= df[df["distance"] < 0.05]

x_data=np.array([0, 3.125e-06, 6.250e-06, 1.250e-05, 2.500e-05, 5.000e-05,
                 1.000e-04, 2.000e-04, 2.000e-01])
y_data = np.array([ 541388.61, 499462.04, 450763.51, 295709.72,
                   141880.68,  69208.92, 39025.37,  29463.05,  12465.16])
y_data = y_data/np.max(y_data) #fix the max to 1

#y_model=p['Finf']+(p['F0']-p['Finf'])/(1+np.power(x/p['xc'],p['n']))

'''
for index, p in df1.iterrows():
    y_model=p['Finf']+(p['F0']-p['Finf'])/(1+np.power(x_data/p['xc'],p['n']))
    plt.plot(x_data, y_model,'--g')
'''
dtp=df

plt.title("")
x=np.copy(x_data)
x[0]=1e-8

for index, p in dtp.iterrows():
    y_model=p['Finf']+(p['F0']-p['Finf'])/(1+np.power(x_data/p['xc'],p['n']))
    plt.plot(x, y_model,'--r')

plt.plot(x, y_model,'--r',label="models")

plt.plot(x,y_data,'-b',label="sg2")
plt.xscale('log')
plt.xlabel("arabinose [%]")
plt.ylabel("Normalized GFP")
plt.legend()
plt.xscale('log')
labels=x_data
plt.xticks(x, labels, rotation='vertical')
#plt.show()
plt.savefig('model_par_abc_smc.pdf', bbox_inches='tight')


sns.pairplot(dtp[['Finf','F0','xc','n']])
plt.savefig('par_plot.pdf', bbox_inches='tight')

'''

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
