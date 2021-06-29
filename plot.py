#ploting parameter

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

number="5" 

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

for index, p in dtp.iterrows():
    y_model=p['Finf']+(p['F0']-p['Finf'])/(1+np.power(x_data/p['xc'],p['n']))
    plt.plot(x_data, y_model,'--r')

plt.plot(x_data,y_data,'-b')
plt.xscale('log')
plt.show()
'''
sns.pairplot(df1[['Finf','F0','xc','n']])
plt.show()
'''
sns.pairplot(dtp[['Finf','F0','xc','n']])
plt.show()
