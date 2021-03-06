#CRISPRi equation + fitting

import matplotlib.pyplot  as plt
import numpy as np


#sg2 data time 100 from javier
x_data=np.array([0, 3.125e-06, 6.250e-06, 1.250e-05, 2.500e-05, 5.000e-05,
                 1.000e-04, 2.000e-04, 2.000e-01])
y_data = np.array([ 541388.61, 499462.04, 450763.51, 295709.72,
                   141880.68,  69208.92, 39025.37,  29463.05,  12465.16])
y_data = y_data/np.max(y_data) #fix the max to 1

#sg1 data time 100 from javier

y2_data= np.array([427075.52, 332918.73, 208716.48,
                   128102.74,  84676.15,  63619.62, 58268.09,  61605.91,  31883.79] )

#sg1t4 data time 100 from javier

y3_data= np.array([482205.10, 389792.31, 283826.17, 174227.41, 113170.48,  76051.62,
                  67850.27,  69234.54,  34065.97])

#sg4 data time 100 from javier

y4_data= np.array([ 363403.53, 291546.26, 215998.83, 145227.34,  99504.62,  70626.14,
                     57525.04,  52571.88,  19049.05])



#params
class params():
    B_gfp=1#1.009157919809281490e+00
    L_gfp=0#4.028399867341930785e-02
    K_crispri =-1# -4.854782966997141180e+00
    n_crispri =2# 1.855743579580714453e+00

    n_X=2
    K_X=-4.0
    B_X=1.0
    L_X=0.1

    
    delta_gfp = 1
    
    K_crispra =1
    n_crispra=1


#chem

ARA=np.array([0.2,	0.0002,	0.0001,	0.00005,	0.000025,	0.0000125,	0.00000625,
              0.000003125,	0])

GFP=0




#other param
maxtime=100
dt=0.1


par=params()

#equation

def model(ARA,par):
    X= par.L_X + (np.power((par.B_X-par.L_X)*ARA,par.n_X))/(np.power(10**par.K_X,par.n_X)+np.power(ARA,par.n_X))
    GFP= par.L_gfp + (par.B_gfp -par.L_gfp) / (1+(X/10**par.K_crispri)**par.n_crispri)
    #GFP = par.L_gfp + (par.B_gfp -par.L_gfp) / (1+(ARA/10**par.K_crispri)**par.n_crispri) 
    return X,GFP


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
i=0
plt.subplot(1,1,1)

for k in [0,1,10,100,1000,10000]:
    par.n_crispri=1
    par.K_crispri=k
    GFP= Integration(maxtime,dt,ARA)
    plt.plot(ARA,GFP[int(maxtime/dt)],"-go")

    par.n_crispri=2
    GFP= Integration(maxtime,dt,ARA)
    plt.plot(ARA,GFP[int(maxtime/dt)],"-bo")
 
    par.n_crispri=8
    GFP= Integration(maxtime,dt,ARA)
    plt.plot(ARA,GFP[int(maxtime/dt)],"-ro")


plt.xscale('log')
plt.ylabel('GFP')
plt.xlabel('ARA')
plt.show()

'''

X,GFP=model(x_data,par)
x=np.copy(x_data)
x[0]=1e-7
#plt.plot(x,y_data,'-bo', label="sg2")
#plt.plot(x,y2_data,'-ro',label="sg1")
#plt.plot(x,y3_data,'--ro',label="sg1t4")
#plt.plot(x,y4_data,'-go',label="sg4")
plt.plot(x,GFP,'-go',label="GFP")
plt.plot(x,X,'-ro',label="X")





plt.title("")
plt.xlabel("arabinose [%]")
plt.ylabel("Normalized GFP")
plt.legend()

plt.xscale('log')
labels=x_data
plt.xticks(x, labels, rotation='vertical')



plt.show()
#plt.savefig('model_par_manuel.pdf', bbox_inches='tight')








