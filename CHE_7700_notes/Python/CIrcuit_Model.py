#### this file creates model for simple dCas9 resource competition between two sgRNAs in regulating different GFP reporters.
####
####
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def model(z,t,para,condition):

	Cas9 = z[0]
	sgRNA_1 = z[1]
	sgRNA_2 = z[2]
	Cas9C_1 = z[3]
	Cas9C_2 = z[4]

	P_1Free = z[5]
	P_2Free = z[6]

	P_1Bound = z[7]
	P_2Bound = z[8]

	GFP_1 = z[9]
	GFP_2 = z[10]


	dCas9dt = para['alpha_cas9']*condition['P_cas9'] - para['gamma_sg1']*Cas9*sgRNA_1 - para['gamma_sg2']*Cas9*sgRNA_2

	dsgRNA_1dt = para['alpha_sg1']*condition['P_sg1'] - para['delta_sg1']*sgRNA_1 - para['gamma_sg1']*Cas9*sgRNA_1

	dsgRNA_2dt = para['alpha_sg2']*condition['P_sg2'] - para['delta_sg2']*sgRNA_2 - para['gamma_sg2']*Cas9*sgRNA_2

	dCas9C_1dt = para['gamma_sg1']*Cas9*sgRNA_1 - para['omega_1']*Cas9C_1*condition['P_1']

	dCas9C_2dt = para['gamma_sg2']*Cas9*sgRNA_2 - para['omega_2']*Cas9C_2*condition['P_2']

	dP_1Freedt = -1* para['omega_1']*Cas9C_1*P_1Free

	dP_2Freedt = -1* para['omega_2']*Cas9C_2*P_2Free

	dP_1Bounddt = para['omega_1']*Cas9C_1*P_1Free

	dP_2Bounddt = para['omega_2']*Cas9C_2*P_2Free

	dGFP_1dt = para['alpha_GFP1']*P_1Free

	dGFP_2dt = para['alpha_GFP2']*P_2Free

	dzdt = [dCas9dt, dsgRNA_1dt, dsgRNA_2dt, dCas9C_1dt, dCas9C_2dt, dP_1Freedt, dP_2Freedt, dP_1Bounddt, dP_2Bounddt, dGFP_1dt, dGFP_2dt]
	return dzdt


para = {'alpha_cas9':0.05,
			'alpha_sg1':0.5,
			'alpha_sg2':0.5,
			'alpha_GFP1':0.05,
			'alpha_GFP2':0.05,
			'delta_sg1':0.5,
			'delta_sg2':0.5,
			'gamma_sg1':10**5,
			'gamma_sg2':10**5,
			'omega_1':10**5,
			'omega_2':10**5
	}






#### number of time points
tf = 7200
#### time points
t = np.linspace(0,7200,tf)
#### store solution


condition = {'P_cas9': 10**(-10),
				'P_sg1': {},
				'P_sg2': {},
				'P_1': 10**(-9),
				'P_2': 10**(-9)
		}



fig,axe = plt.subplots(1,1,sharey = False)

for sub in range(3):
#### initial condition
	

	if sub == 0:
		condition['P_sg1'] = 5*10**(-9)
		condition['P_sg2'] = 0
	elif sub == 1:
		condition['P_sg1'] = 0
		condition['P_sg2'] = 5*10**(-9)
	else:
		condition['P_sg1'] = 5*10**(-9)
		condition['P_sg2'] = 5*10**(-9)

	z0 = [0,0,0,0,0,condition['P_1'],condition['P_2'],0,0,0,0]
	Output = {}
	for i in range(11):
		Output[i] = np.empty_like(t)
		Output[i][0] = z0[i]


	for ij in range(1,tf):
		tspan = [t[ij-1], t[ij]]
		z = odeint(model,z0,tspan,args=(para,condition))
		for j in range(11):
			Output[j][ij] = z[1][j]
		z0 = z[1]

	axe.plot(t,Output[9],'--',linewidth=2)
	# axe[sub].plot(t,Output[7],'r--')
	# axe[sub].plot(t,Output[9],'k--')
	
	
	axe.set_xlabel('time')

#axe[0].set_xlim(0, 1.8*10**(-7))
axe.legend(['sg1','sg2','sg1+sg2'])
axe.set_ylabel('GFP1 Concentration')
plt.show()




