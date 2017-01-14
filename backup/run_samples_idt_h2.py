# evaluate each sample 

import sys
import numpy as np
import cantera as ct

# Public vars
# ===========================================
gas = ct.Solution('H2Reaction_Konnov.xml')	
# p = 1*ct.one_atm
# T = 1000
# T_IGN = 1500
# COMP = 'H2:2,O2:1,AR:3.76'
p = float(sys.argv[1])*ct.one_atm
T = float(sys.argv[2])
COMP = sys.argv[3]
print(p, T, COMP)

gas.TPX = T, p, COMP
gas.equilibrate('HP')
T_equi = gas.T

uq_input = np.loadtxt('data/samples.txt')
index = np.loadtxt('data/samples_index.txt')	
N_sample = uq_input.shape[0]
# ===========================================

def ign_uq(factor):
	gas.set_multiplier(1.0) # reset all multipliers
	for i in range(len(index)):
		gas.set_multiplier(factor[i],index[i]-1) #reaction index start from 1
		# print(gas.reaction_equation(index[i]-1)+' index_reaction:',index[i],'multi_factor:',factor)
	
	gas.TPX = T,p,COMP
	r = ct.IdealGasReactor(gas)
	sim = ct.ReactorNet([r])
	t_end = 10;
	time = []
	temp = []
	# X_OH = []
	while sim.time < t_end and r.T < T_equi - 10.0:
		sim.step()
		time.append(sim.time)
		temp.append(r.T)
		# X_OH.append(r.thermo['OH'].X)

	diff_temp = np.diff(temp)/np.diff(time)
	ign_indice = np.argmax(diff_temp)
	
	return time[ign_indice]

def ign_simple(factor):
	# run this sime ignition delay function fisrt to 
	# determine the ignition break temperature and other parameters
	i = 0
	# gas = ct.Solution('H2Reaction_Konnov.xml')	
	for i in range(len(index)):
		gas.set_multiplier(factor[i],index[i]-1) #reaction index start from 1
		print(gas.reaction_equation(index[i]-1)+' index_reaction:',index[i],'multi_factor:',factor)
		i = i+1
	gas.TPX = T,p, COMP
	r = ct.IdealGasReactor(gas)
	sim = ct.ReactorNet([r])
	t_end = 1;
	time = []
	temp = []
	X_OH = []
	while sim.time < t_end and r.T < T_equi - 10.0:
		sim.step()
		time.append(sim.time)
		temp.append(r.T)
		X_OH.append(r.thermo['OH'].X)

	diff_temp = np.diff(temp)/np.diff(time)
	ign_indice = np.argmax(diff_temp)
	time = np.array(time)
	temp = np.array(temp)
	print(sim.time,r.T-T,temp[ign_indice]-T)

	import matplotlib.pyplot as plt
	plt.plot(time, temp)
	#plt.legend(['Reactor 1','Reactor 2'],2)
	plt.xlabel('Time (s)')
	plt.ylabel('Temperature (K)')
	plt.title('H2:2,O2:1,AR:3.76'+',%.2f [atm] %.0f [K] IDT = %lf [s]'%(p,T,time[ign_indice]) )
	plt.savefig('temp'+'.png',dpi=300)
	plt.show()

	return time[ign_indice]

if __name__ == '__main__':
	# print('\n* INFO: run_sample.py')

	# ==========================DEBUG=======================================================
	# run ign_simple to determine T_IGN, namely the temperature cretition for stop iteration
	# index = [16,29]
	# factor = [1,1]
	# ign = ign_simple(factor)
	# print(ign)
	# exit()
	# ==========================DEBUG=======================================================

	ign_list = np.zeros(N_sample)
	for i in range(0,N_sample):
		factor = uq_input[i,:]
		ign_list[i] = ign_uq(factor)
		# print(ign)
		
	np.savetxt( "data/samples_out_idt.txt", ign_list )