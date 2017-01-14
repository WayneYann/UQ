# evaluate each sample 

import sys,math,os,re,csv
import numpy as np
import cantera as ct

def ign_uq(p,T,T_IGN,index,factor):
	gas = ct.Solution('H2Reaction_Konnov.xml')
	i = 0
	for i in range(len(index)):
		gas.set_multiplier(factor[i],index[i]-1) #reaction index start from 1
		# print(gas.reaction_equation(index[i]-1)+' index_reaction:',index[i],'multi_factor:',factor)
		i = i+1
	
	gas.TPX = T,p*ct.one_atm,'H2:2,O2:1,AR:3.76'
	r = ct.IdealGasReactor(gas)
	sim = ct.ReactorNet([r])
	t_end = 1;
	time = []
	temp = []
	# X_OH = []
	while sim.time < t_end and r.T < T+T_IGN:
		sim.step(t_end)
		time.append(sim.time)
		temp.append(r.T)
		# X_OH.append(r.thermo['OH'].X)

	diff_temp = np.diff(temp)/np.diff(time)
	ign_indice = np.argmax(diff_temp)
	
	return time[ign_indice]

def ign_simple(p,T,T_IGN,index,factor):
	# run this sime ignition delay function fisrt to 
	# determine the ignition break temperature and other parameters
	gas = ct.Solution('H2Reaction_Konnov.xml')
	i = 0
	for i in range(len(index)):
		gas.set_multiplier(factor[i],index[i]-1) #reaction index start from 1
		print(gas.reaction_equation(index[i]-1)+' index_reaction:',index[i],'multi_factor:',factor)
		i = i+1

	gas.TPX = T,p*ct.one_atm,'H2:2,O2:1,AR:3.76'
	r = ct.IdealGasReactor(gas)
	sim = ct.ReactorNet([r])
	t_end = 1;
	time = []
	temp = []
	X_OH = []
	while sim.time < t_end and r.T < T+T_IGN:
		sim.step(t_end)
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
	index = [16,29]
	factor = [1,1]
	uq_input = np.loadtxt('samples.txt')
	index = np.loadtxt('samples_index.txt')
	
	N_sample = uq_input.shape[0]
	 
	T_IGN = 1800
	p = 20
	T = 1000
	# run ign_simple to determine T_IGN, namely the temperature cretition for stop iteration
	# ign = ign_simple(p,T,T_IGN,index,factor)
	# print(ign)
	# exit()
	with open("ignout.txt", "w") as output:
		for i in range(0,N_sample):
			factor = uq_input[i,:]
			#factor[1] = 1
			ign = ign_uq(p,T,T_IGN,index,factor)
			output.write(str(ign)+'\n')
			# print(ign)