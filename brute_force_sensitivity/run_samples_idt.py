# evaluate each sample 

import sys
import numpy as np
import cantera as ct

# Public vars
# ===========================================
# gas = ct.Solution('H2Reaction_Konnov.xml')	
# p = 1*ct.one_atm
# T = 1000
# T_IGN = 1500
# COMP = 'H2:2,O2:1,AR:3.76'
p = float(sys.argv[1])*ct.one_atm
T = float(sys.argv[2])
COMP = sys.argv[3]
mech = sys.argv[4]
n_start = int(sys.argv[5])
n_end = int(sys.argv[6])

gas = ct.Solution(mech)	

print(p, T, COMP, mech, n_start, n_end)

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

if __name__ == '__main__':
	# print('\n* INFO: run_sample.py')

	out_list = np.zeros(n_end-n_start+1)
	for i in range(n_start-1,n_end):
		factor = uq_input[i,:]
		out_list[i-n_start+1] = ign_uq(factor)
		# print(out_list[i-n_start+1])
		
	np.savetxt( "data/samples_out_idt.txt"+"_"+str(n_start)+"_"+str(n_end), out_list )