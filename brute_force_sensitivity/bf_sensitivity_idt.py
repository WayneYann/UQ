# evaluate each sample 

import sys
import numpy as np
import cantera as ct

# Public vars
# ===========================================
# gas = ct.Solution('H2Reaction_Konnov.xml')	
# p = 1*ct.one_atm
# T = 1000
# COMP = 'H2:2,O2:1,AR:3.76'
p = float(sys.argv[1])*ct.one_atm
T = float(sys.argv[2])
COMP = sys.argv[3]
mech = sys.argv[4]
perturbation = float(sys.argv[5])
nrxn = int(sys.argv[6])

gas = ct.Solution(mech)	

print(p, T, COMP, mech)

gas.TPX = T, p, COMP
gas.equilibrate('HP')
T_equi = gas.T

# ===========================================

def ign_uq(factor):
	gas.set_multiplier(1.0) # reset all multipliers
	for i in range(nrxn):
		gas.set_multiplier(factor[i],i)
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

	factor = np.ones(nrxn)
	out0 = ign_uq(factor)

	out_list = np.ones(nrxn)
	for i in range(nrxn):
		factor = np.ones(nrxn)
		factor[i] = perturbation
		out_list[i] = ign_uq(factor)/out0
		# print(out_list[i])
		
	np.savetxt( "data/samples_out_idt.txt", out_list )