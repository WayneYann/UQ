# evaluate each sample 

import sys
import numpy as np
import cantera as ct

# Public vars
# ===========================================
# p = 1*ct.one_atm
# T = 300.0
# COMP = 'H2:2,O2:1,N2:3.76'

p = float(sys.argv[1])*ct.one_atm
T = float(sys.argv[2])
COMP = sys.argv[3]
mech = sys.argv[4]
n_start = int(sys.argv[5])
n_end = int(sys.argv[6])

gas = ct.Solution(mech)	
print(p, T, COMP, mech, n_start, n_end)

uq_input = np.loadtxt('data/samples.txt')
index = np.loadtxt('data/samples_index.txt')	
N_sample = uq_input.shape[0]

initial_grid = np.linspace(0, 0.03, 5)  # m
tol_ss = [1.0e-9, 1.0e-14]  # [rtol atol] for steady-state problem
tol_ts = [1.0e-5, 1.0e-14]  # [rtol atol] for time stepping

gas.TPX = T,p,COMP

# Flame object
f = ct.FreeFlame(gas, initial_grid)
f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)

# Solve with the energy equation disabled
f.energy_enabled = False
f.transport_model = 'Mix'
f.set_max_jac_age(10, 10)
f.set_time_step(1e-5, [2, 5, 10, 20])
f.solve(loglevel=0, refine_grid=False)

# Solve with the energy equation enabled
f.set_refine_criteria(ratio=3, slope=0.07, curve=0.14)
f.energy_enabled = True
f.transport_model = 'Multi'
f.solve(loglevel=0, refine_grid=True)

Su0 = f.u[0]
print('\nmixture-averaged flamespeed = {:7f} m/s\n'.format(f.u[0]))
print('Initial Solution:')
# f.show_stats()

# N_sample = 10
out_list = np.zeros(n_end-n_start+1)
for i in range(n_start-1,n_end):
	gas.set_multiplier(1.0) # reset all multipliers
	factor = uq_input[i,:]
	for k in range(len(index)):
		gas.set_multiplier(factor[k],index[k]-1) #reaction index start from 1
    
	f.solve(loglevel=0, refine_grid=False)
	out_list[i-n_start+1] = f.u[0]
	# f.show_stats()
	print('Sample {:4d} mixture-averaged flamespeed = {:7f} m/s\n'.format(i,f.u[0]))
	
np.savetxt( "data/samples_out_sl.txt"+"_"+str(n_start)+"_"+str(n_end), out_list )