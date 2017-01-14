# evaluate each sample 

import sys
import numpy as np
import cantera as ct
from scipy.integrate import odeint

# Public vars
# ===========================================
gas = ct.Solution('H2Reaction_Konnov.xml')	

p = 1*ct.one_atm
T_in = 300.0
COMP = 'H2:2,O2:1,N2:3.76'
T_ext0 = 1212
Omega_init = 1./(1E-3)

# p = float(sys.argv[1])*ct.one_atm
# T_in = float(sys.argv[2])
# COMP = sys.argv[3]
# T_ext0 = float(sys.argv[4])
# Omega_init = float(sys.argv[5])
# print(p, T_in, COMP)

DT = 300
Ta = T_ext0 - DT
Tb = T_ext0 + DT
gas.TPX = T_in, p, COMP
Y_in = gas.Y

gas.equilibrate('HP')
Y_equi = gas.Y

def PSR_ode(y, t):

	gas.TPY = gas.T, gas.P, y
	HRR = np.dot(gas.net_production_rates, gas.partial_molar_enthalpies)
	omega = -HRR/gas.cp_mass/gas.density/( gas.T - T_in )
	dydt = gas.net_production_rates*gas.molecular_weights/gas.density - omega*( y - Y_in )
	
	return dydt

def gas2ext():

	HRR = np.dot(gas.net_production_rates, gas.partial_molar_enthalpies)
	omega = -HRR/gas.cp_mass/gas.density/( gas.T - T_in )
	ext = 1./omega
	
	return ext

def T2ext(T_ext, y0, omega0):

	t = [0, 20./omega0]	
	gas.TP = T_ext, p
	sol = odeint(PSR_ode, y0, t, rtol=1.e-8)
	ext = gas2ext()

	return ext

def get_ext(TL, TU, y0, omega0):
	# input the interval of temperature for searching the extinction time
	def T2ext_list(T_list, N):
		ext_list = np.zeros( N )
		i = 0
		for T in T_list:
			ext_list[i] = T2ext(T, y0, omega0)
			i = i+1
		return ext_list
	
	N = 200
	rtol = 1.E-4
	T_list = np.linspace(TL, TU, N, endpoint=True)
	ext_list = T2ext_list(T_list, N)
	# print(T_list)
	# print(ext_list)
	ext_index = np.argmin(ext_list)

	# import matplotlib.pyplot as plt
	# plt.semilogx(ext_list, T_list)
	# #plt.legend(['Reactor 1','Reactor 2'],2)
	# plt.xlabel('Time (s)')
	# plt.ylabel('Temperature (K)')
	# plt.show()

	if ext_index == 0 or ext_index == N-1:
		print('please increase the searching range'+' index:'+str(ext_index))
		sys.exit()
		# TL = T_ext0 - DT - 200
		# TU = T_ext0 + DT + 200

	while abs(ext_list[ext_index+1] - ext_list[ext_index]) > rtol*ext_list[ext_index]:
		N = N + 50
		T_list = np.linspace(TL, TU, N, endpoint=True)
		ext_list = T2ext_list(T_list, N)
		ext_index = np.argmin(ext_list)
		print('N_samples\t', 'ext_list\t', 'T_list\t')
		print(N, ext_list[ext_index], T_list[ext_index])
		if N > 400:
			print('reach maxium of N')
			break

	# TL = T_ext0 - DT
	# TU = T_ext0 + DT

	return ext_list[ext_index], T_list[ext_index]

def ext_init():

	ext,T_ext = get_ext(Ta, Tb, Y_equi, Omega_init)

	t = [0, 200*ext]	
	gas.TP = T_ext, p
	sol = odeint(PSR_ode, Y_equi, t, rtol=1.e-8)
	y0 = sol[-1,:]
	ext0 = gas2ext()

	return y0, ext0

# def get_T_ext0():

# 	omega0 = 20
# 	gas.TP = T_ext0, p
# 	sol = odeint(PSR_ode, Y_equi, t, rtol=1.e-14)
# 	HRR = np.dot(gas.net_production_rates, gas.partial_molar_enthalpies)
# 	omega = -HRR/gas.cp_mass/gas.density/( gas.T - T_in )
# 	y0 = sol[-1,:]
# 	print(T_ext0, 1./omega)		

if __name__ == '__main__':

	y0, ext0 = ext_init()
	print('original mech')
	print(ext0, T_ext0)


	# ext_list = T2ext_list(T_list, N)

	# sys.exit()
	uq_input = np.loadtxt('data/samples.txt')
	index = np.loadtxt('data/samples_index.txt')	
	N_sample = uq_input.shape[0]

	# N_sample = 10
	EXT_list = np.zeros([N_sample,2])
	for i in range(0,N_sample):
		gas.set_multiplier(1.0) # reset all multipliers
		factor = uq_input[i,:]

		for k in range(len(index)):
			gas.set_multiplier(factor[k],index[k]-1) #reaction index start from 1
		
		ext,T_ext = get_ext(Ta, Tb, y0, 1./(ext0*1E2))
		EXT_list[i, :] = [ext, T_ext]
		print(i, ext, T_ext)
		
	np.savetxt( "data/samples_out_ext.txt", EXT_list )