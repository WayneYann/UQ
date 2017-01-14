"""
uncertainty uqantification method 
multi-dimensional uniform and normal distribution
"""

import cantera as ct
import numpy as np
import csv as csv
import matplotlib.pyplot as plt

'''
please check these parameters:
size = 4000
sampling_type = 0
files to save 

'''
### 
# inputs:n lists
# output:a txt 
# def append_write_txt

def list_write_txt(name,list_to_write):
    #f=open(name,'w', newline='')
    f=open(name,'w')
    for m in range(len(list_to_write)):        
        f.write(str(list_to_write[m]) + '\n')
    f.close()
    
### Modify the reaction rate coefficients according to the reaction type.
# modify the A factor by A*multiple
    
def modify_specific_reaction(R,multiple):
    
    # Type 1: elementary reaction - Arrhenius form
    #print(R.reaction_type)
    
    if R.reaction_type == 1:
        change_R = ct.ElementaryReaction()
        change_R.reactants = R.reactants
        change_R.products = R.products 
        
        A_change = R.rate.pre_exponential_factor
        b_change = R.rate.temperature_exponent
        E_change = R.rate.activation_energy
        
        change_R.rate = ct.Arrhenius(A_change*multiple[0],b_change,E_change)
        
   # Type 2: ThreeBody reaction - Arrhenius form
    elif R.reaction_type == 2:
        change_R = ct.ThreeBodyReaction()
        change_R.reactants = R.reactants 
        change_R.products = R.products
        change_R.efficiencies = R.efficiencies
        
        A_change=R.rate.pre_exponential_factor
        b_change=R.rate.temperature_exponent
        E_change=R.rate.activation_energy
         
        change_R.rate = ct.Arrhenius(A_change*multiple[0],b_change,E_change)  
        
    # Type 4: Falloff reaction        
    elif R.reaction_type == 4:
        
        # define a new reaction 
        change_R = ct.FalloffReaction()
        change_R.reactants = R.reactants 
        change_R.products = R.products
        change_R.falloff = R.falloff
        change_R.efficiencies = R.efficiencies
        
        # change the reaction rate of high pressure limit
        A_change_high = R.high_rate.pre_exponential_factor
        b_change_high = R.high_rate.temperature_exponent
        E_change_high = R.high_rate.activation_energy
        
        change_R.high_rate = ct.Arrhenius(A_change_high*multiple[1],b_change_high,E_change_high)

        # change the reaction rate of low pressure limit
        A_change_low = R.low_rate.pre_exponential_factor
        b_change_low = R.low_rate.temperature_exponent
        E_change_low = R.low_rate.activation_energy
        
        change_R.low_rate = ct.Arrhenius(A_change_low*multiple[0],b_change_low,E_change_low)
        
    # Type 5: Plog reaction
    elif R.reaction_type == 5:
        change_R = ct.PlogReaction()
        
        change_R.reactants = R.reactants
        change_R.products = R.products 
        
        internal_M=[]
        for index in range(len(R.rates)):

            A_change = R.rates[index][1].pre_exponential_factor
            b_change = R.rates[index][1].temperature_exponent
            E_change = R.rates[index][1].activation_energy
            qwer = (R.rates[index][0],ct.Arrhenius(A_change*multiple[index],b_change,E_change))
            internal_M.append(qwer)

            change_R.rates=internal_M
    
    # other uncommon reaction types
    else:
        raise RaiseTypeError("No reaction type")
    
    return change_R

### reaction type error
class RaiseTypeError(Exception):
    
    def __init__(self, value):
         self.value = value
         
    def __str__(self):
        return repr(self.value)

### the number of all the parameters in a combustion mechanism 
### create a list to record the number of parameters for all the reactions

def parameter_pre_decision(gas):
    
    parameter_list=[]

    for index in range(gas.n_reactions):
        R=gas.reaction(index)
        
        # Type 1: elementary reaction - Arrhenius form
        if R.reaction_type == 1:
           parameter_list.append(1)
            
        # Type 2: ThreeBody reaction - Arrhenius form
        elif R.reaction_type == 2:
            parameter_list.append(1)
            
        # Type 4: Falloff reaction        
        elif R.reaction_type == 4:
            parameter_list.append(2)
            
        # Type 5: Plog reaction
        elif R.reaction_type == 5:
            parameter_list.append(len(R.rates))
            
        # other uncommon reaction types
        else:
            raise RaiseTypeError("No reaction type")
        
    return parameter_list

def sampling_seeds_generator(sampling_type,dimension,size):
    
    # Box-Mueller normal distribution
    if sampling_type == 0:
        
        mean = np.ones(dimension)*0.5
        cov = np.eye(dimension)*(1/6)**2
        
        sampling_seeds = np.random.multivariate_normal(mean, cov, size)
        
        for line_index in range(sampling_seeds.shape[0]):
            for parameter_index in range(sampling_seeds.shape[1]):
                
                if abs(sampling_seeds[line_index,parameter_index]) >= 1:
                    resampling_mark = True
                    break
                else: resampling_mark = False
            
            while resampling_mark == True:
                sampling_seeds[line_index,:] = np.random.normal(0,(1/6)**2,dimension)
                for parameter_index in range(sampling_seeds.shape[1]):
                    
                    if abs(sampling_seeds[line_index,parameter_index]) >= 1:
                        
                        resampling_mark = True
                        break
                    else: resampling_mark = False

    # Box-Mueller uniform distribution
    elif sampling_type == 1:    
        
        sampling_seeds = np.random.uniform(0, 1, size=(size, dimension))
    else:
        raise RaiseTypeError("No sampling type")
        
    return sampling_seeds
        
        
### ignition calculation based on the grad of OH mole fraction
def ignition_dalay_time(tim,OH_fraction):
    
    m=len(tim)
    OH_grad=np.zeros(m)
    
    for n in range(m-1):
        OH_grad[n]=(OH_fraction[n+1]-OH_fraction[n])/(tim[n+1]-tim[n])
        
    max_index=np.argmax(OH_grad)
    
    return tim[max_index]


if __name__ == '__main__':
    
    # gas = ct.Solution('AramcoMech2.0_mech.xml')
    gas = ct.Solution('USC_mech_ver_2.xml')
    #gas = ct.Solution('grimech30.xml')
    
    parameter_list=parameter_pre_decision(gas)       
    
    ### uncertainty factor input 
    #uncertainty_f = open('GRI_uncertainty_my_program.txt', 'r')   
    uncertainty_f = open('USC_sampling_UF_my_program.txt', 'r')   
    
    lines = uncertainty_f.readlines()    
    uncertainty_factor = []
    uncertainty_factor_activity = []
    line_index = 0
    
    for line in lines:
        for single_factor in line.split():            
            
            # inactive reactions --- start with 0
            if single_factor == str(0):
                uncertainty_factor_activity.append(False)
                
            # active reactions --- start with uncertainty factors
            else:                
                uncertainty_factor_activity.append(True)
                uncertainty_factor.append(float(single_factor)) 
            
                if len(line.split()) != parameter_list[line_index]:
                    raise RaiseTypeError("Wrong uncertainty factors")
        
        line_index += 1
    
    active_factors = len(uncertainty_factor)
    all_factors = len(uncertainty_factor_activity)
        
    ### sampling seeds generator. type:
    # 0 for Box-Mueller normal distribution
    # 1 for uniform pseudo-random distribution
    size = 4000
    sampling_type = 0
    
    # only the seeds for those active reactions are generated
    sampling_seeds = sampling_seeds_generator(sampling_type,active_factors,size)
    # sampling_seeds = np.loadtxt('formatUncertainty.txt')
    #np.savetxt('input_free_flame.dat',sampling_seeds)

    ### multiply the seeds by the uncertainty factors
    sampling_multiple=np.ndarray(shape=(size,active_factors))

    for j in range(0,active_factors):
        for i in range(0,size):
            sampling_multiple[i,j]=uncertainty_factor[j]**(sampling_seeds[i,j]*2-1)
    
    pres = ct.one_atm  # pressure [Pa]
    temp = 300.0  # unburned gas temperature [K]
    
    fuel_input = 'iC4H8'
    
    fuel_H_atom = gas.n_atoms(fuel_input,'H')
    fuel_C_atom = gas.n_atoms(fuel_input,'C')
    
    PHIIndex = 1
    fuel_comp = 1
    
    O2_comp = fuel_comp*(fuel_H_atom/4+fuel_C_atom)/PHIIndex
    N2_comp = O2_comp*0.78/0.21
    
    reactants = '{0}:{1},O2:{2},N2:{3}'.format(fuel_input,fuel_comp,O2_comp,N2_comp) 
    
    #initial_grid = [0,0.025,0.05,0.075,0.099,0.1]
    initial_grid = np.linspace(0, 0.01, 10) 
    tol_ss = [1.0e-6, 1.0e-11]  # [rtol atol] for steady-state problem
    tol_ts = [1.0e-6, 1.0e-11]  # [rtol atol] for time stepping
        
    ###  adiabatic flame temperature calculation
    gas.TPX = temp, pres, reactants
    
    # Flame object
    f = ct.FreeFlame(gas, initial_grid)
    f.flame.set_steady_tolerances(default=tol_ss)
    f.flame.set_transient_tolerances(default=tol_ts)
        
    f.set_refine_criteria(ratio=10, slope=0.07, curve=0.14)
    #f.set_refine_criteria(ratio=10, slope=0.05, curve=0.1, prune=0.04)
    f.set_max_jac_age(50, 50)
    f.set_time_step(1.0e-5, [1, 2, 5, 10, 20, 50])

    # Solve with the energy equation enabled
    f.energy_enabled = True
    f.solve(loglevel=1, refine_grid=True)
    
    Su0 = f.u[0]

    ### main-loop 
    for sampling_index in range(size):
    #for sampling_index in range(0,3):   
        #test_reaction_num = 2
        
        ### modify reaction rates
        active_multiple_index = 0
        
        #print(gas.reaction(test_reaction_num).rate.pre_exponential_factor)
        
        for reaction_index in range(gas.n_reactions):
            
            R = gas.reaction(reaction_index)
            
            R_start = sum(parameter_list[0:reaction_index])
            R_end = R_start + parameter_list[reaction_index]

            test_sampling_seeds = []
            multiple = []
            for multiple_index in range(R_start,R_end):
                
                if uncertainty_factor_activity[multiple_index] == True:
                    multiple.append(sampling_multiple[sampling_index,active_multiple_index])
                    test_sampling_seeds.append(sampling_seeds[sampling_index,active_multiple_index])
                    active_multiple_index = active_multiple_index + 1
                    
                # if the reaction is in-active, the multiple factor = 1    
                elif uncertainty_factor_activity[multiple_index] == False:
                    multiple.append(1)
                else:
                    RaiseTypeError('no match')
                    
            change_R = modify_specific_reaction(R,multiple)
            gas.modify_reaction(reaction_index,change_R)
            
            del multiple
            del R
            del change_R
            
        #print(gas.reaction(test_reaction_num).rate.pre_exponential_factor)    
        
        try:
            f.solve(loglevel=0, refine_grid=False)
        except:
            print('refine grid')
            f.solve(loglevel=0, refine_grid=True)
                
        Su = f.u[0]

        free_flame_f=open('output_flame.dat','a')
        free_flame_f.write(str(Su) + '\n')
        free_flame_f.close()   
        
        # reset the value of every A-factor        
        active_multiple_index = 0
        for reaction_index in range(gas.n_reactions):
            
            R = gas.reaction(reaction_index)
            
            R_start = sum(parameter_list[0:reaction_index])
            
            R_end = R_start + parameter_list[reaction_index]

            # if the reaction is in-active, the factor = 1
            multiple = []
            for multiple_index in range(R_start,R_end):
                
                if uncertainty_factor_activity[multiple_index] == True:
                    multiple.append(1/sampling_multiple[sampling_index,active_multiple_index])
                    active_multiple_index += 1
                    
                elif uncertainty_factor_activity[multiple_index] == False:
                    multiple.append(1)
                else:
                    RaiseTypeError('no match')
            
            change_R = modify_specific_reaction(R,multiple)
            gas.modify_reaction(reaction_index,change_R)
           
            del multiple            
            del R
            del change_R
        
        #f.write_csv('free_flame_{}.csv'.format(sampling_index), quiet=False)
            
        #print(gas.reaction(test_reaction_num).rate.pre_exponential_factor)
        
        print('sample:{}   {: 10.3e}'.format(sampling_index, Su))
