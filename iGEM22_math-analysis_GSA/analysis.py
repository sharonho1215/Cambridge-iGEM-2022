from scipy.optimize import fsolve,least_squares,  broyden1, newton_krylov
from model import antithetic_model, model_solve
from math import isclose
import numpy as np
import warnings
from joblib import Parallel, delayed
from tqdm import tqdm


# A function to check whether the input parameters fit the stability criteria or not        
def stability_checker(parameters, get_offset = False):
    repression = parameters[0]
    degradation = parameters[1]
    alpha = parameters[2]
    theta_1 = parameters[3]
    gamma_p =  parameters[4]
    mu_1 = parameters[5]
    eta = parameters[6]
    mu_2 = parameters[7]
    theta_2 = parameters[8]
    psi = mu_1 * eta
    phi = gamma_p * mu_1 * theta_1 * (mu_1 - mu_2) / (alpha * mu_1 * theta_2 + gamma_p* (mu_1 - mu_2))
    if repression and not degradation:
        if mu_2 > +mu_1 and alpha> 1.09* gamma_p * (mu_2 - mu_1)/ (mu_1*theta_2):
            return True
        else:
            return False
    if repression and degradation:
        if not get_offset:
            return True
        if get_offset:
            if mu_2 > 1* mu_1 and alpha> 1.01* gamma_p * (mu_2 - mu_1)/ (mu_1*theta_2):
                return True
            else:
                return False
    if not repression and not degradation:
        return True
    if not repression and degradation:
        return True

# Numerical ODE solver using  different optimal methodologies to avoid asymptotically fragile behavior.
def optimiser_solve(model, init_guess, parameters, get_offset = False, solver = "optimal"):
    stable =  stability_checker(parameters, get_offset = False)
    if stable:
        if solver == "least_squares":
            root = least_squares(model, init_guess, bounds = [(0,0,0),(np.inf,np.inf,np.inf)], args = (0,) + parameters).x
        elif solver == "fsolve":
            root = fsolve(model, init_guess, args = (0,) + parameters)
        elif solver == "optimal":
            try:
                warnings.filterwarnings("error")
                root = fsolve(model, init_guess, args = (0,) + parameters)
                
                if np.min(root) < 0:
                    #print("fsolve returned negative roots")
                    root = least_squares(model, init_guess, bounds = [(0,0,0),(np.inf,np.inf,np.inf)], args = (0,) + parameters).x
                else:
                    #print("fsolve success")
                    pass
            except RuntimeWarning: 
                #print("fsolve failed to converge")
                root = least_squares(model, init_guess, bounds = [(0,0,0),(np.inf,np.inf,np.inf)], args = (0,) + parameters).x


    if not stable:
        root = [np.nan, np.nan, np.nan]
        
    
    return root

#A function to selectively reject parameter combinations that result in instability, only accepts stable ones.
def instability_rejection(model, init_guess, parameter_grid, get_offset = False, parallel = True):
    if parallel:
        root_list = Parallel(n_jobs=-1)(delayed(optimiser_solve)(model, init_guess, parameters, get_offset) for parameters in tqdm(parameter_grid))
    else:
        root_list = [optimiser_solve(model, init_guess, parameters, get_offset) for parameters in parameter_grid]
    unwanted_indices = [i for i, x in enumerate(root_list) if np.nan in x]
    
    #remove elements that result in instability
    if len(unwanted_indices) != 0:
        root_grid = np.array([j for i, j in enumerate(root_list) if i not in unwanted_indices])
        parameter_grid = np.array([j for i, j in enumerate(parameter_grid) if i not in unwanted_indices])
        return root_grid, parameter_grid, len(unwanted_indices)
    else:
        parameter_grid = np.array(parameter_grid)
        root_grid = np.array(root_list)
        return root_grid, parameter_grid, 0

# A function to check bound with delta of choice.
def check_bound(a, reference, delta):
    if abs(a-reference) <= delta* reference:
        return True
    else:
        return False

def check_stability(model, init_condition, final_t, time_step, backward_index, delta, parameters):
    stability_list = []
    Us = model_solve(model, init_condition, final_t, time_step, parameters)[0]
    for i in range(np.shape(Us)[1]):
        truth_array=np.array([check_bound(x, Us[:,i][-1], delta) for x in Us[:,i]])
        if truth_array[-backward_index] == truth_array[-1]:
            stability_list.append(True)
        else:
            stability_list.append(False)

    stable = False
    if len(stability_list) > 0 :
        stable = all(elem == True for elem in stability_list)
    else:
        print("Error")
    return stable

#Improved version of optimiser_solve
def optimiser_solve_2(model, init_guess, parameters, get_offset = False, solver = "optimal"):
    stable =  check_stability(model, [0,0,0], 100000, 10000, 4, 10**-6, parameters)
    if stable:
        if solver == "least_squares":
            root = least_squares(model, init_guess, bounds = [(0,0,0),(np.inf,np.inf,np.inf)], args = (0,) + parameters).x
        elif solver == "fsolve":
            root = fsolve(model, init_guess, args = (0,) + parameters)
        elif solver == "optimal":
            try:
                warnings.filterwarnings("error")
                root = fsolve(model, init_guess, args = (0,) + parameters)
                
                if np.min(root) < 0:
                    #print("fsolve returned negative roots")
                    root = least_squares(model, init_guess, bounds = [(0,0,0),(np.inf,np.inf,np.inf)], args = (0,) + parameters).x
                else:
                    #print("fsolve success")
                    pass
            except RuntimeWarning: 
                #print("fsolve failed to converge")
                root = least_squares(model, init_guess, bounds = [(0,0,0),(np.inf,np.inf,np.inf)], args = (0,) + parameters).x


    if not stable:
        root = [np.nan, np.nan, np.nan]
        
    
    return root

#improved version of instability_rejection
def instability_rejection_2(model, init_guess, parameter_grid, get_offset = False, parallel = True):
    if parallel:
        root_list = Parallel(n_jobs=-1)(delayed(optimiser_solve_2)(model, init_guess, parameters, get_offset) for parameters in tqdm(parameter_grid))
    else:
        root_list = [optimiser_solve_2(model, init_guess, parameters, get_offset) for parameters in parameter_grid]
    unwanted_indices = [i for i, x in enumerate(root_list) if np.nan in x]
    
    #remove elements that result in instability
    if len(unwanted_indices) != 0:
        root_grid = np.array([j for i, j in enumerate(root_list) if i not in unwanted_indices])
        parameter_grid = np.array([j for i, j in enumerate(parameter_grid) if i not in unwanted_indices])
        return root_grid, parameter_grid, len(unwanted_indices)
    else:
        parameter_grid = np.array(parameter_grid)
        root_grid = np.array(root_list)
        return root_grid, parameter_grid, 0


#A function that returns the offset or steady state errors between unperturbed and perturbed system
def get_offset(model, init_guess, parameter_grid, beta, get_offset = True):
    offset = []
    # accept each element of parameter_grid into 'parameters' argument
    root, parameters, rejection = instability_rejection(model, init_guess, parameter_grid, get_offset) 
    param_grid_beta = [list(x).append(beta) for x in parameter_grid ]
    print(param_grid_beta) 
    for i in range(len(root)):
        root_x = root[i][0]
        root_z1 = root[i][1]
        root_z2 = root[i][2]


        offset_x = abs(x_eq - root_x)/x_eq
        offset_z1 = abs(z1_eq - root_z1)/z1_eq
        offset_z2 = abs(z2_eq - root_z2)/z2_eq
        offset.append([offset_x, offset_z1, offset_z2])

    offset_grid = np.array(offset)
    parameter_grid = np.array(parameters)
    
    return offset_grid, parameter_grid, rejection


def max_overshoot(soln, steady_state):
    relative_overshoot = abs(np.max(soln) - steady_state)
    return relative_overshoot