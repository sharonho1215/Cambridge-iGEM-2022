from scipy.optimize import fsolve,least_squares,  broyden1, newton_krylov
from model import antithetic_model, model_solve
from math import isclose
import numpy as np
import warnings
from joblib import Parallel, delayed
from tqdm import tqdm

def convergence_checker(model, parameters, init_guess):
    
    final_t_conv = 1000000
    t_points_conv = 1000000
    Us, ts = model_solve(model, init_guess, final_t_conv, t_points_conv, parameters)
    Us_x = []
    Us_z1 = []
    Us_z2 = []
    converge_x = []
    converge_z1 = []
    converge_z2 = []
    for k in range(len(Us)):
        Us_x.append(Us[k][0])
        Us_z1.append(Us[k][1])
        Us_z2.append(Us[k][2])
    window_x = []
    window_z1 = []
    window_z2 = []
    n = 50
    for i in range(1,100,1):
        window_x.append(np.mean(Us_x[-(n-1)-i:-i]))
        window_z1.append(np.mean(Us_z1[-(n-1)-i:-i]))
        window_z2.append(np.mean(Us_z2[-(n-1)-i:-i]))
    for i in range(len(window_x)-1):
        if window_x[i] - window_x[i+1] < 0.001:
            converge_x.append(True)
        else:
            converge_x.append(False)
        if window_z1[i] - window_z1[i+1] < 0.001:
            converge_z1.append(True)
        else:
            converge_z1.append(False)
        if window_z2[i] - window_z2[i+1] < 0.001:
            converge_z2.append(True)
        else:
            converge_z2.append(False)
    if any(converge_x) == True and any(converge_z1) == True and any(converge_z2) == True:
        return True
    else:
        return False
17:37
def optimiser_solve(model, init_guess, parameters, get_offset = False, solver = "optimal"):
    stable =  convergence_checker(antithetic_model, parameters, init_guess)
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
def instability_rejection(model, init_guess, parameter_grid, get_offset = False, parallel = True):
    if parallel:
        root_list = Parallel(n_jobs=-1)(delayed(optimiser_solve)(model, init_guess, parameters) for parameters in tqdm(parameter_grid))
    else:
        root_list = [optimiser_solve(model, init_guess, parameters) for parameters in parameter_grid]
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
