from cmath import pi
from math import gamma
from random import betavariate
import numpy as np
from scipy.integrate import odeint
from sympy import *


def antithetic_model(U, t, repression=False, degradation=False, alpha = 1, theta_1 = 1, gamma_p = 1, mu_1 = 1, eta = 1, mu_2 = 2, theta_2 = 1, gamma_c = 1, beta =0):
    
    x = U[0]
    z1 = U[1]
    z2 = U[2]
    
    if not degradation:
        gamma_c = 0

    if repression:
        dx = alpha / (theta_1 * z1 + 1) - gamma_p * x + beta
        dz1 = mu_1 - eta*z1*z2 - gamma_c * z1 
        dz2 = mu_2/(theta_2 * x + 1) - eta * z1 * z2 - gamma_c * z2
    if not repression:
        alpha = 0
        mu_2 = 0
        dx = theta_1 * z1  - gamma_p * x 
        dz1 = mu_1 - eta*z1*z2 - gamma_c * z1 
        dz2 = theta_2 * x - eta * z1 * z2 - gamma_c * z2   

    return dx, dz1, dz2


def model_solve(antithetic_model, initial_cond, final_t, t_points, parameters):

    U0 = initial_cond
    ts = np.linspace(0, final_t, t_points)

    Us = odeint(antithetic_model, U0, ts, args = parameters)
    return Us, ts