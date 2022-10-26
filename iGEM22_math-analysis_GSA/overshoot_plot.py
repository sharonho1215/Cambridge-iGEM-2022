
from model import antithetic_model, model_solve
from analysis import instability_rejection, max_overshoot
import numpy as np
from collections import namedtuple
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns

samples = 100
parameter_grid = []

mu_range = np.linspace(0.01,20,samples)

for sample in range(samples):
    for sample_2 in range(samples):
        alpha = 50*1.67/60 #np.random.uniform(0, 2)
        theta_1 = 0.005 #np.random.uniform(0, 2)
        gamma_p =  0.00038 #np.random.uniform(0, 2)
        eta = 0.018/(1.67*60) #np.random.uniform(0, 2)
        mu_2 = mu_range[sample]
        mu_1 = mu_range[sample_2]
        theta_2 = 0.0005 #np.random.uniform(0, 2)
        gamma_c = gamma_p #np.random.uniform(0, 2)
        

        parameters_sample = (
            True, 
            True, 
            alpha, 
            theta_1, 
            gamma_p, 
            mu_1, 
            eta, 
            mu_2, 
            theta_2, 
            gamma_c
            )
        parameter_grid.append(parameters_sample)

root_grid, parameter_grid, rejection = instability_rejection(antithetic_model, [100,10,10], parameter_grid, parallel = True)

tuple_param_grid = [tuple(parameters) for parameters in parameter_grid]
Us = [model_solve(antithetic_model, [0,0,0], 100000, 10000, parameters)[0] for parameters in tuple_param_grid]

overshoot = np.array([[max_overshoot(Us[i][:,n], root_grid[i][n]) for n in range(3)] for i in range(len(parameter_grid))] )

x = parameter_grid[:,5]
y = parameter_grid[:,7]
X,Y = np.meshgrid(x,y)

species = ['x', 'z_1', 'z_2']
f, ax = plt.subplots(1, np.shape(root_grid)[1], figsize=(15, 5), sharex = True, sharey=True)
jet = cm.get_cmap('jet')

for i in range(np.shape(root_grid)[1]):
    z = overshoot[:,i].reshape(samples**2 - rejection)
    #ax[i].plot(mu_range, mu_range, 'k-.')
    points = ax[i].scatter(x, y, c=z, s=15, cmap=jet)
    ax[i].set_title(f'${species[i]}$', fontsize = 16, weight='bold')
    cb = f.colorbar(points, ax=ax[i], location ='bottom')
    cb.ax.tick_params(labelsize = 12)
    cb.set_label('Rel. overshoot',  size= 12)
ax[0].set_xlim(x.min(), y.max())
ax[0].set_ylim(y.min(), y.max())
ax[0].set_ylabel("$\mu_2$", size= 14, weight='bold')
ax[1].set_xlabel("$\mu_1$", size= 14, weight='bold')

plt.show()
