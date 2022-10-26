from model import antithetic_model
from analysis import instability_rejection
import numpy as np
from collections import namedtuple
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm


samples = 150
parameter_grid = []
parameter_grid_beta =[]
mu_range = np.geomspace(0.001,1,samples)

for sample in range(samples):
    for sample_2 in range(samples):
        alpha = 50*1.67/60 #np.random.uniform(0, 2)
        theta_1 = 0.005 #np.random.uniform(0, 2)
        gamma_p =  0.00038 #np.random.uniform(0, 2)
        eta = 0.018/(60*1.67) #np.random.uniform(0, 2)
        mu_2 = mu_range[sample]
        mu_1 = mu_range[sample_2]
        theta_2 = 0.0005 #np.random.uniform(0, 2)
        gamma_c = gamma_p #np.random.uniform(0, 2)
        beta = 0

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
            gamma_c,
            beta
            )
        parameter_grid.append(parameters_sample)


root_grid_nobeta, parameter_grid_after, rejection = instability_rejection(antithetic_model, [10,10,10], parameter_grid, get_offset = False)

for sample in range(samples):
    for sample_2 in range(samples):
        alpha = 50*1.67/60 #np.random.uniform(0, 2)
        theta_1 = 0.005 #np.random.uniform(0, 2)
        gamma_p =  0.00038 #np.random.uniform(0, 2)
        eta = 0.018/(60*1.67) #np.random.uniform(0, 2)
        mu_2 = mu_range[sample]
        mu_1 = mu_range[sample_2]
        theta_2 = 0.0005 #np.random.uniform(0, 2)
        gamma_c = gamma_p #np.random.uniform(0, 2)
        beta = 0.1

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
            gamma_c,
            beta
            )
        parameter_grid_beta.append(parameters_sample)

root_grid_beta, parameter_grid_after, rejection = instability_rejection(antithetic_model, [10,10,10], parameter_grid_beta, get_offset = False)

offset_grid = np.array([abs(root_grid_nobeta[i] - root_grid_beta[i])/root_grid_nobeta[i] for i in range(len(root_grid_nobeta))])
print(offset_grid)

x = parameter_grid_after[:,5]
y= parameter_grid_after[:,7]
X,Y = np.meshgrid(x,y)

species = ['x', 'z_1', 'z_2']
f, ax = plt.subplots()
jet = cm.get_cmap('jet')



z = offset_grid[:,0].reshape(samples**2 - rejection)
points = ax.scatter(x, y, c=z, s=15, cmap=jet)
ax.set_xscale('log')
ax.set_yscale('log')
cb = f.colorbar(points, ax=ax, location ='bottom')
cb.ax.tick_params(labelsize = 12)
cb.set_label(r'Relative offset',  size= 12)
ax.set_ylim(y.min(), y.max())
ax.set_xlim(x.min(), 1)
ax.set_ylabel("$\mu_2$", size= 14, weight='bold')
ax.set_xlabel("$\mu_1$", size= 14, weight='bold')

plt.show()
