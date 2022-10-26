from model import antithetic_model, model_solve
from analysis import instability_rejection
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns


samples = 150
parameter_grid = []

theta_range = np.geomspace(0.001,10,samples)
mu_range = np.geomspace(10,1000,samples)

for sample in range(samples):
    for sample_2 in range(samples):
        alpha = 0
        theta_1 = 0.1 #np.random.uniform(0, 2)
        gamma_p =  0.00038 #np.random.uniform(0, 2)
        eta = 0.018/(60*1.67) #np.random.uniform(0, 2)
        mu_1 = mu_range[sample_2]
        mu_2 = 0
        theta_2 = theta_range[sample] #np.random.uniform(0, 2)
        gamma_c = gamma_p #np.random.uniform(0, 2)
        beta = 10/60 *1.67
        

        parameters_sample = (
            False, 
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

root_grid, parameter_grid_after, rejection = instability_rejection(antithetic_model, [10,10,10], parameter_grid, get_offset = False)

x = parameter_grid_after[:,5]
y = parameter_grid_after[:,8]
X,Y = np.meshgrid(x,y)

species = ['x', 'z_1', 'z_2']
f, ax = plt.subplots()
jet = cm.get_cmap('jet')
z1_z2 = np.array([ root_grid[i,1]/root_grid[i,2] for i in range(np.shape(root_grid)[0]) ])
print(z1_z2)
z = z1_z2.reshape(samples**2 - rejection)
points = ax.scatter(x, y, c=z, s=10, cmap=jet)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title(r'\dfrac: $dfrac{z_1}{z_2}$', fontsize = 16, weight='bold')
cb = f.colorbar(points, ax=ax, location ='bottom')
cb.ax.tick_params(labelsize = 12)
cb.set_label('Log(concentration), ln(nM) ',  size= 12)
ax.set_xlim(x.min(), y.max())
ax.set_ylim(y.min(), y.max())
ax.set_ylabel("$\\theta_2$", size= 14, weight='bold')
ax.set_xlabel("$\\mu_1$", size= 14, weight='bold')


plt.show()