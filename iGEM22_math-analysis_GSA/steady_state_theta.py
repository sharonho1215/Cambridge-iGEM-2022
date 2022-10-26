from model import antithetic_model, model_solve
from analysis import instability_rejection
import numpy as np
from collections import namedtuple
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns


param_tuple = namedtuple('parameters', ["repression", "degradation", 'alpha', 'theta_1', 'gamma_p', 'mu_1', "eta", "mu_2", "theta_2", "gamma_c"])

samples = 100
parameter_grid = []

theta_range = np.linspace(0.00005,0.05,samples)

for sample in range(samples):
    for sample_2 in range(samples):
        alpha = 50*1.67/60 #np.random.uniform(0, 2)
        theta_1 = theta_range[sample]#np.random.uniform(0, 2)
        gamma_p =  0.00038 #np.random.uniform(0, 2)
        eta = 0.018/(60*1.67) #np.random.uniform(0, 2)
        mu_2 = 15
        mu_1 = 0.5
        theta_2 = theta_range[sample_2] #np.random.uniform(0, 2)
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

root_grid, parameter_grid_after, rejection = instability_rejection(antithetic_model, [100,10,10], parameter_grid, get_offset = False)
x = parameter_grid_after[:,3]
y = parameter_grid_after[:,8]
X,Y = np.meshgrid(x,y)

species = ['x', 'z_1', 'z_2']
f, ax = plt.subplots(1, np.shape(root_grid)[1], figsize=(15, 5), sharex = True, sharey=True)
jet = cm.get_cmap('jet')


for i in range(np.shape(root_grid)[1]):
    z = root_grid[:,i].reshape(samples**2 - rejection)
    #ax[i].plot(mu_range, mu_range, 'k-.')
    points = ax[i].scatter(x, y, c=np.log(z), s=15, cmap=jet)
    ax[i].set_title(f'${species[i]}$', fontsize = 16, weight='bold')
    cb = f.colorbar(points, ax=ax[i], location ='bottom')
    cb.ax.tick_params(labelsize = 12)
    cb.set_label('Log(concentration) - $log_{10}$(nM)',  size= 12)
ax[0].set_xlim(x.min(), y.max())
ax[0].set_ylim(y.min(), y.max())
ax[0].set_ylabel("$\\theta_{2}$", size= 14, weight='bold')
ax[1].set_xlabel("$\\theta_{1}$", size= 14, weight='bold')

plt.show()



