from SALib.sample import latin,saltelli
from SALib.analyze import hdmr
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from model import antithetic_model, model_solve
import matplotlib.pyplot as plt
from analysis import instability_rejection_2

problem = {
  'num_vars': 6,
  'names': [ '$\\theta_1$','$\\gamma_p$','$\mu_1$','$\eta$','$\\theta_2$','$\gamma_c$'],
  'bounds':  np.column_stack((np.array([ 10**-2 /60, 0.0003, 10/60, 0.0001/60,10**-2 /60,0.0003])
                              ,np.array([ 100/60, 0.05, 1000 / 60,0.1/60,1000/60,0.005])))
}


# Generate samples
vals = latin.sample(problem, 5000)

# Run model (example)
# numerically soves the ODE
# output is x, z1, and z2 at the end time step
# could save output for all time steps if desired, but requires more memory
Y = np.zeros([len(vals),3])
for i in range(len(vals)):
    parameters = (
        False,
        True,
        0,
        vals[i][0] ,
        vals[i][1] ,
        vals[i][2] ,
        vals[i][3],
        0,
        vals[i][4],
        vals[i][5],
       ) 
    Y[i][:] = model_solve(antithetic_model, [0,0,0], 100000, 100000, parameters)[0][-1]


# completing soboal analysis for each x, z1, and z2
print('\n\n====x Sobol output====\n\n')
Si_x = hdmr.analyze(problem,vals, Y[:,0], print_to_console=True)
print('\n\n====z1 Sobol output====\n\n')
Si_z1 = hdmr.analyze(problem, vals, Y[:,1], print_to_console=True)
print('\n\n====z2 Sobol output====\n\n')
Si_z2 = hdmr.analyze(problem, vals, Y[:,2], print_to_console=True)

Si_x.plot()
Si_z1.plot()
Si_z2.plot()

plt.show()