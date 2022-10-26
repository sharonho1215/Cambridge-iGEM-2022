from SALib.sample import latin
from SALib.analyze import sobol, delta
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from model import antithetic_model, model_solve
import matplotlib.pyplot as plt

problem = {
  'num_vars': 8,
  'names': ['$\\alpha$', '$\\theta_1$','$\\gamma_p$','$\mu_1$','$\eta$','$\mu_2$','$\\theta_2$','$\gamma_c$'],
  'bounds':  np.column_stack((np.array([10/60 *1.67, 10**-2/60, 0.0003, 10/60, 0.001/60, 10/60, 10**-3 /60, 0.0003])
                              ,np.array([1000/60 *1.67, 10**3/60, 0.05, 10**3/60, 10/60, 1000/60, 1000/60, 0.05])))
            }


# Generate samples
vals = latin.sample(problem, 10000)

# initializing matrix to store output
Y = np.zeros([len(vals),1])

# Run model (example)
# numerically soves the ODE
# output is x, z1, and z2 at the end time step
# could save output for all time steps if desired, but requires more memory
Y = np.zeros([len(vals),3])
for i in range(len(vals)):
    parameters = (
        True,
        True,
        vals[i][0] ,
        vals[i][1] ,
        vals[i][2] ,
        vals[i][3],
        vals[i][4],
        vals[i][5],
        vals[i][6],
        vals[i][7]
        )
    Y[i][:] = model_solve(antithetic_model, [0,0,0], 50000, 1000, parameters)[0][-1]


# completing soboal analysis for each x, z1, and z2
print('\n\n====x Sobol output====\n\n')
Si_x = delta.analyze(problem,vals, Y[:,0], print_to_console=True)
print('\n\n====z1 Sobol output====\n\n')
Si_z1 = delta.analyze(problem, vals, Y[:,1], print_to_console=True)
print('\n\n====z2 Sobol output====\n\n')
Si_z2 = delta.analyze(problem, vals, Y[:,2], print_to_console=True)

Si_x.plot()
Si_z1.plot()
Si_z2.plot()
plt.show()