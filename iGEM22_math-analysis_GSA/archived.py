from scipy.interpolate import griddata
#%%
# define grid.
xi = np.linspace(x.min(),x.max(),samples)
yi = np.linspace(y.min(), y.max(),samples)
# grid the data.
zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
# contour the gridded data, plotting dots at the randomly spaced data points.
CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)
plt.colorbar() # draw colorbar
# plot data points.
#plt.scatter(x,y,marker='o',c='b',s=5)
plt.xlim(x.min(),x.max())
plt.ylim(y.min(), y.max())
plt.xlabel("$\mu_2$")
plt.ylabel("$\mu_1$")
plt.legend()
plt.show()



#root = get_ss(antithetic_model, parameters)
#Us, ts = model_solve(antithetic_model, [0,0,0], 300, 10000, parameters)

##labels = ["$x$", "$z_1$", "$z_2$"]

#plt.plot(ts, Us)
#plt.hlines(root, ts[0], ts[-1])


#plt.legend(labels)
#plt.show()

def check_stability(truth_array, backward_index):
    if truth_array[-backward_index] == truth_array[-1]:
        return True
    else:
        return False

def get_ss(model, init_condition, final_t, time_step, delta, parameters):
    ss_species=[]
    Us = model_solve(model, init_condition, final_t, time_step, parameters)[0]
    
    for i in range(np.shape(Us)[1]):
        truth_array=np.array([check_bound(x, Us[:,i][-1], gamma) for x in Us[:,i]])
       
        stable = True
        stable = check_stability(truth_array, 3)
        if not stable:
            ss_species = ['N/A', 'N/A', 'N/A']
            break
        first_false_index = np.where(truth_array==False)[0][-1]
        ss_array = Us[first_false_index+1:,i]
        ss_mean = np.sum(ss_array)/len(ss_array)
        ss_species.append(ss_mean)
    return ss_species

def J_func(U, t, repression=False, degradation=False, alpha = 1, theta_1 = 1, gamma_p = 1, mu_1 = 1, eta = 1, mu_2 = 2, theta_2 = 1, gamma_c = 1):
        x = U[0]
        z1 = U[1]
        z2 = U[2]
        alpha = parameters[2]
        theta_1 = parameters[3]
        gamma_p =  parameters[4]
        mu_1 = parameters[5]
        eta = parameters[6]
        mu_2 = parameters[7]
        theta_2 = parameters[8]
        gamma_c = parameters[9]
        
        return np.array([
            [-gamma_p, -alpha*theta_1/(1+z1*theta_1)**2, 0],
            [0, -z2*eta - gamma_c, -z1*eta - gamma_c] ,
            [-mu_2*theta_2/(x*theta_2+1)**2 , -z2*eta - gamma_c, -z1*eta - gamma_c ]
        ])