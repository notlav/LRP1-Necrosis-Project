# testing matlab exported ODE files

#%% libraries
import numpy as np 
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
import pandas as pd
import seaborn as sns
from IPython.display import display
import matplotlib.colors as mcolors
import matplotlib.cm as cm

# replace with Netflux ODE files 
import LRP1NecrosisMATLABODE as file
import LRP1NecrosisMATLABODE_params as params

# hide warnings
import warnings
warnings.filterwarnings("ignore")

# %% load parameters
[speciesNames, tau, ymax, y0, w, n, EC50] = params.loadParams()

def performSensAnalysis(w1, w2): # args: LRP1ag, ROS

    # Time span
    tspan = [0, 50]

    # Run default simulation to baseline
    sol0 = solve_ivp(file.ODEfunc, tspan,  y0, args=(tau, ymax, w, n, EC50), t_eval=np.linspace(*tspan, 201))
    sol0_SS = sol0.y[:,-1]

    # Run simulation with perturbations to baseline
    # Baseline/control: w_lrp1ag = 0, w_ros = 0.05
    # INDEXES
    idx_lrp1ag = 0
    idx_ros = 9

    w[idx_lrp1ag] = w1
    w[idx_ros] = w2

    sol1 = solve_ivp(file.ODEfunc, tspan,  sol0_SS, args=(tau, ymax, w, n, EC50), t_eval=np.linspace(*tspan, 201))

    # steady state values of baseline simulation
    sol1_SS = sol1.y[:,-1]

    # ODE solutions into dataframe
    sol0_df = pd.DataFrame(sol0.y.T, index = sol0.t, columns = speciesNames).melt(var_name='species', value_name = 'activity', ignore_index= False).reset_index(names='time')
    sol1_df = pd.DataFrame(sol1.y.T, index = sol1.t, columns = speciesNames).melt(var_name='species', value_name = 'activity', ignore_index= False).reset_index(names='time')
                        
    # create new column to label condition
    sol0_df['condition'] = 'base'
    sol1_df['condition'] = f'ROS + LRP1= {w1}'

    # activity levels for each species for each condition
    results = pd.concat([sol0_df, sol1_df])
    display(results) # print df

    # steady state values for each species
    results_ss = results.query('time == @tspan[1]').reset_index().drop(columns=['time', 'index'])
    display(results_ss)

    #create sens matrix
    sens = np.zeros([len(ymax), len(ymax)])

    deltaP = -1 # complete knockdown

    # simulating knockdowns of each species, one by one
    for i in range(0, len(ymax)):
        print('Knockdown #', str(i+1), 'of', str(len(ymax)))
        ymax_kd = ymax.copy()
        # ymax_new = new maximum values of each species AFTER KNOCKDOWN
        # ymax_new should all be 0's since we are doing a complete knockdown
        ymax_kd[i] = (1+deltaP)*ymax[i]

        # KD with perturbation
        w[idx_lrp1ag] = w1
        w[idx_ros] = w2

        sol1_kd = solve_ivp(file.ODEfunc, tspan, sol1_SS, args=(tau, ymax_kd, w, n, EC50))
        sol1_kd_SS = sol1_kd.y[:,-1]

        # change in activity for each species (KD at steady state - baseline at steady state)
        sens[:,i] = sol1_kd_SS-sol1_SS
        print("Calculating change in activity of ", str(i+1))


    # if there are any non-real numbers, change them to 0 
    for i in range(0,len(ymax)):
        for j in range(0, len(ymax)):
            if np.isnan(np.real(sens[i,j]==1)):
                sens[i,j] == 0

    
    sens_df = pd.DataFrame(sens)
    sens_df = sens_df.set_index([speciesNames]).set_axis([speciesNames], axis = 1)

    print(sens_df)
    return sens_df

#%% run sensitivity analysis 
# w_lrp1 = 1, w_ros = 0.4 (simulating MI) - LRP1 + ROS
sens_df = performSensAnalysis(1, 0.4)

#%% control condition: w_lrp1 = 0, w_ros = 0.4 - ROS only
sens_df_control = performSensAnalysis(0, 0.4)

#%%
sums_lrp1 = sensNetworkInfluence(sens_df)
ranked_plot = plotRankedSens(sens_df, sums_lrp1['top10_rows'], sums_lrp1['top10_cols'])

sens_df.to_csv(f'/Volumes/SaucermanLab/Lavie/LRP1/Code/LRP1 Upstream Model/sens_lrp1_{w1}.csv')