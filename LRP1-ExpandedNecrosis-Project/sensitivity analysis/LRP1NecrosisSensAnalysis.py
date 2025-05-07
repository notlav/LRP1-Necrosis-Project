# Lavie Ngo
# May 2025
# Sensitivity analysis for LRP1 CM apoptosis model - paper version
# Figure 4

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
import LRP1_necrosis_ODEs as file
import LRP1_necrosis_params as params

# hide warnings
import warnings
warnings.filterwarnings("ignore")

# %% load parameters
[speciesIDs, y0, ymax, tau, w, n, EC50] = params.loadParams()

def performSensAnalysis(w1, w2): # args: LRP1ag, ROS

    # Time span
    tspan = [0, 50]

    # Run default simulation to baseline
    sol0 = solve_ivp(file.ODEfunc, tspan,  y0, args=(ymax, tau, w, n, EC50), t_eval=np.linspace(*tspan, 201))
    sol0_SS = sol0.y[:,-1]

    # Run simulation with perturbations to baseline
    # Baseline/control: w_lrp1ag = 0, w_ros = 0.05
    # INDEXES
    idx_lrp1ag = 0
    idx_ros = 9

    w[idx_lrp1ag] = w1
    w[idx_ros] = w2

    sol1 = solve_ivp(file.ODEfunc, tspan,  sol0_SS, args=(ymax, tau, w, n, EC50), t_eval=np.linspace(*tspan, 201))

    # steady state values of baseline simulation
    sol1_SS = sol1.y[:,-1]

    # ODE solutions into dataframe
    sol0_df = pd.DataFrame(sol0.y.T, index = sol0.t, columns = speciesIDs).melt(var_name='species', value_name = 'activity', ignore_index= False).reset_index(names='time')
    sol1_df = pd.DataFrame(sol1.y.T, index = sol1.t, columns = speciesIDs).melt(var_name='species', value_name = 'activity', ignore_index= False).reset_index(names='time')
                        
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

        sol1_kd = solve_ivp(file.ODEfunc, tspan, sol1_SS, args=(ymax_kd, tau, w, n, EC50))
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
    sens_df = sens_df.set_index([speciesIDs]).set_axis([speciesIDs], axis = 1)

    print(sens_df)
    return sens_df

#%% Plot sensitivity matrix for all species

def plotSens(sens_df):
    vcenter = 0
    vmin, vmax = -1, 1 #sens_kd.min(), sens_kd.max()
    normalize = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)
    #colormap = cm.RdBu_r
    colormap = 'seismic'

    fig1, ax1 = plt.subplots()
    fig1.set_size_inches(20.5, 14.5)
    ax1 = sns.heatmap(sens_df, norm=normalize, cmap=colormap, xticklabels=speciesIDs, yticklabels=speciesIDs)

    font = {'fontname':'Arial'}
    #plt.title("Sensitivity Analysis", fontsize = 40, **font)
    plt.xlabel("Node Knockdown", fontsize = 30, **font)
    plt.ylabel("Measured Node Activity", fontsize=30, 
    **font)

    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, horizontalalignment='right', fontsize=10, **font)
    ax1.set_yticklabels(ax1.get_yticklabels(), rotation=0, horizontalalignment='right', fontsize=10, **font)

    # change color bar label size and font
    cbar = ax1.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), rotation=0, fontsize=15, **font)

    cbar.ax.yaxis.label.set_size(30)

    return ax1


#%% calculate overall network influence/sensitivity
#Overall network influence analysis 

#Sum columns to calculate total influence of knockdown
# Sum rows to calculate total sensitivity
# Make submatrix with the most sensitive outputs and most influential node knockdowns
def sensNetworkInfluence(sens_df):

    sums = {'rows':[], 'cols':[], 'norm_rows':[], 'norm_cols':[], 'top10_rows':[], 'top10_cols':[]}

    for name, row in sens_df.iterrows():
        sum_row = np.sum(row)
        sum_row_abs = np.abs(sum_row)
        sums['rows'].append(sum_row_abs)


    for name, col in sens_df.items():
        sum_col = np.sum(col)
        sum_col_abs = np.abs(sum_col)
        sums['cols'].append(sum_col_abs)

    # sum of all rows and cols
    sum_all_rows = np.sum(sums['rows'])   
    sum_all_cols = np.sum(sums['cols'])

    # normalize sum of each node by dividing by total sum
    sums['norm_rows'] = sums['rows']/sum_all_rows
    sums['norm_cols'] = sums['cols']/sum_all_cols

    # rank top 10 rows (most sensitive)
    sum_norm_rows_df = pd.DataFrame(sums['norm_rows']).set_index([speciesIDs]).sort_values(by=0, ascending=False)
    # rank top 10 cols (most influential knockdowns)
    sum_norm_cols_df = pd.DataFrame(sums['norm_cols']).set_index([speciesIDs]).sort_values(by=0, ascending=False)
        
    sums['top10_rows']= pd.DataFrame(sum_norm_rows_df[0][0:10]).reset_index(names='Species')
    sums['top10_cols'] = pd.DataFrame(sum_norm_cols_df[0][0:10]).reset_index(names='Species')

    return sums

#%% Plot ranked sensitivity matrix
def plotRankedSens(sens_df, top10_rows, top10_cols):

    # subset matrix
    subset_cols = sens_df[top10_cols['Species']].set_index([speciesIDs])
    top10_matrix = subset_cols.loc[top10_rows['Species']]

    vcenter = 0
    vmin, vmax = -1, 1 #sens_kd.min(), sens_kd.max()
    normalize = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)
    #colormap = cm.RdBu_r
    colormap = 'seismic'

    ax1 = plt.subplots()
    ax1 = sns.heatmap(top10_matrix, norm=normalize, cmap=colormap)
    ax1.tick_params(axis='x', rotation=0, labelsize = 8)

    font = {'fontname':'Arial'}
    plt.xlabel("Node Knockdown", fontsize = 15, **font)
    plt.ylabel("Measured Node Activity", fontsize=15, **font)

    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=0, fontsize=8, **font)
    ax1.set_yticklabels(ax1.get_yticklabels(), rotation=0, horizontalalignment='right', fontsize=10, **font)

    # change color bar label size and font
    cbar = ax1.collections[0].colorbar
    cbar.ax.tick_params(labelsize=20)
    cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), rotation=0, fontsize=8, **font)

    return ax1

#%% run sensitivity analysis 
# w_lrp1 = 1, w_ros = 0.4 (simulating MI) - LRP1 + ROS
sens_df = performSensAnalysis(1, 0.4)

#%% export sensitivity matrix as csv
sens_df.to_csv(f'/Volumes/SaucermanLab/Lavie/LRP1/Code/LRP1 Upstream Model/sens_lrp1.csv')

#%% control condition: w_lrp1 = 0, w_ros = 0.4 - ROS only
sens_df_control = performSensAnalysis(0, 0.4)

#%% plots
sens_plot = plotSens(sens_df)
sens_plot_control = plotSens(sens_df_control)

#%% export plots as svgs
sens_plot.get_figure().savefig(f'/Volumes/SaucermanLab/Lavie/Papers/Paper1_LRP1/Figures/sens_lrp1.svg', format='svg', bbox_inches='tight')
sens_plot_control.get_figure().savefig(f'/Volumes/SaucermanLab/Lavie/Papers/Paper1_LRP1/Figures/sens_lrp1_control.svg', format='svg', bbox_inches='tight')

#%% run overall network influence analysis
sums_lrp1 = sensNetworkInfluence(sens_df)

#%% plot influence
ranked_plot = plotRankedSens(sens_df, sums_lrp1['top10_rows'], sums_lrp1['top10_cols'])

#%% plot control influence, same top 10 rows and cols
ranked_plot_control = plotRankedSens(sens_df_control, sums_lrp1['top10_rows'], sums_lrp1['top10_cols'])

#%% export influence plot as svg
ranked_plot.get_figure().savefig(f'/Volumes/SaucermanLab/Lavie/Papers/Paper1_LRP1/Figures/sens_lrp1_ranked.svg', format='svg', bbox_inches='tight')

#%% export control influence plot as svg
ranked_plot_control.get_figure().savefig(f'/Volumes/SaucermanLab/Lavie/Papers/Paper1_LRP1/Figures/sens_lrp1_control_ranked.svg', format='svg', bbox_inches='tight')

# %% Plot regulators of cell death from sensitivity analysis

sens_celldeath = sens_df[sens_df.index == 'cellDeath'].transpose()

fig1, ax1 = plt.subplots()
fig1.set_size_inches(8,6)

# Filter species that have cellDeath activity > 0.1 or < -0.1
celldeath_nodes = sens_celldeath[(sens_celldeath['cellDeath'] > 0.1) | (sens_celldeath['cellDeath'] < -0.1)]

# sort ranked nodes
ranked_nodes = celldeath_nodes.sort_values(by='cellDeath', ascending=False)

ax1 = sns.barplot(ranked_nodes.transpose(), color = 'black')
# rotate x labels
plt.xticks(rotation=45)

font = {'fontname':'Arial'}
#plt.title("Sensitivity Analysis", fontsize = 40, **font)
plt.xlabel("Node Knockdown", fontsize = 20, **font)
plt.ylabel("Cell Death (Change in Activity)", fontsize=20, **font)

plt.show()

# %% export cell death sensitivity plot as svg
plt.savefig(f'/Volumes/SaucermanLab/Lavie/Papers/Paper1_LRP1/Figures/sensCellDeath.svg', format='svg', dpi=1200, bbox_inches='tight')

