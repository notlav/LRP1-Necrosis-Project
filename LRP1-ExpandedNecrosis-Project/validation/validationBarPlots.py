# Create bar plots from validation results - model vs. experiment
# for LRP1 paper 
# updated May 2025
# author: Lavie
#%% import statements
import numpy as np 
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import LRP1_necrosis_ODEs as model
import LRP1_necrosis_params as params

#%% run model simulation

# LRP1 low and high conditions
LRP1_low = 0.2 # intermediate
LRP1_high = 1

# reaction indices for LRP1 and ROS
idx_LRP1ag = 0
idx_ROS = 9

[speciesIDs, y0, ymax, tau, w, n, EC50] = params.loadParams()

tspan = [0, 50]

# simulate MI
w[idx_ROS] = 0.4

# simulate low LRP1ag
w[idx_LRP1ag] = LRP1_low
sol_low= solve_ivp(model.ODEfunc, tspan, y0, args=(ymax, tau, w, n, EC50), t_eval=np.linspace(*tspan, 201))
lowLRP1_SS = sol_low.y[:,-1]

# simulate high LRP1
w[idx_LRP1ag] = LRP1_high
sol_high = solve_ivp(model.ODEfunc, tspan, lowLRP1_SS, args=(ymax, tau, w, n, EC50), t_eval=np.linspace(*tspan, 201))
highLRP1_SS = sol_high.y[:,-1]

#%% create dataframes for steady state values for each species
lowLRP1_df = pd.DataFrame(sol_low.y.T, index = sol_low.t, columns = speciesIDs).melt(var_name='species', value_name = 'activity', ignore_index= False).reset_index(names='time')

highLRP1_df = pd.DataFrame(sol_high.y.T, index = sol_high.t, columns = speciesIDs).melt(var_name='species', value_name = 'activity', ignore_index= False).reset_index(names='time')

# grab steady state values for each species (time = 50)
# filter time = 50
lowLRP1_ss = lowLRP1_df.query('time == @tspan[1]').reset_index().drop(columns=['time', 'index'])

highLRP1_ss = highLRP1_df.query('time == @tspan[1]').reset_index().drop(columns=['time', 'index'])

# combine low and high LRP1 steady state values and add condition column
lowLRP1_ss['condition'] = 'control'
highLRP1_ss['condition'] = 'LRP1\nagonist'

steady_long = pd.concat([lowLRP1_ss, highLRP1_ss])


#%% plot activity for each output
# plot all species in one figure 
validation_species = ('ERK12', 'PI3K', 'Akt', 'NFkB', 'Bax', 'Bcl2', 'cas3', 'apoptosis')

fig, ax = plt.subplots(2,4, figsize = (20, 10), sharey = True)

palette = sns.color_palette('crest', n_colors=8)

# plot each species in each subplot
for i, species in enumerate(validation_species):
	sns.barplot(
		data = steady_long[steady_long['species'] == species],
		x = 'condition',
		y = 'activity',
		color= palette[i],
		ax = ax[i//4, i%4]
	)
	ax[i//4, i%4].set_title(f'{species} activity', fontsize = 16)
	ax[i//4, i%4].set_xlabel('')
	ax[i//4, i%4].set_ylabel('activity')
	# set font to arial
	plt.rcParams['font.family'] = 'Arial'

#plt.savefig('/Volumes/SaucermanLab/Lavie/LRP1/Figures/' + 'validation_activity_gradient.svg')

#%% plot activity for each output under different ROS conditions
fig, axs = plt.subplots()

# loop through each species and plot activity for each condition
# show each plot separately
# add title
palette = sns.color_palette('crest', n_colors=8)

for i, species in enumerate(validation_species):
	plt.figure(figsize = (3, 3))

	sns.barplot(
		data = steady_long[steady_long['species'] == species],
		x = 'condition',
		y = 'activity',
		color= palette[i],
		width = 0.5,
	)

	plt.xticks(fontsize = 14)
	plt.yticks(fontsize = 14)
	plt.ylabel('activity', fontsize = 16)
	plt.xlabel('', fontsize = 16)

	plt.title(f'{species} activity', fontsize = 16)
	# set font to arial
	plt.rcParams['font.family'] = 'Arial'

	# export to svg
	# change file name each time
	# include file path
	#plt.savefig('/Volumes/SaucermanLab/Lavie/LRP1/Figures/' + f'{species}_activity.svg')

	plt.show()

# %% manually plot experimental cas3 data
exp_cas3 = pd.DataFrame({'condition': ['control', 'LRP1\nagonist'], 'activity': [2.05, 0.57]})
model_cas3 = steady_long[steady_long['species'] == 'cas3']
#
# plot subplot for exp and model cas3
fig, axes = plt.subplots(1, 2, figsize = (8, 3),)
plt.suptitle('caspase 3 activity', fontsize = 16)

#color_palette = sns.color_palette('crest', n_colors=8)

# model prediction
sns.barplot(model_cas3, x = 'condition', y = 'activity', width = 0.5, color = 'gray', ax = axes[0])
axes[0].set_xlabel('Model', fontsize = 16)
axes[0].set_ylabel('activity', fontsize = 12)
axes[0].set_yticks(np.arange(0, 1.2, 0.2))
axes[0].set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

# experimental cas3 data
sns.barplot(exp_cas3, x = 'condition', y = 'activity', yerr = [0.495, 0.4], width = 0.5, color = 'gray', ax = axes[1])
axes[1].set_xlabel('Experiment', fontsize = 16)
axes[1].set_ylabel('-fold increase (vs. sham)', fontsize = 12)
axes[1].set_yticks(np.arange(0, 3, 0.5))
axes[1].set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])


#%%  export cas3 figure
fig.savefig('/Volumes/saucermanlab/Lavie/Papers/Paper1_LRP1/Figures/' + 'cas3ValidationBarplot.svg', bbox_inches = 'tight')


# %% side by side for apoptosis
exp_apoptosis = pd.DataFrame({'condition': ['control', 'LRP1\nagonist'], 'activity': [310.73, 54.83]})
model_apoptosis = steady_long[steady_long['species'] == 'apoptosis']


fig, axes = plt.subplots(1, 2, figsize = (8, 3),)
plt.suptitle('apoptosis', fontsize = 16)

#color_palette = sns.color_palette('crest', n_colors=8)

# model prediction
sns.barplot(model_apoptosis, x = 'condition', y = 'activity', width = 0.5, color = 'gray', ax = axes[0])
axes[0].set_xlabel('Model', fontsize = 16)
axes[0].set_ylabel('activity', fontsize = 12)
axes[0].set_yticks(np.arange(0, 1.2, 0.2))
axes[0].set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

sns.barplot(exp_apoptosis, x = 'condition', y = 'activity', yerr = [36.56, 19.3], width = 0.5, color = 'gray', ax = axes[1])
axes[1].set_xlabel('Experiment', fontsize = 16)
axes[1].set_ylabel('TUNEL+ %cardiomyocytes', fontsize = 12)
axes[1].set_yticks(np.arange(0, 351, 70))  
axes[1].set_yticklabels([f'{x:.1f}' for x in np.arange(0.0, 1.1, 0.2)])  


#%%  export apoptosis figure
fig.savefig('/Volumes/saucermanlab/Lavie/Papers/Paper1_LRP1/Figures/' + 'apoptosisValidationBarplot.svg', bbox_inches = 'tight')


#%% side by side for AKT phosphorylation

exp_akt = pd.DataFrame({'condition': ['control', 'LRP1\nagonist'], 'activity': [64.86, 184.31]})
model_akt = steady_long[steady_long['species'] == 'Akt']

fig, axes = plt.subplots(1, 2, figsize = (8, 3),)
plt.suptitle('Akt phosphorylation', fontsize = 16)

#color_palette = sns.color_palette('crest', n_colors=8)

# model prediction
sns.barplot(model_akt, x = 'condition', y = 'activity', width = 0.5, color = 'gray', ax = axes[0])
axes[0].set_xlabel('Model', fontsize = 16)
axes[0].set_ylabel('activity', fontsize = 12)
axes[0].set_yticks(np.arange(0, 1.2, 0.2))
axes[0].set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])

sns.barplot(exp_akt, x = 'condition', y = 'activity', yerr = [24.53, 69.61], width = 0.5, color = 'gray', ax = axes[1])
axes[1].set_xlabel('Experiment', fontsize = 16)
axes[1].set_ylabel('-fold increase', fontsize = 12)
axes[1].set_yticks(np.arange(0, 351, 70))  
axes[1].set_yticklabels([f'{x:.1f}' for x in np.arange(0.0, 1.1, 0.2)])  


#%%  export AKT figure
fig.savefig('/Volumes/saucermanlab/Lavie/Papers/Paper1_LRP1/Figures/' + 'AktValidationBarplot.svg', bbox_inches = 'tight')

# %%
