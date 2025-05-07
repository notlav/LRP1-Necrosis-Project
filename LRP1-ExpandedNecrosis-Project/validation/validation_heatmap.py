# Lavie Ngo
# May 2025
# Validation heatmap for LRP1 paper

#%% Import statements
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns
from IPython.display import display
import matplotlib.colors as mcolors
import matplotlib.cm as cm

#%% Read results chart from auto validation script 
validation = pd.read_csv('/Volumes/saucermanlab/Lavie/LRP1/LRP1 Project Code/LRP1-ExpandedNecrosis-Project/validation/validationResults.csv')


# %% Take 'Measurement' column and create a new column mapping to numerical values 
lit = np.zeros(len(validation))
model = np.zeros(len(validation))

# Literature results
for i in range(0,len(lit)):
    if validation['Measurement'][i] == 'Increase':
        lit[i] = 1
    elif validation['Measurement'][i] == 'Decrease':
        lit[i] = -1 
    elif validation['Measurement'][i] == 'No Change':
        lit[i] = 0

# Model results (prediction column)
for i in range(0,len(model)):
    if validation['Prediction'][i] == 'Increase':
        model[i] = 1
    elif validation['Prediction'][i] == 'Decrease':
        model[i] = -1 
    elif validation['Prediction'][i] == 'No Change':
        model[i] = 0

# add new column
validation['Literature Number Value'] = lit
validation['Model Number Value'] = model

results = {'Output': validation['Output'], 'Model': model, 'Experiments': lit}

validation_heatmap = pd.DataFrame(results).set_index('Output')
# sort 
validation_heatmap = validation_heatmap.sort_values('Model')

# select outputs 
#validation_heatmap_sub = validation_heatmap.loc[['Akt', 'cas3', 'apoptosis']]     
# %% Plot heatmap

vcenter = 0
vmin, vmax = -1, 1 
normalize = mcolors.TwoSlopeNorm(vcenter=vcenter, vmin=vmin, vmax=vmax)

colors = ["#2171B5", "#D73027"]
colormap = sns.color_palette(colors)

fig1, ax1 = plt.subplots()
fig1.set_size_inches(4, 8)
ax1 = sns.heatmap((validation_heatmap), linewidths=0.5, linecolor = 'black', norm=normalize, cmap=colormap)


# hide colorbar
cbar = ax1.collections[0].colorbar
cbar.remove()

font = {'fontname':'Arial'}
#plt.title("Sensitivity Analysis", fontsize = 40, **font)
plt.xlabel("", fontsize = 30, **font)
plt.ylabel("", fontsize=30, **font)
plt.yticks(rotation = 0, fontsize = 28, **font)
plt.xticks(fontsize = 24, rotation = 45, **font)

plt.show()

# %% export heatmap 
fig1.savefig('/Volumes/Saucermanlab/Lavie/Papers/Paper1_LRP1/Figures/validation_heatmap.svg', dpi=300, bbox_inches='tight')

# %%
