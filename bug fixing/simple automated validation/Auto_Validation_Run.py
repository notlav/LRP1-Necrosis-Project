# Script to run automated validation function
# Last updated by Lavie Ngo 5/9/2024

#%% Import statements
import pandas as pd

# Import validation script
from Auto_Validation_Function import runValidation

# import Netflux ODE function and parameters
import NetfluxODE as model
import NetfluxODE_params as modelparams

# Validation spreadsheet
validation_data = pd.read_excel('/Volumes/SaucermanLab/Lavie/LRP1/LRP1 Project Code/bug fixing/simple automated validation/LRP1_validation_edit.xlsx')

threshold = 0.01 
 
# Call validation function
percentMatch, resultChart = runValidation(model, modelparams, validation_data, threshold)

# print percentMatch and resultChart as dataframe
print('\nPercent Match:', percentMatch)
resultChart

# %%
