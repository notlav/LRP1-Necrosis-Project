# Lavie Ngo
# May 2025
# LRP1 paper validation figures 

#%% Import libraries

import numpy as np
from scipy.integrate import ode
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import pandas as pd
import warnings

# ignore warnings
warnings.filterwarnings('ignore')

data = pd.read_excel('/Volumes/saucermanlab/Lavie/LRP1/LRP1 Project Code/LRP1-ExpandedNecrosis-Project/validation/LRP1NecrosisValidation.xlsx')

# import Netflux ODE function and parameters
import LRP1_necrosis_ODEs as model
import LRP1_necrosis_params as modelparams

[speciesIDs, y0, ymax, tau, w, n, EC50] = modelparams.loadParams()

# Simulation parameters
tspan = [0, 50]

# fill rows that have no data with ''
data = data.fillna('')

# select columns
# make sure column names are the same on excel spreadsheet 
input1 = data['Input']
input2 = data['Input 2']
inputCode = data['Input Code']
measurement = data['Measurement']
output = data['Output']
validationIDs = data['ID']
validationTags = data['in-out']
control = data['Prediction']

# Convert species and reaction names to integer values to map name and reaction

# read reaction IDs from model_data
reactionIDs = []

for i in range(len(w)):
    reactionIDs.append(f'r{i+1}')

# convert IDs to integers 
for i in range(len(speciesIDs)):
    # map each speciesIDs to a number
    exec(speciesIDs[i] + ' = ' + str(i))

for i in range(len(reactionIDs)):
    exec(reactionIDs[i] + ' = ' + str(i))

# Converts parentheses to brackets 
# If there's semicolons at the end, drop them
for i in range(len(inputCode)):
    if '(' in inputCode.iloc[i]:
        inputCode.iloc[i] = inputCode.iloc[i].replace('(', '[')
        
    if ')' in inputCode.iloc[i]:
        inputCode.iloc[i] = inputCode.iloc[i].replace(')', ']')

    if inputCode.iloc[i].endswith(';'):
        inputCode.iloc[i] = inputCode.iloc[i][:-1]
    
    if '^' in inputCode.iloc[i]:
        inputCode.iloc[i] = inputCode.iloc[i].replace('^', '**')
        
#  Set validation threshold change
threshold = 0.01
inc = 'Increase'
dec = 'Decrease'
no_change = 'No change'


# Find indices of output species in model
output_index = np.zeros(len(measurement), dtype=int)
for k in range(len(output)):
    output_index[k] = speciesIDs.index(output.iloc[k])

# Final steady state values for each validation input code
y_end = {}
y_end_all = {}

# Change in activity
activity_change = {}

# store prediction strings
prediction = {}

# predicted change
predicted_change = {}

# vector for matches
matches = np.zeros(len(inputCode))
match_strings = []

# loop over each input
# Number of predictions matching references
numMatch = 0
# number of predictions consistent with literature

# dictionary for different validation classes
# validation categories: in-out, in-int, comb-out, KO-int
byClass = {}
tags = ['In-Int', 'Comb-Out', 'In-Out', 'KO-Int']

for i in range(len(inputCode)):
    print('Validation # ', str(i+1), ' of ', str(len(inputCode)))

    # Reset parameters each time
    [speciesIDs, y0, ymax, tau, w, n, EC50] = modelparams.loadParams()

    tspan = [0, 50]

    # Change to float
    ymax = np.array(ymax, dtype = float)

    # run default simulation: LRP1 off or intermediate
    w[0] = 0.5
    w[9] = 0.4

    sol = solve_ivp(model.ODEfunc, tspan,  y0, args=(ymax, tau, w, n, EC50))

    # grab steady state values
    y_init = sol.y.T[-1,:]
    y0 = sol.y.T[-1,:]

    # Run next step of simulation
    # Evaluate validation conditions
    exec(inputCode.iloc[i])

    sol2 = solve_ivp(model.ODEfunc, tspan,  y0, args=(ymax, tau, w, n, EC50))
    
    # grab final steady state values
    y_end = sol2.y.T[-1,:] 
    
    # add steady state values for each input code in dictionary
    y_end_all[i] = y_end

    # calculate activity change for each validation, for each species in order
    # activity change is calculated by subtracting the steady state value - control of input
    activity_change = np.real(y_end_all[i][output_index[i]])-np.real(y_init[output_index[i]])
    #print(activity_change)

    # Conditions to compare model predictions with threshold for current input code

    # Model prediction is increase
    if activity_change > threshold: # increase
        prediction[i] = 'Increase'
        predicted_change[i] = float(activity_change)
        # check if measurement value is equal to increase
        if measurement.iloc[i].lower() == inc.lower():
            numMatch += 1
            matches[i] = 1 # if simulation matches measurement, add 1 to list
        else:
            matches[i] = 0
    
    # Model prediction is decrease
    elif activity_change < -threshold:
        prediction[i] = 'Decrease'
        predicted_change[i] = float(activity_change)
        if measurement.iloc[i].lower() == dec.lower():
            numMatch += 1
            matches[i] = 1
        else:
            matches[i] = 0

    # No change 
    else: 
        prediction[i] = 'No change'
        predicted_change[i] = float(activity_change)
        if measurement.iloc[i].lower() == no_change.lower():
            numMatch += 1
            matches[i] = 1
        else:
            matches[i] = 0

# copy matches vector
BMatch = matches.copy()

# find index of each validation tag
for i in range(len(tags)):
        
    ind = np.where(np.array(validationTags.str.lower()) == tags[i].lower())[0]

    byClass[tags[i]] = np.sum(BMatch[ind]) / len(BMatch[ind]) * 100 if len(BMatch[ind]) > 0 else 0


# convert matches to yes or no values
for j in range(len(matches)):
    if matches[j] == 1:
        match_strings.append('yes')
    else: 
        match_strings.append('no')

# Output results chart into dataframe
pd.options.display.float_format = '{:.5f}'.format
resultChart = pd.DataFrame({'ID': validationIDs, 'Input': input1, 'Input 2': input2, 'Output': output, 'Control': control, 'Measurement': measurement, 'Prediction': prediction, 'Predicted Change': predicted_change, 'Match': pd.Series(match_strings), 'Tag': validationTags})

# Calculate percent match (% validation)
percentMatch = numMatch/len(measurement)*100

print('Percent match: ', percentMatch)


#%% export results to csv
resultChart.to_csv('/Volumes/saucermanlab/Lavie/LRP1/LRP1 Project Code/LRP1-ExpandedNecrosis-Project/validation/validationResults.csv', index=False)