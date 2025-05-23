# LRP1_necrosis_model_revised_run.py
# Automatically generated by Netflux on 2025-05-06
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import LRP1_necrosis_model_revised_ODEs
import LRP1_necrosis_model_revised_params

speciesNames, y0, ymax, tau, w, n, EC50 = LRP1_necrosis_model_revised_params.loadParams()

# Run single simulation
tspan = [0, 10]
solution = solve_ivp(LRP1_necrosis_model_revised_ODEs.ODEfunc, tspan, y0, rtol=1e-8, args=(ymax, tau, w, n, EC50))

fig, ax = plt.subplots()
ax.plot(solution.t,solution.y.T)
ax.set(xlabel='Time',ylabel='Normalized activity')
ax.legend(speciesNames)
plt.show()
