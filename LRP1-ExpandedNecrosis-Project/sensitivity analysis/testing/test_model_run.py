# test_mode_run.py
# Automatically generated by Netflux on 2025-05-07
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import test_mode_ODEs
import test_mode_params

speciesNames, y0, ymax, tau, w, n, EC50 = test_mode_params.loadParams()

# Run single simulation
tspan = [0, 10]
solution = solve_ivp(test_mode_ODEs.ODEfunc, tspan, y0, rtol=1e-8, args=(ymax, tau, w, n, EC50))

fig, ax = plt.subplots()
ax.plot(solution.t,solution.y.T)
ax.set(xlabel='Time',ylabel='Normalized activity')
ax.legend(speciesNames)
plt.show()