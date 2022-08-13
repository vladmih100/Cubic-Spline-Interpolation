# EXERCISE: Analysis of Net Mass Changes.
#
# Earthquakes are sometimes associated with oil and gas production (taking mass out of the
# ground) and injection (putting it back in) operations.
#
# It has been suggested that injection of water at one particular site, which started midway through
# 1993, has been responsible for a spate of recent earthquakes there. These are the data that were
# plotted in the first exercise of sdlab_practice.py. The operator of the field has claimed they
# cannot be responsible, because injection had been ongoing for almost 10 years before any earthquakes
# occurred.
#
# It has been proposed the earthquakes may be related to NET MASS CHANGES in the field. Therefore,
# it is necessary to understand how this quantity has evolved over time.
#
# Although data from the two production wells (mass extractors) - PW1 and PW2 - are reported regularly,
# data reporting from the injection well, IW1, is more irregular. In addition, the operator only
# reports MASS RATES, not CUMULATIVE production or injection MASS.
#
# TO solve this problem, you will need to use both INTERPOLATION and INTEGRATION.


# TO DO:
# - In sdlab_functions.py, COMPLETE the functions SPLINE_COEFFICIENT_MATRIX, SPLINE_RHS, and
#   SPLINE_INTERPOLATE.
# - Write a Newton-Cotes integration function OR find and use a built-in Python function.
# - Produce a plot of NET MASS CHANGE as a function of time.
# - ANSWER the questions in sdlab_questions.txt

from matplotlib import pyplot as plt  # MATPLOTLIB is THE plotting module for Python
from numpy.linalg import solve
from sdlab_functions import *

# Extracts data from data files
tm, pm1 = np.genfromtxt('PW1.dat', delimiter=',', skip_header=1).T
tq, pm2 = np.genfromtxt('PW2.dat', delimiter=',', skip_header=1).T
ty, iy = np.genfromtxt('IW1.dat', delimiter=',', skip_header=1).T

# Interpolates data for missing data points in PW2 file
A1 = spline_coefficient_matrix(tq)
b1 = spline_rhs(tq, pm2)
ak1 = solve(A1, b1)
pm3 = spline_interpolate(tm, tq, ak1)

# Interpolates data for missing data points in IW1 file
A2 = spline_coefficient_matrix(ty)
b2 = spline_rhs(ty, iy)
ak2 = solve(A2, b2)
iy2 = spline_interpolate(tm, ty, ak2)

#Combines all mass change rates to obtain NET mass change rate
yCom = iy2-pm1-pm3
#Converts time values from year to seconds since start (start = 0 seconds)
tSec = [((t-1980.0417)*31536000) for t in tm]

#Integrates function by using the trapezium rule to obtain net mass change values
cumMass = integrate(tSec, yCom)

#Interpolates the net mass change values corresponding to the times of the earthquakes for plot markers
quakes = [2003.5, 2004.5, 2005]
A3 = spline_coefficient_matrix(tm)
b3 = spline_rhs(tm, cumMass)
ak3 = solve(A3, b3)
iy3 = spline_interpolate(quakes, tm, ak3)

#Creates 1x1 sub plot window
f, ax1 = plt.subplots(nrows=1, ncols=1)

#Plots Net Mass Change over time
ax1.plot(tm, cumMass, 'k')

#Adds markers corresponding to the earthquakes
ax1.plot(quakes[0], iy3[0], 'b*', markersize=9, label='Earthquake 1: M 3.5')
ax1.plot(quakes[1], iy3[1], 'r*', markersize=9, label='Earthquake 2: M 4.0')
ax1.plot(quakes[2], iy3[2], 'g*', markersize=9, label='Earthquake 3: M 4.3')

#Plots horizontal line at y=0
plt.axhline(y=0, color='k', linestyle='--')

#Adds appropriate labels and titles
ax1.set_xlabel('Time [Yr]')
ax1.set_ylabel('Net Mass Change [Kg]')
ax1.set_title('Net Mass Change over Time ')

#Adds legend
ax1.legend()

#Shows plot or saves image of plot depending on save_figure boolean input
save_figure = False
if not save_figure:
    plt.show()
else:
    plt.savefig('sdlab_earthquakes.png', dpi=300)
