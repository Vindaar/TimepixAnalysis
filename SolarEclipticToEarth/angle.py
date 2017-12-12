#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pylab

# set up some LaTeX plotting parameters
# still need to change parameters
fig_width_pt = 478.00812#246.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt  # width in inches
fig_height = fig_width*golden_mean      # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'axes.labelsize':      10,#10,
          'axes.titlesize':      10,
          'font.size':           10,
          'legend.fontsize':     10,#10,
          'xtick.labelsize':     8,#8,
          'ytick.labelsize':     8,#8,
          'text.usetex':         True,
          'text.latex.preamble': [r'\usepackage{siunitx}'],
          'font.family':         'serif',
          'font.serif':          'cm',
          'figure.figsize':      fig_size}
pylab.rcParams.update(params)


daten = open('solar_angle_2014.dat', 'r').readlines()

latitude_2014 = []
days_2014     = []
hours = []

day_iterator  = 0
for i, line in enumerate(daten):
    print i
    if i > 1:
        # split line into columns
        line = line.split()
        # append latitude to list
        latitude_2014.append(float(line[11]))
        
        # now create proper day element
        day = (i-2)*1.0/24.0
        days_2014.append(day)
        hours.append(line[1])

print latitude_2014
print days_2014
print hours

plt.plot(days_2014, latitude_2014)
plt.xlabel('Days from 15/10/2014 to 16/11/2014')
plt.ylabel('Angle to solar equator / deg')
plt.savefig('./angles_2014.pdf')
plt.show()

daten = open('solar_angle_2015.dat', 'r').readlines()

latitude_2015 = []
days_2015     = []
hours = []

day_iterator  = 0
for i, line in enumerate(daten):
    print i
    if i > 1:
        # split line into columns
        line = line.split()
        # append latitude to list
        latitude_2015.append(float(line[11]))
        
        # now create proper day element
        day = (i-2)*1.0/24.0
        days_2015.append(day)
        hours.append(line[1])

print latitude_2015
print days_2015
print hours

plt.plot(days_2015, latitude_2015)
plt.xlabel('Days from 20/06/2015 to 19/11/2015')
plt.ylabel('Angle to solar equator / deg')
plt.savefig('./angles_2015.pdf')
plt.show()
    # plt.plot(order_calc, dist_calc, color='sienna', marker='', linestyle='-')
    # plt.xlabel('Order of Lissajous figure')
    # plt.ylabel('Arbitrary distance $x$ / \si{\cm}')
    # plt.savefig('./Bilder/michelson_plot.eps')

