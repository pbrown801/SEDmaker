# Imports
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import math
import scipy.interpolate
import os
from utilities import sp

# Getting each line in file
lines_list = list(
    map(sp, open('../tablewithmagsandfactors.txt', 'r').readlines()))
# Separating into lists
filename = []
type_ = []
w2_mag = []
m2_mag = []
w1_mag = []
u_mag = []
b_mag = []
v_mag = []
w2_factor = []
m2_factor = []
w1_factor = []
u_factor = []
b_factor = []
v_factor = []

for l in lines_list:
    filename.append(l[0])
    type_.append(l[1])
    w2_mag.append(l[2])
    m2_mag.append(l[3])
    w1_mag.append(l[4])
    u_mag.append(l[5])
    b_mag.append(l[6])
    v_mag.append(l[7])
    w2_factor.append(l[8])
    m2_factor.append(l[9])
    w1_factor.append(l[10])
    u_factor.append(l[11])
    b_factor.append(l[12])
    v_factor.append(l[13])

# Convert lists into a dictionary
dict_vals = {'filename': filename, 'type': type_, 'w2_mag': w2_mag, 'm2_mag': m2_mag, 'w1_mag': w1_mag, 'u_mag': u_mag, 'b_mag': b_mag, 'v_mag': v_mag,
             'w2_factor': w2_factor, 'm2_factor': m2_factor, 'w1_factor': w1_factor, 'u_factor': u_factor, 'b_factor': b_factor, 'v_factor': v_factor}

# Convert dictionary into dataframe
df = pd.DataFrame(dict_vals)

# Convert magnitude and factors from str into float
df['w2_mag'] = df['w2_mag'].astype('float')
df['m2_mag'] = df['m2_mag'].astype('float')
df['w1_mag'] = df['w1_mag'].astype('float')
df['u_mag'] = df['u_mag'].astype('float')
df['b_mag'] = df['b_mag'].astype('float')
df['v_mag'] = df['v_mag'].astype('float')
df['w2_factor'] = df['w2_factor'].astype('float')
df['m2_factor'] = df['m2_factor'].astype('float')
df['w1_factor'] = df['w1_factor'].astype('float')
df['u_factor'] = df['u_factor'].astype('float')
df['b_factor'] = df['b_factor'].astype('float')
df['v_factor'] = df['v_factor'].astype('float')

# w2 - v mag
df_w2_v = df['w2_mag']-df['v_mag']
# m2 - v mag
df_m2_v = df['m2_mag']-df['v_mag']
# m2-w1 mag
df_m2_w1 = df['m2_mag']-df['w1_mag']
# w1 - v mag
df_w1_v = df['w1_mag']-df['v_mag']

p1 = np.poly1d(np.polyfit(df_w2_v.values, df['w2_factor'].values, 2))
t1 = np.linspace(min(df_w2_v.values), max(df_w2_v.values), len(df_w2_v.values))

# PLOT
fig, ax = plt.subplots()
ax.plot(df_w2_v, df['w2_factor'], 'o', label='w2-v vs w2_factor', markersize=3)
ax.plot(df_w2_v, p1(t1) - p1(df['w2_factor']), 'o', label='w2-v vs. Residual', markersize=3)
ax.plot(df_m2_v, p1(t1) - p1(df['w2_factor']), 'o', label='m2-v vs. Residual', markersize=3)
ax.plot(df_m2_w1, p1(t1) - p1(df['w2_factor']), 'o', label='m2-w1 vs. Residual', markersize=3)
ax.plot(df_w1_v, p1(t1) - p1(df['w2_factor']), 'o', label='w1-v vs. Residual', markersize=3)
ax.plot(t1,p1(t1), label = 'poly-fit')

# Plot Settings
ax.set_xlabel("Color")
ax.set_ylabel("Residual (w2_fit_factor - w2_factors")
ax.set_title("Color vs. Residual")
ax.axis = ('equal')
leg=ax.legend()
plt.show()
