# Imports
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import csv
import math
import scipy.interpolate
import os
from utilities import sp

# Using pandas read_csv function I was able to decrease the time of the program.
# It now creates a dataframe with column names and converts the values appropriately in each column directy from the incoming text data.
df = pd.read_csv('../tablewithmagsandfactors.txt', delim_whitespace=True, names=[
                 "filename", "type_", "w2_mag", "m2_mag", "w1_mag", "u_mag", "b_mag", "v_mag", "w2_factor", "m2_factor", "w1_factor", "u_factor", "b_factor", "v_factor"])

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
ax.plot(df_w2_v, p1(t1) - p1(df['w2_factor']),
        'o', label='w2-v vs. Residual', markersize=3)
ax.plot(df_m2_v, p1(t1) - p1(df['w2_factor']),
        'o', label='m2-v vs. Residual', markersize=3)
ax.plot(df_m2_w1, p1(t1) - p1(df['w2_factor']),
        'o', label='m2-w1 vs. Residual', markersize=3)
ax.plot(df_w1_v, p1(t1) - p1(df['w2_factor']),
        'o', label='w1-v vs. Residual', markersize=3)
ax.plot(t1, p1(t1), label='poly-fit')

# Plot Settings
ax.set_xlabel("Color")
ax.set_ylabel("Residual (w2_fit_factor - w2_factors")
ax.set_title("Color vs. Residual")
ax.axis = ('equal')
leg = ax.legend()
plt.show()
