# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 14:02:47 2023

@author: karthik.kumar
"""

import math
import matplotlib.pyplot as plt
import pandas as pd
import tkinter as tk
from tkinter import simpledialog
from scipy.signal import find_peaks
import numpy as np
from tkinter import filedialog as fd
import openpyxl
import sys
from pathlib import Path
sys.path.insert(0, r'\\mg-nas-01\DMS2\Projects\P - Customer Project\P20394 - TechnipFMC – HFP Qualification\DOC61334\3678-006-4-2\04 Analysed Data')
import matplotlib.style as style
style.use('tableau-colorblind10')
from matplotlib.backends.backend_pdf import PdfPages
workbooklocation = fd.askopenfilename()  # Select Excel Sheet
workbook = openpyxl.load_workbook(workbooklocation)

Input_Sheet = workbook["Input Sheet"]
Measurement_Devices = workbook["Measurement Devices"]
Output_Sheet = workbook["Output Sheet"]

Coupon_Number = Input_Sheet['B1'].value
Coupon_ID = Input_Sheet['B3'].value
Coupon_OD = Input_Sheet['B2'].value
Liner_Thickness = Input_Sheet['B16'].value
PEEK_stiffness = Input_Sheet['B17'].value
SP01_gauge_length = Input_Sheet['B5'].value
SP02_gauge_length = Input_Sheet['B5'].value
SP01_starting_length = Input_Sheet['D8'].value
SP02_starting_length = Input_Sheet['D9'].value

#Defining directory for save folder
save_folder = Input_Sheet['B23'].value
dir = Path(save_folder)
pdf_save = Input_Sheet['B24'].value
pdf = PdfPages(pdf_save)

#####################################################################################################################################################################

#Functions register

#Functions for Grid 1, 2, and 3 strains
def Mean_SGC_Grid1(df):
    file_modulus_gauges = df.loc[:,
                          (df.columns.str.contains('SGC')) & (df.columns.str.contains('Grid 1'))]
    Gauge_Average_M_microstrain = file_modulus_gauges.mean(axis=1)
    Gauge_Average_M = Gauge_Average_M_microstrain/10000
    return Gauge_Average_M
def Mean_SGC_Grid2(df):
    file_modulus_gauges = df.loc[:,
                          (df.columns.str.contains('SGC')) & (df.columns.str.contains('Grid 2'))]
    Gauge_Average_M2_microstrain = file_modulus_gauges.mean(axis=1)
    Gauge_Average_M2 = Gauge_Average_M2_microstrain/10000
    return Gauge_Average_M2
def Mean_SGC_Grid3(df):
    file_modulus_gauges = df.loc[:,
                          (df.columns.str.contains('SGC')) & (df.columns.str.contains('Grid 3'))]
    Gauge_Average_M3_microstrain = file_modulus_gauges.mean(axis=1)
    Gauge_Average_M3 = Gauge_Average_M3_microstrain/10000
    return Gauge_Average_M3

#Reading String Pot values from csv file and returning string pot strains
def Stringpot_disp_strain(SP, SP_starting_length, SP_gauge_length):
    sp_diff = SP - SP_starting_length
    sp_strain = (sp_diff/SP_gauge_length)
    return sp_strain

#Calculating cross sectional area
def csa(ID, OD):
    return math.pi * ((OD**2)-(ID**2))/4

#####################################################################################################################################################################

#Reading the Failure data file
df = Input_Sheet['B21'].value
df_failure2 = pd.read_csv(df, thousands=',', header= 0, na_values=["Overflow", "Error", ])   # Reads the csv file
df_failure=df_failure2.iloc[:df_failure2.isna().any(1).idxmax()]

#Averaging the Grid 1, Grid 2 and Grid 3 strains and adding in a new column on the df_failure dataframe.
Gauge_Average_M = Mean_SGC_Grid1(df_failure)
df_failure['Gauge_Average_M'] = Mean_SGC_Grid1(df_failure)
Gauge_Average_M2 = Mean_SGC_Grid1(df_failure)
df_failure['Gauge_Average_M2'] = Mean_SGC_Grid2(df_failure)
Gauge_Average_M3 = Mean_SGC_Grid3(df_failure)
df_failure['Gauge_Average_M3'] = Mean_SGC_Grid3(df_failure)

#Calculating strain from String pots
SP01_raw = df_failure2.loc[:, df_failure2.columns.str.contains('SP01')]
df_failure2['SP01_Strain'] = Stringpot_disp_strain(SP01_raw, SP01_starting_length, SP01_gauge_length)
SP02_raw = df_failure2.loc[:, df_failure2.columns.str.contains('SP02')]
df_failure2['SP02_Strain'] = Stringpot_disp_strain(SP02_raw, SP02_starting_length, SP02_gauge_length)

#Plot between SP strains and time
AverageSPStrainDataframe = df_failure2.loc[:,
                          (df_failure2.columns.str.contains('SP')) & (df_failure2.columns.str.contains('Strain'))]
AverageSPStrainDataframe2 = AverageSPStrainDataframe*100
df_failure2['SPStrainAverage'] = AverageSPStrainDataframe2.mean(axis=1)
plt.figure(1, figsize=(10.17, 6.64), dpi=160)
AverageSPStrainDataframe2.insert(0, 'SP Strain Average', df_failure2['SPStrainAverage'])
AverageSPStrainDataframe2.insert(0, 'Time', df_failure2['Scan #']/(10*60))
AverageSPStrainDataframe2.plot(x='Time')
plt.ylim([-20,1])
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Time (min)', fontsize=14)  # X axis labeling and font size
plt.ylabel('String Pot Strain (%)', fontsize=14)  # Y axis labeling and font size
plt.title('String Pot Axial Strain vs Time', fontsize=14)  # Chart title
plt.legend(fontsize = 12, labels = ['Average SP Strain', 'SP01 Strain', 'SP02 Strain'], loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2)
plt.grid(color = 'lightgrey', ls = '-.', linewidth = 0.25)
plt.savefig(dir/'String Pot Strain vs Time.svg', format = "svg")
pdf.savefig()
plt.show()

#Convert Load from Tonnes to Stress (MPa) i.e., multiplying load cell values by 9.81 and / CSA
CSA = csa(Coupon_ID, Coupon_OD)
Loadintonnes = df_failure.loc[:,(df_failure.columns.str.contains('LC01'))]
df_failure['Stress'] = Loadintonnes*9.81*1000/CSA

#Plot between Load vs SP Strain
fig, ax2 = plt.subplots(figsize=(10.17, 6.64), dpi=160)
Loadintonnes2 = df_failure2.loc[:,(df_failure2.columns.str.contains('LC01'))]
ax2.plot(df_failure2['SPStrainAverage'], Loadintonnes2)
ax2.set_xlabel('Average SP Strain (%)', fontsize=14)  # X axis labeling and font size
ax2.set_ylabel('Load (tonnes)', fontsize=14)  # Y axis labeling and font size
ax2.tick_params(axis='y')
plt.xlim([-20,1])
ax2.set_title('Failure Load vs String Pot Strain', fontsize=14)  # Chart title
ax2.grid(color='lightgrey', ls='-.', linewidth=0.25)
# Set the Y-axis to the right-hand side
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position('right')
# Adjust the plot layout
fig.subplots_adjust(left=0.05, right=0.9, bottom=0.15, top=0.95)
# Save the plot and show
plt.savefig(dir/'Load vs String Pot Strain.svg', format='svg')
pdf.savefig()
plt.show()

# Plot between Stress vs SP Strain
fig, ax2 = plt.subplots(figsize=(10.17, 6.64), dpi=160)
df_failure2['Stress'] = Loadintonnes2 * 9.81 * 1000 / CSA
ax2.plot(df_failure2['SPStrainAverage'], df_failure2['Stress'])
# Find maximum stress point and annotate it with an 'X'
max_stress_idx = df_failure2['Stress'].idxmax()
max_stress_point = (df_failure2.loc[max_stress_idx, 'SPStrainAverage'], df_failure2.loc[max_stress_idx, 'Stress'])
ax2.annotate(f"({max_stress_point[0]:.2f}, {max_stress_point[1]:.2f})", 
             xy=max_stress_point, xycoords='data',
             xytext=(max_stress_point[0], max_stress_point[1] - 20), fontsize=12, ha='center', color = 'red')
ax2.plot(max_stress_point[0], max_stress_point[1], 'x', color='g', markersize = 8)
ax2.axvline(x=max_stress_point[0], color='m', linewidth=0.5, linestyle='--')
ax2.axhline(y=max_stress_point[1], color='m', linewidth=0.5, linestyle='--')
ax2.set_xlabel('Average SP Strain (%)', fontsize=14)
ax2.set_ylabel('Stress (MPa)', fontsize=14)
ax2.tick_params(axis='y')
ax2.set_title('Stress vs String Pot Strain', fontsize=14)
ax2.grid(color='lightgrey', ls='-.', linewidth=0.25)
plt.xlim([-20, 1])
# Set the Y-axis to the right-hand side
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position('right')
# Adjust the plot layout
fig.subplots_adjust(left=0.05, right=0.9, bottom=0.15, top=0.95)
# Save the plot and show
plt.savefig(dir / 'Stress vs Strong Pot Strain.svg', format='svg')
pdf.savefig()
plt.show()

#Plot between Axial strains for the 4 strain gauges and time
df_failure_gauges1 = df_failure.loc[:,
                          (df_failure.columns.str.contains('SGC')) & (df_failure.columns.str.contains('Grid 1')) & (df_failure.columns.str.contains('TRI'))]
df_failure_gauges2 = df_failure_gauges1/10000
plt.figure(1, figsize=(10.17, 6.64), dpi=160)
df_failure['StrainAverage'] = df_failure_gauges2.mean(axis=1)
df_failure_gauges2.insert(0, 'Gauge Strain Average', df_failure['StrainAverage'])
df_failure_gauges2.insert(0, 'Time', df_failure['Scan #']/(10*60))
df_failure_gauges2.plot(x='Time')
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Time (min)', fontsize=14)  # X axis labeling and font size
plt.ylabel('Axial Strain (%)', fontsize=14)  # Y axis labeling and font size
plt.title('Axial Strain vs Time', fontsize=14)  # Chart title
plt.legend(fontsize = 12, labels = ['Average Axial Strain','SGC01-TRI Grid 1', 'SGC02-TRI Grid 1', 'SGC03-TRI Grid 1', 'SGC04-TRI Grid 1'], loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=3)
plt.grid(color = 'lightgrey', ls = '-.', linewidth = 0.25)
plt.savefig(dir/'Axial Strain vs Time.svg', format = "svg")
pdf.savefig()
plt.show()

#Plot between Hoop Strains for the 4 strain gauges and time
df_failure_gauges1 = df_failure.loc[:,
                          (df_failure.columns.str.contains('SGC')) & (df_failure.columns.str.contains('Grid 3')) & (df_failure.columns.str.contains('TRI'))]
df_failure_gauges2 = df_failure_gauges1/10000
df_failure['StrainAverage3'] = df_failure_gauges2.mean(axis=1)
plt.figure(1, figsize=(10.17, 6.64), dpi=160)
df_failure_gauges2.insert(0, 'Gauge Strain Average', df_failure['StrainAverage3'])
df_failure_gauges2.insert(0, 'Time', df_failure['Scan #']/(10*60))
df_failure_gauges2.plot(x='Time')
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Time (min)', fontsize=14)  # X axis labeling and font size
plt.ylabel('Hoop Strain (%)', fontsize=14)  # Y axis labeling and font size
plt.title('Hoop Strain vs Time', fontsize=14)  # Chart title
plt.legend(fontsize = 12, labels = ['Average Hoop Strain', 'SGC01-TRI Grid 3', 'SGC02-TRI Grid 3', 'SGC03-TRI Grid 3', 'SGC04-TRI Grid 3'], loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=3)
plt.grid(color = 'lightgrey', ls = '-.', linewidth = 0.25)
plt.savefig(dir/'Hoop Strain vs Time.svg', format = "svg")
pdf.savefig()
plt.show()

#Plot between Stress and Axial Strains
fig, ax2 = plt.subplots(figsize=(10.17, 6.64), dpi=160)
df_failure_gauges = df_failure.loc[:,
                          (df_failure.columns.str.contains('SGC')) & (df_failure.columns.str.contains('Grid 1')) & (df_failure.columns.str.contains('TRI'))]
df_failure_gauges.insert(0, 'Gauge Strain Average', df_failure['StrainAverage']*10000)
ax2.plot(df_failure_gauges/10000, df_failure['Stress'])
ax2.set_xlabel('Axial Strain (%)', fontsize=14)  # X axis labeling and font size
ax2.set_ylabel('Stress (MPa)', fontsize=14)  # Y axis labeling and font size
ax2.tick_params(axis='y')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax2.set_title('Stress vs Axial Strain', fontsize=14)  # Chart title
ax2.grid(color='lightgrey', ls='-.', linewidth=0.25)
plt.legend(fontsize = 12, labels = ['Average Axial Strain', 'SGC01-TRI Grid 1', 'SGC02-TRI Grid 1', 'SGC03-TRI Grid 1', 'SGC04-TRI Grid 1'], loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=3)
# Set the Y-axis to the right-hand side
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position('right')
# Adjust the plot layout
plt.subplots_adjust(left=0.05, right=0.90, bottom=0.2, top=0.95)
# Save the plot and show
plt.savefig(dir/'Stress vs Axial Strain.svg', format='svg')
pdf.savefig()
plt.show()

#Plot between Stress and Hoop Strains
fig, ax1 = plt.subplots(figsize=(10.17, 6.64), dpi=160)
df_failure_gauges = df_failure.loc[:,
                          (df_failure.columns.str.contains('SGC')) & (df_failure.columns.str.contains('Grid 3')) & (df_failure.columns.str.contains('TRI'))]
df_failure_gauges.insert(0, 'Gauge Strain Average', df_failure['StrainAverage3']*10000)
ax1.plot(df_failure_gauges/10000, df_failure['Stress'])
ax1.set_xlabel('Hoop Strain (%)', fontsize=14)  # X axis labeling and font size
ax1.set_ylabel('Stress (MPa)', fontsize=14)  # Y axis labeling and font size
ax1.tick_params(axis='y')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax1.set_title('Stress vs Hoop Strain', fontsize=14)  # Chart title
ax1.grid(color='lightgrey', ls='-.', linewidth=0.25)
plt.legend(fontsize = 12, labels = ['Average Hoop Strain', 'SGC01-TRI Grid 3', 'SGC02-TRI Grid 3', 'SGC03-TRI Grid 3', 'SGC04-TRI Grid 3'], loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=3)
# Adjust the plot layout
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
# Save the plot and show
plt.savefig(dir/'Stress vs Hoop Strain.svg', format='svg')
pdf.savefig()
plt.show()

#Plot between Stress and Strains overall
fig, ax = plt.subplots(figsize=(10.17, 6.64), dpi=160)
# Defining plot
df_failure_gauges_axial = df_failure.loc[:, (df_failure.columns.str.contains('SGC')) & (df_failure.columns.str.contains('Grid 1')) & (df_failure.columns.str.contains('TRI'))]
df_failure_gauges_hoop = df_failure.loc[:, (df_failure.columns.str.contains('SGC')) & (df_failure.columns.str.contains('Grid 3')) & (df_failure.columns.str.contains('TRI'))]
ax.plot(df_failure_gauges_axial/10000, df_failure['Stress'], label='Axial Strain vs Stress')
ax.plot(df_failure_gauges_hoop/10000, df_failure['Stress'], label = 'Hoop Strain vs Stress')
ax.plot(df_failure2['SPStrainAverage'], df_failure2['Stress'], "-.")
ax.annotate(f"({max_stress_point[0]:.2f}, {max_stress_point[1]:.2f})", 
             xy=max_stress_point, xycoords='data',
             xytext=(max_stress_point[0], max_stress_point[1] - 20), fontsize=12, ha='center', color = 'red')
ax.plot(max_stress_point[0], max_stress_point[1], 'x', color='g', markersize = 8)
ax.axvline(x=max_stress_point[0], color='m', linewidth=0.5, linestyle='--')
ax.axhline(y=max_stress_point[1], color='m', linewidth=0.5, linestyle='--')
plt.ylabel('Stress (MPa)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax.yaxis.set_label_coords(0.7, 0.75)
# Adjust the plot layout
plt.subplots_adjust(left=0.05, right=0.95, bottom=0.2, top=0.95)
# set the x-axis limits based on the data
xmin, xmax = [-10,5]
ax.set_xlim(xmin, xmax)
# set the y-axis to be at the 0 position of the x-axis
ax.spines['left'].set_position('zero')
# Add legend
# add a legend and axis labels
plt.legend(fontsize = 12, labels = ['SGC01 Axial', 'SGC02 Axial', 'SGC03 Axial', 'SGC04 Axial', 'SGC01 Hoop', 'SGC02 Hoop', 'SGC03 Hoop', 'SGC04 Hoop', 'String Pot'], loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=5)
plt.xlabel('Strain (%)', fontsize = 14)
# remove the y-axis from the left position
ax.spines['left'].set_visible(True)
# set the y-ticks to only display positive values
#ax.set_yticks([tick for tick in ax.get_yticks() if tick >= 0])
# add a y-axis line through the x-ticks
ax.axvline(x=0, color='black', linewidth=0.5)
# remove the top border of the plot
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.axvline(x=-10.0, color='black', linewidth=1.25)
# Save the plot and show
plt.grid(color='lightgrey', ls='-.', linewidth=0.25)
plt.title('Stress vs Strain', fontsize=14)
plt.savefig(dir/'Stress vs Strain.svg', format='svg')
pdf.savefig()
plt.show()

#Strain transformation equations:
df_failure['SGC01_Grid1_header']=df_failure.loc[:,(df_failure.columns.str.contains('SGC01')) & (df_failure.columns.str.contains('Grid 1')) & (df_failure.columns.str.contains('TRI'))]
df_failure['SGC01_Grid2_header']=df_failure.loc[:,df_failure.columns.str.contains('SGC01') & df_failure.columns.str.contains('Grid 2') & (df_failure.columns.str.contains('TRI'))]
df_failure['SGC01_Grid3_header']=df_failure.loc[:,df_failure.columns.str.contains('SGC01') & df_failure.columns.str.contains('Grid 3') & (df_failure.columns.str.contains('TRI'))]
df_failure['SGC02_Grid1_header']=df_failure.loc[:,df_failure.columns.str.contains('SGC02') & df_failure.columns.str.contains('Grid 1') & (df_failure.columns.str.contains('TRI'))]
df_failure['SGC02_Grid2_header']=df_failure.loc[:,df_failure.columns.str.contains('SGC02') & df_failure.columns.str.contains('Grid 2') & (df_failure.columns.str.contains('TRI'))]
df_failure['SGC02_Grid3_header']=df_failure.loc[:,df_failure.columns.str.contains('SGC02') & df_failure.columns.str.contains('Grid 3') & (df_failure.columns.str.contains('TRI'))]
df_failure['SGC03_Grid1_header']=df_failure.loc[:,df_failure.columns.str.contains('SGC03') & df_failure.columns.str.contains('Grid 1') & (df_failure.columns.str.contains('TRI'))]
df_failure['SGC03_Grid2_header']=df_failure.loc[:,df_failure.columns.str.contains('SGC03') & df_failure.columns.str.contains('Grid 2') & (df_failure.columns.str.contains('TRI'))]
df_failure['SGC03_Grid3_header']=df_failure.loc[:,df_failure.columns.str.contains('SGC03') & df_failure.columns.str.contains('Grid 3') & (df_failure.columns.str.contains('TRI'))]
df_failure['SGC04_Grid1_header']=df_failure.loc[:,df_failure.columns.str.contains('SGC04') & df_failure.columns.str.contains('Grid 1') & (df_failure.columns.str.contains('TRI'))]
df_failure['SGC04_Grid2_header']=df_failure.loc[:,df_failure.columns.str.contains('SGC04') & df_failure.columns.str.contains('Grid 2') & (df_failure.columns.str.contains('TRI'))]
df_failure['SGC04_Grid3_header']=df_failure.loc[:,df_failure.columns.str.contains('SGC04') & df_failure.columns.str.contains('Grid 3') & (df_failure.columns.str.contains('TRI'))]
fibre_angle_from_axial_direction = Input_Sheet['B13'].value
angle_radians1=math.radians(fibre_angle_from_axial_direction)
angle_radians2=math.radians(fibre_angle_from_axial_direction+90)
angle_radians3 = math.radians(fibre_angle_from_axial_direction+45)
df_failure['shearSGC01'] = df_failure['SGC01_Grid2_header'] - ((df_failure['SGC01_Grid1_header']+df_failure['SGC01_Grid3_header'])/2)
df_failure['shearSGC02'] = df_failure['SGC02_Grid2_header'] - ((df_failure['SGC02_Grid1_header']+df_failure['SGC02_Grid3_header'])/2)
df_failure['shearSGC03'] = df_failure['SGC03_Grid2_header'] - ((df_failure['SGC03_Grid1_header']+df_failure['SGC03_Grid3_header'])/2)
df_failure['shearSGC04'] = df_failure['SGC04_Grid2_header'] - ((df_failure['SGC04_Grid1_header']+df_failure['SGC04_Grid3_header'])/2)
#Calculating the Tape Strain (%)
df_failure['SGC01 Long Fibre Strain']= (((df_failure['SGC01_Grid1_header']+df_failure['SGC01_Grid3_header'])/2) + (((df_failure['SGC01_Grid1_header']-df_failure['SGC01_Grid3_header'])/2)*math.cos(2*angle_radians1))+(df_failure['shearSGC01']*math.sin(2*angle_radians1)))/10000
df_failure['SGC02 Long Fibre Strain']= (((df_failure['SGC02_Grid1_header']+df_failure['SGC02_Grid3_header'])/2) + (((df_failure['SGC02_Grid1_header']-df_failure['SGC02_Grid3_header'])/2)*math.cos(2*angle_radians1))+(df_failure['shearSGC02']*math.sin(2*angle_radians1)))/10000
df_failure['SGC03 Long Fibre Strain']= (((df_failure['SGC03_Grid1_header']+df_failure['SGC03_Grid3_header'])/2) + (((df_failure['SGC03_Grid1_header']-df_failure['SGC03_Grid3_header'])/2)*math.cos(2*angle_radians1))+(df_failure['shearSGC03']*math.sin(2*angle_radians1)))/10000
df_failure['SGC04 Long Fibre Strain']= (((df_failure['SGC04_Grid1_header']+df_failure['SGC04_Grid3_header'])/2) + (((df_failure['SGC04_Grid1_header']-df_failure['SGC04_Grid3_header'])/2)*math.cos(2*angle_radians1))+(df_failure['shearSGC04']*math.sin(2*angle_radians1)))/10000
#Calculating Transverse Tape Strain (%)
df_failure['SGC01 Trans Fibre Strain']= (((df_failure['SGC01_Grid1_header']+df_failure['SGC01_Grid3_header'])/2) + (((df_failure['SGC01_Grid1_header']-df_failure['SGC01_Grid3_header'])/2)*math.cos(2*angle_radians2))+(df_failure['shearSGC01']*math.sin(2*angle_radians2)))/10000
df_failure['SGC02 Trans Fibre Strain']= (((df_failure['SGC02_Grid1_header']+df_failure['SGC02_Grid3_header'])/2) + (((df_failure['SGC02_Grid1_header']-df_failure['SGC02_Grid3_header'])/2)*math.cos(2*angle_radians2))+(df_failure['shearSGC02']*math.sin(2*angle_radians2)))/10000
df_failure['SGC03 Trans Fibre Strain']= (((df_failure['SGC03_Grid1_header']+df_failure['SGC03_Grid3_header'])/2) + (((df_failure['SGC03_Grid1_header']-df_failure['SGC03_Grid3_header'])/2)*math.cos(2*angle_radians2))+(df_failure['shearSGC03']*math.sin(2*angle_radians2)))/10000
df_failure['SGC04 Trans Fibre Strain']= (((df_failure['SGC04_Grid1_header']+df_failure['SGC04_Grid3_header'])/2) + (((df_failure['SGC04_Grid1_header']-df_failure['SGC04_Grid3_header'])/2)*math.cos(2*angle_radians2))+(df_failure['shearSGC04']*math.sin(2*angle_radians2)))/10000
df_failure['eb SGC01'] = (((df_failure['SGC01_Grid1_header']+df_failure['SGC01_Grid3_header'])/2) + (((df_failure['SGC01_Grid1_header']-df_failure['SGC01_Grid3_header'])/2)*math.cos(2*angle_radians3))+(df_failure['shearSGC01']*math.sin(2*angle_radians3)))/10000
df_failure['eb SGC02'] = (((df_failure['SGC02_Grid1_header']+df_failure['SGC02_Grid3_header'])/2) + (((df_failure['SGC02_Grid1_header']-df_failure['SGC02_Grid3_header'])/2)*math.cos(2*angle_radians3))+(df_failure['shearSGC02']*math.sin(2*angle_radians3)))/10000
df_failure['eb SGC03'] = (((df_failure['SGC03_Grid1_header']+df_failure['SGC03_Grid3_header'])/2) + (((df_failure['SGC03_Grid1_header']-df_failure['SGC03_Grid3_header'])/2)*math.cos(2*angle_radians3))+(df_failure['shearSGC03']*math.sin(2*angle_radians3)))/10000
df_failure['eb SGC04'] = (((df_failure['SGC04_Grid1_header']+df_failure['SGC04_Grid3_header'])/2) + (((df_failure['SGC04_Grid1_header']-df_failure['SGC04_Grid3_header'])/2)*math.cos(2*angle_radians3))+(df_failure['shearSGC04']*math.sin(2*angle_radians3)))/10000
#Calculating Tape Shear Strain (%)
df_failure['SGC01 Shear Fibre Strain']= 2*((df_failure['eb SGC01']) - ((df_failure['SGC01 Long Fibre Strain']) + (df_failure['SGC01 Trans Fibre Strain']))/2)
df_failure['SGC02 Shear Fibre Strain']= 2*((df_failure['eb SGC02']) - ((df_failure['SGC02 Long Fibre Strain']) + (df_failure['SGC02 Trans Fibre Strain']))/2)
df_failure['SGC03 Shear Fibre Strain']= 2*((df_failure['eb SGC03']) - ((df_failure['SGC03 Long Fibre Strain']) + (df_failure['SGC03 Trans Fibre Strain']))/2)
df_failure['SGC04 Shear Fibre Strain']= 2*((df_failure['eb SGC04']) - ((df_failure['SGC04 Long Fibre Strain']) + (df_failure['SGC04 Trans Fibre Strain']))/2)
#Calculating the Maximum Principal Strain
df_failure['SGC01 Maximum Principal Strain']= (((df_failure['SGC01_Grid1_header']+df_failure['SGC01_Grid3_header'])/2)+((((((df_failure['SGC01_Grid1_header']-df_failure['SGC01_Grid2_header'])**2)+((df_failure['SGC01_Grid2_header']-df_failure['SGC01_Grid3_header'])**2))/2))**0.5))/10000
df_failure['SGC02 Maximum Principal Strain']= (((df_failure['SGC02_Grid1_header']+df_failure['SGC02_Grid3_header'])/2)+((((((df_failure['SGC02_Grid1_header']-df_failure['SGC02_Grid2_header'])**2)+((df_failure['SGC02_Grid2_header']-df_failure['SGC02_Grid3_header'])**2))/2))**0.5))/10000
df_failure['SGC03 Maximum Principal Strain']= (((df_failure['SGC03_Grid1_header']+df_failure['SGC03_Grid3_header'])/2)+((((((df_failure['SGC03_Grid1_header']-df_failure['SGC03_Grid2_header'])**2)+((df_failure['SGC03_Grid2_header']-df_failure['SGC03_Grid3_header'])**2))/2))**0.5))/10000
df_failure['SGC04 Maximum Principal Strain']= (((df_failure['SGC04_Grid1_header']+df_failure['SGC04_Grid3_header'])/2)+((((((df_failure['SGC04_Grid1_header']-df_failure['SGC04_Grid2_header'])**2)+((df_failure['SGC04_Grid2_header']-df_failure['SGC04_Grid3_header'])**2))/2))**0.5))/10000
#Calculating the Minimum Principal Strain
df_failure['SGC01 Minimum Principal Strain']= (((df_failure['SGC01_Grid1_header']+df_failure['SGC01_Grid3_header'])/2)-((((((df_failure['SGC01_Grid1_header']-df_failure['SGC01_Grid2_header'])**2)+((df_failure['SGC01_Grid2_header']-df_failure['SGC01_Grid3_header'])**2))/2))**0.5))/10000
df_failure['SGC02 Minimum Principal Strain']= (((df_failure['SGC02_Grid1_header']+df_failure['SGC02_Grid3_header'])/2)-((((((df_failure['SGC02_Grid1_header']-df_failure['SGC02_Grid2_header'])**2)+((df_failure['SGC02_Grid2_header']-df_failure['SGC02_Grid3_header'])**2))/2))**0.5))/10000
df_failure['SGC03 Minimum Principal Strain']= (((df_failure['SGC03_Grid1_header']+df_failure['SGC03_Grid3_header'])/2)-((((((df_failure['SGC03_Grid1_header']-df_failure['SGC03_Grid2_header'])**2)+((df_failure['SGC03_Grid2_header']-df_failure['SGC03_Grid3_header'])**2))/2))**0.5))/10000
df_failure['SGC04 Minimum Principal Strain']= (((df_failure['SGC04_Grid1_header']+df_failure['SGC04_Grid3_header'])/2)-((((((df_failure['SGC04_Grid1_header']-df_failure['SGC04_Grid2_header'])**2)+((df_failure['SGC04_Grid2_header']-df_failure['SGC04_Grid3_header'])**2))/2))**0.5))/10000
#Calculating Principal Angle
df_failure['SGC01 Principal Angle']= np.degrees(0.5*np.arctan((2*df_failure['SGC01_Grid2_header']-df_failure['SGC01_Grid1_header']-df_failure['SGC01_Grid3_header'])/(df_failure['SGC01_Grid1_header']-df_failure['SGC01_Grid3_header'])))
df_failure['SGC02 Principal Angle']= np.degrees(0.5*np.arctan((2*df_failure['SGC02_Grid2_header']-df_failure['SGC02_Grid1_header']-df_failure['SGC02_Grid3_header'])/(df_failure['SGC02_Grid1_header']-df_failure['SGC02_Grid3_header'])))
df_failure['SGC03 Principal Angle']= np.degrees(0.5*np.arctan((2*df_failure['SGC03_Grid2_header']-df_failure['SGC03_Grid1_header']-df_failure['SGC03_Grid3_header'])/(df_failure['SGC03_Grid1_header']-df_failure['SGC03_Grid3_header'])))
df_failure['SGC04 Principal Angle']= np.degrees(0.5*np.arctan((2*df_failure['SGC04_Grid2_header']-df_failure['SGC04_Grid1_header']-df_failure['SGC04_Grid3_header'])/(df_failure['SGC04_Grid1_header']-df_failure['SGC04_Grid3_header'])))
df_failure['Avg. SGC Principal Angle']= (df_failure['SGC01 Principal Angle'] + df_failure['SGC02 Principal Angle'] + df_failure['SGC03 Principal Angle'] + df_failure['SGC04 Principal Angle'])/4
#Maximum Shear Strain Calculation
df_failure['SGC01 Maximum Shear Strain']= 2*(((((df_failure['SGC01_Grid1_header']-df_failure['SGC01_Grid3_header'])/2)**2)+((df_failure['shearSGC01'])**2))**0.5)/10000
df_failure['SGC02 Maximum Shear Strain']= 2*(((((df_failure['SGC02_Grid1_header']-df_failure['SGC02_Grid3_header'])/2)**2)+((df_failure['shearSGC02'])**2))**0.5)/10000
df_failure['SGC03 Maximum Shear Strain']= 2*(((((df_failure['SGC03_Grid1_header']-df_failure['SGC03_Grid3_header'])/2)**2)+((df_failure['shearSGC03'])**2))**0.5)/10000
df_failure['SGC04 Maximum Shear Strain']= 2*(((((df_failure['SGC04_Grid1_header']-df_failure['SGC04_Grid3_header'])/2)**2)+((df_failure['shearSGC04'])**2))**0.5)/10000


#Plot between Tape fibre strain (%) and Stress (MPa)
fig, ax1 = plt.subplots(figsize=(10.17, 6.64), dpi=160)
dffailure = df_failure.loc[:,
                          (df_failure.columns.str.contains('SGC')) & (df_failure.columns.str.contains('Long Fibre Strain'))]
ax1.plot(dffailure, df_failure['Stress'])
ax1.set_xlabel('Longitudinal Surface Fibre Strain (%)', fontsize=14)  # X axis labeling and font size
ax1.set_ylabel('Stress (MPa)', fontsize=14)  # Y axis labeling and font size
ax1.tick_params(axis='y')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax1.set_title('Axial Stress vs Strain along fibre direction', fontsize=14)  # Chart title
ax1.grid(color='lightgrey', ls='-.', linewidth=0.25)
plt.legend(fontsize = 12, labels = ['SGC01', 'SGC02', 'SGC03', 'SGC04'], loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2)
# Adjust the plot layout
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
# Save the plot and show
plt.savefig(dir/'Stress vs Tape Strain.svg', format='svg')
pdf.savefig()
plt.show()

#Plot between Transverse to Tape fibre strain (%) and Stress (MPa)
fig, ax1 = plt.subplots(figsize=(10.17, 6.64), dpi=160)
dffailure = df_failure.loc[:,
                          (df_failure.columns.str.contains('SGC')) & (df_failure.columns.str.contains('Trans Fibre Strain'))]
ax1.plot(dffailure, df_failure['Stress'])
ax1.set_xlabel('Transverse Surface Fibre Strain (%)', fontsize=14)  # X axis labeling and font size
ax1.set_ylabel('Stress (MPa)', fontsize=14)  # Y axis labeling and font size
ax1.tick_params(axis='y')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax1.set_title('Axial Stress vs Strain transverse to fibre direction', fontsize=14)  # Chart title
ax1.grid(color='lightgrey', ls='-.', linewidth=0.25)
plt.legend(fontsize = 12, labels = ['SGC01', 'SGC02', 'SGC03', 'SGC04'], loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2)
# Adjust the plot layout
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
# Save the plot and show
plt.savefig(dir/'Stress vs Transverse Tape Strain.svg', format='svg')
pdf.savefig()
plt.show()

#Plot between Tape Shear strain (%) and Stress (MPa)
fig, ax1 = plt.subplots(figsize=(10.17, 6.64), dpi=160)
dffailure = df_failure.loc[:,
                          (df_failure.columns.str.contains('SGC')) & (df_failure.columns.str.contains('Shear Fibre Strain'))]
ax1.plot(dffailure, df_failure['Stress'])
ax1.set_xlabel('Surface Fibre Shear Strain (%)', fontsize=14)  # X axis labeling and font size
ax1.set_ylabel('Stress (MPa)', fontsize=14)  # Y axis labeling and font size
ax1.tick_params(axis='y')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax1.set_title('Axial Stress vs Shear Strain', fontsize=14)  # Chart title
ax1.grid(color='lightgrey', ls='-.', linewidth=0.25)
plt.legend(fontsize = 12, labels = ['SGC01', 'SGC02', 'SGC03', 'SGC04'], loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2)
# Adjust the plot layout
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
# Save the plot and show
plt.savefig(dir/'Stress vs Shear Strain.svg', format='svg')
pdf.savefig()
plt.show()

#Plot between Maximum Principal Strain (%) and Stress
fig, ax1 = plt.subplots(figsize=(10.17, 6.64), dpi=160)
dffailure = df_failure.loc[:,
                          (df_failure.columns.str.contains('SGC')) & (df_failure.columns.str.contains('Maximum Principal Strain'))]
Loadintonnes = df_failure.loc[:,(df_failure.columns.str.contains('LC01'))]
ax1.plot(dffailure, df_failure['Stress'])
ax1.set_xlabel('Maximum Principal Strain (%)', fontsize=14)  # X axis labeling and font size
ax1.set_ylabel('Stress (MPa)', fontsize=14)  # Y axis labeling and font size
ax1.tick_params(axis='y')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax1.set_title('Axial Stress vs Maximum Principal Strain', fontsize=14)  # Chart title
ax1.grid(color='lightgrey', ls='-.', linewidth=0.25)
plt.legend(fontsize = 12, labels = ['SGC01', 'SGC02', 'SGC03', 'SGC04'], loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2)
# Adjust the plot layout
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
# Save the plot and show
plt.savefig(dir/'Stress vs Maximum Principal Strain.svg', format='svg')
pdf.savefig()
plt.show()

#Plot between Minimum Principal Strain (%) and Stress (MPa)
fig, ax1 = plt.subplots(figsize=(10.17, 6.64), dpi=160)
dffailure = df_failure.loc[:,
                          (df_failure.columns.str.contains('SGC')) & (df_failure.columns.str.contains('Minimum Principal Strain'))]
ax1.plot(dffailure, df_failure['Stress'])
ax1.set_xlabel('Minimum Principal Strain (%)', fontsize=14)  # X axis labeling and font size
ax1.set_ylabel('Stress (MPa)', fontsize=14)  # Y axis labeling and font size
ax1.tick_params(axis='y')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax1.set_title('Axial Stress vs Minimum Principal Strain', fontsize=14)  # Chart title
ax1.grid(color='lightgrey', ls='-.', linewidth=0.25)
plt.legend(fontsize = 12, labels = ['SGC01', 'SGC02', 'SGC03', 'SGC04'], loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2)
# Adjust the plot layout
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
# Save the plot and show
plt.savefig(dir/'Stress vs Minimum Principal Strain.svg', format='svg')
pdf.savefig()
plt.show()

#Plot between  Principal Angle (degrees) and Stress (MPa)
fig, ax1 = plt.subplots(figsize=(10.17, 6.64), dpi=160)
dffailure = df_failure.loc[:,
                          (df_failure.columns.str.contains('SGC')) & (df_failure.columns.str.contains('Principal Angle'))]
ax1.plot(dffailure, df_failure['Stress'])
ax1.set_xlabel('Principal Angle (degrees)', fontsize=14)  # X axis labeling and font size
ax1.set_ylabel('Stress (MPa)', fontsize=14)  # Y axis labeling and font size
ax1.tick_params(axis='y')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax1.set_title('Axial Stress vs Principal Angle', fontsize=14)  # Chart title
ax1.grid(color='lightgrey', ls='-.', linewidth=0.25)
plt.legend(fontsize = 12, labels = ['SGC01', 'SGC02', 'SGC03', 'SGC04'], loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2)
# Adjust the plot layout
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
# Save the plot and show
plt.savefig(dir/'Stress vs Principal Angle.svg', format='svg')
pdf.savefig()
plt.show()

#Plot between Stress (MPa) and Maximum Shear strain
fig, ax1 = plt.subplots(figsize=(10.17, 6.64), dpi=160)
dffailure_shearstr = df_failure.loc[:,
                          (df_failure.columns.str.contains('SGC')) & (df_failure.columns.str.contains('Maximum Shear Strain'))]
ax1.plot(dffailure_shearstr, df_failure['Stress'])
ax1.set_xlabel('Maximum Shear strain (%)', fontsize=14)  # X axis labeling and font size
ax1.set_ylabel('Stress (MPa)', fontsize=14)  # Y axis labeling and font size
ax1.tick_params(axis='y')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
ax1.set_title('Stress vs Maximum Shear Strain', fontsize=14)  # Chart title
ax1.grid(color='lightgrey', ls='-.', linewidth=0.25)
plt.legend(fontsize = 12, labels = ['SGC01', 'SGC02', 'SGC03', 'SGC04'], loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2)
# Adjust the plot layout
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
# Save the plot and show
plt.savefig(dir/'Stress vs Maximum Shear Strain.svg', format='svg')
pdf.savefig()
plt.show()
#####################################################################################################################################################################

#Reading the Modulus data file
filename_modulus = Input_Sheet['B20'].value
file_modulus = pd.read_csv(filename_modulus, thousands=',', header= 0, na_values=["Overflow", ])   # Reads the csv file

#Averaging the Grid 1, Grid 2 and Grid 3 strains and adding in a new column on the file_modulus excel sheet.
Gauge_Average_M = Mean_SGC_Grid1(file_modulus)
file_modulus['Gauge_Average_M'] = Mean_SGC_Grid1(file_modulus)
Gauge_Average_M2 = Mean_SGC_Grid1(file_modulus)
file_modulus['Gauge_Average_M2'] = Mean_SGC_Grid2(file_modulus)
Gauge_Average_M3 = Mean_SGC_Grid3(file_modulus)
file_modulus['Gauge_Average_M3'] = Mean_SGC_Grid3(file_modulus)

#Calculating strain from String pots
SP01_raw = file_modulus.loc[:, file_modulus.columns.str.contains('SP01')]
file_modulus['SP01_Strain'] = Stringpot_disp_strain(SP01_raw, SP01_starting_length, SP01_gauge_length)
SP02_raw = file_modulus.loc[:, file_modulus.columns.str.contains('SP02')]
file_modulus['SP02_Strain'] = Stringpot_disp_strain(SP02_raw, SP02_starting_length, SP02_gauge_length)

#Plot between Axial strains for the 4 strain gauges and Scan #
file_modulus_gauges = (file_modulus.loc[:,
                          (file_modulus.columns.str.contains('SGC')) & (file_modulus.columns.str.contains('Grid 1')) & (file_modulus.columns.str.contains('TRI'))])/10000
file_modulus['StrainAverage'] = file_modulus_gauges.mean(axis=1)
plt.figure(1, figsize=(10.17, 6.64), dpi=160)
#file_modulus_gauges.insert(0, 'Gauge Strain Average', file_modulus['StrainAverage'])
file_modulus_gauges.insert(0, 'Time', file_modulus['Scan #']/(10*60))
file_modulus_gauges.plot(x='Time')
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Time (min)', fontsize=14)  # X axis labeling and font size
plt.ylabel('Axial Strain (%)', fontsize=14)  # Y axis labeling and font size
plt.title('Axial Strain vs Time', fontsize=14)  # Chart title
plt.suptitle('(Modulus run)', fontsize = 10, style ='italic', x = 0.85, y = 0.985)
plt.legend(fontsize = 12, labels = ['SGC01-TRI Grid 1', 'SGC02-TRI Grid 1', 'SGC03-TRI Grid 1', 'SGC04-TRI Grid 1'], loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2)
plt.grid(color = 'lightgrey', ls = '-.', linewidth = 0.25)
plt.savefig(dir/'Axial Strain vs Time.svg', format = "svg")
pdf.savefig()
plt.show()

#Plot between Hoop Strains for the 4 strain gauges and Scan #
file_modulus_gauges = (file_modulus.loc[:,
                          (file_modulus.columns.str.contains('SGC')) & (file_modulus.columns.str.contains('Grid 3')) & (file_modulus.columns.str.contains('TRI'))])/10000
#file_modulus['StrainAverage3'] = file_modulus_gauges.mean(axis=1)
plt.figure(1, figsize=(10.17, 6.64), dpi=160)
#file_modulus_gauges.insert(0, 'Gauge Strain Average', file_modulus['StrainAverage3'])
file_modulus_gauges.insert(0, 'Time', file_modulus['Scan #']/(10*60))
file_modulus_gauges.plot(x='Time')
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlabel('Time (min)', fontsize=14)  # X axis labeling and font size
plt.ylabel('Hoop Strain (%)', fontsize=14)  # Y axis labeling and font size
plt.title('Hoop Strain vs Time', fontsize=14)  # Chart title
plt.suptitle('(Modulus run)', fontsize = 10, style ='italic', x = 0.85, y = 0.985)
plt.legend(fontsize = 12, labels = ['SGC01-TRI Grid 3', 'SGC02-TRI Grid 3', 'SGC03-TRI Grid 3', 'SGC04-TRI Grid 3'], loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2)
plt.grid(color = 'lightgrey', ls = '-.', linewidth = 0.25)
plt.savefig(dir/'Hoop Strain vs Time.svg', format = "svg")
pdf.savefig()
plt.show()

#Reading the Modulus data file
filename_modulus = Input_Sheet['B20'].value
file_modulus = pd.read_csv(filename_modulus, thousands=',', header= 0, na_values=["Overflow", "Error", ])   # Reads the csv file

#Calculating Compression modulus from the modulus csv using the peaks and valleys of the Axial Strains
distance = 100
prominence = 0.1
plateau_size = None 

file_modulus_gauges = (file_modulus.loc[:,
                          (file_modulus.columns.str.contains('SGC')) & (file_modulus.columns.str.contains('Grid 1')) & (file_modulus.columns.str.contains('TRI'))])/10000
file_modulus['StrainAverage'] = file_modulus_gauges.mean(axis=1)
plt.plot (file_modulus['StrainAverage'])

peak, height = find_peaks(-1* file_modulus['StrainAverage'], distance=distance, prominence=prominence, plateau_size=plateau_size)
plt.plot(peak, file_modulus['StrainAverage'][peak], "x")
plt.show()

plt.plot (file_modulus['StrainAverage'])
valleys, heights = find_peaks(file_modulus['StrainAverage'], distance=distance, prominence=prominence, plateau_size=plateau_size)
plt.plot(valleys, file_modulus['StrainAverage'][valleys], "o")
plt.show()

def spliting_cycles(df, distance, prominence):
    peak, height = find_peaks(-1 * df['StrainAverage'], distance=distance, prominence=prominence)  # Find peak
    v = df['StrainAverage']  # inverse displacement to calculate minimum points
    valleys, heights = find_peaks(v, distance=distance, prominence=prominence)  # find inverse peaks
    return peak, valleys
peaks_list, valleys_list = spliting_cycles(file_modulus, distance, prominence)

CSA = csa(Coupon_ID, Coupon_OD)
#Convert Load from Tonnes to Stress (MPa) i.e., multiplying load cell values by 9.81 and / CSA
file_modulus_loadcell = file_modulus.loc[:,
                      (file_modulus.columns.str.contains('LC'))]
file_modulus['Stress'] = file_modulus_loadcell*9.81*1000/CSA

upper_limit = Input_Sheet['B27'].value/10000
lower_limit = Input_Sheet['B26'].value/10000

df1 = file_modulus[(file_modulus['StrainAverage'] < -lower_limit) & (file_modulus['StrainAverage'] > -upper_limit)]

plt.plot (df1['StrainAverage'])
cycle1 = df1.loc[0:peaks_list[0]]
cycle2 = df1.loc[valleys_list[0]:peaks_list[1]]
cycle3 = df1.loc[valleys_list[1]:peaks_list[2]]
cycle4 = df1.loc[valleys_list[2]:df1.index[-1]]
modulus_cycle1, intercept1 = np.polyfit(cycle1['StrainAverage'], cycle1['Stress'], 1)
modulus_cycle2, intercept2 = np.polyfit(cycle2['StrainAverage'], cycle2['Stress'], 1)
modulus_cycle3, intercept3 = np.polyfit(cycle3['StrainAverage'], cycle3['Stress'], 1)
modulus_cycle4, intercept4 = np.polyfit(cycle4['StrainAverage'], cycle4['Stress'], 1)
print ('The modulus of cycle 1 is', -modulus_cycle1*100, 'MPa')
print ('The modulus of cycle 2 is', -modulus_cycle2*100, 'MPa')
print ('The modulus of cycle 3 is', -modulus_cycle3*100, 'MPa')
print ('The modulus of cycle 4 is', -modulus_cycle4*100, 'MPa')

modulus_mean_cycles2to4 = (modulus_cycle2 + modulus_cycle3 + modulus_cycle4)/3
print ('\n\nThe average modulus of cycles 2 to 4 is', -modulus_mean_cycles2to4*100)

plt.figure(0, figsize=(10.17, 6.64), dpi=160)
plt.plot(cycle1['StrainAverage'], cycle1['Stress'],
             label='Cycle 1 Modulus = {:.2f} GPa'.format(-modulus_cycle1/10))  # Plot Stress/Strain for Cycle 1
plt.plot(cycle2['StrainAverage'], cycle2['Stress'],
             label='Cycle 2 Modulus = {:.2f} GPa'.format(-modulus_cycle2/10))  # Plot Stress/Strain for Cycle 2
plt.plot(cycle3['StrainAverage'], cycle3['Stress'],
             label='Cycle 3 Modulus = {:.2f} GPa'.format(-modulus_cycle3/10))  # Plot Stress/Strain for Cycle 3
plt.plot(cycle4['StrainAverage'], cycle4['Stress'],
             label='Cycle 4 Modulus = {:.2f} GPa'.format(-modulus_cycle4/10))  # Plot Stress/Strain for Cycle 4
plt.xlabel('Axial Strain (%)', fontsize=14)  # X axis labeling and font size
plt.ylabel('Axial Stress (MPa)', fontsize=14)  # Y axis labeling and font size
plt.title('Modulus Cycles', fontsize=14)  # Chart title
plt.legend(fontsize = 12, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2)
plt.grid(color='lightgrey', ls='-.', linewidth=0.25)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
# Add note on top right of the plot
plt.text(0.75, 0.95, 'Average modulus of cycles 2 to 4 = {:.2f} GPa'.format(-modulus_mean_cycles2to4/10), 
         horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes,
         fontsize=12, color='teal', fontstyle='italic')
# Adjust the plot layout
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
# Save the plot and show
plt.savefig(dir/'Stress vs Strain for Modulus cycles.svg', format='svg')
pdf.savefig()
plt.show()

#Write modulus on the Input Excel Sheet
cycles_modulus = [modulus_cycle1, modulus_cycle2, modulus_cycle3, modulus_cycle4, modulus_mean_cycles2to4]
Output_Sheet['B3'] = modulus_cycle1/10
Output_Sheet['B4'] = modulus_cycle2/10
Output_Sheet['B5'] = modulus_cycle3/10
Output_Sheet['B6'] = modulus_cycle4/10
Output_Sheet['B7'] = modulus_mean_cycles2to4/10

modulus_mean_cycles2to4 = (modulus_cycle2 + modulus_cycle3 + modulus_cycle4)/3
intercept_average = (intercept2+intercept3+intercept4)/3
# Define the slope and intercept of the line
slope = modulus_mean_cycles2to4
intercept = intercept_average
# Generate x values from 2 to -12
x_values = range(2, -13, -1)
# Calculate y values using the slope and intercept
y_values = [slope*x + intercept for x in x_values]
#Plot between Stress and Strains overall
fig, ax = plt.subplots(figsize=(10.17, 6.64), dpi=160)
# Defining plot
df_failure_gauges_axial = df_failure.loc[:, (df_failure.columns.str.contains('SGC')) & (df_failure.columns.str.contains('Grid 1')) & (df_failure.columns.str.contains('TRI'))]
df_failure_gauges_hoop = df_failure.loc[:, (df_failure.columns.str.contains('SGC')) & (df_failure.columns.str.contains('Grid 3')) & (df_failure.columns.str.contains('TRI'))]
AverageGrid1Strain = df_failure_gauges_axial.mean(axis=1)
AverageGrid3Strain = df_failure_gauges_hoop.mean(axis=1)
# Calculate y values using the slope and intercept
y_new = df_failure['Stress'] #MPa
x_new = ((y_new - intercept)/(slope)) #strain % from modulus polynomial
diff1 = (AverageGrid1Strain)/10000 - x_new # Calculate the difference between AverageGrid1Strain and x_new
diff = diff1/x_new
index = np.where(diff[700:] >= 0.2)[0][0]+700 #Find the index value where the difference first reaches 0.2%, but skips the first 700 rows

print(index)
yieldstress_value = df_failure.loc[index, 'Stress']
print('The Yield Stress value is', yieldstress_value) # Find the yield point stress
yieldstrain_value = AverageGrid1Strain.loc[index]
print('The Yield Strain value is', yieldstrain_value/10000) # Find the yield point strain
ax.plot(AverageGrid1Strain/10000, df_failure['Stress'], label='Average Axial Strain')
ax.plot(AverageGrid3Strain/10000, df_failure['Stress'], label = 'Average Hoop Strain')
ax.plot(x_values, y_values, ":", label = 'Calculated Modulus Polynomial Fit')
ax.plot(df_failure2['SPStrainAverage'], df_failure2['Stress'], "-.", label = 'String Pot Strain')
# Add marker at yield point
ax.scatter(yieldstrain_value/10000, yieldstress_value, marker='x', color='g', s=100)
# Add vertical and horizontal dashed lines at yield point
ax.axvline(x=yieldstrain_value/10000, color='m', linewidth=0.5, linestyle='--')
ax.axhline(y=yieldstress_value, color='m', linewidth=0.5, linestyle='--')
ax.annotate(f'({yieldstrain_value/10000:.2f}, {yieldstress_value:.2f})', xy=(yieldstrain_value/10000, yieldstress_value), xytext=(yieldstrain_value/10000-2, yieldstress_value-10), color='red')
plt.ylabel('Stress (MPa)', fontsize = 14)
ax.yaxis.set_label_coords(0.7, 0.75)
# Adjust the plot layout
plt.subplots_adjust(left=0.05, right=0.95, bottom=0.2, top=0.95)
# set the x-axis limits based on the data
xmin, xmax = [-10,5]
ax.set_xlim(xmin, xmax)
ax.set_ylim(0,190)
# set the y-axis to be at the 0 position of the x-axis
ax.spines['left'].set_position('zero')
# Add legend
# add a legend and axis labels
plt.legend(fontsize = 12,loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2)
plt.xlabel('Strain (%)', fontsize = 14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
# remove the y-axis from the left position
ax.spines['left'].set_visible(True)
# set the y-ticks to only display positive values
#ax.set_yticks([tick for tick in ax.get_yticks() if tick >= 0])
# add a y-axis line through the x-ticks
ax.axvline(x=0, color='black', linewidth=0.5)
# remove the top border of the plot
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.axvline(x=-10.0, color='black', linewidth=1.25)
# Save the plot and show
plt.grid(color='lightgrey', ls='-.', linewidth=0.25)
plt.title('Yield Point', fontsize=14)
plt.savefig(dir/'Yield Point.svg', format='svg')
pdf.savefig()
plt.show()

# #Poissons for failure run
df_failure['SPStrainAverage'] = df_failure2['SPStrainAverage']
df_failure['Poisson'] = -df_failure['StrainAverage3']/df_failure['StrainAverage']
AveragePoissonRatio1 = df_failure['Poisson'].mean(axis=0)
# #Plot for Poisson Ratio (all cycles average)
# plt.figure(1, figsize=(10.17, 6.64), dpi=160)
# plt.plot(file_df['Scan #']/(10*60), file_df['Poisson'], label='Poisson ratio')
# plt.xticks(fontsize=12)
# plt.yticks(fontsize=12)
# plt.ylim(0,3)
# #plt.yaxis.set_minor_locator(MultipleLocator(5))
# plt.xlabel('Time (min)', fontsize=14)  # X axis labeling and font size
# plt.ylabel('Poisson Ratio', fontsize=14)  # Y axis labeling and font size
# plt.title('Poisson Ratio', fontsize=14)  # Chart title
# plt.suptitle('(Failure run)', fontsize = 10, style ='italic', x = 0.85, y = 0.985)
# plt.legend(fontsize = 12, loc='lower center', bbox_to_anchor=(0.5, -0.2))
# plt.grid(color='lightgrey', ls='-.', linewidth=0.25)
# plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
# pdf.savefig()
# plt.savefig(dir/'Failure Run Poisson Ratio.svg', format='svg', bbox_inches='tight')
# plt.show()

#plt.plot (file_modulus['StrainAverage3'])
p1 = df_failure[(df_failure['SPStrainAverage'] <= -0.15) & (df_failure['SPStrainAverage'] >= -0.50)] #Selects the dataframe between 50 and 250 bars
p2 = df_failure[(df_failure['SPStrainAverage'] <= -0.5) & (df_failure['SPStrainAverage'] >= -3.5)]
a_1 = df_failure[(df_failure['SPStrainAverage'] <= -0.15) & (df_failure['SPStrainAverage'] >= -0.50)]['Poisson'].mean() #calculating the mean poisson ratio in this range
a_2 = df_failure[(df_failure['SPStrainAverage'] <= -0.5) & (df_failure['SPStrainAverage'] >= -3.5)]['Poisson'].mean()

#Poisson's ratio vs Pressure
plt.figure(1, figsize=(10.17, 6.64), dpi=160)
plt.plot(df_failure['SPStrainAverage'], df_failure['Poisson'], label='Poisson ratio')
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim ([0,0.7])
#plt.ylim(0,3)
#plt.yaxis.set_minor_locator(MultipleLocator(5))
plt.xlabel('Axial Strain (%)', fontsize=14)  # X axis labeling and font size
plt.ylabel('Poisson Ratio', fontsize=14)  # Y axis labeling and font size
plt.title('Strain vs Poisson Ratio', fontsize=14)  # Chart title
plt.suptitle('(Failure run)', fontsize = 10, style ='italic', x = 0.85, y = 0.985)
plt.legend(fontsize = 12, loc='lower center', bbox_to_anchor=(0.5, -0.2))
plt.grid(color='lightgrey', ls='-.', linewidth=0.25)
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.2, top=0.95)
plt.text(0.55, 0.95, 'Average Poisson Ratio between -0.15% and -0.50% axial strain = {:.2f}'.format(a_1), 
         horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes,
         fontsize=12, color='teal', fontstyle='italic')
plt.text(0.55, 0.90, 'Average Poisson Ratio between -0.50% and -3.5% axial strain = {:.2f}'.format(a_2), 
         horizontalalignment='center', verticalalignment='center', transform=plt.gca().transAxes,
         fontsize=12, color='teal', fontstyle='italic')
pdf.savefig()
plt.savefig(dir/'Pressure vs Poisson Ratio.svg', format='svg', bbox_inches='tight')
plt.show()


Output_Sheet['A9'] = 'Yield Stress (MPa)'
Output_Sheet['B9'] = yieldstress_value
Output_Sheet['A10'] = 'Yield Strain (%)'
Output_Sheet['B10'] = yieldstrain_value/10000
Output_Sheet['A11'] = 'Poissons Ratio' 
Output_Sheet['B11'] = a_1

#Save files in the directory indicated on Excel sheet
pdf.close()
new_csv_name = '/Analysed Modulus Data.csv'  # Name of the file
save_folder = Input_Sheet['B23'].value
modulus_path = save_folder + new_csv_name
file_modulus.to_csv(modulus_path)  # Save file to new csv
workbook.save(workbooklocation)