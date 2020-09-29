#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 11:27:36 2019

@author: angelika
"""

specs_for_analysis = r"C:\full\path\to\AnalysisInformation.xlsx"
#%%

import data_import as di
import plotting as pl
#%% 
'''Data Import'''

input_files,save,category_to_plot,variations_in_category,plotting_variations,time_unit,data_exclude_but,data_remove = di.analysis_specification(specs_for_analysis)

data,devices = di.data_import(input_files,time_unit)#,lum_od=False)
data = di.data_excludes_but(data,data_exclude_but,plotting_variations,category_to_plot,variations_in_category)
data = di.data_remove_all(data,data_remove)

'''translation strain-ID to genotype?'''
if True:
    data=di.translateStrainID(data,TranslateTo='combined')

#%% 
'''Plotting'''

if True:# MIC dose response
    for timepoint in ['10 h']:
        pl.plot_IC_interaction(data,plotting_variations,timepoint,save)
        pl.plot_ICs(data,category_to_plot,variations_in_category,plotting_variations,timepoint,save,
                     x_axis=category_to_plot[1],continuous_xaxis=False,IC=0.5,normalize=False,ylog=False,ttest=False)
    pl.plot_dose_response(data,save,category_to_plot,variations_in_category,devices,y='OD',timepoints=['10 h'],SaveInduction=True)

pl.plot_doublingtime(data,save,time_unit,xaxis=category_to_plot[1])
pl.plot_dose_response(data,save,category_to_plot,variations_in_category,devices,y='LUM',timepoints=['30 min'],SaveInduction=True)
pl.plot_dose_response(data,save,category_to_plot,variations_in_category,devices,y='LUM',normalizeX='10 h',timepoints=['30 min'])


for plot_object in ['OD']:#,'LUM']:
    #pl.plot_timedependent_replicates(data,plot_object,category_to_plot,variations_in_category,plotting_variations,time_unit,save,devices)
    pl.plot_timedependent_averaged(data,plot_object,category_to_plot,variations_in_category,plotting_variations,time_unit,save,devices)
 
