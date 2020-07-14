#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 11:27:36 2019

@author: angelika
"""

specs_for_analysis = r"D:\Benutzer\Uni\Documents\CellWallAntibiotics\Lab\Experiments_Australia\200403_CESR_FL\AnalysisInformation.xlsx"

specs_for_analysis = specs_for_analysis.replace('\\','/')

#%%

import data_import as di
import plotting as pl
#%% 
'''Data Import'''

input_files,save,category_to_plot,variations_in_category,plotting_variations,time_unit,data_exclude_but,data_remove = di.analysis_specification(specs_for_analysis)

data,devices = di.data_import(input_files,time_unit,maxtime=2.1)#,lum_od=False)
data = di.data_excludes_but(data,data_exclude_but,plotting_variations,category_to_plot,variations_in_category)
data = di.data_remove_all(data,data_remove)

'''translation strain-ID to genotype?'''
if True:
    data=di.translateStrainID(data,TranslateTo='combined')

#%% 
'''Plotting'''

if True:# MIC dose response
    for timepoint in ['10 h']:
        pl.plot_ICs(data,category_to_plot,variations_in_category,plotting_variations,timepoint,save,
                     xaxis=category_to_plot[1],continuous_xaxis=False,IC=0.5,normalize=False,ylog=False)
    pl.plot_dose_response(data,save,category_to_plot,variations_in_category,devices,y='OD',timepoints=['10 h'])

pl.plot_doublingtime(data,save,time_unit,xaxis=category_to_plot[1])
pl.plot_dose_response(data,save,category_to_plot,variations_in_category,devices,y='LUM',timepoints=['30 min'],SaveInduction=True)
pl.plot_dose_response(data,save,category_to_plot,variations_in_category,devices,y='LUM',normalizeX='10 h',timepoints=['30 min'])

for category1 in set(data[category_to_plot[0]]):
    if len(category_to_plot)>1:
        for category2 in set(data[category_to_plot[1]]):
            
            selection = pl.select_data(data,category_to_plot,variations_in_category,category1,category2,plotting_variations)
            if selection.empty:
                continue
            for y_axis in ['LUM']:
                #pl.plot_timedependent_replicates(selection,category_to_plot,variations_in_category,category1,category2,time_unit,save,plot_object=y_axis)
                pl.plot_timedependent_averaged(selection,category_to_plot,variations_in_category,category1,category2,time_unit,save,devices,plot_object=y_axis)
            
    else:
        category2 =False
        selection =pl.select_data(data,category_to_plot,variations_in_category,category1,category2,plotting_variations)
        
        for y_axis in ['LUM']:
                pl.plot_timedependent_replicates(selection,category_to_plot,variations_in_category,category1,category2,time_unit,save,plot_object=y_axis)
                pl.plot_timedependent_averaged(selection,category_to_plot,variations_in_category,category1,category2,time_unit,save,plot_object=y_axis)#
