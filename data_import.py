#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 11:26:40 2019

@author: angelika
"""
import xlrd
import pandas as pd
import datetime
import numpy as np
import sys
import numbers
from scipy.signal import medfilt
import plotting as pl

def analysis_specification(specs_for_analysis):
    """Imports the specifications from an Excel sheet 
    specs_for_analysis: path to this excel sheet
    """
    if specs_for_analysis[-5:]!='.xlsx':
        return(print('specification of analysis in incorrect file format'))    
    wb = xlrd.open_workbook(filename = specs_for_analysis)
    worksheet = wb.sheet_by_name('Blatt1')
    active=False
    for n_row in range(worksheet.nrows):
        if len(worksheet.cell(n_row, 0).value)<1:   #reset in empty rows
            active=False
        elif not active:
            if worksheet.cell(n_row,0).value=='Input Files':
                active='Input Files'
                input_files = []
            elif worksheet.cell(n_row,0).value=='Saving File':
                active=False
                save = worksheet.cell(n_row+1,0).value
                if save[-1]!='/':
                    save=save+'/'
            elif worksheet.cell(n_row,0).value=='Categories to Plot':
                active=False
                category_to_plot = []
                for n_col in range(2):
                    if worksheet.cell(n_row+1,n_col).value!='':
                        category_to_plot.append(worksheet.cell(n_row+1,n_col).value)
            elif worksheet.cell(n_row,0).value=='Variations in Category':
                active=False
                variations_in_category = []
                for n_col in range(2):
                    if worksheet.cell(n_row+1,n_col).value!='':
                        variations_in_category.append(worksheet.cell(n_row+1,n_col).value)
            elif worksheet.cell(n_row,0).value=='Plotting Variations':
                active='Plotting Variations'
                plotting_variations = dict()
            elif worksheet.cell(n_row,0).value=='time unit to plot':
                active=False
                time_unit = worksheet.cell(n_row+1,0).value
            elif worksheet.cell(n_row,0).value=='Remove all data but with the given specification ':
                active='data_exclude_but'
                data_exclude_but = dict()
            elif worksheet.cell(n_row,0).value=='Remove data with given specification':
                active='data_remove'
                data_remove = dict()
        elif active:
            if active=='Input Files':
                input_files.append(worksheet.cell(n_row, 0).value)
            elif active=='Plotting Variations':
                conc=[]
                for n_col in range(1,worksheet.ncols):
                    if isinstance(worksheet.cell(n_row, n_col).value,numbers.Number):
                        conc.append(worksheet.cell(n_row, n_col).value)
                plotting_variations[change_unit_delimiter(worksheet.cell(n_row, 0).value)]=conc
            elif active=='data_exclude_but':
                options = []
                for n_col in range(1,worksheet.ncols):
                    if worksheet.cell(n_row, n_col).value!='':
                        options.append(worksheet.cell(n_row, n_col).value)             
                if worksheet.cell(n_row, 0).value in data_exclude_but:
                    options.extend(data_exclude_but[worksheet.cell(n_row, 0).value])
                    options = set(options)
                    data_exclude_but[worksheet.cell(n_row, 0).value]=list(options)
                else:
                    data_exclude_but[worksheet.cell(n_row, 0).value]=options
            elif active=='data_remove':
                options = []
                for n_col in range(1,worksheet.ncols):
                    if worksheet.cell(n_row, n_col).value!='':
                        options.append(worksheet.cell(n_row, n_col).value)             
                if worksheet.cell(n_row, 0).value in data_remove:
                    options.extend(data_remove[worksheet.cell(n_row, 0).value])
                    options = set(options)
                    data_remove[worksheet.cell(n_row, 0).value]=list(options)
                else:
                    data_remove[worksheet.cell(n_row, 0).value]=options
    wb.release_resources()
    if len(set(input_files))!=len(input_files):
        raise SystemExit('One experiment is imported twice. Please remove one copy in the AnalysisInformation.xlsx file.')    
    
    return input_files,save,category_to_plot,variations_in_category,plotting_variations,time_unit,data_exclude_but,data_remove


def data_import(input_files,time_unit,lum_od=True,smoothing=True,maxtime=False):
    """imports all plate reader results, given in list with the respective metainfo. 
    Depending of the reader luminescence data is automatically imported as well if it is contained in the same file or is in the same folder as OD data and has the same name with additional '_lum'.
    This is hard coded though. Biotek needs lum data in the same file. Every other reader in another file. 
    
    Keyword arguments:
        input_files: list of absolute paths, pointing towards the OD-readings of experiments
        time_unit: 'h','m','s'
        lum_od: boolean, specifies if luminescence data should be normalized by the OD
        smoothing: boolean, smoothes curves with a median filter, window size=3
        maxtime: num, maximal time in h that is being imported (imports the whole course of the experiment if False)
    """
    start=True
    devices = []
    for file in input_files:
        lum=True
        metainfo,time_between_measurements = metainfo_import(file)
        reader = findreader(file)
        if reader=='biotek':
            dataOD,dataLUM = read_biotek(file,lum_od,metainfo,time_between_measurements,reader,smoothing,time_unit)
            lum=True
        elif reader=='SPECTROstar':
            dataOD,dataLUM,metainfo = read_spectrostar(file,metainfo,time_between_measurements,reader,smoothing,time_unit)
            lum=False
        elif reader=='SPECTROstar TaylorLab':
            dataOD,dataLUM = read_omega(file,lum_od,metainfo,time_between_measurements,reader,time_unit,smoothing)
            lum=False
        elif reader=='CLARIOstar':
            dataOD,dataLUM = read_clariostar(file,lum_od,metainfo,time_between_measurements,reader,time_unit,smoothing)
            lum=True    
        elif reader=='FLUOstar Omega':
            dataOD,dataLUM = read_omega(file,lum_od,metainfo,time_between_measurements,reader,time_unit,smoothing)
            lum=False
        elif reader=='victor':
            dataOD,dataLUM,lum = read_victor(file,lum_od,metainfo,time_between_measurements,reader,smoothing,time_unit)
        else:
            print('ERROR: can\'t idenify plate reader')
            sys.exit()
        tmp_info = determine_doublingtime(dataOD,metainfo)
        tmp_info = join_data_metainfo(metainfo,dataOD,dataLUM,lum_od,lum,maxtime)
        if start:
            Overall_data=tmp_info.copy()
            start=False
        else:
            Overall_data = pd.concat([Overall_data,tmp_info],ignore_index=True,sort=False)
        if not reader in devices:
            devices.append(reader)
    devices=set(devices)
    if len(devices)>1:
        print('ATTENTION: data of more than one reader used for analysis')
    return(Overall_data, devices)
    
def smoothingData(data):
    """Smoothes all OD\LUM data in dataframe with a median filter, window size=3 """
    wells = get_wells()
    for well in wells:
        data[well]=medfilt(data[well], kernel_size=3)
    return(data)    

class InputError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
def metainfo_import(file):
    """Imports metainfo of experiment from excel sheet. Must be located in the same folder as OD data and with the same name +'_metainfo' """
    if file[-5:]=='.xlsx':
        mfile='{}_metainfo{}'.format(file[:-5],file[-5:])
    elif file[-4:]=='.xls':
        mfile='{}_metainfo{}'.format(file[:-4],'.xlsx')
    else:
        return(print('plate reader data in incorrect file format'))
    letters = 'ABCDEFGH'
    metainfo=[['Well']]
    for well in get_wells():
        metainfo.append([well])
    wb = xlrd.open_workbook(filename = mfile)
    try:
        worksheet = wb.sheet_by_name('Blatt1')
    except xlrd.XLRDError:
        raise InputError('Required worksheet of AnalysisInformation not found. Make sure to specify the path of the AnalysisInformation.xlsx sheet as specs_for_analysis.')
    active=False
    conc = False
    time_between_measurements = 0 # in hours
    for i in range(worksheet.nrows):    #Read all rows in metainfo sheet
        if len(worksheet.cell(i, 0).value)<1:   #reset in empty rows
            active=False
            conc=False
        elif len(worksheet.cell(i, 0).value)>1 and not active:    #initiate new column in metainfo when a new metainfotable begins
            if worksheet.cell(i, 0).value=='time between measurements':
                time_between_measurements = worksheet.cell(i, 1).value # in hours
                continue
            if worksheet.cell(i, 0).value=='Date':
                date = worksheet.cell(i, 1).value # in hours
                continue
            active=worksheet.cell(i, 0).value
            metainfo[0].append(active)
            
            if worksheet.cell(i, 14).value=='Concentration':    #allow notation of inducer concentration if nessecary
                metainfo[0].append('{}_Concentration'.format(active))
                conc=True

            else:
                conc=False
        elif worksheet.cell(i, 0).value in letters and active:  #notation of metainfo values
            for c in range(1,13):
                    metainfo_index=letters.index(worksheet.cell(i, 0).value)+1+(c-1)*8  #get index of well in metadata
                    if conc:
                        if isinstance(worksheet.cell(i, c+14).value,float):
                            metainfo[metainfo_index].extend([change_unit_delimiter(worksheet.cell(i, c).value),worksheet.cell(i, c+14).value])
                        else:
                            metainfo[metainfo_index].extend([change_unit_delimiter(worksheet.cell(i, c).value),float('NaN')])
                    else:
                        metainfo[metainfo_index].append(worksheet.cell(i, c).value)                
    metainfo=pd.DataFrame(metainfo[1:],columns=metainfo[0])
    metainfo.insert(0, 'Date', date, True)
    return(metainfo,time_between_measurements)
    
def change_unit_delimiter(inducer):
    if '[' in inducer:
        return inducer
    elif ' /' in inducer:
        inducer=inducer.split('/')
        inducer='{}[{}]'.format(inducer[0],'/'.join(inducer[1:]))
    return inducer

def read_biotek(file,lum_od,metainfo,time_between_measurements,reader,smoothing,time_unit):
    """Reads OD and luminescence data derived from biotek reader. Luminescence data must be contained in the same file"""
    wells=get_wells()
    table=pd.read_excel(file)
    LUMend=-1
    for i_line in range(len(table['Unnamed: 0'])):
        if table['Unnamed: 0'][i_line]=='Read 1:600':
            ODbeg=i_line+2
        elif table['Unnamed: 0'][i_line]=='Read 2:Lum':
            ODend=i_line-1
            LUMbeg=i_line+2
        elif table['Unnamed: 0'][i_line]=='Results':
            LUMend=i_line-1
    dataOD = table.iloc[ODbeg+1:ODend,1:].copy(deep=True)
    dataOD.columns=table.iloc[ODbeg,1:]
    dataOD = dataOD.set_index(pd.Index(range(len(dataOD))))
    try:
        dataOD[wells] = dataOD[wells].apply(pd.to_numeric)
        dataLUM = table.iloc[LUMbeg+1:LUMend,1:].copy(deep=True)
        dataLUM.columns=table.iloc[LUMbeg,1:]
        dataLUM = dataLUM.set_index(pd.Index(range(len(dataLUM))))
        dataLUM[wells] = dataLUM[wells].apply(pd.to_numeric)
    except ValueError:#in case the reader stopped the measurements, cut data to the last full measurement of the plate
        for i_time in range(len(dataOD)):
            if ("OVRFLW" in dataOD.iloc[i_time].unique()) or (dataOD.iloc[i_time].isnull().values.any()):
                dataOD= dataOD[:i_time] # apply cutoff for OD data
                dataOD[wells] = dataOD[wells].apply(pd.to_numeric)
                dataLUM = table.iloc[LUMbeg+1:LUMend,1:].copy(deep=True) #apply cutoff for LUM data
                dataLUM.columns=table.iloc[LUMbeg,1:]
                dataLUM = dataLUM.set_index(pd.Index(range(len(dataLUM))))
                dataLUM = dataLUM[:i_time]
                dataLUM[wells] = dataLUM[wells].apply(pd.to_numeric)
                break
        
    #Background Correction
    dataOD =background_correct(reader,dataOD,metainfo,'OD')
    dataLUM =background_correct(reader,dataLUM,metainfo,'LUM')
    
    #get ideal times and unit convertion
    dataOD=time_conversion(dataOD,time_between_measurements,time_unit)
    dataLUM=time_conversion(dataLUM,time_between_measurements,time_unit)
    if smoothing:
        dataOD=smoothingData(dataOD)
        dataLUM=smoothingData(dataLUM)
    
    # LUMcorrect
    if lum_od:
        dataLUM = lum_correct(dataOD,dataLUM)
    
    return(dataOD,dataLUM)

def find_time_unit(n_row,table):
    time_string=table[6][n_row-1].split('-')[-1].split()
    i=0
    hours=0
    while i<len(time_string):
        hours+=float(time_string[i])*find_timeconversion_factor(time_string[i+1][0],'h')
        i+=2
    if table[6][n_row]==hours:
        return 'h'
    elif table[6][n_row]==hours*find_timeconversion_factor('h', 'm'):
        return 'm'
    elif table[6][n_row]==hours*find_timeconversion_factor('h', 's'):
        return 's'
    else:
        print('Could not find unit of time')

def read_spectrostar(file,metainfo,time_between_measurements,reader,smoothing,time_unit):
    """Reads OD data of the spectrostar reader. This reader cannot read luminescence."""
    wells=get_wells()
    if 'plate run' in metainfo:
        number_of_runs=int(max([float(i) for i in list(set(metainfo['plate run'])) if i!='']))
        file=file.split('.')
        file[0]=file[0][:-1] 
        for i_run in range(number_of_runs):
            table=pd.read_excel('{0}{1}.{2}'.format(file[0],i_run+1,file[1]), header=None)
            begin_of_run=table[10][1].split()
            begin_of_run[1]=begin_of_run[1].split(':')
            if begin_of_run[2]=='PM' and begin_of_run[1][0]!='12':
                pm=12
            else:
                pm=0
            begin_of_run=(pm+float(begin_of_run[1][0]))*find_timeconversion_factor('h','h')+float(begin_of_run[1][1])*find_timeconversion_factor('m','h')+float(begin_of_run[1][2])*find_timeconversion_factor('s','h')
            for n_row in range(len(table)):#find beginning of raw data
                if table[0][n_row]=='Well\nRow' or table[0][n_row]=='Well':
                    n_row+=1
                    break
            data_time_unit=find_time_unit(n_row,table)
            table=table.iloc[n_row:len(table):].copy(deep=True)
            table = table.transpose()
            rename_cols = get_colnames(table,n_row)
            table=table.rename(columns=rename_cols)
            table = table.reset_index(drop=True)
            table= table.drop([0, 1, 2])
            if list(table['Time']).count(0)>1:
                cutoff=[i for i, n in enumerate(list(table['Time'])) if n == 0][1]
                table=table[0:cutoff]
            table = table.reset_index(drop=True)
            table[wells] = table[wells].apply(pd.to_numeric)
            table['Time'] = table['Time'].apply(pd.to_numeric)
            table['Time'] = table['Time']*find_timeconversion_factor(data_time_unit,'h')
            table=table.dropna()
            if i_run==0:
                dataOD=table.copy()
                begin_of_experiment=begin_of_run
                metainfo.loc[metainfo['plate run']==i_run+1,'plate run - start time']=0
            else:
                time_dif=float(begin_of_run-begin_of_experiment)
                table['Time'] = table['Time']+time_dif
                dataOD=pd.concat([dataOD, table],ignore_index=True)
                metainfo.loc[metainfo['plate run']==i_run+1,'plate run - start time']=time_dif
        metainfo['plate run'] = metainfo['plate run'].astype(str)
    else:
        table=pd.read_excel(file, header=None)
        for n_row in range(len(table)):#find beginning of raw data
            if table[0][n_row]=='Well\nRow' or table[0][n_row]=='Well':
                n_row+=1
                break
        data_time_unit=find_time_unit(n_row,table)
        dataOD=table.iloc[n_row:len(table):].copy(deep=True)
        dataOD = dataOD.transpose()
        rename_cols = get_colnames(dataOD,n_row)
        dataOD=dataOD.rename(columns=rename_cols)
        dataOD = dataOD.reset_index(drop=True)
        dataOD= dataOD.drop([0, 1, 2])
        if list(dataOD['Time']).count(0)>1:
            cutoff=[i for i, n in enumerate(list(dataOD['Time'])) if n == 0][1]
            dataOD=dataOD[0:cutoff]
        dataOD = dataOD.reset_index(drop=True)
        dataOD[wells] = dataOD[wells].apply(pd.to_numeric)
        dataOD['Time'] = dataOD['Time'].apply(pd.to_numeric)
        dataOD['Time'] = dataOD['Time']*find_timeconversion_factor(data_time_unit,'h')
#    dataOD= dataOD[0:dataOD.loc[dataOD['Time']==0].index[1]] # temperature or other values might be listed after the ODs
    #Background Correction
    dataOD =background_correct(reader,dataOD,metainfo,'OD')
    #get ideal times and unit convertion
    dataOD=time_conversion(dataOD,time_between_measurements,time_unit)
    if smoothing:
        dataOD=smoothingData(dataOD)
    
    dataLUM = pd.DataFrame()
    return(dataOD,dataLUM,metainfo)
    
def read_clariostar(file,lum_od,metainfo,time_between_measurements,reader,time_unit,smoothing):
    """Reads data of CLARIOstar reader. Must have luminescence data in seperate file."""
    wells=get_wells()
    dataOD=pd.DataFrame()
    dataLUM=pd.DataFrame()
    fileLUM=file.split('.')
    fileLUM='{}_lum.{}'.format(fileLUM[0],fileLUM[1])
    for data,path,dtype in zip([dataOD,dataLUM],[file,fileLUM],['OD','LUM']):
        table=pd.read_excel(path, header=None)
        for n_row in range(len(table)):#find beginning of raw data
            if table[0][n_row]=='Well':
                n_row+=1
                break
        data=table.iloc[n_row:len(table):].copy(deep=True)
        data = data.transpose()
        rename_cols = get_colnames(data,n_row)
        data=data.rename(columns=rename_cols)
        data = data.reset_index(drop=True)
        data = data.drop([0, 1])
        if list(data['Time']).count(0)>1:
            cutoff=[i for i, n in enumerate(list(data['Time'])) if n == 0][1]
            data=data[0:cutoff]
        data = data.reset_index(drop=True)
        data[wells] = data[wells].apply(pd.to_numeric)
        try:
            data['Time'] = data['Time'].apply(pd.to_numeric)
            data['Time'] = data['Time']*find_timeconversion_factor('s','h')
        except ValueError:
            data['Time'] = data['Time'].str.split()
            unit=['h','m','s']
            for i in range(len(data['Time'])):
                ind_unit=0
                tmp=0
                for j in range(0,len(data['Time'].iloc[i]),2):
                    tmp=(tmp+int(data.loc[i,'Time'][j]))*60
                    ind_unit+=1
                data.loc[i,'Time']=tmp*find_timeconversion_factor(unit[ind_unit],'h')
        #Background Correction
        data =background_correct(reader,data,metainfo,dtype)
        if dtype=='OD':
            dataOD=data.copy(deep=True)
        elif dtype=='LUM':
            dataLUM=data.copy(deep=True)
    
    #get ideal times and time convertion to hours
    dataOD=time_conversion(dataOD,time_between_measurements,time_unit)
    dataLUM=time_conversion(dataLUM,time_between_measurements,time_unit)
    
    if smoothing:
        dataOD=smoothingData(dataOD)
        dataLUM=smoothingData(dataLUM)    
    
    if lum_od:
        dataLUM = lum_correct(dataOD,dataLUM)
    return(dataOD,dataLUM)
    
def read_omega(file,lum_od,metainfo,time_between_measurements,reader,time_unit,smoothing):
    """Reads data of FLUOstar Omega reader. """
    wells=get_wells()
    dataOD=pd.DataFrame()
    table=pd.read_excel(file, header=None)
    for n_row in range(len(table)):#find beginning of raw data
        if table[0][n_row]=='Well':
            n_row+=1
            break
    dataOD=table.iloc[n_row:len(table):].copy(deep=True)
    dataOD = dataOD.transpose()
    rename_cols = get_colnames(dataOD,n_row)
    dataOD=dataOD.rename(columns=rename_cols)
    dataOD = dataOD.reset_index(drop=True)
    dataOD = dataOD.drop([0, 1])
    if list(dataOD['Time']).count(0)>1:
        cutoff=[i for i, n in enumerate(list(dataOD['Time'])) if n == 0][1]
        dataOD=dataOD[0:cutoff]
    dataOD = dataOD.reset_index(drop=True)
    dataOD[wells] = dataOD[wells].apply(pd.to_numeric)
    try:
        dataOD['Time'] = dataOD['Time'].apply(pd.to_numeric)
        dataOD['Time'] = dataOD['Time']*find_timeconversion_factor('s','h')
    except ValueError:
        dataOD['Time'] = dataOD['Time'].str.split()
        unit=['h','m','s']
        for i in range(len(dataOD['Time'])):
            ind_unit=0
            tmp=0
            for j in range(0,len(dataOD['Time'].iloc[i]),2):
                tmp=(tmp+int(dataOD.loc[i,'Time'][j]))*60
                ind_unit+=1
            dataOD.loc[i,'Time']=tmp*find_timeconversion_factor(unit[ind_unit],'h')
    #Background Correction
    dataOD =background_correct(reader,dataOD,metainfo,'OD')
    
    #get ideal times and time convertion to hours
    dataOD=time_conversion(dataOD,time_between_measurements,time_unit)
    
    if smoothing:
        dataOD=smoothingData(dataOD)
    dataLUM = pd.DataFrame()
    return(dataOD,dataLUM)

def read_victor(file,lum_od,metainfo,time_between_measurements,reader,smoothing,time_unit):
    """Reads experimental data of victor reader. Automatically recognizes if no luminescene has been measured. 
    Luminescence data must be contained in the same file.
    """
    wells=get_wells()
    table=pd.read_excel(file)
    
    if table.shape[1]<7: #flag in case no luminescence data was measured
        lum=False
    else:
        lum=True
        
    dataOD=pd.DataFrame(columns=wells)
    times=[t*time_between_measurements*find_timeconversion_factor('h',time_unit) for t in range(0,list(table['Well']).count('A01'))]
    dataOD.insert(0,'Time',times)
    if lum:
        dataLUM=pd.DataFrame(columns=wells)
        dataLUM.insert(0,'Time',times)
    else:
        dataLUM=pd.DataFrame()
    t=0
    for i in range(len(table)):
        w=table.loc[i,'Well']
        if w[1]=='0':
            w=w.replace('0','')
        dataOD.loc[t,w]=table.loc[i,'OD600 (A)']
        if lum:
            dataLUM.loc[t,w]=table.loc[i,'CPS (CPS)']
        if w=='H12':
            t+=1
            
    dataOD[wells] = dataOD[wells].apply(pd.to_numeric)
    #Background Correction
    dataOD =background_correct(reader,dataOD,metainfo,'OD')
    
    if smoothing:
        dataOD=smoothingData(dataOD)
    if lum:
        dataLUM[wells] = dataLUM[wells].apply(pd.to_numeric)
        dataLUM =background_correct(reader,dataLUM,metainfo,'LUM')
        if smoothing:
            dataLUM=smoothingData(dataLUM)
        if lum_od:
            dataLUM = lum_correct(dataOD,dataLUM)
    return dataOD,dataLUM,lum

def get_colnames(dataOD,n_row):
    """Returns dictonary to rename columns for BMG readers."""
    rename_cols = {n_row:'Time'}
    i=0
    for col in ['A','B','C','D','E','F','G','H']:
        for j in range(1,13):
            if i >dataOD.shape[1]:
                break
            rename_cols.update({n_row+1+i:'{}{}'.format(col,j)})
            i+=1            
    return rename_cols

def time_conversion(data,time_between_measurements,time_unit):
    """Sets every time point to the ideal measurement according to the given time between measurments. """
    t=0
    ideal_time=[t]
    if isinstance(data['Time'][0], datetime.time) or isinstance(data['Time'][0], str):
        try:
            data.drop('Time',axis=0,inplace=True)
        except KeyError:
            pass
        data['Time']=data['Time'].apply(str)
        for i in range(len(data)):
            tmp=data['Time'][i].split(':')
            if len(tmp[0])>2:
                tmp[0]=24*float(tmp[0][-4])+float(tmp[0][-2:])
            time=float(tmp[0])+float(tmp[1])*find_timeconversion_factor('m','h')+float(tmp[2])*find_timeconversion_factor('s','h')
            data.loc[i,'Time']=time
    while t <max(data['Time']):
        t+= time_between_measurements
        ideal_time.append(t)
    for t in range(len(data['Time'])):
        data.loc[t,'Time']=ideal_time[pl.find_idx_to_closest_value(ideal_time,data.loc[t,'Time'])]*find_timeconversion_factor('h',time_unit)
    return data

def lum_correct(dataOD,dataLUM):
    """Normalizes luminescence data with the OD at the same timepoint.
    Sets all lum values to a low value (10) which are:
        (1) below this value - prevents plotting errors
        (2) whose corresponding OD values are so low that they would artifically show extremely high luminescence/od
    """
    wells = get_wells()
    dataLUM[wells] =dataLUM[wells]/dataOD[wells]
    
    #thresholding
    for well in wells:#avoid artifically high LUM/OD values at very low ODs or negatives
        dataLUM[well].loc[dataLUM[well]<=10**1] = 10**1
        dataLUM[well].loc[dataOD[well]<=3*10**(-3)] = 10**1
    return dataLUM 

def get_wells():
    """Returns list of all wells in 96-well plate."""
    letters = 'ABCDEFGH'
    wells=[]
    for i in range(1,13):
        for l in letters:
            wells.append(l+str(i))
    return(wells)
    
def join_data_metainfo(metainfo,dataOD,dataLUM,lum_od,lum,maxtime):
    """Creates one dataframe that includes all metainfo and od & lum data"""
    ODs=[]
    LUMs=[]
    for well in get_wells():
        if ('plate run' in metainfo) and not ('BG' in str(metainfo.loc[metainfo['Well']==well,'strain'])):
            begin=float(metainfo[metainfo['Well']==well]['plate run - start time'])
            time_idx=pl.find_idx_to_closest_value(dataOD['Time'],begin)
            dataOD['OD']=dataOD[well]
            tmp_data=pd.DataFrame()
            tmp_data['Time']=dataOD['Time'].add(-dataOD.loc[time_idx,'Time'])
            tmp_data['OD']=dataOD['OD']
            if maxtime:
                #shorten all measurements
                maxtime_idx= pl.find_idx_to_closest_value(tmp_data['Time'],maxtime)+1
                ODs.append(tmp_data.iloc[time_idx:maxtime_idx,:].reset_index(drop=True))
            else:
                ODs.append(tmp_data.iloc[time_idx:,:].reset_index(drop=True))
        else:
            if maxtime:
                maxtime_idx= pl.find_idx_to_closest_value(dataOD['Time'],maxtime)+1
                tmpOD=dataOD.iloc[0:maxtime_idx,:].reset_index(drop=True)
                tmpLUM=dataLUM.iloc[0:maxtime_idx,:].reset_index(drop=True)
            else:
                tmpOD=dataOD
                tmpLUM=dataLUM
            tmpOD['OD']=tmpOD[well]
            ODs.append(tmpOD[['Time','OD']])
            if lum_od and lum:
                if lum_od:
                    lum_col='LUM/OD'
                else:
                    lum_col='LUM'
                tmpLUM[lum_col]=tmpLUM[well]
                LUMs.append(tmpLUM[['Time',lum_col]])
        
    metainfo['OD']=ODs
    if lum_od and lum:
        metainfo['LUM/OD']=LUMs
    elif lum:
        metainfo['LUM']=LUMs
    return(metainfo)

    
def background_correct(reader,data,metainfo,dtype):
    """Substracts the measured background from OD values. Luminescence data is not background corrected.
    
    Uses background wells of the specific experiment if available. If not it automatically uses a background standard value for reader & media, but again only if it is available.
    """
    if dtype=='LUM': #don't background correct luminescence data
        return data
    for medium in set(metainfo.media):
        if medium[:4]=='MOPS':
            medium='MOPS'
        elif medium=='':
            continue
        background = metainfo.loc[(metainfo['strain'] == 'BG') & (metainfo['media'].str.match(medium))]
        average= []
        for well in background['Well']:
            average.extend(list(data[well]))#average.append(np.median(list(data[well])))#
        if average:#when there was background measured on the plate take that one
            average= removeOutliers(average)
            average=np.median(average)#average=np.mean(average)#average=min(average)#
            print('OD Background: {}'.format(average))
            if medium=='MOPS' and reader=='CLARIOstar':
                average=0.12
                print('WARNING: OD background fixed to {}'.format(average))
        else:#if no background was measured take a standard one
            if reader=='biotek':
                if dtype=='OD':
                    lookupOD={'MHB':0.0902,'MOPS':0.08366}
                    if medium in lookupOD:
                        average=lookupOD[medium]
                    else:
                        print('No standard {} background saved for {} and {}.'.format(dtype,reader,medium))
                        return data
                elif dtype=='LUM':
                    lookupLUM={'MHB':6.0,'MOPS':5.3}                    
                    if medium in lookupLUM:
                        average=lookupLUM[medium]
                    else:
                        print('No standard {} background saved for {} and {}.'.format(dtype,reader,medium))
                        return data
            elif reader=='CLARIOstar':
                if dtype=='OD':
                    lookupOD={'MOPS':0.12}
                    if medium in lookupOD:
                        average=lookupOD[medium]
                    else:
                        print('No standard {} background saved for {} and {}.'.format(dtype,reader,medium))
                        return data
                else:
                    print('no lum background saved for CLARIOstar')
                    return data
            elif reader=='FLUOstar Omega':
                if dtype=='OD':
                    lookupOD={'MOPS':0.073, 'LB':0.079}
                    if medium in lookupOD:
                        average=lookupOD[medium]
                    else:
                        print('No standard {} background saved for {} and {}.'.format(dtype,reader,medium))
                        return data
                else:
                    print('FLUOstar Omega cannot read luminescence data. But somehow luminescence-background is supposed to be corrected. Recheck all options.')
                    return data
            elif reader=='SPECTROstar':
                if dtype=='OD':
                    lookupOD={'MOPS':0.0768}
                    if medium in lookupOD:
                        average=lookupOD[medium]
                    else:
                        print('No standard {} background saved for {} and {}.'.format(dtype,reader,medium))
                        return data
                else:
                    print('SPECTROstar cannot read luminescence data. But somehow luminescence-background is supposed to be corrected. Recheck all options.')
                    return data
            elif reader=='SPECTROstar TaylorLab':
                if dtype=='OD':
                    lookupOD={'MOPS':0.0768}
                    if medium in lookupOD:
                        average=lookupOD[medium]
                    else:
                        print('No standard {} background saved for {} and {}.'.format(dtype,reader,medium))
                        return data
                else:
                    print('SPECTROstar cannot read luminescence data. But somehow luminescence-background is supposed to be corrected. Recheck all options.')
                    return data
            else:
                print('No standard {} background saved for {} and {}.'.format(dtype,reader,medium))
                return data
        wells=list(metainfo.loc[metainfo['media'].str.match(medium)]['Well'])
        data[wells]=data[wells]-average
    #thresholding
    for well in wells:#avoid negative ODs
        data[well].loc[data[well]<=0] = 0
    return data

def removeOutliers(x, outlierConstant=1.5):
    '''https://gist.github.com/vishalkuo/f4aec300cf6252ed28d3#gistcomment-2847578'''
    a = np.array(x)
    if np.isnan(np.sum(a)):
        a = a[~np.isnan(a)] # just remove nan elements from vector
    upper_quartile = np.percentile(a, 75)
    lower_quartile = np.percentile(a, 25)
    IQR = (upper_quartile - lower_quartile) * outlierConstant
    quartileSet = (lower_quartile - IQR, upper_quartile + IQR)
    

    result = a[np.where((a >= quartileSet[0]) & (a <= quartileSet[1]))]
    
    return result.tolist()

def determine_doublingtime(dataOD,metainfo):
    """Determines doublingtime between OD 0.05 and 0.2.
    Adds new column to metainfo-dataframe.
    """
    metainfo['doublingtime'] = np.full(len(metainfo), np.nan)
    metainfo['unperturbed doublingtime'] = np.full(len(metainfo), np.nan)
    #determine doublingtime between OD0.05-0.2
    for i_wells in range(len(metainfo)):
        well=metainfo['Well'].loc[i_wells]
        beg=False
        end=False
        for t in range(len(dataOD)-1,2,-1):#ignore the first two timepoints in case of high valued artifacts
            if dataOD.loc[t,well]<=0.05:#0.05
                beg=[dataOD.loc[t,well],dataOD.loc[t,'Time']]
                break
            elif dataOD.loc[t,well]>=0.2:#0.2
                end=[dataOD.loc[t,well],dataOD.loc[t,'Time']]
        if beg and end:    
            doubling_t=((end[1]-beg[1])*np.log(2))/(np.log(end[0])-np.log(beg[0]))
            metainfo.loc[i_wells,'doublingtime']=doubling_t
        del beg,end
        
    #insert doubling time of respective growth control
    for i_wells in range(len(metainfo)):
        growthControl = metainfo.loc[i_wells,'control well']
        growthControl = float(metainfo['doublingtime'].loc[metainfo['Well']==growthControl])
        metainfo.loc[i_wells,'unperturbed doublingtime'] = growthControl
    return metainfo

def data_excludes_but(data,data_exclude_but,plotting_variations,category_to_plot,variations_in_category,clean_replicates=True):
    """Removes all data from the set that has other variations in the given categories than the one defined.
        
    plotting_variations: dict.keys: categories; values: lists of variations you want to keep
    data_exclude_but: dict. keys: categories; values: lists of variations you want to keep
    """
    
    if plotting_variations and clean_replicates:
        for cat,var in zip(category_to_plot,variations_in_category):
            data=data.reset_index(drop=True) 
            if any(option in list(data[cat]) for option in plotting_variations.keys()):
                for i in range(len(data)):
                    if data.loc[i,cat] in plotting_variations.keys():
                        if not data.loc[i,var] in plotting_variations[data.loc[i,cat]]:
                            data = data.drop(i, axis=0)
    if data_exclude_but:
        for option in data_exclude_but:
            data=data.reset_index(drop=True)
            for i in range(len(data)):
                if not data.loc[i,option] in data_exclude_but[option]:
                    data = data.drop(i, axis=0)
    data=data.reset_index(drop=True)
    return data

def data_remove_all(data,data_remove):
    """Removes all data from the set that has any of the given variations in the given categories.
    Always removes empty wells from the analysis ('strain'=='BG'/''/'water')
    
    data_remove: dict. keys: categories; values: lists of variations you want to delete
    """    
    if data_remove:
        remainder = pd.DataFrame()
        for option in data_remove:
            remainder = pd.concat([remainder,data.loc[~(data[option].isin(data_remove[option]))]],ignore_index=True,sort=False)
        remainder = pd.concat([remainder,data.loc[~(data['strain'].isin(['BG','','water']))]],ignore_index=True,sort=False)
        return remainder
    else:
        remainder = pd.DataFrame()
        remainder = pd.concat([remainder,data.loc[~(data['strain'].isin(['BG','','water']))]],ignore_index=True,sort=False)
        return remainder

def findreader(file):
    """Automatically identifies the plate reader and with it the format of the given data."""    
    table=pd.read_excel(file, header=None)
    if table[1][8]=='Synergy H1':
        return 'biotek'
    elif any([table[0][i]=='Well\nRow' and table[1][i]=='Well\nCol' for i in range(3,20)]):
        return 'SPECTROstar'
    elif any([table[0][i]=='Well' and table[1][i]=='Content' for i in range(5,20)]):
        if 'Omega' in table[0][1]:
            return 'FLUOstar Omega'
        elif 'CLARIOstar' in table[0][1]:       
            return 'CLARIOstar'
        elif 'SPECTROstar Nano' in table[0][1]:
            return 'SPECTROstar TaylorLab'
    elif table[0][0]=='Plate':
        return 'victor'
    else:
        return(print('Plate reader not recoginzed'))
        
def find_timeconversion_factor(unit_in,unit_out):
    """Returns the factor, with which to muliply a time to get another unit.
    
    unit_in, unit_out: str. options: 's','m','h'
    """
    if unit_out=='doublings':
        unit_out='h'
    if unit_in==unit_out:
        return 1
    elif (unit_in=='s' and unit_out=='m') or (unit_in=='m' and unit_out=='h'):
        return 0.016666666666666666
    elif (unit_in=='s' and unit_out=='h'):
        return 0.0002777777777777778
    elif (unit_in=='m' and unit_out=='s') or (unit_in=='h' and unit_out=='m'):
        return 60
    elif (unit_in=='h' and unit_out=='s'):
        return 3600
    else:
        print('Error in time conversion')
        return

def translateStrainID(data,TranslateTo='combined'):
    """Replaces strain ID numbers with their genotype. Genotype must be listed in excel-table 'GenotypeDictonary.xlsx', which is in the same folder as this script.
    
    TranslateTo: str. options: 'combined','genotype', 'reporter'
    """
    table=pd.read_excel('GenotypeDictonary.xlsx')
    if TranslateTo=='combined':
        values=list()
        v1=list(table['genotype'])
        v1=[' ' if not isinstance(x, str) else x for x in v1]
        v2=list(table['reporter'])
        v2=[' ' if not isinstance(x, str) else x for x in v2]
        for i in range(len(v1)):
            values.append(' '.join([v1[i],v2[i]]).strip())
    else:
        try:
            values=list(table[TranslateTo])
        except KeyError:
            print('strain-ID cannot be translated to {}'.format(TranslateTo))
            return(data)
    keys=list(table['strain ID'])
    for i in range(len(set(data['strain']))):
        dictionary = dict(zip(keys, values))
    data['strain']=data['strain'].replace(dictionary)
    return(data)
    
    
    
    
    
    