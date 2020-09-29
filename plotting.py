#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 17:54:08 2019

@author: angelika
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
new_rc_params = {'text.usetex': False,"svg.fonttype": 'none'}
plt.rcParams.update(new_rc_params)
import numpy as np
import math
import data_import as di
# import timeit
from scipy.stats import ttest_ind



def select_data(data,category_to_plot,variations_in_category,category1,category2,plotting_variations,clean_replicates):
    selection = data.loc[data[category_to_plot[0]]==category1]
    if not variations_in_category[0]==category_to_plot[0]:
        if category1 in plotting_variations and clean_replicates:
            selection = selection.loc[selection[variations_in_category[0]].isin(plotting_variations[category1])]
    if category2:
        selection = selection.loc[selection[category_to_plot[1]]==category2]
        if not variations_in_category[1]==category_to_plot[1] and clean_replicates:
            selection = selection.loc[selection[variations_in_category[1]].isin(plotting_variations[category2])]    
    #check that every replicate line is complete, if some concentrations were also measured together in another dilution series - delete that data
    if clean_replicates:
        if not variations_in_category[0]==category_to_plot[0]:
            selection_clean=pd.DataFrame()
            for cur_date in set(selection['Date']):
                if all(elem in set(selection.loc[selection['Date']==cur_date][variations_in_category[0]])    for elem in plotting_variations[category1]):
                    selection_clean = pd.concat([selection_clean,selection.loc[selection['Date']==cur_date]],ignore_index=True,sort=False)
        else:
            selection_clean=selection                        
        if not category2: #if only one category is used
            return selection_clean

        elif category2: #if more than one category is given
            if not variations_in_category[1]==category_to_plot[1]:
                selection_clean2=pd.DataFrame()
                for cur_date in set(selection_clean['Date']):
                    if all(elem in set(selection_clean.loc[selection_clean['Date']==cur_date][variations_in_category[1]])    for elem in plotting_variations[category2]):
                        selection_clean2 = pd.concat([selection_clean2,selection_clean.loc[selection_clean['Date']==cur_date]],ignore_index=True,sort=False)
                return selection_clean2
            else:
                return selection_clean
    else:
        return selection


def plot_timedependent_replicates(data,plot_object,category_to_plot,variations_in_category,plotting_variations,time_unit,save,devices,clean_replicates=False):
    cmap = sns.color_palette(palette='colorblind')
    translationTable = str.maketrans("μα", "ua")

    if plot_object=='LUM':
        if 'LUM' in data.columns:
            plot_object='LUM'
        elif 'LUM/OD' in data.columns:
            plot_object='LUM/OD'
        else:
            print('Data ({}) not found.'.format(plot_object))
            return
    elif plot_object!='OD':
        print('Data ({}) not found.'.format(plot_object))
        return
            
     
    for category1 in set(data[category_to_plot[0]]):
        for category2 in set(data[category_to_plot[1]]):
    
            #extract revelevant data
            selection = select_data(data,category_to_plot,variations_in_category,category1,category2,plotting_variations,clean_replicates)
            if selection.empty:
                continue
            

            for day in set(selection['Date']):
                for replicate,num in zip(set(selection['unperturbed doublingtime']),range(len(set(selection['unperturbed doublingtime'])))):#if only one replicate was measured that day
                    for variation2 in get_variation_of_category(selection,category_to_plot[1],variations_in_category[1],category2):
                        for variation1 in get_variation_of_category(selection,category_to_plot[0],variations_in_category[0],category1):
                            found_data =False
                            var_data=selection.loc[(selection['Date']==day) & (selection['unperturbed doublingtime']==replicate) & (selection[variations_in_category[0]]==variation1) & (selection[variations_in_category[1]]==variation2)]
                            if var_data.empty:
                                continue 
                            for ind in list(var_data.index.values):
                                found_data=True
                                time_factor = time_modulation(time_unit,var_data.loc[ind,'unperturbed doublingtime'])
                                var_data.loc[ind,plot_object].Time = var_data.loc[ind,plot_object].Time*time_factor
                                sns.lineplot(data=var_data.loc[ind,plot_object],x='Time',y=plot_object,label='{} {}'.format(variation1,get_unit(category1)),palette=cmap)
                        if not found_data:
                            continue
                        if plot_object=='OD':
                            plt.ylim(0.001,2)
                        elif plot_object=='LUM/OD':
                            plt.ylim(10,6*(10**5))
                        elif plot_object=='LUM':
                            plt.ylim(1,3*(10**3))
                        plt.grid(b=True, which = 'major')
                        plt.grid(b=True, which = 'minor',linewidth=0.5)
                        
                        plt.ylabel(find_ylabel(plot_object))
                        plt.xlabel('time [{}]'.format(time_unit))
                        plt.yscale('log')
                        plt.title('{}; {} - {}'.format(get_title(category2,variation2),remove_unit(category1),day))
                        # Shrink current axis by 20%
                        ax=plt.gca()
                        box = ax.get_position()
                        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                        # Put a legend to the right of the current axis
                        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                        plt.savefig('{}timedependent{}_{}_{}_Rep{}.svg'.format(save,plot_object.replace('/','-'),get_save_cat(category2,variation2),remove_unit(category1),num).translate(translationTable), bbox_inches='tight',dpi=300)
                        plt.savefig('{}timedependent{}_{}_{}_Rep{}.png'.format(save,plot_object.replace('/','-'),get_save_cat(category2,variation2),remove_unit(category1),num).translate(translationTable), bbox_inches='tight',dpi=300)
                        plt.close()#show()
    return

def plot_timedependent_averaged(data,plot_object,category_to_plot,variations_in_category,plotting_variations,time_unit,save,devices,clean_replicates=False):
    translationTable = str.maketrans("μα", "ua")

    if plot_object=='LUM':
        if 'LUM' in data.columns:
            plot_object='LUM'
        elif 'LUM/OD' in data.columns:
            plot_object='LUM/OD'
        else:
            print('Data ({}) not found.'.format(plot_object))
            return
    elif plot_object!='OD':
        print('Data ({}) not found.'.format(plot_object))
        return


    for category1 in set(data[category_to_plot[0]]):
        for category2 in set(data[category_to_plot[1]]):
            #extract relevant data
            selection = select_data(data,category_to_plot,variations_in_category,category1,category2,plotting_variations,clean_replicates)
            if selection.empty:
                continue
           
            for variation2 in get_variation_of_category(selection,category_to_plot[1],variations_in_category[1],category2):
                dfplot=pd.DataFrame()
                for variation1 in get_variation_of_category(selection,category_to_plot[0],variations_in_category[0],category1):
                    var_data=selection.loc[(selection[variations_in_category[0]]==variation1) & (selection[variations_in_category[1]]==variation2)]
                    if var_data.empty:
                        continue 
                    
                    time_factor = time_modulation(time_unit,np.mean(var_data['unperturbed doublingtime']))
                    for ind,exp in zip(list(var_data.index.values),range(len(var_data.index.values))):
                        try:
                            time_data = var_data.loc[ind,plot_object].Time*time_factor
                            labelcol = ['{} {}'.format(variation1,get_unit(category1))] * len(time_data)
                            tmp_df = pd.DataFrame(list(zip(time_data,var_data.loc[ind,plot_object][plot_object],labelcol)),columns=['Time /{}'.format(time_unit),plot_object,remove_unit(category1)])
                            dfplot = pd.concat([dfplot,tmp_df],sort=False)
                        except AttributeError:
                            pass
                cmap=sns.color_palette(palette='colorblind',n_colors=len(set(dfplot[remove_unit(category1)])))
                # stop1 = timeit.default_timer()
                sns.lineplot(data=dfplot, x='Time /{}'.format(time_unit) , y=plot_object, hue=remove_unit(category1),palette=cmap)
                # stop2 = timeit.default_timer()
                if plot_object=='OD':
                    plt.ylim(0.001,2)
                elif plot_object=='LUM/OD':
                    if 'CLARIOstar' in devices:
                        ylim_max=2*10**8
                    else:
                        ylim_max=3*(10**3)
                    if not 'CLARIOstar' in devices:
                        ylim_min=5*10**0
                    elif len(devices)>1:
                        ylim_min=5*10**0
                    else:
                        ylim_min=10**3
                    plt.ylim(ylim_min,ylim_max)
                elif plot_object=='LUM':
                    if 'CLARIOstar' in devices:
                        ylim_max=2*10**8
                    else:
                        ylim_max=2*10**5
                    if not 'CLARIOstar' in devices:
                        ylim_min=1
                    elif len(devices)>1:
                        ylim_min=1
                    else:
                        ylim_min=1
                    plt.ylim(ylim_min,ylim_max)
        
                plt.grid(b=True, which = 'major')
                plt.grid(b=True, which = 'minor',linewidth=0.5)
                plt.ylabel(find_ylabel(plot_object))
                plt.xlabel('time [{}]'.format(time_unit))
                plt.yscale('log')
                plt.title('{}; {} (n={})'.format(get_title(category2,variation2),remove_unit(category1),len(var_data)))
                # Shrink current axis by 20%
                ax=plt.gca()
                box = ax.get_position()
                ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
                # Put a legend to the right of the current axis
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                # stop3 = timeit.default_timer()
                plt.savefig('{}AVGtimedependent{}_{}_{}.svg'.format(save,plot_object.replace('/','-'),get_save_cat(category2,variation2),remove_unit(category1)).translate(translationTable), bbox_inches='tight',dpi=300)
                plt.savefig('{}AVGtimedependent{}_{}_{}.png'.format(save,plot_object.replace('/','-'),get_save_cat(category2,variation2),remove_unit(category1)).translate(translationTable), bbox_inches='tight',dpi=300)
                # stop4= timeit.default_timer()
                plt.close()#show()
                
                # print('Preprocessing: {};\n Plotting: {};\n Postprocessing: {};\n Saving: {};\n Total: {}'.format(stop1-start,stop2-stop1,stop3-stop2,stop4-stop3,stop4-start))
              
            print('plotted averaged time dependent behaviour of {} with {}'.format(plot_object,' & '.join([str(category1),str(category2)])))
    return
    
def time_modulation(time_unit,doublingt0):
    if time_unit=='h':
        factor=1
    elif time_unit=='min':
        factor=60
    elif time_unit=='doublings':
        factor=1/doublingt0
    return factor

def find_ylabel(plot_object):
    if plot_object=='OD':
        return 'OD'
    elif plot_object=='LUM':
        return 'luciferase activity [RLU]'
    elif plot_object=='LUM/OD':
        return 'luciferase activity [RLU/OD]'
    
def find_idx_to_closest_value(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def plot_IC_interaction(data,plotting_variations,timepoint,save):
    '''
    Plots a heatmap of the normalized OD with 'Inducer1' and 'Inducer2' as axis

    Parameters
    ----------
    data : TYPE
        DESCRIPTION.
    plotting_variations : TYPE
        DESCRIPTION.
    timepoint : TYPE
        DESCRIPTION.
    save : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    translationTable = str.maketrans("μα", "ua")   

    ODs=[['Inducer1','Inducer1_Concentration','Inducer2','Inducer2_Concentration','OD','OD normalized']]
    for i in range(len(data)):
        tmp=list()
        tmp.extend([data.iloc[i].Inducer1,data.iloc[i].Inducer1_Concentration,data.iloc[i].Inducer2,data.iloc[i].Inducer2_Concentration])

        dt0=data.iloc[i]['unperturbed doublingtime']
        timepoint_h=find_timepoint(timepoint,dt0)
        
        cur_timepoint = find_idx_to_closest_value(data.loc[i,'OD'].Time,timepoint_h)
        od_challenged=float(data.loc[i,'OD'].OD[cur_timepoint])
        control=data.loc[np.logical_and(data['Well']==data.loc[i,'control well'],data['Date']==data.loc[i,'Date'])]
        c_i=list(control.index)
        if len(c_i)>1:
            print('Well has several options as control well. Please check metainfo data.')
        else:
            c_i=int(list(c_i)[0])
        dt0=data.iloc[c_i]['unperturbed doublingtime']
        timepoint_h=find_timepoint(timepoint,dt0)
        cur_timepoint = find_idx_to_closest_value(data.loc[c_i,'OD'].Time,timepoint_h)
        od_control=float(control.loc[c_i,'OD'].OD[cur_timepoint])
        tmp.extend([od_challenged,od_challenged/od_control*100])
        ODs.append(tmp)
        
    ODs=pd.DataFrame(ODs[1:],columns=ODs[0])
    
    for ind1 in set(ODs['Inducer1']):
        for ind2 in set(ODs['Inducer2']):
            specify=ODs.loc[np.logical_and(ODs['Inducer1']==ind1,ODs['Inducer2']==ind2)]
            if len(specify)==0:
                continue
            specify=specify.pivot('Inducer1_Concentration','Inducer2_Concentration','OD normalized')
            ax=sns.heatmap(specify,vmax=100)
            ax.invert_yaxis()
            plt.xlabel(ind2)
            plt.ylabel(ind1)
            plt.savefig('{}{}TimepointIC_{}{}.svg'.format(save,timepoint.replace(' ','_'),remove_unit(ind1),remove_unit(ind2)).translate(translationTable), bbox_inches='tight',dpi=300)
            plt.savefig('{}{}TimepointIC_{}{}.png'.format(save,timepoint.replace(' ','_'),remove_unit(ind1),remove_unit(ind2)).translate(translationTable), bbox_inches='tight',dpi=300)
            plt.close()#show()
    return
        
def find_timepoint(timepoint,dt0):
    '''time in hours at which the MIC is to be measured'''
    if timepoint.split(' ')[-1]=='doublings':
        timepoint_h = dt0*float(timepoint.split(' ')[0])
    elif timepoint.split(' ')[-1]=='h':
        timepoint_h = float(timepoint.split(' ')[0])
    else:
        print('unit of timepoint not recognized')
        return
    return timepoint_h


def plot_ICs(data,category_to_plot,variations_in_category,plotting_variations,timepoint,save,x_axis=False,continuous_xaxis=True,IC=0.3,normalize='min',ylog=False,ttest=False,singles=False):
    '''
    Plots the concentation inhibiting growth to a given percentage after a given timeframe (doublings or hours).

    Parameters
    ----------
    x_axis : TYPE, optional
        DESCRIPTION. The default is False.
    continuous_xaxis : TYPE, optional
        DESCRIPTION. The default is True.
    IC : TYPE, optional
        DESCRIPTION. The default is 0.3.
    normalize : TYPE, optional
        DESCRIPTION. The default is 'min'.
    ylog : TYPE, optional
        DESCRIPTION. The default is False.
    ttest : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    '''
    
    dfpointIC,unit=determine_IC(category_to_plot,variations_in_category,data,timepoint,plotting_variations,IC,estimateIC=True)
    
    if ttest:
        for strain in set(dfpointIC['strain']):
            for ab in set(dfpointIC['Inducer1']):
                sel=dfpointIC.loc[dfpointIC['strain']==strain]
                mic=list(sel.loc[sel['Inducer1']==ab,'Concentration'])
                wt=dfpointIC.loc[dfpointIC['strain']=='WT']
                wtmic=list(wt.loc[wt['Inducer1']==ab,'Concentration'])
                
                stat, p = ttest_ind(mic, wtmic,equal_var=False)
        #        stat, p = f_oneway(mic, wtmic)
                print('stat=%.3f, p=%.3f' % (stat, p))
                if p > 0.05:
                    print('{}, {} has probably the same distribution as WT'.format(strain,ab))
                else:
                    print('{}, {} has probably different distributions as WT'.format(strain,ab))
    
    #IC normalization
    for category1 in set(dfpointIC[category_to_plot[0]]):
        if normalize:
            #find the ic to normalize with
            ic0= np.mean(dfpointIC.loc[(dfpointIC[category_to_plot[0]]==category1) & (dfpointIC[variations_in_category[1]]==norm_with(dfpointIC[variations_in_category[1]], normalize)),'inhibitory concentration'])#
            if math.isnan(ic0):# if the respective condition has no ic value
                x_values=set(dfpointIC[variations_in_category[1]])
                while math.isnan(ic0): # get the ic value of the next best condition
                    x_values.remove(norm_with(x_values, normalize)) 
                    ic0= np.mean(dfpointIC.loc[(dfpointIC[category_to_plot[0]]==category1) & (dfpointIC[variations_in_category[1]]==norm_with(x_values, normalize)),'inhibitory concentration'])
                    if len(x_values)<2: #
                        ic0=100
            dfpointIC.loc[dfpointIC[category_to_plot[0]] == category1, 'IC [%]'] = dfpointIC.loc[dfpointIC[category_to_plot[0]] == category1,'inhibitory concentration']/ic0*100
            dfpointIC.loc[dfpointIC[category_to_plot[0]] == category1, category_to_plot[0]] = '{} ({} {})'.format(category1, np.round(ic0,2), unit[remove_unit(category1)])
            y_axis='IC [%]'
        else:
            y_axis='inhibitory concentration'        
        
    save_ICs(dfpointIC,save,timepoint,category_to_plot,variations_in_category)
    
    if singles:
        plot_singleICs(dfpointIC,category_to_plot,variations_in_category,IC,timepoint,x_axis,continuous_xaxis,y_axis,save)
    else:
        plot_allICs(dfpointIC,category_to_plot,variations_in_category,IC,timepoint,x_axis,continuous_xaxis,y_axis,ylog,normalize,save)
    
    print('IC dose response plottet for timpoint {}'.format(timepoint))
    return

def plot_singleICs(dfpointIC,category_to_plot,variations_in_category,IC,timepoint,x_axis,continuous_xaxis,y_axis,save):
    translationTable = str.maketrans("μα", "ua")
    for category1 in set(dfpointIC[category_to_plot[0]]):
        selection=dfpointIC.loc[dfpointIC[category_to_plot[0]]==category1]
        
        if x_axis=='growth rate [h^-1]':
            # cmap = sns.color_palette(palette='colorblind', n_colors=len(set(dfpointIC[category_to_plot[0]])))
            sns.scatterplot(data=selection,x='growth rate [h^-1]',y=y_axis,hue=category_to_plot[0])
            sns.lineplot(data=selection,x='mean growth rate [h^-1]',y=y_axis,hue=category_to_plot[0],legend=False)
        else:
            if continuous_xaxis:
                # cmap = sns.color_palette(palette='colorblind', n_colors=len(set(dfpointIC[category_to_plot[0]])))
                sns.scatterplot(data=selection,x=x_axis,y=y_axis,hue=category_to_plot[0])
                sns.lineplot(data=selection,x=x_axis,y=y_axis,hue=category_to_plot[0],legend=False)
            else:
                # cmap = sns.color_palette(palette='colorblind', n_colors=len(set(dfpointIC[category_to_plot[0]])))
                sns.barplot(data=selection,x=variations_in_category[1],y=y_axis,hue=category_to_plot[0],ci=None)
                sns.swarmplot(x=variations_in_category[1], y=y_axis, hue=category_to_plot[0], data=selection, dodge=True,color='black')
                
    
    
        plt.ylim(0,max(selection[y_axis])*1.1)
        plt.ylabel('IC{}'.format(IC*100))
                
        plt.grid(which='both',axis='y')
        
        
        actual_times=list(set(selection['IC after']))
        actual_times.sort()
        actual_times = [str(i) for i in actual_times]
        plt.title('IC{} determined after {} ({} h)'.format(int(IC*100),timepoint,' h,'.join(actual_times)))
       
        
        #increase figure width when plotting many categories on xaxis 
        if x_axis!='growth rate [h^-1]' and not continuous_xaxis:
            if len(set(selection[category_to_plot[1]]))>5:
                fig=plt.gcf()
                size=fig.get_size_inches()
                fig.set_size_inches(size[0]/5*len(set(selection[category_to_plot[1]])),size[1])
                plt.xticks(rotation=-90)
            else:
                plt.xticks(rotation=-20)
        # Shrink current axis by 20%
        ax=plt.gca()
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
        # legend to the right of the current axis and remove entries for single replicates(=scatter plot)
        if continuous_xaxis==True or x_axis=='growth rate [h^-1]':        
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        elif len(ax.get_legend_handles_labels()[1])<len(set(selection[x_axis])):
            handles, labels = ax.get_legend_handles_labels()
            handles=handles[:len(set(selection[x_axis]))+1]
            labels=labels[:len(set(selection[x_axis]))+1]
            ax.legend(handles, labels,loc='center left', bbox_to_anchor=(1, 0.5))       
        else:
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[- len(set(selection[category_to_plot[0]])):], labels[- len(set(selection[category_to_plot[0]])):],loc='center left', bbox_to_anchor=(1, 0.5))
    
        plt.savefig('{}{}TimepointIC_{}_{}.svg'.format(save,timepoint.replace(' ','_'),remove_unit(category1).replace(' ',''),x_axis.replace(' ','')).translate(translationTable), bbox_inches='tight',dpi=300)
        plt.savefig('{}{}TimepointIC_{}_{}.png'.format(save,timepoint.replace(' ','_'),remove_unit(category1).replace(' ',''),x_axis.replace(' ','')).translate(translationTable), bbox_inches='tight',dpi=300)
        plt.close()#show()
    
    return

def plot_allICs(dfpointIC,category_to_plot,variations_in_category,IC,timepoint,x_axis,continuous_xaxis,y_axis,ylog,normalize,save):
    translationTable = str.maketrans("μα", "ua")
    if x_axis=='growth rate [h^-1]':
        cmap = sns.color_palette(palette='colorblind', n_colors=len(set(dfpointIC[category_to_plot[0]])))
        sns.scatterplot(data=dfpointIC,x='growth rate [h^-1]',y=y_axis,hue=category_to_plot[0])
        sns.lineplot(data=dfpointIC,x='mean growth rate [h^-1]',y=y_axis,hue=category_to_plot[0],legend=False,ci=None,palette=cmap)
    else:
        if continuous_xaxis:
            cmap = sns.color_palette(palette='colorblind', n_colors=len(set(dfpointIC[category_to_plot[0]])))
            sns.scatterplot(data=dfpointIC,x=x_axis,y=y_axis,hue=category_to_plot[0])
            sns.lineplot(data=dfpointIC,x=x_axis,y=y_axis,hue=category_to_plot[0],legend=False,ci=None,palette=cmap)
        else:
            cmap = sns.color_palette(palette='colorblind', n_colors=len(set(dfpointIC[category_to_plot[0]])))
            sns.barplot(data=dfpointIC,x=variations_in_category[1],y=y_axis,hue=category_to_plot[0],ci=None,palette=cmap)
            sns.swarmplot(x=variations_in_category[1], y=y_axis, hue=category_to_plot[0], data=dfpointIC, dodge=True,color='black')
            
    if ylog:
        plt.yscale('log')
        plt.ylim(min(dfpointIC[y_axis])*0.9,max(dfpointIC[y_axis])*1.1)
    else:
        plt.ylim(0,max(dfpointIC[y_axis])*1.1)
    if not normalize:
        plt.ylabel('IC{}'.format(IC*100))
            
    plt.grid(which='both',axis='y')
    
    
    actual_times=list(set(dfpointIC['IC after']))
    actual_times.sort()
    actual_times = [str(i) for i in actual_times]
    plt.title('IC{} determined after {} ({} h)'.format(int(IC*100),timepoint,' h,'.join(actual_times)))
   
    
    #increase figure width when plotting many categories on xaxis 
    if x_axis!='growth rate [h^-1]' and not continuous_xaxis:
        if len(set(dfpointIC[category_to_plot[1]]))>5:
            fig=plt.gcf()
            size=fig.get_size_inches()
            fig.set_size_inches(size[0]/5*len(set(dfpointIC[category_to_plot[1]])),size[1])
            plt.xticks(rotation=-90)
        else:
            plt.xticks(rotation=-20)
    # Shrink current axis by 20%
    ax=plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # legend to the right of the current axis and remove entries for single replicates(=scatter plot)
    if continuous_xaxis==True or x_axis=='growth rate [h^-1]':        
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    elif len(ax.get_legend_handles_labels()[1])<len(set(dfpointIC[x_axis])):
        handles, labels = ax.get_legend_handles_labels()
        handles=handles[:len(set(dfpointIC[x_axis]))+1]
        labels=labels[:len(set(dfpointIC[x_axis]))+1]
        ax.legend(handles, labels,loc='center left', bbox_to_anchor=(1, 0.5))       
    else:
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[- len(set(dfpointIC[category_to_plot[0]])):], labels[- len(set(dfpointIC[category_to_plot[0]])):],loc='center left', bbox_to_anchor=(1, 0.5))

    plt.savefig('{}{}TimepointIC_{}.svg'.format(save,timepoint.replace(' ','_'),x_axis.replace(' ','')).translate(translationTable), bbox_inches='tight',dpi=300)
    plt.savefig('{}{}TimepointIC_{}.png'.format(save,timepoint.replace(' ','_'),x_axis.replace(' ','')).translate(translationTable), bbox_inches='tight',dpi=300)
    plt.close()#show()
    
    return

def determine_IC(category_to_plot,variations_in_category,data,timepoint,plotting_variations,IC,estimateIC=True):
    columns = ['Date',category_to_plot[0],category_to_plot[1]]
    if category_to_plot[1]!=variations_in_category[1]:
        columns.append(variations_in_category[1])
    columns.extend(['mean growth rate [h^-1]','growth rate [h^-1]','inhibitory concentration','IC after','IC [%]'])
    dfpointIC = [columns]
    unit = dict()
    for category2 in set(data[category_to_plot[1]]):
        for category1 in set(data[category_to_plot[0]]):
            preselection = data.loc[np.logical_and(data[category_to_plot[1]]==category2,data[category_to_plot[0]]==category1)]
            if preselection.empty:
                continue
            for variation2 in set(preselection[variations_in_category[1]]):
                selection = preselection.loc[preselection[variations_in_category[1]]==variation2]
                unit[remove_unit(category1)] = get_unit(category1)
                dt0=np.mean(selection.loc[(selection[category_to_plot[0]]==category1),'unperturbed doublingtime'])
                #time in hours at which the MIC is to be measured
                if timepoint.split(' ')[-1]=='doublings':
                    timepoint_h = dt0*float(timepoint.split(' ')[0])
                elif timepoint.split(' ')[-1]=='h':
                    timepoint_h = float(timepoint.split(' ')[0])
                else:
                    print('unit of timepoint not recognized')
                    return
                
                
                #single out replicates
                for date in set(selection['Date']):
                    for rep in set(selection['control well']):
                        if not ((selection['Date']==date) & (selection['control well']==rep)).any():
                            continue
                        replicate = selection.loc[(selection['Date']==date) & (selection['control well']==rep)]
                        
                        # #remove replicates from analysis, which miss parts of the dataset 
                        # if category1 in plotting_variations:
                        #     if not all(elem in set(replicate[variations_in_category[0]])    for elem in plotting_variations[category1]):
                        #         continue
                        row=[date,category1,category2]
                        if category_to_plot[1]!=variations_in_category[1]:
                            row.append(variation2)
                        
                        replicate =replicate.sort_values(by=variations_in_category[0])
                        concentrations=[]
                        ods = []
                        for idx in list(replicate.index.values):
                            cur_timepoint = find_idx_to_closest_value(replicate.loc[idx,'OD'].Time,timepoint_h)
                            od_challenged=float(replicate.loc[idx,'OD'].OD[cur_timepoint])
                            concentrations.append(replicate.loc[idx,variations_in_category[0]])
                            ods.append(od_challenged)
                        odsCurTimepoint=pd.DataFrame()
                        odsCurTimepoint[variations_in_category[0]]=concentrations
                        odsCurTimepoint['OD']=ods
                        try:
                            od0=float(odsCurTimepoint.loc[odsCurTimepoint[variations_in_category[0]]==0].OD)
                        except TypeError:
                            print('Replicate has no control concentration')
                            od0=1
                        odsCurTimepoint['normalizedODs']=odsCurTimepoint['OD']/od0
                        if estimateIC:
                            m,n =cal_function(variations_in_category,odsCurTimepoint,IC)
                            if m==None:
                                continue
                            else:
                                measuredIC=(IC-n)/m
                        else:
                            measuredIC=min(odsCurTimepoint[odsCurTimepoint['normalizedODs']<IC][variations_in_category[0]])
                        row.extend([np.log(2)/dt0,np.log(2)/replicate.loc[idx,'unperturbed doublingtime'],measuredIC ,round(replicate.loc[idx,'OD'].Time[cur_timepoint],2), float('NaN')])
                        dfpointIC.append(row)
    dfpointIC=pd.DataFrame(dfpointIC[1:],columns=dfpointIC[0])
    return(dfpointIC,unit)

def func(x, m, n):
    """Standard linear function"""
    return m*x+n

def cal_function(variations_in_category,odsCurTimepoint,IC):
    lower_than_IC=max(odsCurTimepoint[odsCurTimepoint['normalizedODs']>IC][variations_in_category[0]])
    lower_than_IC=odsCurTimepoint[odsCurTimepoint[variations_in_category[0]]==lower_than_IC]
    try:
        higher_than_IC=min(odsCurTimepoint[odsCurTimepoint['normalizedODs']<IC][variations_in_category[0]])
        higher_than_IC=odsCurTimepoint[odsCurTimepoint[variations_in_category[0]]==higher_than_IC]
    except ValueError:
        return(None,None)
    x=float(higher_than_IC[variations_in_category[0]])-float(lower_than_IC[variations_in_category[0]])
    y=float(higher_than_IC['normalizedODs'])-float(lower_than_IC['normalizedODs'])
    m=y/x
    n=-m*float(lower_than_IC[variations_in_category[0]])+float(lower_than_IC['normalizedODs'])
    return(m,n)

def norm_with(column,min_or_max):
    """Returns the value to normalize with.
    
    min_or_max: str. option: 'min','max'
        no normalization  if min_or_max is neither 'min' nor 'max'
    """
    if len(column)==0:
        return 1
    elif min_or_max=='min':
        return min(column)
    elif min_or_max=='max':
        return max(column)
    else:
        return 1 # normalization with 1 == no normalization
        
def save_ICs(dfpointIC,save,timepoint,category_to_plot,variations_in_category):
    translationTable = str.maketrans("μα", "ua")
    if variations_in_category[1]==category_to_plot[1]:
        ICs=[[category_to_plot[1],category_to_plot[0],'IC mean','IC std','IC values']]
        for category2 in set(dfpointIC[category_to_plot[1]]):
            sel=dfpointIC.loc[dfpointIC[category_to_plot[1]]==category2]
            for category1 in set(sel[category_to_plot[0]]):
                ic=list(sel.loc[sel[category_to_plot[0]]==category1,'Concentration'])
                mean=np.mean(ic)
                std=np.std(ic)
                ICs.append([category2,category1,mean,std,ic])
    else:
        ICs=[[category_to_plot[1],variations_in_category[1],category_to_plot[0],'IC mean','IC std','IC values']]
        for category2 in set(dfpointIC[category_to_plot[1]]):
                pre_sel=dfpointIC.loc[dfpointIC[category_to_plot[1]]==category2]
                for variation2 in set(pre_sel[variations_in_category[1]]):
                    sel=pre_sel.loc[pre_sel[variations_in_category[1]]==variation2]
                    for category1 in set(sel[category_to_plot[0]]):
                        ic=list(sel.loc[sel[category_to_plot[0]]==category1,'inhibitory concentration'])
                        mean=np.mean(ic)
                        std=np.std(ic)
                        ICs.append([category2,variation2,']'.join(category1.split(']')[:-1])+']',mean,std,ic])
    ICs=pd.DataFrame(ICs[1:],columns=ICs[0])
    ICs.to_csv('{}{}TimepointIC.csv'.format(save,timepoint.replace(' ','_')).translate(translationTable),index=False)
    return

def remove_unit(string):
    string = str(string).split('[')[0]
    string = string.strip()
    return string
        
def get_unit(string):
    string=str(string).split('[')[1:]
    string='['.join(string)
    string=str(string).split(']')[:-1]
    string=']'.join(string)
    return string

def get_variation_of_category(data,category_to_plot,variations_in_category,category):
    selection = data.loc[data[category_to_plot]==category]
    concentrations =list(set(selection[variations_in_category]))
    concentrations.sort(reverse=True)
    return concentrations
    
def get_title(category,variation):
    cat = remove_unit(category)
    unit = get_unit(category)
    if cat==str(variation):
        return cat
    elif unit:
        if unit[0]=='%':
            return '{} {}{}'.format(cat,variation,unit)
        else:
            return '{} {} {}'.format(cat,variation,unit)
    else:
        return '{} {} {}'.format(cat,variation,unit)
    
def get_save_cat(category,variation):
    cat = remove_unit(category)
    unit = get_unit(category)
    cat = str(cat).replace('.',',').replace(' ','')
    variation = str(variation).replace('.',',').replace(' ','')
    if cat==variation:
        return cat
    elif unit and unit[0]=='%':
        return '{}{}'.format(cat,variation)
    else:
        return '{}{}{}'.format(cat.replace(' ',''),variation,unit.replace(' ','').replace('/',''))
    
def plot_doublingtime(data,save,time_unit,xaxis,csv=False):
    control_only=data.copy()
    for i in range(len(control_only)):
        if control_only.loc[i,'Well']!=control_only.loc[i,'control well']:
            control_only = control_only.drop(i, axis=0)
    sns.barplot(data=control_only,x=xaxis,y='unperturbed doublingtime',ci=None)
    sns.swarmplot(data=control_only,x=xaxis,y='unperturbed doublingtime', dodge=True,color='black')
    plt.ylabel('unperturbed doublingtime [1/{}]'.format(time_unit))
    if len(set(control_only[xaxis]))>5:
        plt.xticks(rotation=-90)
    plt.savefig('{}Doublingtime-{}.svg'.format(save,xaxis), bbox_inches='tight',dpi=300)
    plt.savefig('{}Doublingtime-{}.png'.format(save,xaxis), bbox_inches='tight',dpi=300)
    plt.close()#show()
    if csv:
        table=[[xaxis,'doublingtime mean [h]','doublingtime std','doublingtime replicates [h]']]
        for cat in set(control_only[xaxis]):
            dts=list(control_only.loc[control_only[xaxis]==cat,'unperturbed doublingtime'])
            table.append([cat,np.mean(dts),np.std(dts),dts])
        table=pd.DataFrame(data=table[1:],columns=table[0])
        table.to_csv('{}Doublingtimes_{}.csv'.format(save,xaxis),index=False,sep=';')
    print('Doubling time plotted')
    return

def doublingtime_growthrate_conversion(one):
    other = np.log(2)/one
    return other

def plot_dose_response(data,save,category_to_plot,variations_in_category,devices,y='LUM',normalizeX=False,timepoints=['gr'],SaveInduction=False):
    """Dose response plots of either luciferase activity (lum) or OD (normalized to the control well) at the specified timepoint.
    
    Keyword arguments
    y: Specifies the y-axis; can be 'LUM' or 'OD'
    normalizeX: specifies if and with what the X axis should be normalized. Either False, or timepoint of MIC - plot_MICs of that timepoint has to be executed earlier
    timepoints: timepoints for analysis, list of strings with time and unit. e.g '10 min', '5 doublings' 
    SaveInduction: boolean, saving fold-induction of every concentration
    """
    translationTable = str.maketrans("μα", "ua")
    dfDoseResponse=[[]]
    

    
    #Data analysis
    for i in range(2): # only add variation col if its not the same as the category (otherwise error)
        if category_to_plot[i]==variations_in_category[i]:
            dfDoseResponse[0].append(category_to_plot[i])
        else:
            dfDoseResponse[0].extend([category_to_plot[i],variations_in_category[i]])
            
    if y!='growth rate':
        for t in timepoints:
            dfDoseResponse[0].append('Response after {}'.format(t))
    else:
        dfDoseResponse[0].append('Response')
    for i in range(len(data)):
        lum=True
        tmp=[]
        for k in range(2):
            if category_to_plot[k]==variations_in_category[k]:
                tmp.append(data.loc[i,category_to_plot[k]])
            else:
                tmp.extend([data.loc[i,category_to_plot[k]],data.loc[i,variations_in_category[k]]])
        if y=='growth rate':
            if math.isnan(data.loc[i,'doublingtime']):
                relative_gr=0
            else:                
                relative_gr=100*(doublingtime_growthrate_conversion(data.loc[i,'doublingtime'])/doublingtime_growthrate_conversion(data.loc[i,'unperturbed doublingtime']))
            tmp.append(relative_gr)
            dfDoseResponse.append(tmp)
        else:
            for t in timepoints:#give time for measurement in minutes
                dt0=data.loc[i,'unperturbed doublingtime']
                #time in hours at which the MIC is to be measured
                if t.split(' ')[-1]=='doublings':
                    t_h = dt0*float(t.split(' ')[0])
                elif t.split(' ')[-1]=='h':
                    t_h = float(t.split(' ')[0])
                elif t.split(' ')[-1]=='min':
                    t_h = float(t.split(' ')[0])*di.find_timeconversion_factor('m', 'h')
                else:
                    print('unit of timepoint not recognized')
                    return
                
                
                
                if y=='LUM':
                    try:
                        idx=find_idx_to_closest_value(data.loc[i,'LUM/OD'].Time,t_h)
                        tmp.append(data.loc[i,'LUM/OD'].loc[idx,'LUM/OD'])
                        label_y = 'luciferase acitvity [RLU/OD]'
                    except KeyError:
                        try:
                            idx=find_idx_to_closest_value(data.loc[i,'LUM'].Time,t_h)
                            tmp.append(data.loc[i,'LUM'].loc[idx,'LUM'])     
                            label_y= 'LUM [RLU]'
                        except KeyError:
                            print('{} dose response NOT plotted due to missing luminesence data'.format(y))
                            return
                    except AttributeError:
                        print('no luminescence data found')
                        lum=False
                elif y=='OD':
                    idx=find_idx_to_closest_value(data.loc[i,'OD'].Time,t_h)
                    control=data.loc[(data['Date']==data.loc[i,'Date']) & (data['Well']==data.loc[i,'control well'])]
                    control.reset_index(drop=True,inplace=True)
                    if len(control)>1:
                        print('Error, several instances possible as control well. Code is not equipped to intercept. Dose response of OD cannot be performed')
    #                    return
                    control = control.loc[0,'OD'].loc[idx,'OD']
                    tmp.append(100/control*data.loc[i,'OD'].loc[idx,'OD'])
                else:
                    print('No adequate response for dose response given. Try y=\'LUM\' or y=\'OD\'')
    #                return
            if y=='OD':
                dfDoseResponse.append(tmp)
            elif lum: # don't add to Dataset if no luminescent data has been measured
                dfDoseResponse.append(tmp)
    
    
    dfDoseResponse=pd.DataFrame(dfDoseResponse[1:],columns=dfDoseResponse[0])
    
    
    for t in timepoints:
        #some LUM/OD entries are 0, to plot on log and not have pixels in infity remove those values
        tmpDR = dfDoseResponse.copy()
        if y=='LUM':
            for i in range(len(tmpDR)):
                if tmpDR.loc[i,'Response after {}'.format(t)]==0:
                    tmpDR =tmpDR.drop(i, axis=0)
            
        if normalizeX:
            ICs=pd.read_csv('{}{}TimepointIC.csv'.format(save,normalizeX.replace(' ','_')).translate(translationTable))
            for category1 in set(ICs[category_to_plot[0]]):
                tmpIC=ICs.loc[ICs[category_to_plot[0]]==category1]
                for category2 in set(tmpIC[category_to_plot[1]]):
                    IC=float(tmpIC.loc[ICs[category_to_plot[1]]==category2,'IC mean'])
                    tmpDR.loc[np.logical_and((tmpDR[category_to_plot[0]]==category1), (tmpDR[category_to_plot[1]]==category2)),variations_in_category[0]]=tmpDR.loc[np.logical_and((tmpDR[category_to_plot[0]]==category1), (tmpDR[category_to_plot[1]]==category2)),variations_in_category[0]]/IC
        
        #save data as csv
        for inducer in set(data[category_to_plot[0]]):
            if SaveInduction:
                table='{} \t {} \t {} \t fold-change \n'.format(category_to_plot[0],category_to_plot[1],variations_in_category[0])
                for strain in set(data[category_to_plot[1]]):
                    indDR= tmpDR.loc[tmpDR[category_to_plot[0]]==inducer]
                    indDR= indDR.loc[indDR[category_to_plot[1]]==strain]
                    if y=='growth rate':
                        control=np.mean(indDR.loc[indDR[variations_in_category[0]]==0,'Response'])
                    else:
                        control=np.mean(indDR.loc[indDR[variations_in_category[0]]==0,'Response after {}'.format(t)])
                    for conc in set(indDR[variations_in_category[0]]):
                        if y=='growth rate':
                            mean = np.mean(indDR.loc[indDR[variations_in_category[0]]==conc,'Response'])
                        else:
                            mean = np.mean(indDR.loc[indDR[variations_in_category[0]]==conc,'Response after {}'.format(t)])
                        foldchange='{:.2f}'.format(mean/control)
                        table = '{}{}\t {} \t {} \t {}\n'.format(table,inducer,strain,conc,foldchange)
                if y=='growth rate':
                    with open('{}FoldChange{}_{}.csv'.format(save,y,remove_unit(inducer).translate(translationTable)),'w', encoding="utf-8") as writer:
                        writer.write(table)
                else:                    
                    with open('{}FoldChange{}{}_{}.csv'.format(save,y,t.replace(' ',''),remove_unit(inducer).translate(translationTable)),'w', encoding="utf-8") as writer:
                        writer.write(table)
            #Plotting
            cmap = sns.color_palette(palette='colorblind',n_colors=len(set(tmpDR.loc[tmpDR[category_to_plot[0]]==inducer,variations_in_category[1]])))
            if y=='growth rate':
                sns.lineplot(data=tmpDR.loc[tmpDR[category_to_plot[0]]==inducer],x=variations_in_category[0],y='Response',hue=variations_in_category[1],style=category_to_plot[0],markers=True,palette=cmap)
            else:
                sns.lineplot(data=tmpDR.loc[tmpDR[category_to_plot[0]]==inducer],x=variations_in_category[0],y='Response after {}'.format(t),hue=variations_in_category[1],style=category_to_plot[0],markers=True,palette=cmap)
            if y=='LUM':
                plt.yscale('log')
                if 'CLARIOstar' in devices:
                    ylim_max=3*10**7
                    if label_y== 'LUM [RLU]':
                        ylim_max=2*10**5
                else:
                    ylim_max=2*10**5
                if not 'CLARIOstar' in devices:
                    ylim_min=5*10**0
                    if label_y== 'LUM [RLU]':
                        ylim_min=5*10**2
                elif len(devices)>1:
                    ylim_min=5*10**0
                else:
                    ylim_min=10**4
                if ('PbcrC-lux' in list(tmpDR[category_to_plot[1]])) and (len(set(tmpDR[category_to_plot[1]]))==1):
                    ylim_min=10**5
                    ylim_max=10**7
                
                
            elif y=='OD':
                label_y='growth relative to unperturbed [%]'
                ylim_min=0
                ylim_max=150
            elif y=='growth rate':
                label_y='growth rate relative to unperturbed [%]'
                ylim_min=0
                ylim_max=150
            plt.ylim(ylim_min,ylim_max)

            if normalizeX:
                plt.xlabel('{} concentration, normalized by the MIC'.format(remove_unit(inducer)))
            else:
                plt.xlabel(inducer)
            # #manipulate x axis range
            # xpl = np.array([5])
            # ypl = np.array([10**6])
            # # plot the data
            # plt.plot(xpl,ypl)
            
            plt.ylabel(label_y)
            if y=='growth rate':
                plt.title('Growth rate response')
            else:
                plt.title('Response after {}'.format(t))
            plt.grid(b=True, which = 'major')
            plt.grid(b=True, which = 'minor',linewidth=0.5)
            # Shrink current axis by 20%
            ax=plt.gca()
            handles, labels = ax.get_legend_handles_labels()
            if not (len(ax.get_legend_handles_labels()[1])<=len(set(tmpDR[category_to_plot[1]]))):
                handles=handles[:len(set(tmpDR.loc[tmpDR[category_to_plot[0]]==inducer,variations_in_category[1]]))+1]
                labels=labels[:len(set(tmpDR.loc[tmpDR[category_to_plot[0]]==inducer,variations_in_category[1]]))+1]
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
            # Put a legend to the right of the current axis
            ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))
            if y=='growth rate':
                if normalizeX:
                    plt.savefig('{}DoseResponse_{}_{}_normed.svg'.format(save,y,remove_unit(inducer)).translate(translationTable), bbox_inches='tight',dpi=300)
                    plt.savefig('{}DoseResponse_{}_{}_normed.png'.format(save,y,remove_unit(inducer)).translate(translationTable), bbox_inches='tight',dpi=300)
                else:
                    plt.savefig('{}DoseResponse_{}_{}.svg'.format(save,y,remove_unit(inducer)).translate(translationTable), bbox_inches='tight',dpi=300)
                    plt.savefig('{}DoseResponse_{}_{}.png'.format(save,y,remove_unit(inducer)).translate(translationTable), bbox_inches='tight',dpi=300)
            else:
                if normalizeX:
                    plt.savefig('{}DoseResponse{}{}_{}_normed.svg'.format(save,y,t.replace(' ',''),remove_unit(inducer)).translate(translationTable), bbox_inches='tight',dpi=300)
                    plt.savefig('{}DoseResponse{}{}_{}_normed.png'.format(save,y,t.replace(' ',''),remove_unit(inducer)).translate(translationTable), bbox_inches='tight',dpi=300)
                else:
                    plt.savefig('{}DoseResponse{}{}_{}.svg'.format(save,y,t.replace(' ',''),remove_unit(inducer)).translate(translationTable), bbox_inches='tight',dpi=300)
                    plt.savefig('{}DoseResponse{}{}_{}.png'.format(save,y,t.replace(' ',''),remove_unit(inducer)).translate(translationTable), bbox_inches='tight',dpi=300)
            plt.close()#show()
    print('{} dose response plotted at timepoints {}'.format(y,str(timepoints)))
    return
    
def time_of_growth(data,category_to_plot,variations_in_category,save,threshold_OD=10**(-1)):
    df=[[category_to_plot[0],variations_in_category[0],category_to_plot[1],variations_in_category[1],threshold_OD,'time till reaching OD threshold']]
    for ind in range(len(data)):
        tmp=[data.loc[ind,category_to_plot[0]],data.loc[ind,variations_in_category[0]],data.loc[ind,category_to_plot[1]],data.loc[ind,variations_in_category[1]],threshold_OD]
        tmp_df=data.loc[ind,'OD']
        t=float('NaN')
        for time in range(2,len(tmp_df)):
            if tmp_df.loc[time,'OD']>=threshold_OD:
                t=tmp_df.loc[time,'Time']
                break
        tmp.append(t)
        df.append(tmp)
    df=pd.DataFrame(df[1:],columns=df[0])
    df.to_csv('{}TimeOfGrowth_{}.csv'.format(save,threshold_OD),index=False,sep=';')
    return