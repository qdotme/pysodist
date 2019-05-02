import os
import pandas as pd
import glob
import pylab
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from pyteomics.mass.mass import calculate_mass


from davis_lab_quant_report_defs import*

def extract_skyline_sub(full_dataframe, sample, proteins=None, isotope='light'):
    '''
    extract_skyline_sub takes a full skyline report as input and outputs a filtered dataframe with only those columns 
    relevant to the proteins and sample of interest

    :param full_dataframe: pandas dataframe directly uploaded from skyline
    :param sample: string of injection sample name
    :param proteins: string of uniprot identifier for particular proteins of interest (default None includes all proteins)
    :param isotope: string for filtering output dataframe such that columns from a particular isotope type are included (default 'light')

    :return: a pandas dataframe with the only the sample and proteins of interest and only those columns needed for subsequent analyses
    '''
    required_cols = [isotope+' '+sample+' '+RT_START_FIELD,isotope+' '+sample+' '+RT_END_FIELD,PROTEIN_FULL_NAME_FIELD,PEPTIDE_MOD_SEQ_FIELD,isotope+' '+ISOTOPE_FIND_FIELD,PEPTIDE_CHARGE_FIELD,isotope+' '+sample+' '+DETECT_Q_VALUE_FIELD]#fill this in with the columns you need - use the isotope info to get the right ones
    new_data_frame = full_dataframe.loc[:,required_cols]
    if proteins is None:
        return new_data_frame #get all the proteins
    return new_data_frame[new_data_frame[PROTEIN_FULL_NAME_FIELD].isin(proteins)]

def parse_sub_skyline(skyline_sub_df, sample, q_value=0.05, isotope='light', IO=True, output_file=None):
    '''
    parse_sub_skyline is a helper function to parse a sample and protein filtered skyline dataframe such that a properly 
    sorted and labeled dataframe is output
    
    :param skyline_sub_df: pandas dataframe with only one sample injection and only those proteins of interest
    :param sample: string of injection sample name
    :param q_value: float with a optional q_value cutoff to filter your data (default 0.05)
    :param isotope: string for filtering output dataframe such that columns from a particular isotope type are included (default 'light')
    :param IO: bool determining if you want to display output (number of peptides, plots of filtering)
    :param output_file: string determining if you want to write the dataframe to disk (default None won't write to disk)
    
    :return: a pandas dataframe with no index, columns of peptide_modified_sequence, rt_start (seconds), rt_end (seconds), charge, 
    mz [note that this may only be approxiamate for modified peptides], protein_IDs
    '''
    if IO:
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('reading peptides from skyline dataframe ' + sample)
    
    #generate pandas file, filter, parse peptides, save encyc_rts values, save charge values, calculate mz, save proteins
    filt_entries = skyline_sub_df.loc[skyline_sub_df[isotope+' '+sample+' '+DETECT_Q_VALUE_FIELD] < q_value]
    peptides = [i for i in filt_entries[PEPTIDE_MOD_SEQ_FIELD]]
    rts_start = [round(t,2) for t in filt_entries[isotope+' '+sample+' '+RT_START_FIELD]]
    #print(isotope+' '+sample+' '+RT_START_FIELD)
    rts_end = [round(t,2) for t in filt_entries[isotope+' '+sample+' '+RT_END_FIELD]]
    zs = [i for i in filt_entries[PEPTIDE_CHARGE_FIELD]]
    prot_names = [i for i in filt_entries[PROTEIN_FULL_NAME_FIELD]]
    mzs = [mz for mz in filt_entries[isotope+' '+ISOTOPE_FIND_FIELD]]

    #calc number of used peptides, and make sure all the arrays are the same size
    used_peptides = len(peptides)
    assert used_peptides == len(rts_start)
    assert used_peptides == len(rts_end)
    assert used_peptides == len(zs)
    assert used_peptides == len(prot_names)
    assert used_peptides == len(mzs)
    
    #make a dataframe with your newly calculated entries
    pandas_dataframe = pd.DataFrame({peptides[i]:
    {'peptide_modified_sequence':peptides[i], 'rt_start':rts_start[i], 'rt_end':rts_end[i], 
    'charge':zs[i], 'mz':mzs[i], 'protein_IDs':prot_names[i]} for i in range(used_peptides)}).T
    pandas_dataframe = pandas_dataframe.sort_values(by=['rt_start'])
    # show plots for q-value criteria
    if IO:
        skyline_sub_df.rename(columns={isotope+' '+sample+' '+DETECT_Q_VALUE_FIELD : 'q-value'},inplace=True)
        
        print(str(skyline_sub_df.shape[0]) + ' peptides identified')
        print(str(len(peptides)) + ' peptides meet q-value criteria')
        print(str(skyline_sub_df.shape[0]-len(peptides)) + ' peptides fail q-value criteria')

        fig, ax = pylab.subplots() # create a new figure with a default 111 subplot
        ax.hist(skyline_sub_df['q-value'].dropna(), bins=100)
        ax.vlines([q_value], 0, ax.get_ylim()[1], colors='red', linestyle='dashed', label='q-cuttoff')
        #axins = zoomed_inset_axes(ax, 2, loc=1) # zoom-factor: 2.5, location: upper-left
        #axins.hist(skyline_sub_df['q-value'].dropna(), bins=100)
        #axins.set_xlim(0, q_value*2) # apply the x-limits
        #axins.set_ylim(0, ax.get_ylim()[1]/4) # apply the x-limits
        #axins.vlines([q_value], 0, ax.get_ylim()[1]/4, colors='red', linestyle='dashed', label='q-cuttoff')
        #mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
        ax.set_xlabel('q-value')
        ax.set_ylabel('peptides')
        pylab.show()
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    
    #return the sorted dataframe
    return pandas_dataframe

def parse_skyline(path_to_skyline_csv, output_file_prefix = None, samples=None, proteins=None, isotope='light', q_value=0.05, IO=True):
    '''
    Wrapper function to execute extract_skyline_sub and parse_sub_skyline on all proteins of interest

    :param path_to_skyline_csv: string of directory path to location of skyline report
    :param output_file_prefix: string user defined prefix identifier for output dataframe as tab separated table (default: None results in no csv outputs)
    :param samples:  list of strings of individual sample identifiers in skyline report (default None results in analysis of all samples)
    :param proteins: string of uniprot identifier for particular proteins of interest (default None includes all proteins)
    :param isotope: string for filtering output dataframe such that columns from a particular isotope type are included (default 'light')
    :param q_value: float with a optional q_value cutoff to filter your data (default 0.05)
    :param IO: bool determining if you want to display output (number of peptides, plots of filtering)
    
    :return: a list of pandas dataframes with no indices, columns of peptide_modified_sequence, rt_start (seconds), rt_end (seconds), charge, 
    mz [note that this may only be approxiamate for modified peptides], protein_IDs
    '''
    skyline_complete = pd.read_csv(path_to_skyline_csv+'.csv')
    if samples is None:
        samp = list(set([i.split(' ')[1] for i in skyline_complete.columns if 'Min Start Time' in i])) #NOTE that this requires no spaces in the column names, perhaps revise
    else:
        samp = samples
    output_list = []
    for sam in samp:
        extracted_sub = extract_skyline_sub(skyline_complete, sam, proteins=proteins, isotope=isotope)
        output_list.append(parse_sub_skyline(extracted_sub, sam, q_value=q_value, isotope=isotope, IO=IO, output_file=output_file_prefix+sam))
        if not(output_file_prefix is None):
            output_list[-1].to_csv(output_file_prefix+'_'+sam+'.tsv', sep='\t', columns = ['rt_start', 'rt_end', 'peptide_modified_sequence', 'charge', 'mz', 'protein_IDs'], index=False)
    return output_list

if __name__ == "__main__":
    skyline_df = parse_skyline('C:\\Users\\Andrew\\Dropbox (Personal)\\Davis Lab\\Skyline\\20190212\\bt_timecourse_pysodist_test_report', isotope='light', samples = ['BT_control_7hdia'],output_file_prefix='test',IO = True)