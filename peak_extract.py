import os
import pandas as pd
import glob
import scipy
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pyteomics as pyt
from pyteomics import mzml
import math

os.chdir('C:\\Users\\Andrew\\Dropbox (Personal)\\Davis Lab\\Skyline\\20190212')

from davis_lab_quant_report_defs import*

def mzml_retrieve(mzml_iterator, scans = 500, old_slice = None):
    '''
    retrieves the scans from a portion of an mzml file and generates arrays of mz versus intensity at each RT

    :param mzml_iter: string of directory path to location of mzml file
    :param scans: number of scans to retrieve from mzml file
    :param old_slice: a portion of the previous mzml slice that is wider than the widest expected peak 

    :return: lists of mz versus intensity dictionaries and retention times
    '''
    counter = 0
    list_of_mz_and_intensity_dict = []
    rt_list = []

    if not(old_slice is None):
        rt_list = old_slice[0]
        list_of_mz_and_intensity_dict = old_slice[1]
    
    for scan in mzml_iterator:
        mz_and_intensity_dict = {'mz_data':scan['m/z array'], 'intensity_data':scan['intensity array']}
        list_of_mz_and_intensity_dict.append(mz_and_intensity_dict)
        new_rt = float(scan['scanList']['scan'][0]['scan start time'])
        rt_list.append(new_rt)
        counter += 1
        if counter >= scans:
            break

    return [rt_list, list_of_mz_and_intensity_dict]

def mz_window(peptide_modified_sequence, peptide_charge, peptide_mz, extra_mz, labeling = 'None'):
    '''
    calculates the mz window required to accomodate the range of mzs possible for all isotopes present in a given peptidde

    :param peptide_modified_sequence: string of peptide amino acid sequence including modifications
    :param peptide_charge: int of peptide charge state
    :param peptide_mz: float of mz for the unlabeled monoisotopic peptide
    :param extra_mz: additional mass to include in the window
    :param labeling: isotope labeleing method ie. "K8R10" or "N15" (default None: unlabeled)

    :return: list of lower and upper range of possible mz values +/- 8 mass units for each peptide
    '''
    if labeling == 'N15':
        labeled_mz = ((sum([AANITROGENS[i]*(N15MASS-N14MASS) for i in peptide_modified_sequence if i in AANITROGENS.keys()]))/peptide_charge)+peptide_mz
    elif labeling == 'k8R10':
        labeled_mz = ((sum([AANITROGENS[i]*(N15MASS-N14MASS)+AACARBONS[i]*(C13MASS-C12MASS) for i in peptide_modified_sequence if i in ['K','R']]))/peptide_charge)+peptide_mz
    
    elif labeling == 'C13':
        labeled_mz = ((sum([AACARBONS[i]*(C13MASS-C12MASS) for i in peptide_modified_sequence if i in AACARBONS.keys()]))/peptide_charge)+peptide_mz

    else:
        labeled_mz = peptide_mz

    return [peptide_mz-(extra_mz/peptide_charge),labeled_mz+(extra_mz/peptide_charge)]

def peaks_find(scan, mz_range):
    '''
    from a single scan record all intensities at each mz that are within the mz range of interest

    :param scan: array of mz versus intensity at a particular retention time
    :param mz_range: lower and upper ranges of possible mz values to consider in the peak

    :return: an np array of mz values and intensities that fall within the range of interest
    '''
    mz = []
    intensity = []

    for i in range(len(scan['mz_data'])):
        if (scan['mz_data'][i] >= mz_range[0]) and (scan['mz_data'][i] <= mz_range[1]):
            mz.append(scan['mz_data'][i])
            intensity.append(scan['intensity_data'][i])
    mz_int = [mz,intensity]
    
    return(np.array(mz_int))

def interpolate_mz_array(mz_int_array, mz_interp):
    '''
    interpolates each mz_intensity_array to enable combination at discrete mz values for each peak
    across the retention time window

    :param mz_int_array: 2D np array of observed mass-to-charge versus intensities extracted from a given scan
    :param mz_interp: np array of evenly spaced mz values upon which to generate interpolated intensities

    :return: numpy array of interpolated intensities for mz_interp
    '''
    f = interp1d(mz_int_array[0],mz_int_array[1])
    interp_intensity = f(mz_interp)
    
    return np.array(interp_intensity)

def peaks_extract(path_to_peptide_rt_tsv, path_to_mzml, labeling = None,output_prefix = None, batch_size=500, min_batch_size=125, interp_space_mz = 0.005):
    '''
    Wrapper function to extract peaks from each peptide within the peptide rt csv

    :param path_to_peptide_rt_tsv: string of directory path to peptide retention time tsv
    :param path_to_mzml: string of directory path to mzml file
    :param labeling: isotope labeleing method ie. "SILAC" or "N15" (default None: unlabeled)
    :param output_prefix: string user defined prefix identifier for output text files (default: None results in no text file)

    :return text files of mz versus interpolated and summed intensities across the mz and rt ranges specific to each peptide
    '''
    peptide_rt_df = pd.read_csv(path_to_peptide_rt_tsv+'.tsv',sep='\t')
    mzml_iter = pyt.mzml.MzML(path_to_mzml+'.mzML')

    mzml_slice = mzml_retrieve(mzml_iter, scans = batch_size, old_slice = None) #get first/next 1000 scans in the rt indexed array format
    while not len(mzml_slice[0]) == min_batch_size:
        for i in range(len(peptide_rt_df)):
            if (peptide_rt_df['rt_end'][i] <=  max(mzml_slice[0])) and (peptide_rt_df['rt_start'][i] >= min(mzml_slice[0])): #if the beginning and end rts of the ith peptide in the rt tsv falls within the slice of the mzml that we've taken:
                mz_range_wide = mz_window(peptide_rt_df['peptide_modified_sequence'][i], peptide_rt_df['charge'][i], peptide_rt_df['mz'][i], 16, labeling = labeling) # calculate mz_range to include as inputs for interpolation function
                mz_range_narrow = mz_window(peptide_rt_df['peptide_modified_sequence'][i], peptide_rt_df['charge'][i], peptide_rt_df['mz'][i], 8, labeling = labeling) # calculate the mz_range to include as outputs for interpolation function
                #further slice the rt indexed mz/intensity array such that it only includes the rt window relevant to this peptide
                pep_slice_indices = [j for j in range(len(mzml_slice[0])) if (peptide_rt_df['rt_start'][i] <= mzml_slice[0][j]) and (peptide_rt_df['rt_end'][i] >= mzml_slice[0][j])] 
                mzml_peptide_slice = [mzml_slice[0][min(pep_slice_indices):max(pep_slice_indices)],mzml_slice[1][min(pep_slice_indices):max(pep_slice_indices)]]
                
                np.set_printoptions(suppress = True)

                while True:
                    try:
                        peaks = [peaks_find(scn,mz_range_wide) for scn in mzml_peptide_slice[1]] #generate peaks list with mz values extending out to the wide range
                        interp_mz = np.arange(mz_range_narrow[0],mz_range_narrow[1], interp_space_mz) #generate mz array upon which to interpolate intensities, ensure that these values are within the mz range of interest
                        interp_intens = [interpolate_mz_array(pks, interp_mz) for pks in peaks] #generate interpolated peaks list
                    
                    except ValueError:
                        mz_range_wide[0] -= 1
                        mz_range_wide[1] += 1
                        #print(min(peaks[0][0]))
                        #print(max(peaks[0][0]))
                        #print(mz_range_wide)
                        continue
                    break


                sum_intens = sum(interp_intens)
                mz_intens_sums = np.array([interp_mz,sum_intens])
                mz_intens_sums_out = np.transpose(mz_intens_sums)

                if not(output_prefix is None):
                    np.savetxt(output_prefix+peptide_rt_df['peptide_modified_sequence'][i]+'+'+str(peptide_rt_df['charge'][i])+'.tsv', mz_intens_sums_out,fmt = '%0.5f', delimiter='\t') #Create a 'PEPTIDEK+2' unique ID to accompany the text file in the output
    #currently there's no return statement, but I could store all of the summed peaks in their own list
        mzml_slice = mzml_retrieve(mzml_iter, scans=batch_size, old_slice=[sublist[-min_batch_size:] for sublist in mzml_slice]) #get first/next 500 scans in the rt indexed array format
       
                
if __name__ == "__main__":
    peaks_extract('test_BT_control_7hdia','BT_control_7hdia',labeling= 'N15', output_prefix= 'test_')