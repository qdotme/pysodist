import pandas
import glob
import pylab
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from pyteomics.mass.mass import calculate_mass

def calc_mz(seq, charge):
    mod_mass=0.0
    if '[' in seq:
        split_open = seq.split('[')
        seq = split_open[0]
        for i in split_open[1:]:
            seq = seq + ''.join(i.split(']')[1])
            mod_mass = mod_mass + float(i.split(']')[0])
    return calculate_mass(seq, charge=charge) + mod_mass/(charge*1.0)

def parse_encyclopedia(path, filename=None, q_value=0.05, IO=True):
    """Helper to generate a list of peptides an analyze from an encyclopedia output.
    Filters the input to the q-value specified by the user (default of 0.05)
    Typical use:
       peptides = gen_peptide_list('./')
       **assumes that the encyclopedia file you want to work with is the only file in the directory you are using
    Optional arguments:
        **  q_value=[float] | specify the q-value cutoff to use to filter peptides 
            (0.0 will take all peptides in the input)
        **  filename | directly specify the full path and filename to analyze instead 
            of assuming you will use the first one (alphabetically) in the directory
        ** IO | True will show outputs, false will supress
    Return:
        Returns a pandas dataframe indexed by the peptide sequence and bearing columns of 'mz' and 'z' (mass/charge, charge)
    """
    #check for optional arguments then read proper file
    if filename is None:
        filename = glob.glob(path+'/*.encyclopedia.txt')[0]
    if IO:
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('reading peptides from encyclopedia file ' + filename)
    
    #generate pandas file, filter, parse peptides, save encyc_rts values, save charge values, calculate mz, save proteins
    encyc = pandas.read_csv(filename, sep='\t')
    all_entries = encyc[encyc['q-value']<q_value]
    peptides = [i[2:-2] for i in all_entries['peptide']]
    encyc_rts = [float(i.split(':')[1])/60.0 for i in all_entries['PSMId']]
    zs = [i.split('+')[-1] for i in all_entries['PSMId']]
    prot_names = [i for i in all_entries['proteinIds']]
    mzs = []
    for index, pepseq in enumerate(peptides):
        mzs.append(calc_mz(pepseq, int(zs[index])))
    #mzs = [calc_mz(pair[0], int(zs[pair[1]])) for pair in enumerate(peptides))]

    assert len(peptides) == len(encyc_rts)
    assert len(peptides) == len(zs)
    assert len(peptides) == len(prot_names)
    assert len(peptides) == len(mzs)
    
    #convert to a dictionary
    pdl = {}
    for i, p in enumerate(peptides):
        pdl[p] = [encyc_rts[i], zs[i], prot_names[i], mzs[i]]
    
    # show plots for q-value criteria
    if IO:
        print(str(encyc.shape[0]) + ' peptides identified')
        print(str(len(peptides)) + ' peptides meet q-value criteria')
        print(str(encyc.shape[0]-len(peptides)) + ' peptides fail q-value criteria')

        fig, ax = pylab.subplots() # create a new figure with a default 111 subplot
        ax.hist(encyc['q-value'].dropna(), bins=100)
        axins = zoomed_inset_axes(ax, 2, loc=1) # zoom-factor: 2.5, location: upper-left
        axins.hist(encyc['q-value'].dropna(), bins=100)
        axins.set_xlim(0, q_value*2) # apply the x-limits
        axins.set_ylim(0, ax.get_ylim()[1]/4) # apply the x-limits
        axins.vlines([q_value], 0, ax.get_ylim()[1]/4, colors='red', linestyle='dashed', label='q-cuttoff')
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
        ax.set_xlabel('q-value')
        ax.set_ylabel('peptides')
        pylab.show()
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    
    # convert to a pandas dataframe
    mz_df = pandas.DataFrame.from_dict(pdl, orient='index')
    mz_df = mz_df.rename(index=str, columns={0:'encyc_rt', 1:'z', 2:'proteinIds', 3:'mz'})
    return mz_df
    
def parse_rt_fit(path, mz_dataframe, filename = None, IO=True):
    """Helper to get peptides from the rt_ft file from encyclopedia output.
    Typical use:
       full_df = parse_rt_fit('./', peptides)
       **assumes that the encyclopedia rt_fit file you want to work with is the only file in the directory you are using
    Optional arguments:
        **  filename | directly specify the full path and filename to analyze instead 
            of assuming you will use the first one (alphabetically) in the directory
        ** IO | True will show outputs, false will supress
    """
    #read rt_fit file
    if filename is None:
        filename = glob.glob('./*.rt_fit.txt')[0]
    
    rt_fit = pandas.read_csv(filename, sep='\t')
    if IO:
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
        print('reading rt_fit file '  + filename)
        print(str(mz_dataframe.shape[0]) + ' peptides in filtered .encyclopedia file')
        print(str(rt_fit.shape[0]) + ' peptides in rt_fit file')
        print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    
    rt_fit = rt_fit.set_index('sequence')
    rt_fit['z'] = mz_dataframe['z']
    rt_fit['mz'] = mz_dataframe['mz']
    rt_fit['proteinIds'] = mz_dataframe['proteinIds']
    rt_fit['encyc_rt'] = mz_dataframe['encyc_rt']
    rt_fit = rt_fit.rename(index=str, columns={'actual':'observed_rt'})

    return rt_fit
    
def output_gMassacre(full_df, outfile='./gMassacre', sort_field='observed_rt'):
    to_output = full_df[['observed_rt', 'encyc_rt', 'mz', 'z', 'proteinIds']]
    to_output = to_output.sort_values(by=[sort_field])
    to_output.to_csv(outfile)

