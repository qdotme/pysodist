import dists_and_params
import parsers
import sys
import scipy.optimize as optimize
import time
import pylab
#############################
#helper functions for fitting
#############################
def chisq(array1,array2):
    n=len(array1)
    sq=sum((array1-array2)**2)
    return sq/n

def residual(a1,a2):
    return sum((a1-a2)**2)

def print_param_labels():
    print(list(map(lambda par:par.symbol,params)))
    
def extract_param_vector():
    return list(map(lambda par:par.value,params))

def apply_param_vector(new_vals):
    for i in range(len(params)):
        params[i].change_value(new_vals[i])
     
def eval_param_vector(param_vec,spectrum,peak_data):
    apply_param_vector(param_vec)
    calc_spec=spectrum.get_mz()
    #scale up to roughly match size (experimental)
    #calc_spec*=max(peak_data)/max(calc_spec)
    return calc_spec-peak_data

if __name__ == "__main__":
    #Write out Liscence Information
    print("""\npysodist Copyright (C) 2018
    This program comes with ABSOLUTELY NO WARRANTY
    This is free software, and you are welcome to redistribute it under 
    the terms of the GNU General Public Licence\n""")

    print("When running pysodist, please provide the following unless you want to run with default values:")
    print("python pysodist.py [input1] [input2] [input3] [input4_OPTIONAL]")
    print("input1 || xml_input_file.xml(default)")
    print("input2 || results_file.tsv(default)")
    print("input3 || display fits and spectra [True/False(default)]")
    print("input4 || path_to_fit_directory [optional]\n")

    #parse the initial files to get oriented
    print("++++++++++++++++++++++++++++++++")
    try:
        input_file = sys.argv[1]
    except IndexError:
        print("Info file missing. Please provide a .in file simliar to the example xml_inputs.xml")
        print("Using default file of './xml_inputs.xml'")
        input_file='./xml_inputs.xml'
    print("Input file : " + input_file)
    print("++++++++++++++++++++++++++++++++")

    try:
        results_tsv = sys.argv[2]
    except IndexError:
        print("No results file specified. Please provide a results file to save your results")
        print("Using default results file of ./results.tsv")
        results_tsv = './results.tsv'
    print("Results file : " + results_tsv)
    print("++++++++++++++++++++++++++++++++")

    try:
        DISPLAY_RESULTS = sys.argv[3]
    except IndexError:
        print("No fits or spectra will be displayed. To see spectra and fits, pass 'True' as argument 3.")
        DISPLAY_RESULTS = False
    
    try:
        fit_directory = sys.argv[4]
        SAVE_FITS=True
        print("Fit output directory : " + fit_directory)
    except IndexError:
        print("No output directory specified. Please provide a directory to store fits")
        SAVE_FITS=False
    print("++++++++++++++++++++++++++++++++")

    MainInfo=parsers.MainInfoParser(input_file)
    AtomLib=parsers.AtomInfoParser(MainInfo.atomfile)
    ResInfo=parsers.ResInfoParser(MainInfo.resfile,AtomLib)
    PeakInfo=parsers.PeakInfoParser(MainInfo.peakfile,ResInfo)

    #parameters
    b=MainInfo.b_param
    m_off=MainInfo.m_off_param
    gw=MainInfo.gw_param
    amps=ResInfo.amp_params
    thetas=ResInfo.theta_params
    phis=ResInfo.phi_params

    params=[b,m_off,gw] #make a list of all params
    params.extend(amps)
    params.extend(thetas)
    params.extend(phis)

    #prepare variable atoms and residues
    var_atom_syms=map(lambda par:par.label,thetas)
    var_atoms=list(map(lambda sym:AtomLib.atom_dict[sym],var_atom_syms))
    for i,theta in enumerate(thetas):
        theta.link_to(var_atoms[i])

    var_res_syms=map(lambda par:par.label,phis)
    var_residues=list(map(lambda sym:ResInfo.res_dict[sym][1],var_res_syms))
    for i,phi in enumerate(phis):
        phi.link_to(var_residues[i])

    #main loop
    results_dict={}
    for i in range(PeakInfo.num_peaks):
        m_hd=PeakInfo.get_m_hd(i)
        if b.tag=='auto':
            b.value=PeakInfo.auto_baseline(i)
        peak_data=PeakInfo.get_peak_data(i)
        peak_data/=max(peak_data)/3
        
        #prepare shift
        shift=dists_and_params.Distribution()
        shift.build_as_shift(m_hd,m_off.value)
        m_off.link_to(shift)
        
        #prepare gaussian
        gaussian=dists_and_params.Distribution()
        gaussian.build_as_gaussian(gw.value)
        gw.link_to(gaussian)
        
        #prepare stick spectrum
        stick_spectrum=dists_and_params.Distribution()
        stick_spec_amps=list(map(lambda param:param.value,amps))
        stick_spectrum.build_as_dist_sum_with_amps(PeakInfo.peak_list[i],stick_spec_amps)
        for amp in amps:
            amp.link_to(stick_spectrum)
        
        #prepare overall spectrum
        spectrum=dists_and_params.Distribution()
        spectrum.build_as_dist_sum_with_mults([stick_spectrum,gaussian,shift],[1,1,1])
        
        #do fitting

        t1 = time.time()
        optimize.leastsq(lambda param_vec: eval_param_vector(param_vec,spectrum,peak_data),
                        extract_param_vector(),maxfev=500)
        param_vector = extract_param_vector()
        t2 = time.time()
        resid_value = str(round(residual(spectrum.get_mz(), peak_data),2))
        peptide_UID = PeakInfo.get_peptide_UID(i)
        peptide_seq = PeakInfo.get_peptide_seq(i)
        charge = PeakInfo.get_peptide_charge(i)
        protein = PeakInfo.get_protein_name(i)
        retention_time = PeakInfo.get_retention_time(i)

        fit_time = str(round(t2-t1,2))
        print('++Fit ' + peptide_UID+'['+retention_time+']'+ ' in ' + fit_time + ' seconds | res = ' + resid_value + '++')
        print_param_labels()
        print(param_vector)
        
        if DISPLAY_RESULTS:
            PeakInfo.display_spectra(i)
            spectrum.display_fit()

        if SAVE_FITS:
            fit_file = fit_directory+peptide_UID+'_'+str(retention_time)+'.res'

            fit_mz = [str(i) for i in spectrum.get_mz_axis()]
            fit_int = [str(i) for i in dists_and_params.ifft(spectrum.get_ft())]
            spectra_int = [str(i) for i in PeakInfo.get_peak_data(i)]

            peaks_file = open(fit_file,'w')
            
            peaks_file.write('mz\texp_int\tfit_int\n')
            for i in range(len(fit_mz)):
                peaks_file.write(fit_mz[i]+'\t'+spectra_int[i]+'\t'+fit_int[i]+'\n')
            peaks_file.close()
        else:
            fit_file='not_saved'

        results_dict[peptide_UID+retention_time] = {'peptide_seq':peptide_seq, 'charge':charge, 'retention_time':retention_time, 'protein':protein, 
                                                'AMP_U':round(param_vector[3],6), 'AMP_L':round(param_vector[4],6), 'FRAC_NX':round(param_vector[5], 6),
                                                'B':param_vector[0], 'M_OFF':round(param_vector[1],6), 'GW':round(param_vector[2],6), 
                                                'residual':resid_value, 'fit_file':fit_file}
    #output results
    r_file = open(results_tsv,'w')
    values_to_print=list(results_dict[peptide_UID+retention_time].keys())

    header = '\t'.join(values_to_print)
    r_file.write('peptide_UID\t'+header+'\n')

    for pep in results_dict.keys():
        out_line = [pep]
        for val in values_to_print:
            out_line.append(str(results_dict[pep][val]))
        r_file.write('\t'.join(out_line))
        r_file.write('\n')
    r_file.close()