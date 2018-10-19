import dists_and_params
import xmltodict
import numpy
import matplotlib.pyplot as plt

class MainInfoParser(): #parses the .in file to hold relevant variables
    def __init__(self,infile):
        self.file=infile
        f=open(infile,'r')
        
        lines = f.readlines()
        message = ''.join(lines)
        d = xmltodict.parse(message)['pysodist_settings']
        program_settings = d['program_settings']
        input_files = d['input_files']
        fit_settings = d['fit_settings']
                
        self.opt=program_settings['mode']
        self.peakfile=input_files['batch']
        self.atomfile=input_files['atom']
        self.resfile=input_files['res']
        self.fit_iter=int(fit_settings['niter'])
        self.sig_global=float(fit_settings['sigma'])
        
        b_init=fit_settings['baseline_offset']
        b_tag=fit_settings['baseline_mode']
        self.b_param=dists_and_params.Parameter()
        self.b_param.build_as_baseline_offset(b_init,b_tag)
                         
        m_off_init=float(fit_settings['mass_offset'])
        m_off_tag=fit_settings['mass_offset_mode']
        self.m_off_param=dists_and_params.Parameter()
        self.m_off_param.build_as_machine_offset(m_off_init,m_off_tag)
                        
        gw_init=float(fit_settings['gaussian_width'])
        gw_tag=fit_settings['gaussian_width_mode']
        self.gw_param=dists_and_params.Parameter()
        self.gw_param.build_as_gaussian_width(gw_init,gw_tag)
        
        f.close()

class AtomInfoParser():
    def __init__(self,atomfile):
        self.file=atomfile
        print('AtomFile:',atomfile)
        
        self.atom_dict=dict() #dict to map atomic symbols to their distributions
                   
        f=open(atomfile,'r')
        
        reading=True
        while reading:
            newline=f.readline()
            if newline=='':
                reading=False
                break
            
            atom_info=newline.split()
            
            sym=atom_info[1] #symbol of new atom
            niso=int(atom_info[2]) #number of isotopes of new atom
            
            masses=[] #list of isotope masses for new atom
            freqs=[] #list of isotope freqeuncies for new atom
            for i in range(0,niso):
                iso_info=f.readline().split()
                masses.append(float(iso_info[0]))
                freqs.append(float(iso_info[1]))
            
            atom=dists_and_params.Distribution()
            atom.build_as_vals_freqs(masses,freqs)
            self.atom_dict[sym]=atom
        
        f.close()

class ResInfoParser():
    def __init__(self,resfile,AtomLib):
        self.file=resfile
        print('ResFile:',resfile)
        
        f=open(resfile,'r')
        f.readline() #skips over comment line
        
        #learn about the species 
        self.num_species=int(f.readline().split()[0])
        self.amp_params=[] #for species
        
        for i in range(0,self.num_species):
            spec_info=f.readline().split()
            label=spec_info[0]
            amp_init=float(spec_info[1])
            new_amp_param=dists_and_params.Parameter()
            new_amp_param.build_as_amplitude(label,amp_init,i)
            self.amp_params.append(new_amp_param)
            
        #learn about the atoms used
        
        num_atoms=int(f.readline().split()[0]) #number of active atoms (the ones used in residues)
        active_atoms=[]
        self.theta_params=[] #for variable atoms
        
        for i in range(0,num_atoms):
            atom_info=f.readline().split()
            sym=atom_info[0]
            
            if len(atom_info)<2:
                atom_info.append('default')
                
            atom_type=atom_info[1] #default, fixed, or variable
            
            current_atom=AtomLib.atom_dict[sym]
            active_atoms.append(current_atom)
            
            if atom_type=='default':
                current_atom.assert_fixed()
                #pass
            elif atom_type=='fixed':
                current_atom.change_param(float(atom_info[2]))
            elif atom_type=='variable':
                new_theta_param=dists_and_params.Parameter()
                new_theta_param.build_as_atom_frac(atom_info[0],float(atom_info[2]))
                current_atom.change_param(float(atom_info[2]))
                self.theta_params.append(new_theta_param)                
        
        #learn about the residues
        self.res_dict=dict() #dict to map a residue symbol to its species distributions
        self.phi_params=[] #for variable residues
        
        reading=True
        while reading:
            newline=f.readline()
            if newline=='':
                reading=False
                break
            
            res_header=newline.split()
            sym=res_header[0]
            res_type=res_header[1] #default, fixed, or variable            
            
            res_dists_init=[]
            for i in range(0,self.num_species):
                atom_comp=list(map(int,f.readline().split()[:num_atoms]))
                new_dist=dists_and_params.Distribution()
                new_dist.build_as_dist_sum_with_mults(active_atoms,atom_comp)
                res_dists_init.append(new_dist)
            
            if res_type=='default':
                self.res_dict[sym]=res_dists_init
            if res_type=='fixed':
                res_dists_final=[res_dists_init[0]]
                second_dist=dists_and_params.Distribution()
                second_dist.build_as_frac_dist_sum([res_dists_init[0],res_dists_init[1]],[1-float(res_header[2]),float(res_header[2])])
                res_dists_final.append(second_dist)
                self.res_dict[sym]=res_dists_final
            if res_type=='variable':
                res_dists_final=[res_dists_init[0]]
                second_dist=dists_and_params.Distribution()
                second_dist.build_as_frac_dist_sum([res_dists_init[0],res_dists_init[1]],[1-float(res_header[2]),float(res_header[2])])
                res_dists_final.append(second_dist)
                self.res_dict[sym]=res_dists_final
                new_phi_param=dists_and_params.Parameter()
                new_phi_param.build_as_residue_frac(sym,float(res_header[2]))
                self.phi_params.append(new_phi_param)
                
        f.close()

class PeakInfoParser(): #parses the .batch file
    def __init__(self,batchfile,ResLib):
        print('PeakFile:',batchfile)
        
        self.num_peaks=0
        self.peak_list=[] #holds list of species dists for each peak
        self.datfiles=[] #holds address of experimental data for each peak
        self.peptide_list=[] #holds the peptide sequence
        self.protein_list=[]
        self.retention_time_list=[]
        self.charges=[]
                      
        f=open(batchfile,'r')
        lines = f.readlines()
        f.close()
        fields = lines[0].split()
        assert fields[0] == 'pep_mod_seq'
        assert fields[1] == 'charge'
        assert fields[2] == 'retention_time'
        assert fields[3] == 'spectra_file'
        #assert fields[4] == 'protein[optional]'
        for newline in lines[1:]:
            peak_info=newline.split()
            charge=int(peak_info[1])
            pepseq=peak_info[0]+'Z'+charge*'X' #adds on a water and hydrogens to account for charge
            retention_time=peak_info[2]
            datfile=peak_info[3]
            try:
                prot_name=peak_info[4]
            except IndexError:
                prot_name='NOT_GIVEN'
        
            res_syms=list(set(pepseq))
            res_mults=list(map(lambda sym:str.count(pepseq,sym),res_syms))  
            res_dists=list(map(lambda sym:ResLib.res_dict[sym],res_syms))
            
            peak_dists=[]
            for i in range(ResLib.num_species):
                new_dist=dists_and_params.Distribution()
                parent_dists=list(map(lambda res:res[i],res_dists))
                new_dist.build_as_dist_sum_with_mults(parent_dists,res_mults)
                peak_dists.append(new_dist)
            
            self.peak_list.append(list(peak_dists))
            self.datfiles.append(datfile)
            self.charges.append(charge)
            self.num_peaks+=1
            self.peptide_list.append(pepseq)
            self.retention_time_list.append(retention_time)
            self.protein_list.append(prot_name)
        
    def get_m_hd(self,peak_idx):
        datfile=self.datfiles[peak_idx]
        
        f=open(datfile,'r')
        first_line=f.readline().split()
        f.close()
        
        first_mz=float(first_line[0])
        first_m=self.charges[peak_idx]*first_mz
        return first_m
    
    def auto_baseline(self,peak_idx):
        datfile=self.datfiles[peak_idx]
        intensities=[]
        
        f=open(datfile,'r')
        for line in f.readlines():
            new_intensity=float(line[1])
            intensities.append(new_intensity)
            
        f.close()
        
        intensities.sort()
        lowest_val=intensities[0]
        first_quartile_val=intensities[len(intensities)//4]
        return (lowest_val+first_quartile_val)/2
        
    def get_peak_data(self,peak_idx):
        datfile=self.datfiles[peak_idx]
        charge=self.charges[peak_idx]
        intensity=dists_and_params.new_mz_0_array()
        m_hd=None
        
        f=open(datfile,'r')
        for line in f.readlines():
            point_data=list(map(float,line.split()))   
            if m_hd is None:
                m_hd=point_data[0]*charge    
                               
            intensity[int((point_data[0]*charge-m_hd)*dists_and_params.scale_mz)]+=point_data[1]*charge
        f.close()
        
        return numpy.array(intensity)
    
    def get_peptide_seq(self,peak_idx):
        return self.peptide_list[peak_idx].split('Z')[0]
    
    def get_peptide_UID(self,peak_idx):
        return self.peptide_list[peak_idx].split('Z')[0]+'[+'+str(self.charges[peak_idx])+']'
    
    def get_peptide_charge(self, peak_idx):
        return self.charges[peak_idx]

    def get_protein_name(self,peak_idx):
        return self.protein_list[peak_idx]
    
    def get_retention_time(self,peak_idx):
        return self.retention_time_list[peak_idx]

    def display_spectra(self,peak_idx):
        intensity=self.get_peak_data(peak_idx)
        plt.plot(dists_and_params.mz_axis,intensity)
        plt.show()