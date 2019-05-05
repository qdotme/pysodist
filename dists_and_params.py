import numpy
import math
import matplotlib.pyplot as plt

#np=65536
np=65536
ncp=np//2+1
scale_mz=100

fft=lambda x:numpy.fft.rfft(x)
ifft=lambda x:numpy.fft.irfft(x)

mz_axis=list(map(lambda i:i/scale_mz,range(0,np)))
mu_axis=list(map(lambda i:i/np,range(0,ncp)))
     
new_mz_0_array=lambda: numpy.zeros(np,dtype=float)
new_mu_0_array=lambda: numpy.zeros(ncp,dtype=complex)
new_mu_1_array=lambda: numpy.ones(ncp,dtype=complex)

class Distribution(object):
    #############
    #constructors
    #############
    def __init__(self):
        self.mode='unset'
        self.fixed_free='free' #if a distribution will not change (ie a non-variable atom), it may be designated as 'fixed'. This makes Fourier computations faster for it's children.
        self.parents=[]
        self.children=[]
        self.fourier_rep=None
    def build_as_fourier_rep(self,fourier_rep):
        self.fourier_rep=fourier_rep
        self.mode='fourier_rep'
    def build_as_gaussian(self,gw):
        self.gw=gw
        self.mode='gaussian'
    def build_as_shift(self,m_hd,m_off):
        self.m_hd=m_hd
        self.m_off=m_off
        self.mode='shift'
    def build_as_vals_freqs(self,vals,freqs): #used for atoms
        self.values=vals
        self.freqs=freqs
        self.mode='vals_and_freqs'
    def build_as_dist_sum_with_mults(self,dists,mults): #used for residue species, peak species
        self.parents_fixed=[]
        self.parents_free=[]
        self.parent_fixed_mults=[]
        self.parent_free_mults=[]
        self.fourier_fixed=None
        self.fourier_free=None
        self.fixed_free='fixed'
        for i,dist in enumerate(dists):
            if mults[i]!=0:
                if dist.fixed_free=='fixed':
                    self.parents_fixed.append(dist)
                    self.parent_fixed_mults.append(mults[i])
                elif dist.fixed_free=='free':
                    self.parents_free.append(dist)
                    self.parent_free_mults.append(mults[i])
                    self.fixed_free='free'
                dist.children.append(self)
        self.mode='dist_sum_with_mults'
    def build_as_frac_dist_sum(self,dists,freqs): #used for fractional residues
        self.parent_freqs=[]
        for i,dist in enumerate(dists):
            self.parents.append(dist)
            self.parent_freqs.append(freqs[i])
            dist.children.append(self)
        self.mode='frac_dist_sum'
    def build_as_dist_sum_with_amps(self,dists,amps): #used for building the final peak
        self.parent_amps=[]
        for i,dist in enumerate(dists):
            self.parents.append(dist)
            self.parent_amps.append(amps[i])
            dist.children.append(self)
        self.mode='dist_sum_with_amps'
    ##################
    #internal routines
    ##################
    def ft_fixed_internal(self):
        if self.mode=='dist_sum_with_mults':
            ft=new_mu_1_array()
            for i,dist in enumerate(self.parents_fixed):
                ft*=dist.get_ft()**self.parent_fixed_mults[i]
            self.fourier_fixed=ft
    def ft_internal(self):
        if self.mode=='gaussian':
            ft=[]
            _2s = 2*self.gw**2 
            _2a = 2**.5*math.pi*self.gw 
            for mu in mu_axis:
                ft.append(math.e**((-mu**2)/_2s)/_2a)
            self.fourier_rep=numpy.array(ft)
        elif self.mode=='shift':
            ft=[]
            for mu in mu_axis:
                ft.append(math.e**((1j)*2*math.pi*mu*(-self.m_off+self.m_hd)*scale_mz))
            self.fourier_rep=numpy.array(ft)
        elif self.mode=='vals_and_freqs':
            pre_ft=new_mz_0_array()
            for i,val in enumerate(self.values):
                pre_ft[int(val*scale_mz)]=self.freqs[i]
            self.fourier_rep=fft(pre_ft)
        elif self.mode=='dist_sum_with_mults':
            if self.fourier_fixed is None:
                self.ft_fixed_internal()
            ft=new_mu_1_array()
            for i,dist in enumerate(self.parents_free):
                ft*=dist.get_ft()**self.parent_free_mults[i]
            self.fourier_free=ft
            self.fourier_rep=ft*self.fourier_fixed
        elif self.mode=='frac_dist_sum':
            ft=new_mu_0_array()
            for i,dist in enumerate(self.parents):
                ft+=self.parent_freqs[i]*dist.get_ft()
            self.fourier_rep=ft
        elif self.mode=='dist_sum_with_amps':
            ft=new_mu_0_array()
            for i,dist in enumerate(self.parents):
                ft+=self.parent_amps[i]*dist.get_ft()
            self.fourier_rep=ft
        else:
            raise Exception('FT Unknown for Mode: '+self.mode)
    def invalidate_ft(self):
        if self.fourier_rep is not None:
            self.fourier_free=None
            self.fourier_rep=None
            for child in self.children:
                child.invalidate_ft()
    ###################
    #interface routines
    ###################
    def assert_fixed(self):
        self.fixed_free='fixed'
    def change_param(self,new_val,i=None):
        if self.mode=='gaussian':
            self.gw=new_val
            self.invalidate_ft()
        elif self.mode=='shift':
            self.m_off=new_val
            self.invalidate_ft()
        elif self.mode=='vals_and_freqs':
            self.freqs=[1-new_val,new_val]
            self.invalidate_ft()
        elif self.mode=='frac_dist_sum':
            self.parent_freqs=[1-new_val,new_val]
            self.invalidate_ft()
        elif self.mode=='dist_sum_with_amps':
            self.parent_amps[i]=new_val
            self.invalidate_ft()
        else:
            raise Exception('Change Parameter Unknown for Mode: '+self.mode)
    def get_ft(self):
        if self.fourier_rep is None:
            self.ft_internal()
        return self.fourier_rep
    def get_mz(self):
        return ifft(self.get_ft())
    def display_fit(self):
        intensity_array=ifft(self.get_ft())
        plt.plot(list(mz_axis),intensity_array)
        plt.show()
    def save_fit(self, path):
        mz=ifft(self.get_ft())
        f = open(path, 'w')
        f.write('m/z\tintensity\n')
        for i in len(mz_axis):
            f.write(mz_axis[i]+'\t'+mz[i]+'\n')
        f.close()
    def get_mz_axis(self):
        return mz_axis
    

class Parameter(object):
    #############
    #constructors
    #############
    def __init__(self):
        self.mode=None
        self.symbol=None
        self.value=None
        self.modifies=None
    def build_as_baseline_offset(self,init_value,tag):
        self.mode='baseline_offset'
        self.symbol='B'
        self.value=init_value
        self.tag=tag
    def build_as_gaussian_width(self,init_value,tag):
        self.mode='gaussian_width'
        self.symbol='GW'
        self.value=init_value
        self.tag='fixed' if tag=='fixed' else 'free'
    def build_as_machine_offset(self,init_value,tag):
        self.mode='machine_offset'
        self.symbol='M_OFF'
        self.value=init_value
        self.tag='fixed' if tag=='fixed' else 'free'
    def build_as_amplitude(self,label,init_value,spec_idx):
        self.mode='amplitude'
        self.symbol='A_'+label
        self.value=init_value
        self.spec_idx=spec_idx
    def build_as_atom_frac(self,label,init_value):
        self.mode='atom_frac'
        self.symbol='FRAC_'+label
        self.value=init_value
        self.label=label
    def build_as_residue_frac(self,label,init_value):
        self.mode='residue_frac'
        self.symbol='PHI_'+label
        self.value=init_value
        self.label=label
    #########
    #routines
    #########
    def update_object(self):
        if self.modifies is not None:
            if self.mode in ['gaussian_width','machine_offset','atom_frac','residue_frac']:
                self.modifies.change_param(self.value)
            elif self.mode=='amplitude':
                self.modifies.change_param(self.value,self.spec_idx)
    ##########
    #interface
    ##########
    def link_to(self,obj):
        self.modifies=obj
    def change_value(self,new_val):
        self.value=new_val
        self.update_object()
