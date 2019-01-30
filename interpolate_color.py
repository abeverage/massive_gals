import numpy as np
from astropy.table import Table, vstack, Column, hstack,unique
import matplotlib.pyplot as plt

def mag(flux):
    mag = 25-2.5*np.log10(flux)
    return(mag)

def color(fluxb,fluxr):
    color = -2.5*np.log10(fluxb/fluxr)
    return(color)

def interpolate_color(cat):
    '''
    Input:
    3dhst catalog
    
    Functionality:
    Take 3dhst full catalog (used to generate models) and 
    calculated the upper and lower bound of color (f125w-f160w)
    for an assortment of redshift bins
    
    Return:
    upper and lower color bound to be interpolated and applied to 
    grism sample
    
    '''
    mag_140 = Column(25 - 2.5*np.log10(np.maximum(cat['f_F140W'], 1e-4)),name = 'mag_F140W')
    mag_160 = Column(25 - 2.5*np.log10(np.maximum(cat['f_F160W'], 1e-4)),name = 'mag_F160W')
    mag_125 = Column(25 - 2.5*np.log10(np.maximum(cat['f_F125W'], 1e-4)),name = 'mag_F125W')
    color = Column((mag_125 - mag_160),name = 'color')
    
    clip = (cat['star_flag'] != 1) & (cat['use_phot'] == 1)
    clip &=  (mag_160 > 0) & (mag_160 < 28) & (mag_125 > 0) & (mag_125 < 28) & np.isfinite(cat['lmass']) & (cat['lmass']>9)
    clip &= (cat['z_peak'] > 0.02) & (cat['star_flag'] != 1) & (cat['use_phot'] == 1) #& (cat['z_peak']<1.875)

    cat_interp = cat[clip]
    cat_interp.add_column(mag_160[clip], index=0)
    cat_interp.add_column(mag_140[clip], index=0)
    cat_interp.add_column(mag_125[clip], index=0)
    cat_interp.add_column(color[clip], index=0)

    #print(len(cat_interp))

    cat_interp['colorf125'] = cat_interp['mag_F125W']
    # Bin uvista data and calculate hmag cut at each bin
    bins = np.arange(0.25,2,0.25)
    inds = np.digitize(cat_interp['z_peak'],bins+.125)
    
    color_up = np.zeros(7)
    color_low = np.zeros(7)
    for i in (np.arange(7)):
        zslice = cat_interp[inds==i]
        #print(np.median(zslice['z_peak']))
        color_up[i],color_low[i] = np.percentile(zslice['color'],[97,3])

    # Check out the hmag cuts

    plt.scatter(bins,color_up)
    plt.scatter(bins,color_low)
    plt.xlabel('z_peak')
    plt.ylabel('color')
    #plt.savefig('../figures/uvista_hmag_cut.png')
    plt.show()
    
    
    return(bins,color_up,color_low)