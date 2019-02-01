from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.visualization import LogStretch

def normalize_stack(image):
        m, M = np.min(image), np.max(image)
        m, M = np.percentile(image,[1,99])
        M *= 4
        m = -0.1*M
        return (image-m) / (M-m)
    
        
def select_stack_gals(row,final):
    if row['mass'] == 10.5:
        m = final['mass'] < 11  #select mass range
    
    elif row['mass'] == 11:
        m = final['mass'] > 11  #select mass range
    
    z = (final['redshift']> row['zmin']) & (final['redshift']<= row['zmax']) # select z range
    cat = final['cat'] == row['cat'] # select catalog
    c = m&z&cat
    
    return c
    
def convert_mags(table,flux_convert):
    apcorr = table['flux_auto']/table['flux_aper_1']
    filt = table['hmag_filt']

    # Column 'stack_mag' has the converted magnitudes
    table['stack_mag'] = table['hmag_aliza']
    table['stack_mag'][filt =='F160W'] = 23.9-2.5*np.log10(table[filt=='F160W']['flux_aper_1']*flux_convert*apcorr[filt=='F160W'])

    # Column 'flux_convert' is 1 for F140W and != 1 for F160W
    table['flux_convert'] = np.ones(len(table))
    table['flux_convert'][filt=='F160W'] = flux_convert
    
    m0 = stack_median_mag(table)
    
    return table,m0

def stack_median_mag(table):
    
    m0 = np.median(table['stack_mag']) # This column is defined in 'convert_mags'
    
    return m0

def Norm(row,m0):
    
    norm = 10**(0.4*(row['stack_mag']-m0))
        
    return norm

def get_data(row,flux_convert):
    full = fits.open('../final_data/full_new/{0}_{1:05d}.full.fits'.format(row['root'],row['id']))
    seg = full['SEG'].data
    data = full['DSCI',row['hmag_filt']].data * flux_convert
    dwht = full['DWHT',row['hmag_filt']].data / flux_convert**2
    dwht_new = dwht/np.sum(dwht)
    
    data_line = full['LINE','Ha'].data
    dwht_line = full['LINEWHT','Ha'].data
    dwht_line_new = dwht_line/np.sum(dwht_line)
    
    psf_fit = fits.open('../complete/psfs/{0}-f{1}w_psf.fits'.format(row['root'],row['hmag_filt'][1:4]))
    try:
        psf = psf_fit['PSF','DRIZ1'].data * flux_convert #that's the ratio F140w/F160w
    except:
        psf = psf_fit['PSF','DRIZ'].data * flux_convert
    
    return seg,data,dwht,dwht_new,data_line,dwht_line,dwht_line_new,psf

def show_stacks(row,stack_n,im1,im2,imtype,save):
    log_stretch = LogStretch(a=1)
    
    fig = plt.figure(figsize=(8,4))
    fig.suptitle('{0}: {1:.1f} < z < {2:.1f}, N = {3}, M = {4}'.format(row['cat'],row['zmin'],row['zmax'],stack_n,row['mass']),fontsize=15)
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax.imshow(log_stretch(normalize_stack(im1)),origin='lower', cmap='gray')
    ax.text(5,140,'F140W',fontsize = 22,color='white')
    
    ax2.imshow(log_stretch(normalize_stack(im2)),origin='lower', cmap='gray',)
    ax2.text(5,140,r'H$\alpha$',fontsize = 22,color='white')
    
    if save == True:
        fig.savefig('../final_data/stack/stacks/z{0:.1f}-{1:.1f}_M{2}_{3}_{4}.png'.format(row['zmin'],row['zmax'],row['mass'],row['cat'],imtype),dpi = 300)


    
    
 
    