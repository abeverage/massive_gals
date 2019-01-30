import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

def stack(samp_stack,figname):
    flux_convert = np.median(cat_3dhst_clean['f_F140W']/cat_3dhst_clean['f_F160W'])
    
    cont,line,weights_line,weights,psfs,psfs_weight,masks_line,masks = [],[],[],[],[],[],[],[]
    print(len(samp_stack))
    
    # convert f160w global fluxes to f140w. Use this to find median value of stack.
    f160w = samp_stack['hmag_filt'] == 'F160W'
    f140w = samp_stack['hmag_filt'] == 'F140W'
    
    apcorr = samp_stack['flux_auto'][f160w]/samp_stack['flux_aper_1'][f160w]
    mag160= np.array(23.9-2.5*np.log10(samp_stack[f160w]['flux_aper_1']*flux_convert*apcorr))
    mag140 = np.array(samp_stack['hmag_aliza'][f140w])
    
    mag = np.hstack([mag160,mag140])
    m0 = np.median(mag) #median magnitude of stack *in F140W*
    print(m0)
    
    xmin = np.min(samp_stack['flux_Ha'])
    xmax = np.max(samp_stack['flux_Ha'])
    print(np.median(samp_stack['flux_Ha']))

    for i,tab in enumerate(samp_stack):
        #print(tab['hmag_filt'],tab['root'],tab['id'])
        full = fits.open('../final_data/full_new/{0}_{1:05d}.full.fits'.format(tab['root'],tab['id']))
        filt = tab['hmag_filt']
        
        if filt == 'F160W':
            apcorr = tab['flux_auto']/tab['flux_aper_1']
            m_i = 23.9-2.5*np.log10(tab['f160w_flux_aper_1']*flux_convert*apcorr) # new mag is flux*flux_convert
            norm = 10**(0.4*(tab['hmag_aliza']-m0))   # use this to scale up faint galaxies
            
            ## Get psf
            psf_fit = fits.open('../complete/psfs/{0}-f{1}w_psf.fits'.format(tab['root'],filt[1:4]))
            try:
                psf = psf_fit['PSF','DRIZ1'].data * flux_convert #that's the ratio F140w/F160w
            except:
                psf = psf_fit['PSF','DRIZ'].data * flux_convert
            
            
            data = full['DSCI',filt].data * flux_convert #that's the ratio F140w/F160w
            dwht = full['DWHT',filt].data / flux_convert**2
            dwht_new = dwht/np.sum(dwht)
            
        elif filt =='F140W':
            norm = 10**(0.4*(tab['hmag_aliza']-m0))   # use this to scale up faint galaxies
            ## Get psf
            psf_fit = fits.open('../complete/psfs/{0}-f{1}w_psf.fits'.format(tab['root'],filt[1:4]))
            try:
                psf = psf_fit['PSF','DRIZ1'].data
            except:
                psf = psf_fit['PSF','DRIZ'].data
                
            data = full['DSCI',filt].data
            dwht = full['DWHT',filt].data
            dwht_new = dwht/np.sum(dwht)
        
        data_line = full['LINE','Ha'].data
        dwht_line = full['LINEWHT','Ha'].data
        dwht_line_new = dwht_line/np.sum(dwht_line)
        
        segmap = full['SEG'].data

        mask = ((segmap == tab['id']) | (segmap == 0)) & (dwht!=0)
        mask_line = (dwht_line != 0)
        
        weight = mask*dwht_new/(norm**2)
        weights.append(weight)
        weight_line = mask_line*dwht_line_new/(norm**2)
        weights_line.append(weight_line)
        print(tab['root'],tab['id'],tab['hmag_aliza'],norm)#,norm,tab['flux_Ha'])
        #plt.imshow(weight_line)
        #plt.show()
                
#         if ((tab['root'] == 'j100025+021706') & (tab['id'] == 1674)):
#             print(norm,tab['hmag_aliza'],m0)
#             plt.imshow(log_stretch(normalize_stack(data)),origin='lower')
#             plt.show()
#             plt.imshow(log_stretch(normalize_stack(data_line)),origin='lower')
#             plt.show()
#             plt.imshow(weight,origin='lower')
#             plt.show()
        
        cont.append(data*norm*weight)  # This is top of sum - im*mask/(sig^2*norm)
        line.append(data_line*norm*weight_line)
        
        psfs.append(psf)
        psfs_weight.append(1.)
    
    cont_total_weight = np.sum(weights,axis=0)
    print('cont_weight')
    plt.imshow(cont_total_weight,origin='lower')
    plt.show()
    line_total_weight = np.sum(weights_line,axis=0)
    print('line_weight')
    plt.imshow(line_total_weight,origin='lower')
    plt.show()
    avg_cont = np.sum(cont,axis=0)/cont_total_weight
    avg_cont[~np.isfinite(avg_cont)] = 0
    avg_line = np.sum(line,axis=0)/line_total_weight
    avg_line[~np.isfinite(avg_line)] = 0
    
    psf_stack = sum(psfs)/len(psfs_weight)
    
    hdu_cont = fits.PrimaryHDU(avg_cont)
    hdu_line = fits.ImageHDU(avg_line)
    hdu_cont_wht = fits.ImageHDU(cont_total_weight) #1/sqrt(wht)
    hdu_line_wht = fits.ImageHDU(line_total_weight)
    hdu_psf = fits.ImageHDU(psf_stack)
    hdu_total = fits.HDUList([hdu_cont,hdu_cont_wht,hdu_line,hdu_line_wht,hdu_psf])
#     0 - continuum
#     1 - continuum weight
#     2 - line
#     3 - line weight
#     4 - psf

    hdu_total.writeto('../final_data/stack/fits_aas/stack_{0}.fits'.format(figname),overwrite=True)
#    hdu_total.writeto('../final_data/stack/fits/{0}.fits'.format(figname),overwrite=True)
    print('cont')
    plt.imshow(log_stretch(normalize_stack(avg_cont)),origin='lower')
    plt.show()
    print('line')
    plt.imshow(log_stretch(normalize_stack(avg_line)),origin='lower')  
    plt.show()
