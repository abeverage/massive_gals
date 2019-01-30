import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import time
import statmorph
from collections import OrderedDict
import math

def run_everything(file,seg_l,seg_c,tab):
    
    full = fits.open(file)
    # 0 - continuum
    # 1 - continuum weight
    # 2 - line
    # 3 - line weight
    # 4 - psf
    cont = full[0].data
    cont_wht = full[1].data
    line = full[2].data
    line_wht = full[3].data
    psf = full[4].data
    
    seg_sci = (seg_c ==seg_c[tab[0],tab[1]])
    seg_line = (seg_c ==seg_c[tab[0],tab[1]])
    
    mask_sci = (seg_c != seg_c[tab[0],tab[1]]) & (seg_c != 0)
    mask_line = (seg_c != seg_c[tab[0],tab[1]]) & (seg_l != 0)

    plt.imshow(mask_line)
    plt.show()
    
    start = time.time()
    source_morphs = statmorph.source_morphology(cont, segmap = seg_sci, mask = mask_sci, weightmap=1./np.sqrt(cont_wht), psf=psf, sersic_maxiter=5000)
    morph_cont = source_morphs[0]
    
    start = time.time()
    source_morphs = statmorph.source_morphology(line, segmap = seg_line, mask = mask_line, weightmap=1./np.sqrt(line_wht), psf=psf, sersic_maxiter=5000)
    morph_line = source_morphs[0]
    
    #write to output dictionary
    output = OrderedDict() # like dict but better
    for k in morph_cont.__dict__:
        if not k.startswith('_'):
            output['cont_'+k] = morph_cont.__dict__[k]
            output['line_'+k] = morph_line.__dict__[k]
    
    
    ms, vm = make_statmorph_fig(line, morph_line, savename=file[31:],psf=psf)
    ms, vm = make_statmorph_fig(cont, morph_cont, savename=file[31:],psf=psf,ms=ms,vm=vm)
    
    return output

def make_statmorph_fig(data,morph,savename,psf,ms = None, vm = None):
    snp = 25
    ny, nx = data.shape
    y, x = np.mgrid[0:ny, 0:nx] + 0.5
    fitted_model = statmorph.ConvolvedSersic2D(
      amplitude=morph.sersic_amplitude,
      r_eff=morph.sersic_rhalf,
      n=morph.sersic_n,
      x_0=morph.sersic_xc,
      y_0=morph.sersic_yc,
      ellip=morph.sersic_ellip,
      theta=morph.sersic_theta)
    fitted_model.set_psf(psf)  # always required when using ConvolvedSersic2D
    image_model = fitted_model(x, y)
    bg_noise = (1.0 / snp) * np.random.standard_normal(size=(ny, nx))
    fig = plt.figure(figsize=(15,5))
    ax = fig.add_subplot(131)
    
    if ms==None:
        normed, m, M = normalize(data, ms)
        vm = log_stretch(normalize(data,ms=(m,M))).min(),log_stretch(normalize(data,ms=(m,M))).max()
    
    plt.imshow(log_stretch(normalize(data,ms=(m,M))),origin='lower', cmap='gray',vmin=vm[0],vmax=vm[1])
    ax.set_title('Original image')
    ax = fig.add_subplot(132)
    ax.imshow(log_stretch(normalize((image_model + bg_noise),ms=(m,M))), origin='lower', cmap='gray',vmin=vm[0],vmax=vm[1])
    ax.set_title('Fitted model')
    ax = fig.add_subplot(133)
    residual = data - image_model
    ax.imshow(log_stretch(normalize(residual,ms=(m,M))),origin='lower', cmap='gray',vmin=vm[0],vmax=vm[1])
    ax.set_title('Residual')

    for ax in fig.axes:
        ax.set_xticklabels([]); ax.set_yticklabels([])

    fig.tight_layout(pad = 2)
    plt.show()
    #plt.savefig('../final_data/stack/statmorph/statmorph_{0}.png'.format(savename.replace('.fits','')))

    return ms, vm