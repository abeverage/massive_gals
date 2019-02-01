from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

def normalize_stack(image):
        m, M = np.min(image), np.max(image)
        m, M = np.percentile(image,[1,99])
        M *= 4
        m = -0.1*M
        return (image-m) / (M-m)


def plot_stacks(row, img):
    im = fits.open(img)
    cont = im[0].data
    line = im[2].data
    
    xc = row['xc']
    yc = row['yc']
    
    fig = plt.figure(figsize=(8,4))
    fig.suptitle('{0}: {1:.1f} < z < {2:.1f} stacks:  ({3})'.format(row['cat'],row['zmin'],row['zmax'],row['n']),fontsize=15)
    ax = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    tokpc = row['arc_per_pix']*row['kpcperarc']

    
    ax.imshow(log_stretch(normalize_stack(cont)),origin='lower', cmap='gray',
              extent=[-xc*tokpc,xc*tokpc,-xc*tokpc,xc*tokpc])
    ax.set_xlim([-15,15])
    ax.set_ylim([-15,15])
    ax.set_title('F140W',fontsize = 15)
    ax.set_xlabel('[kpc]',fontsize = 15)
    
    ax2.imshow(log_stretch(normalize_stack(line)),origin='lower', cmap='gray',
               extent=[-xc*tokpc,xc*tokpc,-xc*tokpc,xc*tokpc])
    ax2.set_xlim([-25,25])
    ax2.set_ylim([-25,25])
    ax2.set_title(r'H$\alpha$',fontsize = 15)
    ax2.set_xlabel('[kpc]',fontsize = 15)
    plt.show()
    
#    fig.savefig('../final_data/stack/{0}_mass{1}_z{2:.1f}_z{3:.1f}.png'.format(cat,mass,0.7+dz*i,0.7+dz+dz*i))
    print('../final_data/stack/{0}_mass{1}_z{2:.1f}_z{3:.1f}.png'.format(cat,mass,0.7+dz*i,0.7+dz+dz*i))