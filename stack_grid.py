import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

def stack_grid(table,line,filename):
    log_stretch = LogStretch(a=1)
    fig = plt.figure(figsize=(12, 12))
    columns = math.ceil(np.sqrt(len(table)))
    rows = columns
    
    for i in range(1,len(table)):
        data = fits.open('../final_data/full/{0}_{1:05d}.full.fits'.format(table['root'][i],table['id'][i]))
        seg = fits.open('../final_data/segmap/{0}_{1:05d}.seg.fits'.format(table['root'][i],table['id'][i]))[0].data
        if table['cat'][i] == 'wisp':
            ngrow = 5
            segmap = nd.maximum_filter((seg == table['id'][i]),ngrow)
        else:
            segmap = (seg == table['id'][i])
        
        if line == 1:
            norm,m,M = normalize(data['Line','Ha'].data * segmap,ms=None)
            img = data['Line','Ha'].data
        elif line == 0:
            norm,m,M = normalize(data['DSCI'].data * segmap,ms=None)
            img = data['DSCI'].data
        fig.add_subplot(rows, columns, i)
        plt.imshow(log_stretch(normalize(img,ms=(m,M))),origin='lower', cmap='gray')
        plt.text()
        plt.axis('off')
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.savefig('../final_data/stack/grids/grid_{0}.png'.format(filename),dpi=400)    
    plt.show()
