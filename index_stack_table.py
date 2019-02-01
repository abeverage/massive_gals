import math
import numpy as np
from astropy.table import Table

samp_stack = Table.read('../final_data/clean_sample.fits')

zmin = np.zeros(12)
zmax = np.zeros(12)
mass = np.zeros(12)
cat = np.zeros(12).astype('str')
n= np.zeros(12)

# Make the stacks:
for i in range(4):
    zmin[2*i] = math.floor((0.7+.2*i)*10)/10
    zmin[2*i+1] = math.floor((0.7+.2*i)*10)/10
    zmax[2*i] = math.floor((0.9+.2*i)*10)/10
    zmax[2*i+1] = math.floor((0.9+.2*i)*10)/10
    mass[2*i] = 10.5
    mass[2*i+1] = 10.5
    
    z = (samp_stack['redshift']> math.floor((0.7+.2*i)*10)/10) & (samp_stack['redshift']<= math.floor((0.9+.2*i)*10)/10)
    m = samp_stack['mass'] < 11
    
    cat[2*i] = 'WISP'
    n[2*i] = ((samp_stack['cat']=='wisp') & z & m).sum()
    
    cat[2*i+1] = '3DHST'
    n[2*i+1] = ((samp_stack['cat']!='wisp') & z & m).sum()

    
    
    
for i in range(2):
    j = i+4
    
    zmin[2*j] = math.floor((0.7+.4*i)*10)/10 
    zmax[2*j] = math.floor((1.1+.4*i)*10)/10
    mass[2*j] = 11
    zmin[2*j+1] = math.floor((0.7+.4*i)*10)/10 
    zmax[2*j+1] = math.floor((1.1+.4*i)*10)/10
    mass[2*j+1] = 11
    
    z = (samp_stack['redshift']> math.floor((0.7+.4*i)*10)/10) & (samp_stack['redshift']<= math.floor((1.1+.4*i)*10)/10)
    m = samp_stack['mass'] >= 11
    
    cat[2*j] = 'WISP'
    n[2*j] = ((samp_stack['cat']=='wisp') & z & m).sum()
    
    cat[2*j+1] = '3DHST'
    n[2*j+1] = ((samp_stack['cat']!='wisp') & z & m).sum()



table = Table([zmin,zmax,mass,cat,n],names=('zmin','zmax','mass','cat','n'))
print(table)

table.write('../final_data/stack/index_stack.fits',overwrite=True)
