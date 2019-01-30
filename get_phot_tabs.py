# Generate script that downloads photcats for all fields
x = utils.GTable.read('../complete/grismcats/wisps-aug10.fits')
fields = np.unique(x['root'])
f = open('get_psf.sh','w+')
print(len(fields))
for i,field in enumerate(fields):
    #f.write('curl -L -o ../complete/photcats/3dhst/{0}_phot.fits https://s3.amazonaws.com/aws-grivam/Pipeline/{1}/Extractions/{1}_phot.fits\n'.format(field,field.replace('+','%2B')))
    #f.write('curl -L -o ../complete/psfs/{0}-f105w_psf.fits https://s3.amazonaws.com/aws-grivam/Pipeline/{1}/Extractions/{1}-f105w_psf.fits\n'.format(field,field.replace('+','%2B')))
    f.write('curl -L -o ../complete/psfs/{0}-f098m_psf.fits https://s3.amazonaws.com/aws-grivam/Pipeline/{1}/Extractions/{1}-f098m_psf.fits\n'.format(field,field.replace('+','%2B')))
f.close()

#!bash get_phots.sh