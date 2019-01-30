# Read in the wisp tables you have. the full wisp table contains all objects from 
# early wisp table. Match and remove duplicates.

cat0 = utils.GTable.read('../complete/grismcats/wisps-aug10.fits')
cat0['cat'] = 'wisp'
cat1 = utils.GTable.read('../complete/grismcats/grizli-18.05.17-full.fits')
cat1['cat'] = 'other'
cat2 = utils.GTable.read('../complete/grismcats/j021726-051246.info.fits')
cat2['cat'] = '3dhst'
cat3 = utils.GTable.read('../complete/grismcats/j100025+021706.info.fits')
cat3['cat'] = '3dhst'
cat4 = utils.GTable.read('../complete/grismcats/j141923+525013.info.fits')
cat4['cat'] = '3dhst'
cat = vstack([cat0,cat1,cat2,cat3,cat4])
cat.write('../final_data/cat.fits')

phot1 = utils.GTable.read('../complete/photcats/3dhst_full_phot.fits')
phot2 = utils.GTable.read('../complete/photcats/wisp_full_phot.fits')
phot3 = utils.GTable.read('../complete/photcats/other_full_phot.fits')
phot = vstack([phot1,phot2,phot3])
phot.write('../final_data/phot.fits',overwrite = True)

# Match cats with photcat
phot = utils.GTable.read('../final_data/phot.fits')
cat = utils.GTable.read('../final_data/cat.fits')

idx, dr = phot.match_to_catalog_sky(cat)
print(len(phot),len(cat),len(idx))
new_idx = idx[dr.value<0.4]
photcat = phot[new_idx]; print(len(photcat))
cats = cat[dr.value<0.4]; print(len(photcat))
print(cats['ra'][15000],photcat['ra'][15000])

cats['f160w_flux_aper_2'] = photcat['f160w_flux_aper_2']
cats['f140w_flux_aper_2'] = photcat['f140w_flux_aper_2']
cats['f160w_flux_aper_1'] = photcat['f160w_flux_aper_1']
cats['f140w_flux_aper_1'] = photcat['f140w_flux_aper_1']
cats['f160w_mask_aper_1'] = photcat['f160w_mask_aper_1']
cats['f140w_mask_aper_1'] = photcat['f140w_mask_aper_1']
cats['f160w_flux_aper_0'] = photcat['f160w_flux_aper_0']
cats['f140w_flux_aper_0'] = photcat['f140w_flux_aper_0']
cats['flux_auto'] = photcat['flux_auto']
cats['flux_aper_1'] = photcat['flux_aper_1']
cats['flux_aper_0'] = photcat['flux_aper_0']
cats['flux_aper_2'] = photcat['flux_aper_2']
cats['flux_aper_3'] = photcat['flux_aper_3']
cats['npix'] = photcat['npix']
cats['flag'] = photcat['flag']

cats.write('../full_with_phot.fits',overwrite=True)