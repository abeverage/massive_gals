# Compile photometric catalog for 3dhst
version = 'v4.1'

cats = [ utils.read_catalog('../tables/{0}_3dhst.{1}.cats/Catalog/{0}_3dhst.{1}.cat'.format(field, version)) for field in ['uds','aegis','cosmos','goodss','goodsn']]
cat = vstack(cats)
cat['cat'] = '3dhst'

zouts = [ utils.read_catalog('../tables/{0}_3dhst.{1}.cats/Eazy/{0}_3dhst.{1}.zout'.format(field, version)) for field in ['uds','aegis','cosmos','goodss','goodsn']]
zout = vstack(zouts)

fouts = [ utils.read_catalog('../tables/{0}_3dhst.{1}.cats/Fast/{0}_3dhst.{1}.fout'.format(field, version)) for field in ['uds','aegis','cosmos','goodss','goodsn']]
fout = vstack(fouts)

colors = [ utils.read_catalog('../tables/{0}_3dhst.{1}.cats/RF_colors/{0}_3dhst.{1}.master.RF'.format(field, version)) for field in ['uds','aegis','cosmos','goodss','goodsn']]
color = vstack(colors)

cat_3dhst = hstack([cat,zout,fout,color])

ok = (cat_3dhst['use_phot'] == 1)
ok &= cat_3dhst['star_flag'] !=1
ok &= (cat_3dhst['z_peak'] > 0.1) & (cat_3dhst['z_peak'] < 1.85) & np.isfinite(cat_3dhst['z_peak'])

cat_3dhst_clean = cat_3dhst[ok]
cat_3dhst_clean.write('../final_data/cat_3dhst_clean.fits')

len(cat_3dhst)