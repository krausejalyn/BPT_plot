# CODE WRITTEN BY JALYN KRAUSE ON SEPT. 21, 2020

# PLOTTING BPT DIAGRAM FOR eBOSS DATA FOLLOWING L. KEWLEY ET AL. 2006

import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import pandas as pd
from matplotlib import colors as mcolors
colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)


hdu = fits.open('/run/media/krause/Hypatia/spLine_trim_dr16_eboss.fits')
hdu2 = fits.open('/run/media/krause/Hypatia/eboss_simplefit_v1_part1.fits')

# DISPLAY HEADER
#print(hdu2[1].header)

z = hdu[1].data['Z'].tolist()
ha = hdu[1].data['H_ALPHA_FLUX'].tolist()
ha_err = hdu[1].data['H_ALPHA_FLUX_ERR'].tolist()
hb = hdu[1].data['H_BETA_FLUX'].tolist()
hb_err = hdu[1].data['H_BETA_FLUX_ERR'].tolist()
#oi = hdu[1].data['OI_6300_FLUX'].tolist()
oiii = hdu[1].data['OIII_5007_FLUX'].tolist()
oiii_err = hdu[1].data['OIII_5007_FLUX_ERR'].tolist()
sii = hdu[1].data['SII_6717_FLUX'].tolist()
sii_err = hdu[1].data['SII_6717_FLUX_ERR'].tolist()
nii = hdu[1].data['NII_6584_FLUX'].tolist()
nii_err = hdu[1].data['NII_6584_FLUX_ERR'].tolist()
nev = hdu2[1].data['NEV_3426_FLUX'].tolist()
nev_err = hdu2[1].data['NEV_3426_EW'].tolist()


# CREATE DICTIONARY AND DATA FRAME FOR ORGANIZATION

data = ['z', 'ha', 'ha_err', 'hb', 'hb_err', 'oiii', 'oiii_err', 'sii', 'sii_err', 'nii', 'nii_err']
data_dict = {'z':z, 'ha':ha, 'ha_err':ha_err, 'hb':hb, 'hb_err':hb_err, 'oiii':oiii, 'oiii_err':oiii_err, 'sii':sii, 'sii_err':sii_err,  'nii':nii,  'nii_err':nii_err}
df = pd.DataFrame(data_dict)

for name in data: 
    df = df.loc[df[name].notnull()]

# S/N CUTS FOLLOWING L. KEWLEY ET AL.
df = df.loc[(df['ha']/df['ha_err']) >= 3]
df = df.loc[(df['hb']/df['hb_err']) >= 3]
df = df.loc[(df['nii']/df['nii_err']) >= 3]
df = df.loc[(df['sii']/df['sii_err']) >= 3]
df = df.loc[(df['oiii']/df['oiii_err']) >= 3]

 
df1 = df.query('z > 0.1 & z <=0.2')

def nev(df1): 
    x_values = df1.x.values
    y_values = df1.y.values.tolist()
    
    z_1 = df1['z'].values
    ha_1 = df1['ha'].values
    ha_err_1 = df1['ha_err'].values
    hb_1 = df1['hb'].values
    hb_err_1 = df1['hb_err'].values
    oiii_1 = df1['oiii'].values
    oiii_err_1 = df1['oiii_err'].values
    sii_1 = df1['sii'].values
    sii_err_1 = df1['sii_err'].values
    nii_1 = df1['nii'].values
    nii_err_1 = df1['nii_err'].values
    
    kewl_agn_class_line_mod = (0.61/(x_values-0.02-0.1833*(df1.loc[:,"z"].median())))+1.4+0.03*(df1.loc[:,"z"].median())
    
    df1 = []
    
    for x, y, agn_line in zip(x_values, y_values, kewl_agn_class_line_mod):

        if x >=-0.8 and y > agn_line:
            df1.append('agn')
        elif x > 0.0:
            df1.append('agn')
        else: 
            df1.append('starform')
            
    return df1, ha_1, ha_err_1, hb_1, hb_err_1, nii_1, nii_err_1, oiii_1, oiii_err_1, z_1



df1['x'] = np.log10(df1['nii']/df1['ha'])
df1['y'] = np.log10(df1['oiii']/df1['hb'])

df2 = nev(df1) 

df2 = pd.DataFrame(df2)
df2 = df2.transpose()
df2.columns = ['Type', 'ha', 'ha_err', 'hb', 'hb_err', 'nii', 'nii_err', 'oiii', 'oiii_err', 'z']

dfagn = df2[df2.Type != 'starform']
dfsf = df2[df2.Type != 'agn']
dflowmetagn = df2[df2.Type != 'agn']

dfagn = dfagn[dfagn.Type != 'lowmetagn']
dfsf = dfsf[dfsf.Type != 'lowmetagn']
dflowmetagn = dflowmetagn[dflowmetagn.Type != 'starform']

def agn_classification_line_mod():
    x = np.linspace(-2, 0, 500)
    y = (0.61/(x-0.02-0.1833*(df2.loc[:,"z"].median())))+1.4+0.03*(df2.loc[:,"z"].median())
    return x, y

def z_zero_kewl_line():
    x = np.linspace(-2, 0, 500)
    y = 0.61/(x-0.02)+1.2
    return x, y

def plot_bpt(nii,halpha,oiii,hbeta,nii2,halpha2,oiii2,hbeta2,nii3,halpha3,oiii3,hbeta3,fig=None,**kwargs):
    halpha=np.array(halpha,dtype=float)
    halpha2=np.array(halpha2,dtype=float)
    halpha3=np.array(halpha3, dtype=float)
    nii=np.array(nii,dtype=float)
    nii2=np.array(nii2,dtype=float)
    nii3=np.array(nii3,dtype=float)
    hbeta=np.array(hbeta,dtype=float)
    hbeta2=np.array(hbeta2,dtype=float)
    hbeta3=np.array(hbeta3,dtype=float)
    oiii=np.array(oiii,dtype=float)
    oiii2=np.array(oiii2,dtype=float)
    oiii3=np.array(oiii3,dtype=float)
    n2ha=np.log10(nii/halpha)
    o3hb=np.log10(oiii/hbeta)
    n2ha2=np.log10(nii2/halpha2)
    o3hb2=np.log10(oiii2/hbeta2)
    n2ha3=np.log10(nii3/halpha3)
    o3hb3=np.log10(oiii3/hbeta3)
    
    plt.style.use('grayscale')
    
    fig = plt.figure(figsize=(8,5))
    ax = fig.add_subplot()
    
    plt.xlim(-2.0,0.75)
    plt.ylim(-1.25,1.75)
    plt.suptitle('BPT eBOSS Data, z = 0.1-0.2 (PRACTICE)')
    
    #plt.subplot(1,2,1)
    plt.ylabel(r'$LOG([OIII](\lambda 5007)/H\beta)$')
    plt.xlabel(r'$LOG([NII](\lambda 6584)/H\alpha)$')
    
    starform = ax.hexbin(n2ha, o3hb, gridsize = 100, cmap = "viridis_r", bins = "log", mincnt = 1) # Kewley et al. 2001a extreme starburst line
    agn = ax.hexbin(n2ha3, o3hb3, gridsize = 100, cmap = "viridis_r", bins = "log", mincnt = 1) # Kauffman et al. 2003a classification line
    
    cb = fig.colorbar(agn, ax=ax)
    cb.set_label('LOG(number of obj.)')
    plt.text(-0.5, 1.35, 'AGN', fontsize=15)
    plt.text(-1.2, -0.1, 'HII', fontsize=15)
    plt.text(-0.2, -0.9, 'Comp', fontsize=15)
    plt.text(0.2, 0, 'LINER', fontsize=15)
    
    agn_class = agn_classification_line_mod()
    z_zero = z_zero_kewl_line()
    
    #plotting the mixing line
    #mix = plt.plot([-0.466, 0.003], [-0.378, 0.979], c='b', ms=3)
    agn_0 = plt.plot(agn_class[0], agn_class[1], c='r', ms=3)
    agn_1 = plt.plot(z_zero[0], z_zero[1], c='b', ms=3, linestyle='--')
    
    
    

    #plt.subplot(1,2,2)
    #plt.xlabel(r'$LOG([SII](\lambda\lambda 6717,31)/H\alpha)$')
    plt.show()

halpha = dfsf['ha'].values
halpha2 = dflowmetagn['ha'].values
halpha3 = dfagn['ha'].values
hbeta = dfsf['hb'].values
hbeta2 = dflowmetagn['hb'].values
hbeta3 = dfagn['hb'].values
nii = dfsf['nii'].values
nii2 = dflowmetagn['nii'].values
nii3 = dfagn['nii'].values
oiii = dfsf['oiii'].values
oiii2 = dflowmetagn['oiii'].values
oiii3 = dfagn['oiii'].values

plot_bpt(nii, halpha, oiii, hbeta, nii2, halpha2, oiii2, hbeta2, nii3, halpha3, oiii3, hbeta3, alpha = .5)
