
from astropy.io import fits as pf
import numpy as np
from astropy import units as u
from astropy import cosmology
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

import sys
from scipy.interpolate import interp1d
from scipy import integrate
import treecorr

# setting up config file and correlator for treecorr
config_file = 'default.params'
config = treecorr.config.read_config(config_file)

zmin = float(sys.argv[1])
zmax = float(sys.argv[2])
Nz = int(sys.argv[3])
Maglim1 = float(sys.argv[4])
Maglim2 = float(sys.argv[5])
lambmin = float(sys.argv[6])
lambmax = float(sys.argv[7])

clusters = pf.open(sys.argv[8])[1].data
randoms = pf.open(sys.argv[9])[1].data

galaxies = pf.open(sys.argv[10])[1].data
galaxy_randoms = pf.open(sys.argv[11])[1].data
jkid = int(sys.argv[12])
tot_area = float(sys.argv[13])
outfile = sys.argv[14]
nR = int(sys.argv[15])

# read and mask JK sub-sample
JK = clusters['JK']
mask = (JK!=jkid)*(clusters['Z']>=zmin)*(clusters['Z']<zmax)*(clusters['LAMBDA']>=lambmin)*(clusters['LAMBDA']<lambmax)
RA = clusters['RA'][mask]
DEC = clusters['DEC'][mask]
Z = clusters['Z'][mask]
LAMB = clusters['LAMBDA'][mask]
#print len(Z)

JK_ran = randoms['JK']
mask = (JK_ran!=jkid)*(randoms['Z']>=zmin)*(randoms['Z']<zmax)*(randoms['LAMBDA']>=lambmin)*(randoms['LAMBDA']<lambmax)
RA_ran = randoms['RA'][mask]
DEC_ran = randoms['DEC'][mask]
Z_ran = randoms['Z'][mask]
W_ran = randoms['W'][mask]
#print len(Z_ran)

JK_gal = galaxies['JK']
#type_gal = pf.open('/project/kicp/chihway/brutus/splashback/data/sdss/sdss_gal_galtype_wcolor.fits')[1].data['class']
mask = (JK_gal!=jkid)
ra = galaxies['RA'][mask]
dec = galaxies['DEC'][mask]
mag = galaxies['MAG_I'][mask]

JK_gal_ran = galaxy_randoms['JK']
N_allgal = len(JK_gal_ran)
mask = (JK_gal_ran!=jkid)
ra_ran = galaxy_randoms['RA'][mask]
dec_ran = galaxy_randoms['DEC'][mask]
N_jkgal = len(ra_ran)

n1 = np.histogram(Z, range=(zmin,zmax), bins=Nz)
zmid = (n1[1][1:]+n1[1][:-1])/2

# calculate area
#tot_area = 1403.96471333 # this is for i<21.5
#tot_area = 1410.29568398 # this is for i<19.0
area = tot_area*(N_jkgal*1.0/N_allgal)

# Measurement parameters
h = 0.7
#nR = 30
Rmin = 0.1/h   #Mpc
Rmax = 30.0/h
lnrperp_bins = np.linspace(np.log(Rmin), np.log(Rmax), num = nR+1)
R_edge = np.exp(lnrperp_bins)
R_mid = np.sqrt(R_edge[:-1] * R_edge[1:])
bslop = 0.03 

# treecorr
Xi = []
Var = []
N = []
Nran = []
n = []
nran = []
Ave_dens = []
W = []

for i in range(Nz):
    
    print("bin:", i, "mean z:", zmid[i])    

    mask = (Z>=n1[1][i])*(Z<n1[1][i+1])
    mask_ran = (Z_ran>=n1[1][i])*(Z_ran<n1[1][i+1])
    
    z_cen_bin = Z[mask]
    lamb_cen_bin = LAMB[mask]
    ra_cen_bin = RA[mask]
    dec_cen_bin = DEC[mask] 
    
    z_ran_bin = Z_ran[mask_ran]
    ra_ran_bin = RA_ran[mask_ran]
    dec_ran_bin = DEC_ran[mask_ran]  
    weight_ran_bin = W_ran[mask_ran] 
    
    M = mag - 5*(np.log10(cosmo.luminosity_distance(zmid[i]).value*1e6) - 1)
    M = M - 5.0*np.log10(0.7)
    mask_gal = (M>Maglim1)*(M<Maglim2)
    
    ra_gal_bin = ra[mask_gal]
    dec_gal_bin = dec[mask_gal]
    
    ra_gal_ran_bin = ra_ran.copy()
    dec_gal_ran_bin = dec_ran.copy()

    area_Mpch = area*(np.pi/180.)**2*(cosmo.comoving_distance(zmid[i]).value)**2
    #ave_density = len(z_gal_bin)*1.0/(area/factor2) # number per (Mpc/h)^2
    ave_density = len(ra_gal_bin)*1.0/area_Mpch
    #Ave_dens.append(ave_density)

    if len(z_cen_bin)>1 and len(ra_gal_bin)>10:

        Ave_dens.append(ave_density)

        # Convert physical to angular distance at this zl
        D_l = cosmo.comoving_distance(zmid[i]).value # comiving distance
        #D_l = cosmo.comoving_distance(zmid[i]) / (1.+ zmid[i]) # physical distance
      
        thmin = np.arctan(Rmin / D_l) * (180./np.pi) * 60.      #arcmin
        thmax = np.arctan(Rmax / D_l) * (180./np.pi) * 60.

        #thmin = (Rmin / lmean**(4./9) / D_l) * (180./pi) * 60.*lmeanall**(4./9)        
        #thmax = (Rmax / lmean**(4./9) / D_l) * (180./pi) * 60.*lmeanall**(4./9)


        #thmin = (Rmin / D_l) * (180./np.pi) * 60.                     # arcmin
        #thmax = (Rmax / D_l) * (180./np.pi) * 60.
        #thmin = Rmin*(cosmo.arcsec_per_kpc_comoving(zmid[i])*1000./60)
        #thmax = Rmax*(cosmo.arcsec_per_kpc_comoving(zmid[i])*1000./60)

        #thmin = Rmin/(cosmo.arcsec_per_kpc_comoving(zmid[i])*1000./60)
        #thmax = Rmax/(cosmo.arcsec_per_kpc_comoving(zmid[i])*1000./60)

        # Define catalogs
        cen_cat = treecorr.Catalog(ra=ra_cen_bin, dec=dec_cen_bin,
                                    ra_units='degrees', dec_units='degrees')
        ran_cat = treecorr.Catalog(ra=ra_ran_bin, dec=dec_ran_bin,
                                    ra_units='degrees', dec_units='degrees', w=weight_ran_bin)
        gal_cat = treecorr.Catalog(ra=ra_gal_bin, dec=dec_gal_bin,
                                   ra_units='degrees', dec_units='degrees')
        gal_ran_cat = treecorr.Catalog(ra=ra_gal_ran_bin, dec=dec_gal_ran_bin,
                                   ra_units='degrees', dec_units='degrees')

        # Numerator
        dd = treecorr.NNCorrelation(nbins = nR, min_sep = thmin, max_sep = thmax,
                                        bin_slop = bslop, sep_units = 'arcmin', verbose=2)

        rd = treecorr.NNCorrelation(nbins = nR, min_sep = thmin, max_sep = thmax,
                                        bin_slop = bslop, sep_units = 'arcmin', verbose=2) 

        dr = treecorr.NNCorrelation(nbins = nR, min_sep = thmin, max_sep = thmax,
                                        bin_slop = bslop, sep_units = 'arcmin', verbose=2)

        rr = treecorr.NNCorrelation(nbins = nR, min_sep = thmin, max_sep = thmax,
                                        bin_slop = bslop, sep_units = 'arcmin', verbose=2) 

        dd.process(cen_cat, gal_cat)    # For cross-correlation.
        rd.process(ran_cat, gal_cat)
        dr.process(cen_cat, gal_ran_cat)
        rr.process(ran_cat, gal_ran_cat)

        W.append(rr.npairs*(len(ra_cen_bin)*1.0/len(ra_ran_bin)*len(ra_gal_bin)/len(ra_gal_ran_bin))) # RR
        #W.append(len(ra_cen_bin))  # D      
        xi,varxi = dd.calculateXi(rr, dr=dr, rd=rd)
        Xi.append(xi)
        Var.append(varxi)
        #N.append(dd.npairs)
        N.append(len(ra_cen_bin))
        Nran.append(len(ra_ran_bin))
        n.append(len(ra_gal_bin))
        nran.append(len(ra_gal_ran_bin))

Xi = np.array(Xi)
Var = np.array(Var)
N = np.array(N)
Nran = np.array(Nran)
n = np.array(n)
nran = np.array(nran)
Ave_dens = np.array(Ave_dens)
W = np.array(W)


dens = np.sum(Ave_dens*N)/np.sum(N)
mean = np.sum(Xi*W, axis=0)/np.sum(W, axis=0)*dens
std = (np.sum(Var*W, axis=0)/np.sum(W, axis=0)**2)**-0.5

np.savez(outfile, R=R_mid, mean=mean, std=std, N=N, Nran=Nran, n=n, nran=nran, Ave_dens=Ave_dens)


