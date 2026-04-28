from astropy.io import fits
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u


"""
AGN Analysis for the SDSS DR7 MPA/JHU Catalogs

all the files:
1. galinfo: includes redshifts and coords (https://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/SDSS_info.html)
2. galline: includes emission line fluxes and errors (https://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/SDSS_line.html)
3. stellarmasses: includes log scale stellar masses (https://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/Data/stellarmass.html)
the files contain 927,552 galaxies

goals:
1. classify galaxies as star forming, agn, etc (this will help with the agn fraction!)
2. emission line properties
3. host galaxy properties comparison between classes of galaxies

focusing on:
1. agn fraction using kauffman 2003 and kewley 2001 standards for BPT diagram
2. stellar mass comparisons
3. OIII 5007 luminosity comparison as outlined in kauffman paper

"""

info = Table.read(r"/Users/isabel-macpro/agnlab/research practice/data analysis/galinfo")
lines = Table.read(r"/Users/isabel-macpro/agnlab/research practice/data analysis/galline")
mass = Table.read(r"/Users/isabel-macpro/agnlab/research practice/data analysis/stellarmasses")

print(f"galaxies in catalog: {len(info):,}")

#make this ONE dataframe

g = Table()

#from galinfo
g["ra"] = info["RA"]
g["dec"] = info["DEC"]
g["z"] = info["Z"]
g["z_warning"] = info["Z_WARNING"]
g["signalnoise_median"] = info["SN_MEDIAN"]

#from gallines
g["oiii5007_flux"] = lines["OIII_5007_FLUX"]
g["oiii5007_flux_err"] = lines["OIII_5007_FLUX_ERR"]

g["hbeta_flux"] =  lines["H_BETA_FLUX"]
g["hbeta_flux_err"] = lines["H_BETA_FLUX_ERR"]

g["nii6584_flux"] = lines["NII_6584_FLUX"]
g["nii6584_flux_err"] = lines["NII_6584_FLUX_ERR"]

g["halpha_flux"] = lines["H_ALPHA_FLUX"]
g["halpha_flux_err"] = lines["H_ALPHA_FLUX_ERR"]

#this is all from the veilleux & osterbrock paper

#from stellarmasses
g["log_mass"] = mass["MEDIAN"]
g["log_mass_low_uncertainty"] = mass["P16"]
g["log_mass_high_uncertainty"] = mass["P84"]

print(f"  galaxy table built: {len(g):,} rows")

#checking signal to noise ratios
def snr(flux,err):
    return np.divide(flux, err, out = np.full_like(flux, np.nan, dtype = float), where = (err > 0))

oiii5007_snr = snr(g["oiii5007_flux"], g["oiii5007_flux_err"])
hbeta_snr = snr(g["hbeta_flux"], g["hbeta_flux_err"])
nii6584_snr = snr(g["nii6584_flux"], g["nii6584_flux_err"])
halpha_snr = snr(g["halpha_flux"], g["halpha_flux_err"])

snr_threshold = 3 #snr threshold from kauffman 2003

mask = (
    (g["z_warning"] == 0) &

     #going to use kaufmann et al 2003 for low redshift bounds with SDSS DR7
     #also check passbands for sdss spectrograph (3800 - 6150 Å (blue) and 5800 - 9200 Å (red))
    (g["z"] > 0.04) &
    (g["z"] < 0.2)&

    #now gonna check signal to noise ratio for lines!
    (oiii5007_snr > snr_threshold) &
    (hbeta_snr > snr_threshold) &
    (nii6584_snr > snr_threshold) &
    (halpha_snr > snr_threshold) &

    #now gonna just check for errors
    (g["oiii5007_flux"] > 0) &
    (g["hbeta_flux"] > 0) &
    (g["nii6584_flux"] > 0) &
    (g["halpha_flux"] > 0)
)

g_allclean = g[mask].copy()
print(f"\n  after quality cutting: {len(g_allclean):,} galaxies")
print(f"  ({len(g_allclean)/len(g)*100:.1f}% of full catalog)")


#line ratios, finally
g_allclean["log_nii_halpha"] = np.log10(g_allclean["nii6584_flux"]/g_allclean["halpha_flux"])
g_allclean["log_oiii_hbeta"] = np.log10(g_allclean["oiii5007_flux"]/g_allclean["hbeta_flux"])

def kewley(log_nii_halpha):
    #this is the definition for extreme starburst classification lines
    #above this line (> function below), there is no normal starforming hii region
    #CANNOT be starforming --> definite agn (except in cases of shocks/noon-stellar ionization?)
    return 0.61/(log_nii_halpha - 0.47) + 1.19

def kauffman(log_nii_halpha):
    #the demarcation between starburst galaxies and AGN
    #above this line (> than the function below) CANNOT be starforming -> must be agn
    return 0.61/(log_nii_halpha - 0.05) + 1.3

#x axis should be log_nii_halpha and y should be log_oiii_hbeta
x = g_allclean["log_nii_halpha"]
y = g_allclean["log_oiii_hbeta"]

#kewley for upper limit, kauffman for lower
is_starforming  = y < kauffman(x)
is_agn          = y > kewley(x)
is_questionable = (~is_starforming) & (~is_agn)

classification = np.full(len(g_allclean), "questionable", dtype=object)

classification[is_starforming] = "starforming"
classification[is_agn] = "agn"

g_allclean["classification"] = classification

#oiii luminosity
#need to use FlatLamdaCDM to compute distances using redshift

cosmo = FlatLambdaCDM(H0 = 70, Om0 = 0.3) # from astropy website
luminositydistance = cosmo.luminosity_distance(np.array(g_allclean["z"])).to(u.cm).value #cause flux in CGS units according to galline website (techinically 1e-17 erg/s/cm^2)
#use the standard formula F = L / ( 4 * pi * distance^2)
g_allclean["log_oiii_lum"] = np.log10(4 * np.pi * luminositydistance**2 * np.array(g_allclean["oiii5007_flux"]) * 1e-17)

#results

colors = {
    "agn": "#963F3F",
    "questionable": "#673F96",
    "starforming": "#340D6B",
}

classification_array = np.array(g_allclean["classification"])
log_mass_array = np.array(g_allclean["log_mass"])
z_array = np.array(g_allclean["z"])


total = len(g_allclean)
labels, counts = np.unique(g_allclean["classification"], return_counts = True)
print("classification results:")
for label, n in zip(labels, counts):
    print(f"  {label:<15}: {n:>6,}  ({n/total*100:.1f}%)")

agn_index = np.where(labels == "agn")[0][0]
AGN_FRACTION = counts[agn_index]/total * 100

print(f"\n  agn fraction of full catalogs: {AGN_FRACTION:.2f}%")

print("\nstellar mass by classification (log M*/Msun):")
for label in ["starforming", "questionable", "agn"]:
    selection = classification_array == label
    masses = log_mass_array[selection & (log_mass_array > 0)]
    print(f"  {label:<15}: median = {np.median(masses):.2f}  mean = {np.mean(masses):.2f}  std = {np.std(masses):.2f}")

L_array = np.array(g_allclean["log_oiii_lum"])
classification_array = np.array(g_allclean["classification"])

print("\n[OIII] luminosity by classification (log erg/s):")
for label in ["starforming", "questionable", "agn"]:
    selection = (classification_array == label) & np.isfinite(L_array)
    print(f"  {label:<15}: median = {np.median(L_array[selection]):.2f}  mean = {np.mean(L_array[selection]):.2f}  std = {np.std(L_array[selection]):.2f}")

#plotting and visualizations :p
"""
going to do plots for
1. BPT diagram
2. stellar mass distribution by the classification
3. oiii luminosity

"""

#plot 1: bpt
fig1, ax1 = plt.subplots(figsize=(7, 6))

for label, color in colors.items():
    selection = classification_array == label
    ax1.scatter(x[selection], y[selection], s = 0.5, alpha = 0.25, color = color, label=label, rasterized = True)

x_kauf = np.linspace(-1.5,  0.0,  300)
x_kewl = np.linspace(-1.5,  0.35, 300)
ax1.plot(x_kauf, kauffman(x_kauf), color = "#6B0D0D", lw = 1.5, ls = "--", label = "Kauffman2003")
ax1.plot(x_kewl, kewley(x_kewl), color = "#6B0D0D", lw = 1.5, ls = "-", label = "Kewley2001")

ax1.set_xlabel(r"log([NII] 6584 / H$\alpha$)", fontsize=11)
ax1.set_ylabel(r"log([OIII] 5007 / H$\beta$)",  fontsize=11)
ax1.set_title("BPT Diagram", fontsize=13)
ax1.set_xlim(-1.5, 0.8)
ax1.set_ylim(-1.2, 1.5)
ax1.set_facecolor("#F1EDF7")
ax1.legend(markerscale=10, fontsize=8, loc="upper left")
plt.tight_layout()
plt.savefig("bpt_diagram.png", dpi=150, bbox_inches="tight")
#plt.show()

#plot 2: stellar masses
fig2, ax2 = plt.subplots(figsize = (7, 6))

for label, color in colors.items():
    selection = (classification_array == label) & (log_mass_array > 0)
    ax2.hist(log_mass_array[selection], bins = 40, alpha = 0.6, color = color, label = label, density = True) #density normalizes histogram

ax2.set_xlabel(r"log($M_*$ / $M_\odot$)", fontsize=11)
ax2.set_ylabel("Normalized Count", fontsize=11)
ax2.set_title("Stellar Mass by Classification", fontsize=13)
ax2.legend(fontsize=8)
ax2.set_facecolor("#F1EDF7")
plt.tight_layout()
plt.savefig("stellar_mass.png", dpi=150, bbox_inches="tight")
#plt.show()

#plot 3: oiii luminosity
fig3, ax3 = plt.subplots(figsize=(7, 6))

L_array = np.array(g_allclean["log_oiii_lum"])
for label, color in colors.items():
    selection = (classification == label) & np.isfinite(L_array)
    ax3.hist(L_array[selection], bins = 40, alpha = 0.6, color = color, label = label, density = True)

ax3.set_xlabel(r"log($L_{\rm [OIII]}$ / erg s$^{-1}$)", fontsize=11)
ax3.set_ylabel("Normalized Count", fontsize=11)
ax3.set_title("[OIII] 5007 Luminosity by Class", fontsize=13)
ax3.set_xlim(38, 44)
ax3.legend(fontsize=8)
ax3.set_facecolor("#F1EDF7")
plt.tight_layout()
plt.savefig("oiii_luminosity.png", dpi=150, bbox_inches="tight")
plt.show()

print("plots saved!")
