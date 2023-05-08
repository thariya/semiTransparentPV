# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import pandas
from numpy import pi, linspace, inf, array, arange
from scipy.interpolate import interp1d
import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.constants import h, c, k

from colour_system import ColourSystem, xyz_from_xy

def get_spectra(lambdas):
    rdf = pd.read_excel(r"Spectra.xlsx",sheet_name="Reflectance",skiprows=1, header=None)
    tdf = pd.read_excel(r"Spectra.xlsx", sheet_name="Transmittance", skiprows=1, header=None)
    adf = pd.read_excel(r"Spectra.xlsx", sheet_name="Absorptance", skiprows=1, header=None)

    wavelengths=pd.read_excel(r"Spectra.xlsx",sheet_name="Reflectance",nrows=1, header=None)
    wavelengths=wavelengths.iloc[:,1:].loc[0]
    wavelengths=wavelengths.to_numpy()
    widths=rdf[0].to_numpy()

    del rdf[0]
    del tdf[0]
    del adf[0]

    rdatalist=[]
    tdatalist = []
    adatalist = []
    for i in range(len(widths)):
        valsr = rdf.loc[i].to_numpy()
        valst = tdf.loc[i].to_numpy()
        valsa = adf.loc[i].to_numpy()

        fnr = interp1d(wavelengths, valsr, kind='linear', fill_value='extrapolate')
        fnt = interp1d(wavelengths, valst, kind='linear', fill_value='extrapolate')
        fna = interp1d(wavelengths, valsa, kind='linear', fill_value='extrapolate')

        rdatalist.append(fnr(lambdas))
        tdatalist.append(fnt(lambdas))
        adatalist.append(fna(lambdas))

    rdf=pandas.DataFrame(rdatalist)
    tdf = pandas.DataFrame(tdatalist)
    adf = pandas.DataFrame(adatalist)

    return widths, rdf, tdf, adf

def get_csv_columns_interp(file,xcolumn, ycolumn):
    df = pd.read_csv(file)
    wavelengths=list(df[xcolumn].to_numpy())
    vals=list(df[ycolumn].transpose().to_numpy())
    fn = interp1d(wavelengths, vals, kind='linear', fill_value='extrapolate')
    return fn

def planck(lam, T):
    """ Returns the spectral radiance of a black body at temperature T.

    Returns the spectral radiance, B(lam, T), in W.sr-1.m-2 of a black body
    at temperature T (in K) at a wavelength lam (in nm), using Planck's law.

    """

    lam_m = lam / 1.e6
    fac = h*c/lam_m/k/T
    B = 2*h*c**2/lam_m**5 / (np.exp(fac) - 1)
    return B

def main():

    solar_spectrum_file="PVL spectrum.csv"

    wavelengths=np.array([x/1000.0 for x in range(780,379,-1)])
    widths, rdf, tdf, adf = get_spectra(wavelengths)

    #solar_spectrum_fn=get_csv_columns_interp("si_solar_cell_1d_AM15_approx.csv","Wavelength (nm)","AM 1.5 approx")
    solar_spectrum_fn = get_csv_columns_interp(solar_spectrum_file, "Wavelength (nm)", "Spectral irradiance (W m-2 nm-1)")
    solar_spectrum=solar_spectrum_fn(wavelengths*1000)*1e15
    #solar_spectrum=planck(wavelengths,6500)
    print(solar_spectrum)

    dg = 1000
    ng_fn = get_csv_columns_interp("Glass.csv", "Wavelength(nm)", "Glass_n")
    ng = np.array(ng_fn(wavelengths * 1000))
    kg_fn = get_csv_columns_interp("Glass.csv", "Wavelength(nm)", "Glass_k")
    kg = np.array(kg_fn(wavelengths * 1000))
    alpha=4*math.pi*kg/wavelengths
    indexg=ng+1j*kg

    Rg=abs((1-indexg)/(1+indexg))**2
    Tg = abs((4*indexg) / (1 + indexg)**2)

    newrdf = rdf.copy()
    newtdf = tdf.copy()
    newadf = adf.copy()

    numwidths=len(widths)

    illuminant_D65 = xyz_from_xy(0.3127, 0.3291)
    cs = ColourSystem(red=xyz_from_xy(0.64, 0.33),
                      green=xyz_from_xy(0.30, 0.60),
                      blue=xyz_from_xy(0.15, 0.06),
                      white=illuminant_D65, wavelengths=wavelengths*1000)

    plt.figure("Reflection Spectrum")
    plt.title("Reflection Spectrum")
    plt.xlabel("Wavelength(um)")
    plt.ylabel("Reflectance")
    plt.figure("Transmission Spectrum")
    plt.title("Transmission Spectrum")
    plt.xlabel("Wavelength(um)")
    plt.ylabel("Transmittance")
    plt.figure("EQE Spectrum")
    plt.title("EQE Spectrum")
    plt.xlabel("Wavelength(um)")
    plt.ylabel("EQE(W/m^2/nm)")
    plt.plot(wavelengths, solar_spectrum, label='Sunlight')

    plt.figure("Reflectance colors")
    plt.title("Reflectance colors")
    fig = plt.gcf()
    ax = fig.gca()
    html_rgb = cs.spec_to_rgb(solar_spectrum, out_fmt='html')
    x, y = 0 % 6, -(0 // 6)
    circle = Circle(xy=(x, y * 1.2), radius=0.4, fc=html_rgb)
    ax.add_patch(circle)
    ax.annotate('Sunlight', xy=(x, y * 1.2 - 0.5), va='center', ha='center', color=html_rgb)


    plt.figure("Transmittance colors")
    plt.title("Transmittance colors")
    fig = plt.gcf()
    ax = fig.gca()
    html_rgb = cs.spec_to_rgb(solar_spectrum, out_fmt='html')
    x, y = 0 % 6, -(0 // 6)
    circle = Circle(xy=(x, y * 1.2), radius=0.4, fc=html_rgb)
    ax.add_patch(circle)
    ax.annotate('Sunlight', xy=(x, y * 1.2 - 0.5), va='center', ha='center', color=html_rgb)

    solar_power=np.trapz(solar_spectrum,wavelengths*1000)


    absorbed_power=[]
    transmitted_power=[]
    reflected_power=[]

    for i in range(numwidths):
        width=widths[i]
        R0= rdf.loc[i].to_numpy()
        T0 = tdf.loc[i].to_numpy()
        A0 = adf.loc[i].to_numpy()

        ones=np.ones(np.shape(R0))
        R=Rg+Tg*Tg*R0*np.exp(-2*alpha*dg)/(ones-R0*Rg*np.exp(-2*alpha*dg))
        T=T0*Tg*np.exp(-2*alpha*dg)/(ones-R0*Rg*np.exp(-2*alpha*dg))
        A=Tg*A0*np.exp(-2*alpha*dg)/(ones-R0*Rg*np.exp(-2*alpha*dg))
        R = R * solar_spectrum
        T = T * solar_spectrum
        A = A * solar_spectrum

        newrdf.loc[i]=R
        newtdf.loc[i] = T
        newadf.loc[i] = A

        plt.figure("Reflection Spectrum")
        plt.plot(wavelengths,R, label='{0:.2f} um'.format(width))

        plt.figure("Transmission Spectrum")
        plt.plot(wavelengths, T, label='{0:.2f} um'.format(width))

        plt.figure("EQE Spectrum")
        plt.plot(wavelengths, A, label='{0:.2f} um'.format(width))

        plt.figure("Reflectance colors")
        fig=plt.gcf()
        ax=fig.gca()
        html_rgb = cs.spec_to_rgb(R, out_fmt='html')
        x, y = (i+1) % 6, -((i+1) // 6)
        circle = Circle(xy=(x, y * 1.2), radius=0.4, fc=html_rgb)
        ax.add_patch(circle)
        ax.annotate('{0:.2f} um'.format(width), xy=(x, y * 1.2 - 0.5), va='center',
                    ha='center', color=html_rgb)

        plt.figure("Transmittance colors")
        fig = plt.gcf()
        ax = fig.gca()
        html_rgb = cs.spec_to_rgb(T, out_fmt='html')
        x, y = (i+1) % 6, -((i+1) // 6)
        circle = Circle(xy=(x, y * 1.2), radius=0.4, fc=html_rgb)
        ax.add_patch(circle)
        ax.annotate('{0:.2f} um'.format(width), xy=(x, y * 1.2 - 0.5), va='center',
                    ha='center', color=html_rgb)

        absorbed_power.append(np.trapz(A,wavelengths*1000)/solar_power)
        reflected_power.append(np.trapz(R, wavelengths * 1000)/solar_power)
        transmitted_power.append(np.trapz(T, wavelengths * 1000)/solar_power)

    plt.figure("Reflectance colors")
    fig = plt.gcf()
    ax = fig.gca()
    ax.set_xlim(-0.5, 5.5)
    ax.set_ylim(-5.55, 0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor('k')
    ax.set_aspect("equal")

    plt.figure("Transmittance colors")
    fig = plt.gcf()
    ax = fig.gca()
    ax.set_xlim(-0.5, 5.5)
    ax.set_ylim(-5.55, 0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor('k')
    ax.set_aspect("equal")

    plt.figure("Reflection Spectrum")
    plt.legend(loc="upper left")

    plt.figure("Transmission Spectrum")
    plt.legend(loc="upper left")

    plt.figure("EQE Spectrum")
    plt.legend(loc="upper left")

    plt.figure("Power Distribution")
    plt.title("Solar Power distribution")
    plt.xlabel("Width of the window")
    plt.ylabel("Fraction of power")
    plt.plot(widths,absorbed_power,label="Absorption")
    plt.plot(widths, reflected_power, label="Reflection")
    plt.plot(widths, transmitted_power, label="Transmission")
    plt.legend(loc="upper left")


    plt.show()

if __name__ == '__main__':
    main()

