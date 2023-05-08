# colour_system.py
import numpy as np
from scipy.interpolate import interp1d
import pandas as pd

def xyz_from_xy(x, y):
    """Return the vector (x, y, 1-x-y)."""
    return np.array((x, y, 1-x-y))

def get_csv_columns_interp(file,xcolumn, ycolumn):
    df = pd.read_csv(file)
    wavelengths=list(df[xcolumn].to_numpy())
    vals=list(df[ycolumn].transpose().to_numpy())
    fn = interp1d(wavelengths, vals, kind='cubic', fill_value='extrapolate')
    return fn

class ColourSystem():
    """A class representing a colour system.

    A colour system defined by the CIE x, y and z=1-x-y coordinates of
    its three primary illuminants and its "white point".

    TODO: Implement gamma correction

    """

    # The CIE colour matching function for 380 - 780 nm in 5 nm intervals


    def __init__(self, red, green, blue, white, wavelengths):
        """Initialise the ColourSystem object.

        Pass vectors (ie NumPy arrays of shape (3,)) for each of the
        red, green, blue  chromaticities and the white illuminant
        defining the colour system.

        """
        color_file = "CIEcoords.csv"
        xx=get_csv_columns_interp(color_file,'Wavelength', 'x')
        yy = get_csv_columns_interp(color_file, 'Wavelength', 'y')
        zz = get_csv_columns_interp(color_file, 'Wavelength', 'z')
        self.wavelengths=wavelengths
        self.cmf = np.stack([xx(wavelengths), yy(wavelengths),zz(wavelengths)], axis=1)

        # Chromaticities
        self.red, self.green, self.blue = red, green, blue
        self.white = white
        # The chromaticity matrix (rgb -> xyz) and its inverse
        self.M = np.vstack((self.red, self.green, self.blue)).T
        self.MI = np.linalg.inv(self.M)
        # White scaling array
        self.wscale = self.MI.dot(self.white)
        # xyz -> rgb transformation matrix
        self.T = self.MI / self.wscale[:, np.newaxis]

    def xyz_to_rgb(self, xyz, out_fmt=None):
        """Transform from xyz to rgb representation of colour.

        The output rgb components are normalized on their maximum
        value. If xyz is out the rgb gamut, it is desaturated until it
        comes into gamut.

        By default, fractional rgb components are returned; if
        out_fmt='html', the HTML hex string '#rrggbb' is returned.

        """

        rgb = self.T.dot(xyz)
        if np.any(rgb < 0):
            # We're not in the RGB gamut: approximate by desaturating
            print("Outside RGB")
            w = - np.min(rgb)
            rgb += w
        if not np.all(rgb==0):
            # Normalize the rgb vector
            rgb /= np.max(rgb)

        if out_fmt == 'html':
            return self.rgb_to_hex(rgb)
        return rgb

    def xyz_to_rgb_mat(self,xyz):
        mat=np.array([[3.240479,-1.537150,-0.498535],[-0.969256 ,1.875992 ,0.041556],[0.055648,-0.204043 ,1.057311]])
        rgb=mat.dot(xyz)

        return self.postprocess_rawrgb_to_hex(rgb)


    def postprocess_rawrgb_to_hex(self,rawrgb):
        if np.any(rawrgb < 0):
            # We're not in the RGB gamut: approximate by desaturating
            print("Outside RGB")
            w = - np.min(rawrgb)
            rawrgb += w
        if not np.all(rawrgb==0):
            # Normalize the rgb vector
            rawrgb /= np.max(rawrgb)
        return self.rgb_to_hex(rawrgb)


    def rgb_to_hex(self, rgb):
        """Convert from fractional rgb values to HTML-style hex string."""

        hex_rgb = (255 * rgb).astype(int)
        return '#{:02x}{:02x}{:02x}'.format(*hex_rgb)

    def spec_to_xyz(self, spec):
        """Convert a spectrum to an xyz point.

        The spectrum must be on the same grid of points as the colour-matching
        function, self.cmf: 380-780 nm in 5 nm steps.

        """

        XYZ = np.sum(spec[:, np.newaxis] * self.cmf, axis=0)
        den = np.sum(XYZ)
        if den == 0.:
            return XYZ
        return XYZ / den

    def cri_test_xyz(self):
        fname="CIE-test.csv"
        ncolors=14
        output=np.zeros((ncolors,3))
        for i in range(ncolors):
            cfn=get_csv_columns_interp(fname, "Wavelength(nm)", "TCS"+str(i+1))
            cspectrum=cfn(self.wavelengths)
            output[i,:]=self.spec_to_xyz(cspectrum)
        return output

    def cri_test_uvY(self,spec,wspec):
        fname="CIE-test.csv"
        ncolors=14
        output=np.zeros((ncolors,3))
        for i in range(ncolors):
            cfn=get_csv_columns_interp(fname, "Wavelength(nm)", "TCS"+str(i+1))
            cspectrum=np.array(cfn(self.wavelengths))*spec
            inter=self.spec_to_xyz(cspectrum)
            output[i, :]=self.xyz_to_uvY(inter)

        woutput = np.zeros((ncolors, 3))
        for i in range(ncolors):
            cfn = get_csv_columns_interp(fname, "Wavelength(nm)", "TCS" + str(i + 1))
            cspectrum = np.array(cfn(self.wavelengths)) * wspec
            inter = self.spec_to_xyz(cspectrum)
            woutput[i, :] = self.xyz_to_uvY(inter)
        return output, woutput

    def xyz_to_uvY(self,xyz):
        temp = np.zeros((3))
        temp[0] = 4 * xyz[0] / (xyz[0] + xyz[1] * 15.0 + 3 * xyz[2])
        temp[1] = 6 * xyz[1] / (xyz[0] + xyz[1] * 15.0 + 3 * xyz[2])
        temp[2] = xyz[1]
        return temp

    def cri(self,spec,white_spectrum):
        uvytest, uvyref=self.cri_test_uvY(spec,white_spectrum)
        ref=self.xyz_to_uvY(self.spec_to_xyz(white_spectrum))

        uref=ref[0]
        vref=ref[1]
        yref=ref[2]

        ncolors=14
        val=0.0
        for i in range(ncolors):
            yrefi= uvyref[i][2]
            urefi= uvyref[i][0]
            vrefi= uvyref[i][1]
            ytesti = uvytest[i][2]
            utesti = uvytest[i][0]
            vtesti = uvytest[i][1]


            Lref=116.0*np.cbrt(yrefi/yref)-16.0
            Ltest=116.0*np.cbrt(ytesti/yref)-16.0
            delL=Lref-Ltest
            delu=13.0*Lref*(urefi-uref)-13.0*Ltest*(utesti-uref)
            delv = 13.0 * Lref * (vrefi - vref) - 13.0 * Ltest * (vtesti - vref)

            val=val+np.sqrt(delL**2+delu**2+delv**2)


        output=100-4.6*val/ncolors
        return output

    def cri_min(self,spec,white_spectrum):
        uvytest, uvyref=self.cri_test_uvY(spec,white_spectrum)
        ref=self.xyz_to_uvY(self.spec_to_xyz(white_spectrum))

        uref=ref[0]
        vref=ref[1]
        yref=ref[2]

        ncolors=14
        val=0.0
        for i in range(ncolors):
            yrefi= uvyref[i][2]
            urefi= uvyref[i][0]
            vrefi= uvyref[i][1]
            ytesti = uvytest[i][2]
            utesti = uvytest[i][0]
            vtesti = uvytest[i][1]


            Lref=116.0*np.cbrt(yrefi/yref)-16.0
            Ltest=116.0*np.cbrt(ytesti/yref)-16.0
            delL=Lref-Ltest
            delu=13.0*Lref*(urefi-uref)-13.0*Ltest*(utesti-uref)
            delv = 13.0 * Lref * (vrefi - vref) - 13.0 * Ltest * (vtesti - vref)

            val=max(val,np.sqrt(delL**2+delu**2+delv**2))


        output=100-4.6*val
        return output

    def spec_to_rgb(self, spec, out_fmt=None):
        """Convert a spectrum to an rgb value."""

        xyz = self.spec_to_xyz(spec)
        #return self.xyz_to_rgb(xyz, out_fmt)
        return self.xyz_to_rgb_mat(xyz)

# illuminant_D65 = xyz_from_xy(0.3127, 0.3291)
# cs_hdtv = ColourSystem(red=xyz_from_xy(0.67, 0.33),
#                        green=xyz_from_xy(0.21, 0.71),
#                        blue=xyz_from_xy(0.15, 0.06),
#                        white=illuminant_D65)
#
# cs_smpte = ColourSystem(red=xyz_from_xy(0.63, 0.34),
#                         green=xyz_from_xy(0.31, 0.595),
#                         blue=xyz_from_xy(0.155, 0.070),
#                         white=illuminant_D65)
#
# cs_srgb = ColourSystem(red=xyz_from_xy(0.64, 0.33),
#                        green=xyz_from_xy(0.30, 0.60),
#                        blue=xyz_from_xy(0.15, 0.06),
#                        white=illuminant_D65)
