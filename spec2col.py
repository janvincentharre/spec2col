# by Jan-Vincent Harre (janvincent.harre@stud.uni-goettingen.de) and Rene Heller (heller@mps.mpg.de)
# January 2021
# for details, see Harre and Heller (2021): https://arxiv.org/abs/2101.06254

import numpy as np
import os
from astropy.io import fits
import statistics


#### Start of code section by courtesy of Christian Hill ####
#### https://scipython.com/blog/converting-a-spectrum-to-a-colour
def xyz_from_xy(x, y):
    """Return the vector (x, y, 1-x-y)."""
    return np.array((x, y, 1-x-y))

class ColourSystem:
    """A class representing a colour system.

    A colour system defined by the CIE x, y and z=1-x-y coordinates of
    its three primary illuminants and its "white point".

    TODO: Implement gamma correction

    """

    # The CIE colour matching function for 380 - 780 nm in 5 nm intervals
    #cmf = np.loadtxt('cie-cmf.txt', usecols=(1,2,3))
    cmf = np.loadtxt('cmf_precise.txt', usecols=(1,2,3))	

    def __init__(self, red, green, blue, white):
        """Initialise the ColourSystem object.

        Pass vectors (ie NumPy arrays of shape (3,)) for each of the
        red, green, blue  chromaticities and the white illuminant
        defining the colour system.

        """

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
            w = - np.min(rgb)
            rgb += w
        if not np.all(rgb==0):
            # Normalize the rgb vector
            rgb /= np.max(rgb)

        if out_fmt == 'html':
            return self.rgb_to_hex(rgb)
        return rgb

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

    def spec_to_rgb(self, spec, out_fmt=None):
        """Convert a spectrum to an rgb value."""

        xyz = self.spec_to_xyz(spec)
        return self.xyz_to_rgb(xyz, out_fmt)

illuminant_D65 = xyz_from_xy(0.3127, 0.3291)
cs_hdtv = ColourSystem(red=xyz_from_xy(0.67, 0.33),
                       green=xyz_from_xy(0.21, 0.71),
                       blue=xyz_from_xy(0.15, 0.06),
                       white=illuminant_D65)

cs_smpte = ColourSystem(red=xyz_from_xy(0.63, 0.34),
                        green=xyz_from_xy(0.31, 0.595),
                        blue=xyz_from_xy(0.155, 0.070),
                        white=illuminant_D65)

cs_srgb = ColourSystem(red=xyz_from_xy(0.64, 0.33),
                       green=xyz_from_xy(0.30, 0.60),
                       blue=xyz_from_xy(0.15, 0.06),
                       white=illuminant_D65)
#### End of code section by courtesy of Christian Hill ####


def shape_spectrum(twod_spec):
    """Match the input spectrum to the CIE CMF table with lambda from 3900 Angstrom to 8300 Angstrom, in steps of 1 Angstrom,
    if the input spectrum is not in that shape
    """
    indices = []
    spec_result = []
    wl_range = np.arange(3900, 8301, 1) #reference to how it needs to be
    comparison = np.intersect1d(twod_spec[0], wl_range, return_indices=True) #check where the arrays are the same
    indices.append(comparison[1]) #extract the indices
    spec_result = twod_spec[1][indices[0]] #resulting spectrum value, without the wavelength values
    return np.array(spec_result)

def spec_to_color(spec):
    """small function using the functions of the module color_system in order to convert the spectrum to rgb and hex colors
        returns an array with the rgb and hex values
    """
    spec = shape_spectrum(spec) #bring the input spectrum into the right shape
    rgb = cs_hdtv.spec_to_rgb(spec) #get rgb color values
    hex_ = cs_hdtv.spec_to_rgb(spec, out_fmt='html') #get hex color value
    return np.array([rgb,[hex_,0,0]])

def cut_data(wave):
    """ takes out the wavelength and flux values out of the data, which cannot be used for comparing the data with the CMF values (hence all floats) and converts the data to integers
    """
    eps = 0.49
    ind = []
    wl_int = []
    for i in range(len(wave)):
        if wave[i] % 1 < eps and int(wave[i-1]) != int(wave[i]): # making sure that no wavelength is taken multiple times
            if int(wave[i]) >= 3900 and int(wave[i]) <= 8300:
                wl_int.append(int(wave[i]))
                ind.append(i)
    return ind, wl_int

def extract_color_fits(wavelength_data):
    path = "spectra/"
    wavelengths = get_wavelength_fits(wavelength_data)
    datalist = os.listdir(path)
    indices, wavelengths_int = cut_data(wavelengths)
    result = open('color_result.txt','w')
    result.write('# spectrum name \t R G B (Space) \t Hex (Space) \n')
    for a in range(len(datalist)):
        flux_fits = path + datalist[a]
        hdul_flux_fits = fits.open(flux_fits)
        flux = hdul_flux_fits[0].data # flux array
        model_name = flux_fits.split(".fits")[0].split("/")[-1] #get the name and from that temperature, log_g and Z
        print(model_name)
        flux_int = []
        for i in indices:
            flux_int.append(int(flux[i])) # integer flux values, needed for the next line
        spectrum_fits = np.vstack((wavelengths_int, flux_int))    
        color = spec_to_color(spectrum_fits) # Color computation
        result.write("%s, %s,%s,%s, %s \n" % (datalist[a], color[0][0], color[0][1], color[0][2], color[1][0]))
    result.close()
    print('Done!')

def get_wavelength_fits(wavelength_file):
    wave_fits = wavelength_file # always the same, so it can be loaded beforehand
    hdul_wave_fits = fits.open(wave_fits)
    return hdul_wave_fits[0].data # wavelength array

def extract_color_table():
    path = "spectra/"
    datalist = os.listdir(path)
    result = open('color_result.txt','w')
    result.write('# spectrum name \t R G B (Space) \t Hex (Space) \n')
    for i in range(len(datalist)):
        print("Determining color of %s now." % (datalist[i]))
        data = np.loadtxt(path + datalist[i], delimiter=",") # Load the data values, change delimiter here
        data = np.transpose(data)
        wavelengths = data[0]
        flux = data[1]
        indices, wavelengths_int = cut_data(wavelengths)
        flux_int = []
        for j in indices:
            flux_int.append(int(flux[j])) # integer flux values, needed for the next line
        spectrum = np.vstack((wavelengths_int, flux_int))    
        color = spec_to_color(spectrum) # Color computation
        result.write("%s, %s,%s,%s, %s \n" % (datalist[i], color[0][0], color[0][1], color[0][2], color[1][0]))
    result.close()
    print('Done!')

##################################################    MAIN     #############################################################
datatype = os.listdir("spectra/")[0].split(".")[-1]
if datatype == "fits":
    wavelength_data = "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
    extract_color_fits(wavelength_data)
elif datatype == "csv" or datatype == "txt":
    extract_color_table()

