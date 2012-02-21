"""Random collection of utility functions"""
import numpy as np

def fluxScale(flux_value):
    return 10 ** np.floor(np.log10(flux_value) + 0.5)


def meanEnergy(emin, emax, index_value):
    x = emax / emin
    if index_value == -2.0:
        eflux = emax * np.log(x) / (x - 1)
    elif index_value == -1.0:
        eflux = emin * (x - 1) / np.log(x)
    else:
        factor1 = emin * (index_value + 1) / (index_value + 2)
        factor2 = (x ** (index_value + 2) - 1) / (x ** (index_value + 1) - 1)
        eflux = factor1 * factor2
    return eflux


def calcAngSepDeg(ra0, dec0, ra1, dec1):
    '''Return the angular separation between two objects. Use the
    special case of the Vincenty formula that is accurate for all
    distances'''
    C = np.pi / 180
    d0 = C * dec0
    d1 = C * dec1
    r12 = C * (ra0 - ra1)
    cd0 = np.cos(d0)
    sd0 = np.sin(d0)
    cd1 = np.cos(d1)
    sd1 = np.sin(d1)
    cr12 = np.cos(r12)
    sr12 = np.sin(r12)
    num = np.sqrt((cd0 * sr12) ** 2 + (cd1 * sd0 - sd1 * cd0 * cr12) ** 2)
    den = sd0 * sd1 + cd0 * cd1 * cr12
    return np.arctan2(num, den) / C


def etag(emin, emax, fmt='%07d'):
    return ('emin_%s_emax_%s' % (fmt, fmt)) % (emin, emax)


def cube_to_image(cube, slicepos=None, mean=False):
    """ Make an image out of a cube.
    Both in- and output shoud by pyfits.HDUs"""
    from pyfits import PrimaryHDU
    header = cube.header.copy()
    header['NAXIS'] = 2
    del header['NAXIS3']
    del header['CRVAL3']
    del header['CDELT3']
    if slicepos:
        data = cube.data[slicepos]
    else:
        if mean:
            data = cube.data.mean(0).astype(cube.data.dtype)
        else:
            data = cube.data.sum(0).astype(cube.data.dtype)

    return PrimaryHDU(data, header)


def make_fluxmap(counts_cube_file, exposure_cube_file,
                 counts_map_file, exposure_map_file,
                 flux_map_file,
                 theta = 1):
    import pyfits
    from image.utils import tophat
    print 'Making flux map'

    # Get header and data
    header = pyfits.getheader(counts_cube_file)
    header['NAXIS'] = 2
    del(header['NAXIS3'])

    counts_cube = pyfits.getdata(counts_cube_file)
    exposure_cube = pyfits.getdata(exposure_cube_file)

    theta_pix = int(theta / header['CDELT2'])
    print 'theta_pix', theta_pix

    # Make maps and save to file
    counts_map = counts_cube.sum(0)
    exposure_map = exposure_cube.mean(0)

    # Compute flux and flux error maps
    flux_cube = counts_cube / exposure_cube
    flux_map = flux_cube.sum(0)
    
    # Save all the maps to file    
    pyfits.writeto(counts_map_file, counts_map, 
                   header, clobber=True)
    pyfits.writeto(exposure_map_file, exposure_map, 
                   header, clobber=True)
    pyfits.writeto(flux_map_file, flux_map, 
                   header, clobber=True)

    # Save correlated versions of the files
    # TODO: replace by ftool ftfilter?
    counts_map_corr_file = counts_map_file.replace('.fits', '_corr.fits')
    flux_map_corr_file = flux_map_file.replace('.fits', '_corr.fits')

    counts_map_corr = tophat_correlate(counts_map, theta_pix)
    flux_map_corr = tophat_correlate(flux_map, theta_pix)

    pyfits.writeto(counts_map_corr_file, counts_map_corr, 
                   header, clobber=True)
    pyfits.writeto(flux_map_corr_file, flux_map_corr, 
                   header, clobber=True)
