"""Random collection of utility functions"""

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
