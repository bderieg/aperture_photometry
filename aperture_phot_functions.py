from astropy.io import fits
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from photutils import RectangularAperture
from photutils import SkyEllipticalAperture
from photutils import EllipticalAperture
import spreadsheet_interface as spin
import matplotlib.pyplot as plt

# Constants
background_width = 10
fc_K2V_W1 = 1.0038
fc_K2V_W2 = 1.0512
fc_K2V_W3 = 1.0030
fc_K2V_W4 = 1.0013
W1_ZMFD = 306.682
W2_ZMFD = 170.663
W3_ZMFD = 29.045
W4_ZMFD = 8.284


# TODO for whole script: This version does not have image splicing functionality

# FOLLOWING: Static functions

def to_mag(zp, flux):
    return zp - 2.5*np.log10(flux)


def mag_to_jy(mag, zmfd, fc):
    return (zmfd/fc)*(10**(-mag/2.5))


def correction_type(band, fc):
    if fc == 'K2V':
        if band == 'W1':
            return fc_K2V_W1, W1_ZMFD
        elif band == 'W2':
            return fc_K2V_W2, W2_ZMFD
        elif band == 'W3':
            return fc_K2V_W3, W3_ZMFD
        elif band == 'W4':
            return fc_K2V_W4, W4_ZMFD


def list_rms(in_list):
    rms = 0
    for i in range(len(in_list)):
        rms += in_list[i]**2
    rms /= len(in_list)
    rms = np.sqrt(rms)
    return rms


# Converts flux from integrated image values to Jy/px, depending on the telescope (and, thus, the input unit)
def flux_conversion(inp_flux, band, header):
    fc = 'K2V'  # Color correction for WISE

    if "ALMA" in band:
        return inp_flux * (np.pi / 4*np.log(2)) * (header['CDELT2']/header['BMAJ']) * \
        (header['CDELT2']/header['BMIN'])
    elif ("Wide" in band) or ("N" in band):
        return ((inp_flux * 1000000) / (4.25 * (10 ** 10))) * ((header['CD2_2'] * 3600) ** 2)
    elif "W" in band:  # FIXME: This assumes one file
        return mag_to_jy(to_mag(header['MAGZP'], inp_flux), correction_type(band, fc)[1], correction_type(band, fc)[0])
    elif ("IRAC" in band) or ("MIPS" in band):
        return ((inp_flux * 1000000) / (4.25 * (10 ** 10))) * (header['PXSCAL2'] ** 2)


# FOLLOWING: File read functions and auxiliaries

# For use with the aperture_file_read function.  Gets the next number in the file, and converts to degrees if necessary.
def get_num(line, index):
    num = ""
    while (line[index] != ",") and (line[index] != "\"") and (line[index] != ")"):
        num += line[index]
        index += 1
    num = float(num)
    if line[index] == "\"":
        index += 1
        num /= 3600  # If in arc seconds, convert to degrees
    index += 1
    return num, index


# Determines the number of elliptical apertures contained in an aperture file.
def aperture_count(file_name):
    file = open(file_name)
    count = 1
    cur_line = file.readline()
    while "ellipse" not in cur_line:
        cur_line = file.readline()
    cur_line = file.readline()
    while "ellipse" in cur_line:
        cur_line = file.readline()
        count += 1
    return count


# Reads an elliptical aperture from an aperture file (file_name).  It reads the (n)th aperture from a file of multiple.
def aperture_file_read(file_name, n):
    # Open file and define variables
    file = open(file_name)
    index = 0
    cur_line = file.readline()
    ra = ""
    dec = ""
    major_ax = ""
    minor_ax = ""
    tilt = ""

    # Find line with aperture data
    while "ellipse" not in cur_line:
        cur_line = file.readline()

    # If num>1, iterate until next ellipse line
    for i in range(n):
        cur_line = file.readline()

    # Find character where first piece of data starts
    while cur_line[index] != "(":
        index += 1
    index += 1

    # Find values
    ra, index = get_num(cur_line, index)
    dec, index = get_num(cur_line, index)
    major_ax, index = get_num(cur_line, index)
    minor_ax, index = get_num(cur_line, index)
    tilt, index = get_num(cur_line, index)

    file.close()

    return ra, dec, major_ax, minor_ax, tilt

# FOLLOWING: Main aperture photometry function(s)

# This function takes an image_data array (the entire image), and aperture_data array (the aperture over which to
# perform aperture photometry, and the pixel coordinates of the background center (tuple: x, y), and returns the total
# flux through the aperture.
# WARNING: The aperture_data mask must have been made with method='center' or else negative flux values might result
def aperture_phot_manual(image_data, aperture_data, background_center=(0, 0)):
    # Find median background flux
    background_aperture = RectangularAperture(background_center, background_width, background_width)
    mask = background_aperture.to_mask(method='center')
    mask_data = mask.multiply(image_data)
    mask_data = mask_data[1:-1]
    temp_data = np.zeros((background_width-1, background_width-1))
    for i in range(background_width-1):
        temp_data[i] = mask_data[i][1:-1]
    mask_data = temp_data
    background_flux_med = np.median(mask_data)

    # Find total flux by summing aperture data
    total_flux = 0
    num_pix = 0
    for i in range(len(aperture_data)):
        for j in range(len(aperture_data[0])):
            if aperture_data[i, j] > 0:
                total_flux += aperture_data[i, j]
                num_pix += 1

    # Subtract out background flux
    total_flux -= num_pix * background_flux_med

    return total_flux


# This function takes some image, defines the aperture, and [uses other functions to] integrate the flux over the
# aperture.  It returns the flux density over some aperture, as well as the statistical uncertainty gained from the
# 'calc_error' function below.  NOTE: This does not use the unc file to calculate statistical uncertainty.
def image_aperture_phot(file_name, aperture_file_name, background_x, background_y, band):
    # Open .fits file
    print(file_name)
    cur_file = fits.open(file_name)

    # Define some things
    main_x, main_y, main_major, main_minor, main_tilt = aperture_file_read(aperture_file_name, 0)
    sub_aps = []
    for i in range(1, aperture_count(aperture_file_name)):
        sub_aps.append(aperture_file_read(aperture_file_name, i))
    data = cur_file[0].data
    if "ALMA" in band:  # The data is 'thicker' with ALMA
        data = cur_file[0].data[0][0]
    background_center = (background_x, background_y)
    main_center = (main_x, main_y)
    if ("W" in band) and ("Wide" not in band):  # The rotation is weird for WISE apertures
        main_tilt = main_tilt - 90 + cur_file[0].header['CROTA2']

    # Get total flux for subtraction apertures
    sub_total = 0
    for i in range(len(sub_aps)):
        sub_x, sub_y, sub_major, sub_minor, sub_tilt = sub_aps[i]
        sub_coord = SkyCoord(sub_x, sub_y, frame='icrs', unit='deg')
        sub_ap = SkyEllipticalAperture(sub_coord, sub_major*u.deg, sub_minor*u.deg, theta=sub_tilt*u.deg)
        sub_ap = sub_ap.to_pixel(WCS(cur_file[0].header))
        sub_ap_mask = sub_ap.to_mask(method='center')
        sub_ap_mask = sub_ap_mask.multiply(data)
        if aperture_phot_manual(data, sub_ap_mask, background_center) > 0:
            sub_total += aperture_phot_manual(data, sub_ap_mask, background_center)
        else:  # On the off chance the flux value is negative
            pass

    cur_wcs = WCS(cur_file[0].header, naxis=[1, 2])
    # Define aperture, create masks, and pass to manual aperture photometry function
    coord = SkyCoord(main_center[0], main_center[1], frame='icrs', unit='deg')
    ap = SkyEllipticalAperture(coord, main_major*u.deg, main_minor*u.deg, theta=main_tilt*u.deg)
    ap = ap.to_pixel(cur_wcs)
    ap_mask = ap.to_mask(method='center')
    ap_mask = ap_mask.multiply(data)

    # Pass mask to aperture_phot_manual function to integrate flux over aperture
    if aperture_phot_manual(data, ap_mask, background_center) > 0:
        total_flux = aperture_phot_manual(data, ap_mask, background_center)
    else:  # On the off chance the flux value is negative
        return 0.0, 0.0
    total_flux -= sub_total

    # If this function is to evaluate the noise in the aperture, do that instead, overwriting the previous block
    if "Noise" in file_name:
        total_flux = 4*aperture_rms(ap_mask)*len(np.array(ap_mask).ravel())

    # Close .fits file
    cur_file.close()

    # Calculate error (does not yet work for Spitzer bands, so just return)
    try:
        error = calc_error(ap, data)
    except:
        return flux_conversion(total_flux, band, cur_file[0].header), 0

    # Convert units and return flux, error
    return flux_conversion(total_flux, band, cur_file[0].header), flux_conversion(error, band, cur_file[0].header)


def aperture_rms(aperture_data):
    flat_list = []
    for i in range(len(aperture_data)):
        for j in range(len(aperture_data[0])):
            flat_list.append(aperture_data[i][j])
    i = 0
    while i < len(flat_list):
        if flat_list[i] == 0:
            flat_list.pop(i)
            i -= 1
        i += 1

    return list_rms(flat_list)


# FOLLOWING: Image splicing functions for WISE.  These aren't currently called by anything, but rather they're here for
# FOLLOWING: archival purposes, in case they need to be used in the future.

# An aperture photometry function (for WISE images) for when the target is spread over two images vertically.
# Takes an upper and a lower fits file, the coordinates for background subtraction for both files (in pixel), and the
# ICRS coordinates of the aperture to perform photometry over.
def splice_WISE_vertical(file_name_top, file_name_bottom, aperture_file_name, background_x_top, background_y_top,
                         background_x_bottom, background_y_bottom, band):
    # Define some things
    fc = "K2V"
    main_x, main_y, main_major, main_minor, main_tilt = aperture_file_read(aperture_file_name, 0)
    background_center_top = (background_x_top, background_y_top)
    background_center_bottom = (background_x_bottom, background_y_bottom)
    top_file = fits.open(file_name_top)
    top_data = top_file[0].data
    bottom_file = fits.open(file_name_bottom)
    bottom_data = bottom_file[0].data
    magzp = top_file[0].header['MAGZP']
    main_center = (main_x, main_y)
    coord = SkyCoord(main_center[0], main_center[1], frame='icrs', unit='deg')
    full_ap = SkyEllipticalAperture(coord, main_major*u.deg, main_minor*u.deg, theta=main_tilt*u.deg)

    # Define top, then bottom apertures
    top_ap = full_ap.to_pixel(WCS(top_file[0].header))
    top_mask = top_ap.to_mask(method='center')
    top_mask_data = top_mask.multiply(top_data)

    bottom_ap = full_ap.to_pixel(WCS(bottom_file[0].header))
    bottom_mask = bottom_ap.to_mask(method='center')
    bottom_mask_data = bottom_mask.multiply(bottom_data)

    start = (0, 0)
    max_val = 0

    # Find center (in pixel coordinates on the top image) by finding the max value, then zero values below the center
    for i in range(len(top_mask_data)):
        for j in range(len(top_mask_data[0])):
            if top_mask_data[i, j] > max_val:
                max_val = top_mask_data[i, j]
                start = (i, j)
    for i in range(start[0]):
        for j in range(len(top_mask_data[0])):
            top_mask_data[i, j] = 0

    # Find center (in pixel coordinates on the bottom image) by finding the max value,
    # then zero values above and including the center
    for i in range(len(bottom_mask_data)):
        for j in range(len(bottom_mask_data[0])):
            if bottom_mask_data[i, j] > max_val:
                max_val = bottom_mask_data[i, j]
                start = (i, j)
    for i in range(start[0], len(bottom_mask_data)):
        for j in range(len(bottom_mask_data[0])):
            bottom_mask_data[i, j] = 0

    # Perform manual aperture photometry on top and bottom
    top_flux = aperture_phot_manual(top_data, top_mask_data, background_center_top)
    bottom_flux = aperture_phot_manual(bottom_data, bottom_mask_data, background_center_bottom)

    # Convert units to Jy
    top_flux = mag_to_jy(to_mag(magzp, top_flux), correction_type(band, fc)[1], correction_type(band, fc)[0])
    bottom_flux = mag_to_jy(to_mag(magzp, bottom_flux), correction_type(band, fc)[1], correction_type(band, fc)[0])

    # Close files
    top_file.close()
    bottom_file.close()

    return top_flux + bottom_flux


# Almost exactly the same as the splice_WISE_vertical function, but for when a target is spread over two WISE images
# horizontally.
def splice_WISE_horizontal(file_name_left, file_name_right, aperture_file_name, background_x_left, background_y_left,
                         background_x_right, background_y_right, band):
    # Define some things
    fc = "K2V"
    main_x, main_y, main_major, main_minor, main_tilt = aperture_file_read(aperture_file_name, 0)
    background_center_left = (background_x_left, background_y_left)
    background_center_right = (background_x_right, background_y_right)
    left_file = fits.open(file_name_left)
    left_data = left_file[0].data
    right_file = fits.open(file_name_right)
    right_data = right_file[0].data
    magzp = left_file[0].header['MAGZP']
    main_center = (main_x, main_y)
    coord = SkyCoord(main_center[0], main_center[1], frame='icrs', unit='deg')
    full_ap = SkyEllipticalAperture(coord, main_major * u.deg, main_minor * u.deg, theta=main_tilt * u.deg)

    # Define left, then right apertures
    left_ap = full_ap.to_pixel(WCS(left_file[0].header))
    left_mask = left_ap.to_mask(method='center')
    left_mask_data = left_mask.multiply(left_data)

    right_ap = full_ap.to_pixel(WCS(right_file[0].header))
    right_mask = right_ap.to_mask(method='center')
    right_mask_data = right_mask.multiply(right_data)

    start = (0, 0)
    max_val = 0

    # Find center (in pixel coordinates on the left image) by finding the max value,
    # then zero values right of and including the center
    for col in range(len(left_mask_data[0])):
        for row in range(len(left_mask_data)):
            if left_mask_data[row, col] > max_val:
                max_val = left_mask_data[row, col]
                start = (row, col)
    for col in range(start[1], len(left_mask_data[0])):
        for row in range(len(left_mask_data)):
            left_mask_data[row, col] = 0

    # Find center (in pixel coordinates on the right image) by finding the max value,
    # then zero values left of the center
    for col in range(len(right_mask_data[0])):
        for row in range(len(right_mask_data)):
            if right_mask_data[row, col] > max_val:
                max_val = right_mask_data[row, col]
                start = (row, col)
    for col in range(start[1]):
        for row in range(len(right_mask_data)):
            right_mask_data[row, col] = 0

    # Perform manual aperture photometry on left and right
    left_flux = aperture_phot_manual(left_data, left_mask_data, background_center_left)
    right_flux = aperture_phot_manual(right_data, right_mask_data, background_center_right)

    # Convert units to Jy
    left_flux = mag_to_jy(to_mag(magzp, left_flux), correction_type(band, fc)[1], correction_type(band, fc)[0])
    right_flux = mag_to_jy(to_mag(magzp, right_flux), correction_type(band, fc)[1], correction_type(band, fc)[0])

    # Close files
    left_file.close()
    right_file.close()

    return left_flux + right_flux


# FOLLOWING: Error estimation functions

def unc_error(aperture_data):
    # Sum over aperture
    total = 0
    num_pix = 0
    for i in range(len(aperture_data)):
        for j in range(len(aperture_data[0])):
            if aperture_data[i, j] > 0:
                total += aperture_data[i, j]
                num_pix += 1

    # Find error from sum
    mean = total / num_pix
    error = mean * np.sqrt(num_pix)

    return error

# TODO: Make 'factor' automatically detect.  It should be 1 unless the major axis of the aperture is more than half the
# TODO: image size
# FIXME: This doesn't work for Spitzer images
# This function takes the main aperture, copies it 20 random times around the image that don't overlap with the main
# aperture, and returns the standard deviation of the flux through these
def calc_error(main_ap, image_data):
    # Define some things
    factor = .7
    majorax = main_ap.a
    minorax = main_ap.b
    tilt = main_ap.theta
    main_center = main_ap.positions
    flux_values = []

    # Iterate through 20 random apertures
    for i in range(20):
        # Randomly produce a coordinate for the new aperture.  We first randomly find an x component, then choose the y
        # component so that the center of the new aperture is outside of the main aperture
        new_center = main_center.copy()
        new_center[0] = np.random.randint(main_center[0], len(image_data[0]) - .5*majorax)
        allowance = main_center[0]  # This variable is the min y value for one side
        if (1-((np.abs(new_center[0] - main_center[0]))/majorax)) < 1:
            allowance = main_center[0]+factor*majorax*(1-((np.abs(new_center[0] - main_center[0]))/majorax))
        new_center[1] = np.random.randint(allowance, len(image_data[0]) - .5*majorax)
        if np.random.randint(0, 2) == 0:  # Account for the other side
            new_center[1] -= 2 * (new_center[1] - main_center[1])

        # Create new aperture and add the flux through it to the list
        temp_ap = EllipticalAperture(new_center, majorax, minorax, theta=tilt)
        temp_ap_mask = temp_ap.to_mask(method='center')
        temp_ap_mask = temp_ap_mask.multiply(image_data)
        flux_values.append(aperture_phot_manual(image_data, temp_ap_mask, new_center))

    # If for some reason there are nan values in the list, get rid of them
    i = 0
    while i < len(flux_values):
        if flux_values[i] != flux_values[i]:
            flux_values.pop(i)
            i -= 1
        i += 1

    return np.std(flux_values)  # Return the standard deviation of the flux through each aperture
