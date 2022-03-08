import openpyxl as opxl
import aperture_phot_functions as apf

# Directory names
main_path = "C:\\Users\\bderi\\Box\\School\\Research Boizelle\\"
aperture_folder = "Saved Apertures\\"
fits_folder = "FITS\\"
spreadsheet_loc = main_path + "Flux Measurements.xlsx"

# Constants
start_index = 4
end_index = 54
target_name_col = "A"
band_col = "B"
flux_col = "C"
err_col = "D"
ap_name_col = "E"
background_x_col = "F"
background_y_col = "G"
background_x_col_2 = "H"
background_y_col_2 = "I"
num_files_col = "J"

def get_parameters(t):
    # Open sheet, define some things
    wb = opxl.load_workbook(spreadsheet_loc)
    sheet = wb.active
    n = str(t)
    target_name_index = target_name_col + n
    ap_name_index = ap_name_col + n
    i = 1
    while sheet[target_name_index].value is None:
        target_name_index = target_name_col + str(t-i)
        i += 1
    file_num = ""
    two_files = 'None'
    if sheet[num_files_col + n].value == 'H':
        file_num = "_2"
        two_files = 'H'
    elif sheet[num_files_col + n].value == 'V':
        file_num = "_2"
        two_files = 'V'

    band = sheet[band_col + n].value
    sci_name = main_path + fits_folder + sheet[target_name_index].value + band + ".fits"
    sci_name_2 = main_path + fits_folder + sheet[target_name_index].value + band + file_num + ".fits"
    unc_name = main_path + fits_folder + sheet[target_name_index].value + band + "unc.fits"
    unc_name_2 = main_path + fits_folder + sheet[target_name_index].value + band + file_num + "unc.fits"
    ap_name = main_path + aperture_folder + sheet[ap_name_index].value
    background_x = sheet[background_x_col + n].value
    background_y = sheet[background_y_col + n].value
    background_x_2 = sheet[background_x_col_2 + n].value
    background_y_2 = sheet[background_y_col_2 + n].value

    return sci_name, unc_name, ap_name, background_x, background_y, background_x_2, background_y_2, band, sci_name_2, \
        unc_name_2, two_files

def update_spreadsheet(start=start_index, end=end_index, curve_of_growth=False, show_image=False):
    wb = opxl.load_workbook(spreadsheet_loc)
    sheet = wb.active

    # Update flux value column
    for i in range(start, end):
        sci_name = get_parameters(i)[0]
        unc_name = get_parameters(i)[1]
        ap_name = get_parameters(i)[2]
        background_x = get_parameters(i)[3]
        background_y = get_parameters(i)[4]
        background_x_2 = get_parameters(i)[5]
        background_y_2 = get_parameters(i)[6]
        band = get_parameters(i)[7]
        sci_name_2 = get_parameters(i)[8]
        unc_name_2 = get_parameters(i)[9]
        two_files = get_parameters(i)[10]

        if curve_of_growth:
            if "ALMA" in band:
                apf.curve_of_growth(sci_name, ap_name, background_x, background_y, band, background_sub=False)
            else:
                apf.curve_of_growth(sci_name, ap_name, background_x, background_y, band, background_sub=True)
        update_values = apf.image_aperture_phot(sci_name, ap_name, background_x, background_y, band,
                                                show_image=show_image)

        sheet[flux_col + str(i)] = round(update_values[0], 8)
        sheet[err_col + str(i)] = round(update_values[1], 8)

    wb.save(spreadsheet_loc)
