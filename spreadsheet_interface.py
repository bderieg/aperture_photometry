import openpyxl as opxl
import aperture_phot_functions as apf

# Directory names
main_path = "C:\\Users\\bderi\\Desktop\\Research Boizelle\\"
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

def update_spreadsheet(start=start_index, end=end_index):
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

        if "ALMA" in band:
            sheet[flux_col + str(i)] = round(apf.aperture_phot_ALMA(sci_name, unc_name, ap_name, background_x,
                                                                    background_y)[0], 4)
        elif("Wide" in band) or ("N" in band):
            sheet[flux_col + str(i)] = round(apf.aperture_phot_Akari(sci_name, unc_name, ap_name, background_x,
                                                                     background_y)[0], 4)
        elif "W" in band:
            if two_files == 'V':
                sheet[flux_col + str(i)] = round(apf.splice_WISE_vertical(sci_name, sci_name_2, ap_name, background_x,
                                                                          background_y, background_x_2, background_y_2,
                                                                          band), 4)
            elif two_files == 'H':
                sheet[flux_col + str(i)] = round(apf.splice_WISE_horizontal(sci_name, sci_name_2, ap_name, background_x,
                                                                          background_y, background_x_2, background_y_2,
                                                                          band), 4)
            else:
                sheet[flux_col + str(i)] = round(apf.aperture_phot_WISE(sci_name, unc_name, ap_name, background_x,
                                                                        background_y, band)[0], 4)
        elif "IRAC" or "MIPS" in band:
            sheet[flux_col + str(i)] = round(apf.aperture_phot_Spitzer(sci_name, unc_name, ap_name, background_x,
                                                                       background_y)[0], 4)

    # Update error value column
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

        if "ALMA" in band:
            sheet[err_col + str(i)] = apf.aperture_phot_ALMA(sci_name, unc_name, ap_name, background_x,
                                                             background_y)[1]
        elif ("Wide" in band) or ("N" in band):
            sheet[err_col + str(i)] = apf.aperture_phot_Akari(sci_name, unc_name, ap_name, background_x,
                                                              background_y)[1]
        elif "W" in band:
            sheet[err_col + str(i)] = apf.aperture_phot_WISE(sci_name, unc_name, ap_name, background_x,
                                                             background_y, band)[1]
        elif "IRAC" or "MIPS" in band:
            sheet[err_col + str(i)] = apf.aperture_phot_Spitzer(sci_name, unc_name, ap_name, background_x,
                                                                background_y)[1]

    wb.save(spreadsheet_loc)
