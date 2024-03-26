from astropy.io import fits as fits
import numpy as np
from matplotlib import pyplot as plt
import os
import subprocess
import streamlit as st
import json
from io import BytesIO
from streamlit_ace import st_ace
import sys
from code_editor import code_editor
import importlib
from datetime import date
from datetime import datetime
import shutil
        


def create_molec_input(wav,flux,error,sciencefile, minwav, maxwav, relatives, crop = True, directory = ""):
    """Creates files necessary to run molecfit
    
    Positional Arguments, in order:
    wav -- list, containing wavelength data
    flux -- list, containing flux data
    error -- list, containing errors in flux data
    sciencefile -- string; filename of original science data, for grabbing the header
    minwav -- float; contains minimum wavelength for region wanted for running molecfit by chip
    maxwav -- float; contains maximum wavelength for region wanted for running molecfit by chip
    relatives -- list; relative chemical abundances in order [O2, CO2, H2O, CO, CH4]    
    
    """
    numchips = len(wav)

    if len(directory) == 0:
        directory = "./MOLECFIT_Outputs"
    
    if len(sciencefile) > 0:
        hd = fits.open(sciencefile)
        primary_hdu = hd[0]
        hdul_output = fits.HDUList([primary_hdu])
        hdu_win = fits.HDUList([primary_hdu])
        hdu_mol = fits.HDUList([primary_hdu])
        hdu_atm = fits.HDUList([primary_hdu])
        hdu_con = fits.HDUList([primary_hdu])
        hdu_cor = fits.HDUList([primary_hdu])

    else: 
        hdul_output = fits.HDUList()
        hdu_win = fits.HDUList()
        hdu_mol = fits.HDUList()
        hdu_atm = fits.HDUList()
        hdu_con = fits.HDUList()
        hdu_cor = fits.HDUList()
    
    if crop:
        for i in range(len(wav)-1):
            if wav[i] < minwav and wav[i+1] > minwav:
                mindex = i
            if wav[i] < maxwav and wav[i+1] > maxwav:
                maxdex = i + 2
        
        cropwav = wav[mindex:maxdex]
        cropflux = flux[mindex:maxdex]
        croperr = error[mindex:maxdex]
            
        col1 = fits.Column(name='WAVE', format='D', array=cropwav) # create fits column for WAVE
        col2 = fits.Column(name='SPEC', format='D', array=cropflux) # create fits column for SPEC
        col3 = fits.Column(name='ERR', format='D', array=croperr) # create fits column for ERR

    else:
        col1 = fits.Column(name='WAVE', format='D', array=wav) # create fits column for WAVE
        col2 = fits.Column(name='SPEC', format='D', array=flux) # create fits column for SPEC
        col3 = fits.Column(name='ERR', format='D', array=error) # create fits column for ERR

    table_hdu = fits.BinTableHDU.from_columns([col1, col2, col3]) # create fits table HDU with WAVE SPEC ERR
    hdul_output.append(table_hdu) # append table HDU to output HDU list
    hdul_output.writeto(f"{directory}/SCIENCE.fits", overwrite=True)
    
    # Create FITS file with WAVE_INCLUDE
    if isinstance(minwav, float) and isinstance(maxwav, float):
        col_wmin = fits.Column(name="LOWER_LIMIT", format="D", array=[minwav])
        col_wmax = fits.Column(name="UPPER_LIMIT", format="D", array=[maxwav])
    else:
        col_wmin = fits.Column(name="LOWER_LIMIT", format="D", array=minwav)
        col_wmax = fits.Column(name="UPPER_LIMIT", format="D", array=maxwav)
    col_map = fits.Column(name="MAPPED_TO_CHIP", format="I", array=[1])
    col_wlc = fits.Column(name="WLC_FIT_FLAG", format="I", array=[1])
    col_cont = fits.Column(name="CONT_FIT_FLAG", format="I", array=[1])
    columns = [col_wmin, col_wmax, col_map, col_wlc, col_cont]
    #columns = [col_wmin, col_wmax]
    table_hdu = fits.BinTableHDU.from_columns(columns)
    hdu_win.append(table_hdu)
    hdu_win.writeto(f"{directory}/WAVE_INCLUDE.fits", overwrite=True)
    
    # Create FITS file with MOLECULES
    molecules = ["O2", "CO2", "H2O", "CO", "CH4", "N2O", "O3"]
    boolflags = [1,1,1,1,1,1,1]
    col_molecules = fits.Column(name="LIST_MOLEC", format="10A", array=molecules)
    col_boolflags = fits.Column(name="FIT_MOLEC", format="I", array=boolflags)
    col_relatives = fits.Column(name="REL_COL", format="D", array=relatives)
    columns = [col_molecules, col_boolflags, col_relatives]
    table_hdu = fits.BinTableHDU.from_columns(columns)
    hdu_mol.append(table_hdu)
    hdu_mol.writeto(f"{directory}/MOLECULES.fits", overwrite=True)
    
    # Create FITS file with MAPPING_ATMOSPHERIC
    name = "ATM_PARAMETERS_EXT"
    col_atm = fits.Column(name=name, format="K", array=[0,1])
    #col_atm = fits.Column(name=name, format="K", array=[0])
    table_hdu = fits.BinTableHDU.from_columns([col_atm])
    hdu_atm.append(table_hdu)
    fits_file = f"{directory}/MAPPING_ATMOSPHERIC.fits"
    hdu_atm.writeto(fits_file, overwrite=True)

    # Create FITS file with MAPPING_CONVOLVE
    name = "LBLRTM_RESULTS_EXT"
    col_conv = fits.Column(name=name, format="K", array=[0,1])
    #col_conv = fits.Column(name=name, format="K", array=[0])
    table_hdu = fits.BinTableHDU.from_columns([col_conv])
    hdu_con.append(table_hdu)
    fits_file = f"{directory}/MAPPING_CONVOLVE.fits"
    hdu_con.writeto(fits_file, overwrite=True)

    # Create FITS file with MAPPING_CORRECT
    name = "TELLURIC_CORR_EXT"
    col_corr = fits.Column(name=name, format="K", array=[0,1])
    #col_corr = fits.Column(name=name, format="K", array=[0])
    table_hdu = fits.BinTableHDU.from_columns([col_corr])
    hdu_cor.append(table_hdu)
    fits_file = f"{directory}/MAPPING_CORRECT.fits"
    hdu_cor.writeto(fits_file, overwrite=True)

def create_sofs(location):

    if not os.path.exists("./configs"):
        os.mkdir("configs")

    base = os.getcwd()    

    base = base + "/" + location + "/"

    sciencename = base + f"SCIENCE.fits SCIENCE \n"
    wavename = base + f"WAVE_INCLUDE.fits WAVE_INCLUDE \n"
    moleculename = base + f"MOLECULES.fits MOLECULES \n"
    mapcorrname = base + f"MAPPING_CORRECT.fits MAPPING_CORRECT \n"
    mapatmosname = base + f"MAPPING_ATMOSPHERIC.fits MAPPING_ATMOSPHERIC \n"
    mapconvname = base + f"MAPPING_CONVOLVE.fits MAPPING_CONVOLVE \n"
    
    atmparamname = base + f"ATM_PARAMETERS.fits ATM_PARAMETERS \n"
    modelmolname = base + f"MODEL_MOLECULES.fits MODEL_MOLECULES \n"
    bestfitname = base + f"BEST_FIT_PARAMETERS.fits BEST_FIT_PARAMETERS \n"
    tellurname = base + f"TELLURIC_CORR.fits TELLURIC_CORR \n"
    
    
    modelsof = open(f"configs/model.sof", "w")
    modelsof.write(sciencename)
    modelsof.write(wavename)
    modelsof.write(moleculename)
    modelsof.close()
    
    calctransof = open(f"configs/calctrans.sof", "w")
    calctransof.write(sciencename)
    calctransof.write(mapatmosname)
    calctransof.write(mapconvname)
    calctransof.write(atmparamname)
    calctransof.write(modelmolname)
    calctransof.write(bestfitname)
    calctransof.close()
    
    correctsof = open(f"configs/correct.sof","w")
    correctsof.write(sciencename)
    correctsof.write(mapcorrname)
    correctsof.write(tellurname)
    correctsof.close()
    


def save_corrected_data(wavs, fluxes, errors, sciencefile=""):
    """Saves corrected data as a fits file, segmented by chip
    
    Positional Arguments, in order:
    wavs -- list of lists; wavelength data, segmented by chip
    fluxes -- list of lists; flux data, segmented by chip
    errors -- list of lists; errors in flux data, segmented by chip
    sciencefile -- string; name of original data file, for grabbing initial header
    
    Keyword Arguments:
    onlytrim -- boolean; True: adds _TRIMMED to filename
                         False (default): adds _CORRECTED to filename
    """
    
    hd = fits.open(sciencefile)
    primary_hdu = hd[0]

    hdul_output = fits.HDUList([primary_hdu])
    if len(sciencefile) == 0:
        filename = "corrected_data.fits"
    
    filename = sciencefile[:-5] + "_CORRECTED.fits"
    
    for i,j,k in zip(wavs,fluxes,errors):
        col1 = fits.Column(name='WAVE', format='D', array=i) # create fits column for WAVE
        col2 = fits.Column(name='SPEC', format='D', array=j) # create fits column for SPEC
        col3 = fits.Column(name='ERR', format='D', array=k) # create fits column for ERR
        table_hdu = fits.BinTableHDU.from_columns([col1, col2, col3]) # create fits table HDU with WAVE SPEC ERR
        hdul_output.append(table_hdu) # append table HDU to output HDU list
    hdul_output.writeto(filename, overwrite=True)
    #input(f"Successfully Saved Data to {filename}, input anything to continue: ")


def load_css(filename):
    with open(filename) as f:
        st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)


def makeplot(x, y, telx, tely, boxranges=[]):
    fig, ax1 = plt.subplots(figsize=[10, 5.1])
    
    #fig.set_size_inches(240.5, 140.5)
    #fig.set_size_inches(120.25, 70.25)
    #fig, ax1 = plt.subplots(figsize=[40, 20.4])
    ax1.plot(x, y)
    
    #ax1.set_xlabel(f"Wavelength ({NAME_SHORT[NAME_CONVERSION[st.session_state.wavunit]]})",fontsize=24)
    ax1.set_xlabel(f"Wavelength ({NAME_SHORT[NAME_CONVERSION[st.session_state.wavunit]]})",fontsize=12)
    #ax1.set_ylabel("Flux", color= 'tab:blue', fontsize=24)
    ax1.set_ylabel("Flux", color= 'tab:blue', fontsize=12)
    color = 'tab:blue'
    color2 = 'tab:orange'
    #ax1.tick_params(axis='y', labelcolor=color, labelsize=18)
    ax1.tick_params(axis='y', labelcolor=color, labelsize=12)
    #ax1.tick_params(axis='x', labelsize=18)
    ax1.tick_params(axis='x', labelsize=12)
    #ax1.set_xlim([x[0],x[-1]])
    ax1.set_xlim([st.session_state.dispmin, st.session_state.dispmax])
    chunk = [y[z] for z in range(len(x)) if (x[z] > st.session_state.dispmin and x[z] < st.session_state.dispmax)]
    ax1.set_ylim(min(chunk)*0.98, max(chunk)*1.02)
    if len(st.session_state.telluricwavelengths) > 0 and len(st.session_state.tellurictransmission) > 0:
        ax2 = ax1.twinx()
        #ax2.set_ylabel('Transmission', color=color2, fontsize=24)
        #ax2.tick_params(axis='y', labelcolor=color2, labelsize=18)
        ax2.set_ylabel('Transmission', color=color2, fontsize=12)
        ax2.tick_params(axis='y', labelcolor=color2, labelsize=12)
        ax2.plot(telx,tely,color=color2)
        #ax2.set_xlim([x[0],x[-1]])
        ax2.set_ylim([0,1])
    #if len(boxranges) > 0:
        #try:
         #   for i in boxranges:
          #      ax1.axvspan(boxranges[i][0], boxranges[i][1], color='green', alpha=0.2)
        #except:
         #   pass
    for i in range(st.session_state.winnum):
        try:
            ax1.axvspan(float(st.session_state.wavmin[i]), float(st.session_state.wavmax[i]), color='green', alpha=0.2)
        except:
            pass

    return fig


def update_vars(code):
    _locals = locals()
    wavelengths = []
    fluxes = []
    errors = []
    telluricwavelengths = []
    tellurictransmission = []
    sciencefile = ""
    wavunit = ""
    ldict = {}
    if os.path.exists("userscript.py"):
        os.remove("userscript.py")
    with open("userscript.py", "w") as scriptfile:
        newfunction = "def usercode():\n\twavelengths = []\n\tfluxes = []\n\terrors = []\n\ttelluricwavelengths = []\n\ttellurictransmission = []\n\tsciencefile = \"\"\n\twavunit = \"micron\"\n"
        scriptfile.write(newfunction)
        for i in code.split("\n"):
            newline = '\t' + i + '\n'
            scriptfile.write(newline)
        scriptfile.write("\treturn wavelengths, fluxes, errors, telluricwavelengths, tellurictransmission, sciencefile, wavunit")

    if 'usercode' not in sys.modules:
        import userscript as uscr
    importlib.reload(sys.modules["userscript"])
    wavelengths, fluxes, errors, telluricwavelengths, tellurictransmission, sciencefile, wavunit = uscr.usercode()

    if os.path.exists("userscript.py"):
        os.remove("userscript.py")

    return wavelengths, fluxes, errors, telluricwavelengths, tellurictransmission, sciencefile, wavunit


def final_button(save = True):
    if save:
        st.session_state.fluxes[st.session_state.chipnum] = st.session_state.tempfluxes
        st.session_state.errors[st.session_state.chipnum] = st.session_state.temperrors
        st.session_state.disp_save = True
        st.session_state.save_directory = st.session_state.directory_list[st.session_state.appliedit]
        save_corrected_data(st.session_state.wavelengths, st.session_state.fluxes, st.session_state.errors, sciencefile = st.session_state.sciencefile)
        save_log(st.session_state.sciencefile, st.session_state.wavmin, st.session_state.wavmax, st.session_state.loglines, st.session_state.correctselect)
    for i in st.session_state.directory_list:
        if save and i == st.session_state.directory_list[st.session_state.appliedit]:
            with open(st.session_state.directory_list[st.session_state.appliedit] + "/WAVELENGTHS.txt", 'w') as f:
                f.write("WAVELENGTH RANGES FOR THIS MOLECFIT RUN:\n")
                for i,j in zip(st.session_state.wavmin, st.session_state.wavmax):
                    f.write(f"[{i}, {j}]\n")
        else:
            shutil.rmtree(i)
    st.session_state.itmode = False
    st.session_state.selectmode = True
    st.session_state.wavmax = [""]
    st.session_state.wavmin = [""]
    st.session_state.winnum = 1
    st.session_state.tried_abundances = []
    st.session_state.chis = []
   # st.session_state.dip_ratios = []
    st.session_state.mingood = False
    st.session_state.maxgood = False
    st.session_state.pastflux = []
    st.session_state.pasterr = []
    st.session_state.appliedit = 0
    st.session_state.loglines = []
    st.session_state.directory_list = []
    st.session_state.okaymax = [0]
    st.session_state.okaymin = [0]
    st.session_state.resetmax = True
    st.session_state.resetmin = True
    st.rerun()

def streamlit_run_molecfit_together():
    """Primary function for running all of molecfit
    
    Positional Arguments, in order:
    wavs -- list of lists; wavelength data by chip
    fluxes -- list of lists; flux data by chip
    errors -- list of lists; errors in flux data by chip
    relatives -- list; relative abundances in order [O2, CO2, H2O, CO, CH4]
    min -- float; wavelength minimum
    max -- float; wavelength maximum
    """
    
    #Starting Abundances for O2, CO2, H2O, CO, CH4 respectively
    #Derived from: http://eodg.atm.ox.ac.uk/RFM/atm/ngt.atm
    #relatives = [9.98E-1, 1.38E-3, 5.07E-4, 7.74E-05, 2.40E-06]
    #relatives = [0.9, 0.05, 0.025, 0.015, 0.01]
    #relatives = [1, 1, 1, 1, 1]

    chip = st.session_state.chipnum
    wav = [z * CONVERSION[st.session_state.wavunit] for z in st.session_state.wavelengths[chip]]
    flux = st.session_state.fluxes[chip]
    error = st.session_state.errors[chip]
    smallest_flux = "nada"
    
    wavmin = [float(z) * CONVERSION[st.session_state.wavunit] for z in st.session_state.wavmin]
    wavmax = [float(z) * CONVERSION[st.session_state.wavunit] for z in st.session_state.wavmax]

    st.session_state.tried_abundances.append(st.session_state.abundances.copy())
    print("Running MOLECFIT_Model...")
    tday = str(date.today())
    ahora = str(datetime.now())[11:19].replace(":", ".")
    output_directory = f"MOLECFIT_Outputs/{tday}/{ahora}"
    if not os.path.exists("MOLECFIT_Outputs"):
        os.mkdir("MOLECFIT_Outputs")
    if not os.path.exists(f"MOLECFIT_Outputs/{tday}"):
        os.mkdir(f"MOLECFIT_Outputs/{tday}")
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    st.session_state.directory_list.append(output_directory)

    create_sofs(output_directory)
    create_molec_input(wav,flux,error,st.session_state.sciencefile,wavmin,wavmax,st.session_state.abundances, crop=False, directory = output_directory)


    with open(os.devnull, 'wb') as devnull:
        if not st.session_state.verbose:
            subprocess.check_call(['esorex', '--recipe-config=configs/molecfit_model.rc', 
                                   f'--output-dir={output_directory}', 'molecfit_model', 'configs/model.sof'], 
                                   stdout=devnull, stderr=subprocess.STDOUT)
        else:
            subprocess.check_call(['esorex', '--recipe-config=configs/molecfit_model.rc', 
                                   f'--output-dir={output_directory}', 'molecfit_model', 'configs/model.sof'], 
                                   stderr=subprocess.STDOUT)
        #try:
            #subprocess.check_output(...)
        #except subprocess.CalledProcessError as e:
            #print(e.output)
        
    print("Running Molecfit_Calctrans...")
    with open(os.devnull, 'wb') as devnull:
        if not st.session_state.verbose:
            subprocess.check_call(['esorex', '--recipe-config=configs/molecfit_calctrans.rc', 
                                   f'--output-dir={output_directory}', 'molecfit_calctrans', 'configs/calctrans.sof'], 
                                   stdout=devnull, stderr=subprocess.STDOUT)
        else:
            subprocess.check_call(['esorex', '--recipe-config=configs/configs/molecfit_calctrans.rc', 
                                   f'--output-dir={output_directory}', 'molecfit_calctrans', 'configs/calctrans.sof'], 
                                   stderr=subprocess.STDOUT)
        
    print("Running Molecfit_Correct...")
    with open(os.devnull, 'wb') as devnull:
        if not st.session_state.verbose:
            subprocess.check_call(['esorex', '--recipe-config=configs/molecfit_correct.rc', 
                                   f'--output-dir={output_directory}', 'molecfit_correct', 'configs/correct.sof'], 
                                   stdout=devnull, stderr=subprocess.STDOUT)
        else:
            subprocess.check_call(['esorex', '--recipe-config=configs/molecfit_correct.rc', 
                                   f'--output-dir={output_directory}', 'molecfit_correct', 'configs/correct.sof'], 
                                    stderr=subprocess.STDOUT)
        
    new_fits = fits.open(f"{output_directory}/SCIENCE_TELLURIC_CORR_SCIENCE.fits")
    chi_file = fits.open(f"{output_directory}/BEST_FIT_PARAMETERS.fits")
    chi = [i[1] for i in chi_file[1].data if i[0] == 'reduced_chi2']
    st.session_state.chis.append(f"{chi[0]:.3f}")

    for i in range(1,len(new_fits)):
        new_data = new_fits[i].data
        corrw = []
        corrf = []
        corre = []
        for ii in new_data:
            corrw.append(ii[0])
            corrf.append(ii[1])
            corre.append(ii[2])
    combined_flux = []
    combined_error = []
    corrit = 0
    for it in range(len(wav)):
        if wav[it] == corrw[corrit]:
            combined_flux.append(corrf[corrit])
            combined_error.append(corre[corrit])
            if not (corrit == len(corrw) - 1):
                corrit += 1
        else:
            combined_error.append(error[it])
            combined_flux.append(flux[it])
            
    #dip_val = np.mean(combined_flux) - combined_flux[smallest_flux_index]
    #dip_rat = 1.0 - abs(dip_val) / base_dip
    #st.session_state.dip_ratios.append(f"{dip_rat:.3f}")

    st.session_state.tempfluxes = combined_flux
    st.session_state.temperrors = combined_error
    st.session_state.pastflux.append(st.session_state.tempfluxes.copy())
    st.session_state.pasterr.append(st.session_state.temperrors.copy())

    st.session_state.appliedit = len(st.session_state.pasterr) - 1

    it = 0

    if not os.path.exists("logs"):
        os.mkdir("logs")
    if not os.path.exists(f"logs/{tday}"):
        os.mkdir(f"logs/{tday}")
    
    shutil.copy("esorex.log", f"logs/{tday}/esorex{ahora}.log")



def streamlit_run_molecfit_seperate():
    """Primary function for running all of molecfit
    
    Positional Arguments, in order:
    wavs -- list of lists; wavelength data by chip
    fluxes -- list of lists; flux data by chip
    errors -- list of lists; errors in flux data by chip
    relatives -- list; relative abundances in order [O2, CO2, H2O, CO, CH4, N2O, O3]
    min -- float; wavelength minimum
    max -- float; wavelength maximum
    """
    
    #Starting Abundances for O2, CO2, H2O, CO, CH4 respectively
    #Derived from: http://eodg.atm.ox.ac.uk/RFM/atm/ngt.atm
    #relatives = [9.98E-1, 1.38E-3, 5.07E-4, 7.74E-05, 2.40E-06]
    #relatives = [0.9, 0.05, 0.025, 0.015, 0.01]
    #relatives = [1, 1, 1, 1, 1]

    chip = st.session_state.chipnum
    wav = [z * CONVERSION[st.session_state.wavunit] for z in st.session_state.wavelengths[chip]]
    flux = st.session_state.fluxes[chip]
    error = st.session_state.errors[chip]

    st.session_state.tempfluxes = flux
    st.session_state.temperrors = error

    
    tempchi = []
    tempdip = []

    tday = str(date.today())
    ahora = str(datetime.now())[11:19].replace(":", ".")
    output_directory = f"MOLECFIT_Outputs/{tday}/{ahora}"
    if not os.path.exists("MOLECFIT_Outputs"):
        os.mkdir("MOLECFIT_Outputs")
    if not os.path.exists(f"MOLECFIT_Outputs/{tday}"):
        os.mkdir(f"MOLECFIT_Outputs/{tday}")
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    st.session_state.directory_list.append(output_directory)

    create_sofs(output_directory)
    
    st.session_state.tried_abundances.append(st.session_state.abundances.copy())
    for it,(mini,maxi) in enumerate(zip(st.session_state.wavmin, st.session_state.wavmax)):
        mini = float(mini) * CONVERSION[st.session_state.wavunit]
        maxi = float(maxi) * CONVERSION[st.session_state.wavunit]
        create_molec_input(wav,flux,error,st.session_state.sciencefile, mini, maxi,st.session_state.abundances, directory=output_directory)
        winflux = [flux[i] for i in range(len(flux)) if (wav[i] > mini and wav[i] < maxi)]

        print(f"Running MOLECFIT_Model for Window {it+1}...")
        with open(os.devnull, 'wb') as devnull:
            if not st.session_state.verbose:
                subprocess.check_call(['esorex', '--recipe-config=configs/molecfit_model.rc', 
                                       f'--output-dir={output_directory}', 'molecfit_model', 'configs/model.sof'], 
                                       stdout=devnull, stderr=subprocess.STDOUT)
            else:
                subprocess.check_call(['esorex', '--recipe-config=configs/molecfit_model.rc', 
                                       f'--output-dir={output_directory}', 'molecfit_model', 'configs/model.sof'], 
                                       stderr=subprocess.STDOUT)
            #try:
                #subprocess.check_output(...)
            #except subprocess.CalledProcessError as e:
                #print(e.output)
            
        print(f"Running Molecfit_Calctrans for Window {it+1}...")
        with open(os.devnull, 'wb') as devnull:
            if not st.session_state.verbose:
                subprocess.check_call(['esorex', '--recipe-config=configs/molecfit_calctrans.rc', 
                                       f'--output-dir={output_directory}', 'molecfit_calctrans', 'configs/calctrans.sof'], 
                                       stdout=devnull, stderr=subprocess.STDOUT)
            else:
                subprocess.check_call(['esorex', '--recipe-config=configs/molecfit_calctrans.rc', 
                                       f'--output-dir={output_directory}', 'molecfit_calctrans', 'configs/calctrans.sof'], 
                                       stderr=subprocess.STDOUT)
            
        print(f"Running Molecfit_Correct for Window {it+1}...")
        with open(os.devnull, 'wb') as devnull:
            if not st.session_state.verbose:
                subprocess.check_call(['esorex', '--recipe-config=configs/molecfit_correct.rc', 
                                       f'--output-dir={output_directory}', 'molecfit_correct', 'configs/correct.sof'], 
                                       stdout=devnull, stderr=subprocess.STDOUT)
            else:
                subprocess.check_call(['esorex', '--recipe-config=configs/molecfit_correct.rc', 
                                       f'--output-dir={output_directory}', 'molecfit_correct', 'configs/correct.sof'], 
                                       stderr=subprocess.STDOUT)
            
        new_fits = fits.open(f"{output_directory}/SCIENCE_TELLURIC_CORR_SCIENCE.fits")
        chi_file = fits.open(f"{output_directory}/BEST_FIT_PARAMETERS.fits")
        chi = [i[1] for i in chi_file[1].data if i[0] == 'reduced_chi2']
        tempchi.append(f"{chi[0]:.3f}")
        
        for i in range(1,len(new_fits)):
            new_data = new_fits[i].data
            corrw = []
            corrf = []
            corre = []
            for ii in new_data:
                corrw.append(ii[0])
                corrf.append(ii[1])
                corre.append(ii[2])
        combined_flux = []
        combined_error = []
        corrit = 0
        for jit in range(len(wav)):
            if wav[jit] == corrw[corrit]:
                combined_flux.append(corrf[corrit])
                combined_error.append(corre[corrit])
                if not (corrit == len(corrw) - 1):
                    corrit += 1
            else:
                combined_error.append(st.session_state.temperrors[jit])
                combined_flux.append(st.session_state.tempfluxes[jit])
        
        folderit = 0
        while(os.path.exists(f"logs/{tday}/esorex{folderit}")):
            folderit += 1
        if not os.path.exists("logs"):
            os.mkdir("logs")
        if not os.path.exists(f"logs/{tday}/esorex{ahora}"):
            os.mkdir(f"logs/{tday}/esorex{ahora}")
        shutil.copy("esorex.log", f"logs/{tday}/esorex{ahora}/esorex{it}.log")


    
        st.session_state.tempfluxes = combined_flux
        st.session_state.temperrors = combined_error

    st.session_state.pastflux.append(st.session_state.tempfluxes.copy())
    st.session_state.pasterr.append(st.session_state.temperrors.copy())
    st.session_state.appliedit = len(st.session_state.pasterr) - 1
    #st.session_state.dip_ratios.append(tempdip)
    st.session_state.chis.append(tempchi)


def save_log(sciencefile, mins, maxes, lines, only_correct_segment):
    shortname = ""
    if "/" in sciencefile or "\\" in sciencefile:
        noext = True
        for i in sciencefile[::-1]:
            if noext:
                if i == '.':
                    noext = False 
                else:
                    pass
            elif ~noext and i != "/" and i != "\\":
                shortname = shortname + i
            elif ~noext:
                break
        shortname = shortname[::-1]
    else:
        shortname = sciencefile
    if not os.path.exists("./logs"):
        os.mkdir("logs")
    tday = str(date.today())
    if not os.path.exists(f"logs/{tday}"):
        os.mkdir(f"logs/{tday}")
    fulllogpath = f"logs/{tday}/{shortname}.txt"
    if not os.path.exists(fulllogpath):
        with open(fulllogpath, 'w') as f:
            f.write(f"Correction Log for {shortname} on {tday}\n")
            f.write("------------------------------------------------------------------------------------------\n")
    with open(fulllogpath, 'a') as f:
        if isinstance(mins, list) or isinstance(mins, np.ndarray):
            wavstring = "Wavelength Ranges: "
            for i,j in zip(mins, maxes):
                if i == mins[-1]:
                    wavstring = wavstring + f"[{i}, {j}]\n"
                else:
                    wavstring = wavstring + f"[{i}, {j}], "
        else:
            wavstring = f"Wavelength Range: [{mins}, {maxes}]\n"
        f.write(wavstring)
        if only_correct_segment:
            f.write("Corrected Only Selected Region\n")
        else:
            f.write("Corrected Entire Chip\n")
        for it,i in enumerate(lines):
            if it == st.session_state.appliedit:
                f.write(i[:-2] + " | APPLIED\n")
            else:
                f.write(i)
        f.write("------------------------------------------------------------------------------------------\n")

    os.remove("esorex.log")



if 'wavelengths' not in st.session_state:
    st.session_state.wavelengths = []
if 'fluxes' not in st.session_state:
    st.session_state.fluxes = []
if 'errors' not in st.session_state:
    st.session_state.errors = []
if 'telluricwavelengths' not in st.session_state:
    st.session_state.telluricwavelengths = []
if 'tellurictransmission' not in st.session_state:
    st.session_state.tellurictransmission = []
if 'sciencefile' not in st.session_state:
    st.session_state.sciencefile = []
if 'selectmode' not in st.session_state:
    st.session_state.selectmode = True
if 'chipnum' not in st.session_state:
    st.session_state.chipnum = None
if 'currentplot' not in st.session_state:
    st.session_state.currentplot = None
if 'disablechip' not in st.session_state:
    st.session_state.disablechip = False
if 'wavmin' not in st.session_state:
    st.session_state.wavmin = [""]
if 'wavmax' not in st.session_state:
    st.session_state.wavmax = [""]
if 'abundances' not in st.session_state:
    st.session_state.abundances = [0.9, 0.05, 0.025, 0.015, 0.01, 0.01, 0.01]
if 'mingood' not in st.session_state:
    st.session_state.mingood = False
if 'maxgood' not in st.session_state:
    st.session_state.maxgood = False
if 'overlap' not in st.session_state:
    st.session_state.overlap = False
if 'itmode' not in st.session_state:
    st.session_state.itmode = False
if 'tempfluxes' not in st.session_state:
    st.session_state.tempfluxes = []
if 'temperrors' not in st.session_state:
    st.session_state.temperrors = []
if 'tried_abundances' not in st.session_state:
    st.session_state.tried_abundances = []
if 'chis' not in st.session_state:
    st.session_state.chis = []
#if 'dip_ratios' not in st.session_state:
    #st.session_state.dip_ratios = []
if 'disp_save' not in st.session_state:
    st.session_state.disp_save = False
if 'first_load' not in st.session_state:
    st.session_state.first_load = True
if 'code_input' not in st.session_state:
    st.session_state.code_input = ""
if 'wavmsg' not in st.session_state:
    st.session_state.wavmsg = ":red[wavelengths not set.]"
if 'fluxmsg' not in st.session_state:
    st.session_state.fluxmsg = ":red[fluxes not set.]"
if 'errmsg' not in st.session_state:
    st.session_state.errmsg = ":red[errors not set.]"
if 'telwavmsg' not in st.session_state:
    st.session_state.telwavmsg = ":orange[telluricwavelengths not set.]"
if 'teltransmsg' not in st.session_state:
    st.session_state.teltransmsg = ":orange[tellurictransmission not set.]"
if 'scimsg' not in st.session_state:
    st.session_state.scimsg = f":red[sciencefile not set.]"
if 'unitmsg' not in st.session_state:
    st.session_state.unitmsg = f":red[wavunit not set.]"
if 'showtel' not in st.session_state:
    st.session_state.showtel = False
if 'wavunit' not in st.session_state:
    st.session_state.wavunit = ""
if 'showtelselect' not in st.session_state:
    st.session_state.showtelselect = True
if 'dispmin' not in st.session_state:
    st.session_state.dispmin = 0
if 'dispmax' not in st.session_state:
    st.session_state.dispmax = 0
if 'correctselect' not in st.session_state:
    st.session_state.correctselect = False
if 'winnum' not in st.session_state:
    st.session_state.winnum = 1
if 'pastflux' not in st.session_state:
    st.session_state.pastflux = []
if 'pasterr' not in st.session_state:
    st.session_state.pasterr = []
if 'appliedit' not in st.session_state:
    st.session_state.appliedit = 0
if 'tagnums' not in st.session_state:
    st.session_state.tagnums = [0]
if 'skipit' not in st.session_state:
    st.session_state.skipit = False
if 'loglines' not in st.session_state:
    st.session_state.loglines = []
if 'directory_list' not in st.session_state:
    st.session_state.directory_list = []
if 'save_directory' not in st.session_state:
    st.session_state.save_directory = ""
if 'okaymax' not in st.session_state:
    st.session_state.okaymax = [0]
if 'okaymin' not in st.session_state:
    st.session_state.okaymin = [0]
if 'resetmax' not in st.session_state:
    st.session_state.resetmax = True
if 'resetmin' not in st.session_state:
    st.session_state.resetmin = True
if 'abunderrs' not in st.session_state:
    st.session_state.abunderrs = [0, 0, 0, 0, 0, 0, 0]
if 'verbose' not in st.session_state:
    st.session_state.verbose = False


CONVERSION = {'micron':1.0, 'microns':1.0, 'mum':1.0, 'µm':1.0, 'm':1E6, 'meter':1E6, 'meters':1E6, 'cm':1E4, 'centimeter':1E4, 'centimeters':1E4, 'mm':1E3, 'millimeter':1E3, 'millimeters':1E3, 'nm':1E-3, 'nanometer':1E-3, 'nanometers':1E-3, 'angstrom':1E-4, 'angstroms':1E-4, 'Å':1E-4}
NAME_CONVERSION = {
                   'micron':'micron', 'microns':'micron', 'mum':'micron', 'µm':'micron', 
                   'm':'meter', 'meter':'meter', 'meters':'meter', 
                   'cm':'centimeter', 'centimeters':'centimeter', 'centimeter':'centimeter', 
                   'mm':'millimeter', 'millimeters':'millimeter', 'millimeter':'millimeter', 
                   'nm':'nanometer', 'nanometer':'nanometer', 'nanometers':'nanometer', 
                   'angstrom':'angstrom', 'angstroms':'angstrom', 'Å':'angstrom'
                   }
NAME_SHORT = {'micron':'µm', 'meter':'m', 'centimeter':'cm', 'millimeter':'mm', 'angstrom':'Å', 'nanometer':'nm'}
SHOW_LABEL = {True: "visible", False: "hidden"}
ACTIVE_BUTTON_LABEL = {True: "Current", False: "Reapply"}
ACTIVE_BUTTON_COLOR = {True: "primary", False: "secondary"}

def streamlit_molecfit():

    st.set_page_config(layout="wide")
    load_css("css/streamlit.css")

    inputtab, molectab = st.tabs(["Data Input", "MOLECFIT"])
    savejson = {}
    goodtogo = True
    if 'code_default' not in st.session_state:
        try:
            with open("./save.json") as f:
                savejson = json.load(f)
                st.session_state.code_default = savejson["code_default"]
        except:
            st.session_state.code_default = ""

    with inputtab:
        st.header("Input Python Code Defining Data")
        st.text("You must define variables for correction, as a single series or as lists of lists segmented however is desired for correction (both lists and numpy arrays work):\n'wavelengths' \n'fluxes' \n'errors'")
        st.markdown("You must also define the variables: \n\'wavunit\' (defines the wavelength unit the data is in, ex: \'micron\', \'Å\', \'angstrom\', \'nm\', etc)\n\'sciencefile\' (original data fits file, for observation data in header)  \n:grey[TIP: to continue work from a previous session, load in data from {sciencefile}_CORRECTED.fits, setting sciencefile to the real original]")
        st.text("You can also define optional variables: \n\'telluricwavelengths' (Wavelengths for telluric visualization, segmented/unsegmented makes no difference) \n\'tellurictransmission' (Transmission values for telluric visualization, segmented/unsegmented makes no difference)")
        col1, col2 = st.columns([2, 1])
        #col1, col2 = st.columns(2)
        with col1:
            #form = st.form(key="user_form")
            #code_input = form.text_area("Python Code Defining Variables (saves on variable update):", value=st.session_state.code_default, height=300, placeholder="Input Python Here")
            if st.session_state.first_load:
                st.session_state.code_input = st.session_state.code_default
                st.session_state.first_load = False
            st.session_state.code_input = st_ace(value = st.session_state.code_input, language = 'python', height = 400)
            #print(st.session_state.code_input)
            #text_code_input = code_editor(st.session_state.code_default, lang="python", height=300)
            #code_input = text_code_input['text']
            #print(text_code_input)
            #if form.form_submit_button("Update Variables"):
            #if code_input != st.session_state.code_default:
                #st.session_state.code_default = code_input
            

        with col2:
            st.text("")
            #st.text("")
            #st.text("")
            if st.button("Update Variables", type="primary"):
                savejson["code_default"] = st.session_state.code_input
                with open("save.json", "w") as outfile:
                    json.dump(savejson, outfile)
                st.session_state.wavelengths, st.session_state.fluxes, st.session_state.errors, st.session_state.telluricwavelengths, st.session_state.tellurictransmission, st.session_state.sciencefile, st.session_state.wavunit = update_vars(st.session_state.code_input)

                if not isinstance(st.session_state.wavelengths, (list, np.ndarray)):
                    st.session_state.fluxes = f":red[fluxes not set, non-list type {type(st.session_state.wavelengths)} detected.]"
                if len(st.session_state.wavelengths) > 0:
                    if all([isinstance(x, (list, np.ndarray)) for x in st.session_state.wavelengths]) and all([isinstance(xx, (int, float)) for x in st.session_state.wavelengths for xx in x]):
                        st.session_state.wavmsg = f":green[wavelengths set, Segment 0: \[{st.session_state.wavelengths[0][0]:.3f} ... {st.session_state.wavelengths[0][-1]:.3f}\]]"
                    elif all(isinstance(x, (int, float)) for x in st.session_state.wavelengths):
                        temporary = st.session_state.wavelengths.copy()
                        st.session_state.wavelengths = []
                        st.session_state.wavelengths.append(temporary)
                        st.session_state.wavmsg = f":green[wavelengths set, single series: \[{st.session_state.wavelengths[0][0]:.3f} ... {st.session_state.wavelengths[0][-1]:.3f}\]]"
                    else:
                        #for i in st.session_state.wavelengths:
                           #print(type(i))
                        st.session_state.wavmsg = f":red[wavelengths not set: unsupported format loaded]"
                        goodtogo = False
                else:
                    st.session_state.wavmsg = ":red[wavelengths not set.]"
                    goodtogo = False
    
                if not isinstance(st.session_state.fluxes, (list,np.ndarray)):
                    st.session_state.fluxes = f":red[fluxes not set, non-list type {type(st.session_state.fluxes)} detected.]"
                elif len(st.session_state.fluxes) > 0:
                    if all(isinstance(x, (list,np.ndarray)) for x in st.session_state.fluxes) and all([isinstance(xx, (int, float)) for x in st.session_state.fluxes for xx in x]):
                        if all((len(x) == len(y) for x,y in zip(st.session_state.wavelengths,st.session_state.fluxes))):
                            st.session_state.fluxmsg = f":green[fluxes set, Segment 0: \[{st.session_state.fluxes[0][0]:.3f} ... {st.session_state.fluxes[0][-1]:.3f}\]]"
                        else:
                            st.session_state.fluxmsg = f":red[fluxes not set, not all segment lengths match wavelength segment lengths]"
                            goodtogo = False
                    elif all(isinstance(x, (int, float)) for x in st.session_state.fluxes):
                        temporary = st.session_state.fluxes.copy()
                        st.session_state.fluxes = []
                        st.session_state.fluxes.append(temporary)
                        if len(st.session_state.fluxes[0]) == len(st.session_state.wavelengths[0]):
                            #print("green light")
                            #print(st.session_state.fluxes[0], st.session_state.wavelengths[0])
                            st.session_state.fluxmsg = f":green[fluxes set, single series: \[{st.session_state.fluxes[0][0]:.3f} ... {st.session_state.fluxes[0][-1]:.3f}\]]"
                            #print(st.session_state.fluxmsg)
                        else:
                            #print("red light")
                            #print(st.session_state.fluxes[0], st.session_state.wavelengths[0])
                            st.session_state.fluxmsg = f":red[fluxes not set, single input series length doesn't match wavelength series length]"
                            #print(st.session_state.fluxmsg)
                            goodtogo = False
                    else:
                        st.session_state.fluxmsg = f":red[fluxes not set: unsupported format loaded]"
                        goodtogo = False
                else:
                    st.session_state.fluxmsg = ":red[fluxes not set.]"
                    goodtogo = False

                if not isinstance(st.session_state.errors, (list, np.ndarray)):
                    st.session_state.errmsg = f":red[telluricwavelengths not set, non-list type {type(st.session_state.errors)} detected.]"
                elif len(st.session_state.errors) > 0:
                    if all(isinstance(x, (list,np.ndarray)) for x in st.session_state.errors) and all([isinstance(xx, (int, float)) for x in st.session_state.errors for xx in x]):
                        if all((len(x) == len(y) for x,y in zip(st.session_state.wavelengths,st.session_state.errors))):
                            st.session_state.errmsg = f":green[errors set, Segment 0: \[{st.session_state.errors[0][0]:.3f} ... {st.session_state.errors[0][-1]:.3f}\]]"
                        else:
                            st.session_state.errmsg = f":red[errors not set, not all segment lengths match wavelength segment lengths]"
                            goodtogo = False
                    elif all(isinstance(x, (int, float)) for x in st.session_state.errors):
                        temporary = st.session_state.errors.copy()
                        st.session_state.errors = []
                        st.session_state.errors.append(temporary)
                        if len(st.session_state.errors[0]) == len(st.session_state.wavelengths[0]):
                            st.session_state.errmsg = f":green[errors set, single series: \[{st.session_state.errors[0][0]:.3f} ... {st.session_state.errors[0][-1]:.3f}\]]"
                        else:
                            st.session_state.errmsg = f":red[errors not set, single input series length doesn't match wavelength series length]"
                            goodtogo = False
                    else:
                        st.session_state.errmsg = f":red[errors not set: unsupported format loaded]"
                        goodtogo = False
                else:
                    st.session_state.errmsg = ":red[errors not set.]"
                    goodtogo = False

                if not isinstance(st.session_state.wavunit, str):
                    st.session_state.unitmsg = f":red[wavunit not set, non-string type {type(st.session_state.wavunit)} detected.]"
                    st.session_state.wavunit = ""
                elif len(st.session_state.wavunit) > 0:
                    if st.session_state.wavunit not in CONVERSION.keys():
                        st.session_state.unitmsg = f":red[wavunit not set, given name not in unit dictionary]"
                    else:
                        st.session_state.unitmsg = f":green[wavunit set: \"{NAME_CONVERSION[st.session_state.wavunit]}\"]"
                else:
                    st.session_state.wavmsg = f":red[wavunit not set.]"

                if not isinstance(st.session_state.sciencefile, str):
                    st.session_state.scimsg = f":red[sciencefile not set, non-string type {type(st.session_state.sciencefile)} detected.]"
                elif len(st.session_state.sciencefile) > 0:
                    if os.path.exists(st.session_state.sciencefile):
                        st.session_state.scimsg = f":green[sciencefile set: \"{st.session_state.sciencefile}\"]"
                    else:
                        st.session_state.scimsg = f":red[sciencefile not set, file \"{st.session_state.sciencefile}\" not found]"
                        st.session_state.sciencefile = ""
                else:
                    st.session_state.scimsg = f":red[sciencefile not set.]"

                if not isinstance(st.session_state.telluricwavelengths, (list, np.ndarray)):
                    st.session_state.telwavmsg = f":red[telluricwavelengths not set, non-list type {type(st.session_state.telluricwavelengths)} detected.]"
                elif len(st.session_state.telluricwavelengths) > 0:
                    temporary = []
                    bad = False
                    for i in st.session_state.telluricwavelengths:
                        if isinstance(i, (list, np.ndarray)):
                            for ii in i:
                                temporary.append(ii)
                        elif isinstance(i, (float,int)):
                            temporary.append(i)
                        else:
                            badtype = type(i)
                            bad = True
                            break
                    st.session_state.telluricwavelengths = temporary
                    if bad:
                        st.session_state.teltransmsg = f":red[telluricwvelengths not set; unsupported data type {badtype} detected in list]"
                        st.session_state.telluricwavelengths = []
                    else:
                        st.session_state.telwavmsg = f":green[telluricwavelengths set: \[{st.session_state.telluricwavelengths[0]:.3f} ... {st.session_state.telluricwavelengths[-1]:.3f}\]]"
                else:
                    st.session_state.telwavmsg = ":orange[telluricwavelengths not set.]"

                if not isinstance(st.session_state.tellurictransmission, (list, np.ndarray)):
                    st.session_state.teltransmsg = f":red[tellurictransmission not set, non-list type {type(st.session_state.tellurictransmission)} detected.]"
                elif len(st.session_state.tellurictransmission) > 0:
                    temporary = []
                    bad = False
                    for i in st.session_state.tellurictransmission:
                        if isinstance(i, (list, np.ndarray)):
                            for ii in i:
                                temporary.append(ii)
                        elif isinstance(i, (float,int)):
                            temporary.append(i)
                        else:
                            badtype = type(i)
                            bad = True
                            break
                    st.session_state.tellurictransmission = temporary
                    if bad:
                        st.session_state.teltransmsg = f":red[tellurictransmission not set; unsupported data type {badtype} detected in list]"
                        st.session_state.tellurictransmission = []
                    else:
                        st.session_state.teltransmsg = f":green[tellurictransmission set: \[{st.session_state.tellurictransmission[0]:.3f} ... {st.session_state.tellurictransmission[-1]:.3f}\]]"
                else:
                    st.session_state.teltransmsg = ":orange[tellurictransmission not set.]"

                


                st.rerun()
            st.markdown(st.session_state.wavmsg)
            st.markdown(st.session_state.fluxmsg)
            st.markdown(st.session_state.errmsg)
            st.markdown(st.session_state.unitmsg)
            st.markdown(st.session_state.scimsg)
            st.markdown(st.session_state.telwavmsg)
            st.markdown(st.session_state.teltransmsg)
            


        #if sciencefile != "":
            #print(sciencefile) 
                

    with molectab:
        preminerror, premaxerror = False, False
        updateplot = False
        st.header("MOLECFIT Correction")
        molectopcol1, molectopcol2, molectopcol3 = st.columns([0.4,2,0.4])
        if any(':red[' in z for z in [st.session_state.wavmsg, st.session_state.fluxmsg, st.session_state.errmsg, st.session_state.scimsg, st.session_state.unitmsg]):
            goodtogo = False
        if not goodtogo:
            st.markdown(":red[MUST SUCCESSFULLY IMPORT WAVELENGTHS, FLUXES, ERRORS, UNITS, AND SCIENCE FILE BEFORE CORRECTION]")
        else:
            with molectopcol1:
                for i in range(10):
                    st.text("")
                
                selectedchipnum = st.selectbox("Select Segment to Correct", [z for z in range(len(st.session_state.wavelengths))], disabled=(not st.session_state.selectmode))
                if selectedchipnum != st.session_state.chipnum:
                    st.session_state.chipnum = selectedchipnum
                    st.session_state.dispmin = st.session_state.wavelengths[selectedchipnum][0]
                    st.session_state.dispmax = st.session_state.wavelengths[selectedchipnum][-1]
                dispminstr = st.text_input(f'Display Wavelength Minimum (MIN: {st.session_state.wavelengths[st.session_state.chipnum][0]:.1f})', value=f"{st.session_state.dispmin:.2f}")
                if float(dispminstr) != st.session_state.dispmin:
                    st.session_state.dispmin = float(dispminstr)
                    updateplot = True
                    #st.rerun()
                dispmaxstr = st.text_input(f'Display Wavelength Maximum (MAX: {st.session_state.wavelengths[st.session_state.chipnum][-1]:.1f})', value=f"{st.session_state.dispmax:.2f}")
                if float(dispmaxstr) != st.session_state.dispmax:
                    st.session_state.dispmax = float(dispmaxstr)
                    updateplot = True
                    #st.rerun()
            with molectopcol2:
                if isinstance(st.session_state.chipnum, int):
                    if st.session_state.selectmode == True:
                        if st.session_state.showtelselect == True:
                            plot = makeplot(st.session_state.wavelengths[st.session_state.chipnum], st.session_state.fluxes[st.session_state.chipnum], st.session_state.telluricwavelengths, st.session_state.tellurictransmission)
                        else:
                            plot = makeplot(st.session_state.wavelengths[st.session_state.chipnum], st.session_state.fluxes[st.session_state.chipnum], [], [])
                        buf = BytesIO()
                        plot.savefig(buf, format="png", facecolor='#DDDDDD', bbox_inches='tight')
                        st.session_state.currentplot = buf
                    if st.session_state.itmode == True:
                        fig, ax = plt.subplots(figsize=[10, 5.1])
                        ax.plot(st.session_state.wavelengths[st.session_state.chipnum], st.session_state.fluxes[st.session_state.chipnum])
                        if not st.session_state.overlap:
                            ax.plot(st.session_state.wavelengths[st.session_state.chipnum], [w + np.mean(st.session_state.tempfluxes)/3 for w in st.session_state.tempfluxes])
                        else:
                            ax.plot(st.session_state.wavelengths[st.session_state.chipnum], st.session_state.tempfluxes)
                        #ax.axvspan(st.session_state.wavmin, st.session_state.wavmax, color='green', alpha=0.2)
                        for i in range(st.session_state.winnum):
                            try:
                                ax.axvspan(float(st.session_state.wavmin[i]), float(st.session_state.wavmax[i]), color='green', alpha=0.2)
                            except:
                                pass
                        oldfluxfull = st.session_state.fluxes[st.session_state.chipnum]
                        oldwavfull = st.session_state.wavelengths[st.session_state.chipnum]
                        oldchunk = [oldfluxfull[z] for z in range(len(oldwavfull)) if (oldwavfull[z] > st.session_state.dispmin and oldwavfull[z] < st.session_state.dispmax)]
                        if st.session_state.overlap:
                            newfluxfull = st.session_state.tempfluxes
                        else:
                            newfluxfull = [w + np.mean(st.session_state.tempfluxes)/3 for w in st.session_state.tempfluxes]
                        newchunk = [newfluxfull[z] for z in range(len(oldwavfull)) if (oldwavfull[z] > st.session_state.dispmin and oldwavfull[z] < st.session_state.dispmax)]
                        ax.set_ylim(min(min(oldchunk), min(newchunk))*0.95, max(max(oldchunk), max(newchunk))*1.02)
                        ax.set_xlabel(f"Wavelength ({NAME_SHORT[NAME_CONVERSION[st.session_state.wavunit]]})")
                        ax.set_ylabel("Flux")
                        ax.set_title("Uncorrected and Corrected Spectra")
                        #ax.set_xlim([0.9999 * st.session_state.wavelengths[st.session_state.chipnum][0], 1.0001 * st.session_state.wavelengths[st.session_state.chipnum][-1]])
                        ax.set_xlim([0.9999 * st.session_state.dispmin, 1.0001 * st.session_state.dispmax])
                        if st.session_state.overlap:
                            if isinstance(st.session_state.chis[-1],list):
                                ax.legend(["Uncorrected Spectrum", "Corrected Spectrum", f"Reduced X^2: {[float(z) for z in st.session_state.chis[-1]]}"])
                            else:
                                ax.legend(["Uncorrected Spectrum", "Corrected Spectrum", f"Reduced X^2: {float(st.session_state.chis[-1])}"])
                        else:
                            if isinstance(st.session_state.chis[-1],list):
                                ax.legend(["Uncorrected Spectrum", "Corrected Spectrum + Offset", f"Reduced X^2: {[float(z) for z in st.session_state.chis[-1]]}"])
                            else:
                                ax.legend(["Uncorrected Spectrum", "Corrected Spectrum + Offset", f"Reduced X^2: {float(st.session_state.chis[-1])}"])
                        if st.session_state.showtel:
                            ax2 = ax.twinx()
                            color2 = "green"
                            ax2.set_ylabel('Transmission', color=color2, fontsize=24)
                            ax2.tick_params(axis='y', labelcolor=color2)
                            ax2.plot(st.session_state.telluricwavelengths,st.session_state.tellurictransmission,color=color2)
                            ax2.set_xlim([0.9999 * st.session_state.dispmin, 1.0001 * st.session_state.dispmax])
                            ax2.set_ylim([0,1])
                        buf = BytesIO()
                        fig.savefig(buf, format="png", facecolor='#DDDDDD', bbox_inches='tight')
                        st.session_state.currentplot = buf
    
                    st.image(st.session_state.currentplot, use_column_width="always")
                    
                    if st.session_state.selectmode:
                        if st.button("CORRECT THIS SEGMENT"):
                            st.session_state.selectmode = False
                            st.session_state.disablechip = True
                            st.rerun()
    
            with molectopcol3:
                if isinstance(st.session_state.chipnum, int) and st.session_state.selectmode == True:
                    for i in range(12):
                        st.text("\n")
                    showtelselectbutton = st.checkbox("Show Tellurics", value=st.session_state.showtelselect)
                    if showtelselectbutton != st.session_state.showtelselect:
                        st.session_state.showtelselect = showtelselectbutton
                        st.rerun()
                if st.session_state.disp_save:
                    st.session_state.disp_save = False
                    if len(st.session_state.sciencefile) > 0:
                        newname = ""
                        for letter,z in enumerate(st.session_state.sciencefile):
                            if z == "." and letter != 0:
                                break
                            else:
                                newname = newname + z
                        st.markdown(f":green[Corrected data saved to: \"{newname}_CORRECTED.fits\"]")
                    else:
                        st.markdown(f":green[Corrected data saved to: \"corrected_data.fits\"]")
                    st.markdown(f":green[MOLECFIT files stored in: {st.session_state.save_directory}]")
                if isinstance(st.session_state.chipnum, int) and st.session_state.selectmode == False:
                    for i in range(10):
                        st.text("")
                    if st.button("Run MOLECFIT",type="primary"):
                        if (not all([m == 0 for m in st.session_state.okaymax]) or not all([m==0 for m in st.session_state.okaymin]) or not all([m==0 for m in st.session_state.abunderrs])):
                            st.markdown(":red[Must fix all invalid entries]")
                        else:
                            if st.session_state.itmode == False:
                                st.session_state.itmode = True
                            with st.spinner("Running MOLECFIT..."):
                                if st.session_state.correctselect:
                                    streamlit_run_molecfit_seperate()
                                else:
                                    streamlit_run_molecfit_together()
                            st.rerun()
                    correctselectradio = st.radio("Correct All or Selection", ["Correct entire segment\n (fit only selected region)", "Correct only selected region"], disabled = st.session_state.itmode, label_visibility = "hidden")
                    if not st.session_state.itmode and correctselectradio == "Correct only selected region":
                        st.session_state.correctselect = True
                    elif not st.session_state.itmode:
                        st.session_state.correctselect = False
                    rerunchange = False
                    if st.session_state.itmode:
                        rerunchange = False
                        overlapbutton = st.checkbox("Overlap Old/New", value=st.session_state.overlap)
                        if overlapbutton != st.session_state.overlap:
                            st.session_state.overlap = overlapbutton
                            rerunchange = True
                        showtelbutton = st.checkbox("Show Tellurics", value=st.session_state.showtel)
                        if showtelbutton != st.session_state.showtel:
                            st.session_state.showtel = showtelbutton
                            rerunchange = True
                        st.session_state.verbose = st.checkbox("Verbose Terminal")
                        if rerunchange:
                            st.rerun()

                        if st.button("Finalize Correction",type="secondary"):
                            final_button()
                    else:
                        showtelselectbutton = st.checkbox("Show Tellurics", value=st.session_state.showtelselect)
                        if showtelselectbutton != st.session_state.showtelselect:
                            st.session_state.showtelselect = showtelselectbutton
                            if st.session_state.showtelselect == True:
                                plot = makeplot(st.session_state.wavelengths[st.session_state.chipnum], st.session_state.fluxes[st.session_state.chipnum], st.session_state.telluricwavelengths, st.session_state.tellurictransmission)
                            else:
                                plot = makeplot(st.session_state.wavelengths[st.session_state.chipnum], st.session_state.fluxes[st.session_state.chipnum], [], [])
                            buf = BytesIO()
                            plot.savefig(buf, format="png", facecolor='#DDDDDD', bbox_inches='tight')
                            st.session_state.currentplot = buf
                            st.rerun()
                        st.session_state.verbose = st.checkbox("Verbose Terminal")
                        
                    if st.button("Cancel Correction", type="secondary"):
                        final_button(save = False)
    
    
            if st.session_state.selectmode == False:
                initialwinnum = st.session_state.winnum
                for i in range(st.session_state.winnum):
                    if i == initialwinnum - 1 and st.session_state.skipit:
                        st.session_state.skipit = False
                        continue
                    wavcol1, wavcol2, wavcol3, wavcol4 = st.columns(4)
                    with wavcol2:
                        textmin = st.text_input("Wavelength Minimum:", label_visibility = SHOW_LABEL[i == 0], disabled = st.session_state.itmode, key=f"min{st.session_state.tagnums[i]}")
                        if textmin != "" and len(st.session_state.wavmin) == st.session_state.winnum and i < len(st.session_state.wavmin) and textmin != st.session_state.wavmin[i]:
                            st.session_state.wavmin[i] = textmin
                            st.session_state.resetmin = True
                            try:
                                float(textmin)
                                if float(textmin) > st.session_state.wavelengths[st.session_state.chipnum][-1] or float(textmin) < st.session_state.wavelengths[st.session_state.chipnum][0]:
                                    gotoexcept()
                            except:
                                st.session_state.okaymin[i] = 1
                                st.session_state.resetmin = False
                                #st.markdown(f":red[MUST INPUT NUMBER WITHIN RANGE ABOVE: ({st.session_state.wavelengths[st.session_state.chipnum][0]:.3f}, {st.session_state.wavelengths[st.session_state.chipnum][-1]:.3f})]")
                            else:
                                try:
                                    if len(st.session_state.wavmax[i]) > 0 and float(st.session_state.wavmax[i]) <= float(textmin):
                                        gotoexcept()
                                except:
                                    st.session_state.resetmin = False
                                    st.session_state.okaymin[i] = 2

                                
                                updateplot = True


                        if st.session_state.resetmin:
                            st.session_state.okaymin[i] = 0

                        if st.session_state.okaymin[i] == 1:
                            st.markdown(f":red[MUST INPUT NUMBER WITHIN RANGE ABOVE: ({st.session_state.wavelengths[st.session_state.chipnum][0]:.3f}, {st.session_state.wavelengths[st.session_state.chipnum][-1]:.3f})]")
                        elif st.session_state.okaymin[i] == 2:
                            st.markdown(":red[ABOVE MIN MUST BE LESS THAN MAX]")

                    with wavcol3:
                        #textmax = st.text_input("Wavelength Maximum:", label_visibility = SHOW_LABEL[i == 0], value=st.session_state.wavmax[i], key=f"max{i}")
                        textmax = st.text_input("Wavelength Maximum:", label_visibility = SHOW_LABEL[i == 0], disabled = st.session_state.itmode, key=f"max{st.session_state.tagnums[i]}")
                        if textmax != "" and len(st.session_state.wavmax) == st.session_state.winnum and i < len(st.session_state.wavmax) and textmax != st.session_state.wavmax[i]:
                            st.session_state.resetmax = True
                            st.session_state.wavmax[i] = textmax
                            try:
                                float(textmax)
                                if float(textmax) > st.session_state.wavelengths[st.session_state.chipnum][-1] or float(textmax) < st.session_state.wavelengths[st.session_state.chipnum][0]:
                                    gotoexcept()
                            except:
                                st.session_state.okaymax[i] = 1
                                st.session_state.resetmax = False
                                #st.markdown(f":red[MUST INPUT NUMBER WITHIN RANGE ABOVE: ({st.session_state.wavelengths[st.session_state.chipnum][0]:.3f}, {st.session_state.wavelengths[st.session_state.chipnum][-1]:.3f})]")
                            else:
                                try:
                                    if len(st.session_state.wavmax[i]) > 0 and float(st.session_state.wavmin[i]) >= float(textmax):
                                        gotoexcept()
                                except:
                                    st.session_state.okaymax[i] = 2
                                    st.session_state.resetmax = False
                                    #st.markdown(":red[ABOVE MAX MUST BE GREATER THAN MIN]")
                                
                                updateplot = True
                        if st.session_state.resetmax:
                            st.session_state.okaymax[i] = 0

                        if st.session_state.okaymax[i] == 1:
                            st.markdown(f":red[MUST INPUT NUMBER WITHIN RANGE ABOVE: ({st.session_state.wavelengths[st.session_state.chipnum][0]:.3f}, {st.session_state.wavelengths[st.session_state.chipnum][-1]:.3f})]")
                        elif st.session_state.okaymax[i] == 2:
                            st.markdown(":red[ABOVE MAX MUST BE GREATER THAN MIN]")

                    if st.session_state.winnum > 1:
                        with wavcol4:
                            st.markdown('<p style="font-size:9px; opacity:0"> #</p>', unsafe_allow_html=True)
                            if st.button("X", disabled = st.session_state.itmode, key=f"del{i}"):
                                updateplot = True
                                st.session_state.wavmin.pop(i)
                                st.session_state.wavmax.pop(i)
                                st.session_state.tagnums.pop(i)
                                st.session_state.okaymax.pop(i)
                                st.session_state.okaymin.pop(i)
                                st.session_state.winnum -= 1
                                if i != initialwinnum - 1:
                                    st.session_state.skipit = True

                addcol1, addcol2, addcol3, addcol4 = st.columns(4)
                with addcol2:
                    if ~st.session_state.itmode:
                        if st.button("Add Window", disabled = st.session_state.itmode, type="secondary"):
                            st.session_state.winnum += 1
                            st.session_state.tagnums.append(st.session_state.tagnums[-1] + 1)
                            st.session_state.wavmin.append("")
                            st.session_state.wavmax.append("")
                            st.session_state.okaymax.append(0)
                            st.session_state.okaymin.append(0)
                            st.rerun()

                
                st.session_state.mingood = True
                st.session_state.maxgood = True
                if st.session_state.mingood and st.session_state.maxgood and ~st.session_state.itmode and updateplot:
                    updateplot = False
                    plot = makeplot(st.session_state.wavelengths[st.session_state.chipnum], st.session_state.fluxes[st.session_state.chipnum], st.session_state.telluricwavelengths, st.session_state.tellurictransmission, boxranges=[st.session_state.wavmin,st.session_state.wavmax])
                    buf = BytesIO()
                    plot.savefig(buf, format="png", facecolor='#DDDDDD', bbox_inches='tight')
                    st.session_state.currentplot = buf
                    st.rerun()
    
    
                thirdscol1, thirdscol2, thirdscol3 = st.columns([1, 6, 1])
                #with thirdscol2:
                    #st.text("Elemental Abundances")
                #st.markdown("<h1 style='text-align: center; font-size: 24px; position: relative; top: -20px;'>Elemental Abundances</h1>", unsafe_allow_html=True)
                with thirdscol2:
                    elecol1, elecol2, elecol3, elecol4, elecol5, elecol6, elecol7 = st.columns(7)
                    with elecol1:
                        texto2 = st.text_input("O2", value=st.session_state.abundances[0])
                        try:
                            st.session_state.abundances[0] = float(texto2)
                            st.session_state.abunderrs[0] = 0
                        except:
                            if len(texto2) > 0:
                                st.markdown(f":red[MUST INPUT NUMBER]")
                                st.session_state.abunderrs[0] = 1
                    with elecol2:
                        textco2 = st.text_input("CO2", value=st.session_state.abundances[1])
                        try:
                            st.session_state.abundances[1] = float(textco2)
                            st.session_state.abunderrs[1] = 0
                        except:
                            if len(textco2) > 0:
                                st.markdown(f":red[MUST INPUT NUMBER]")
                                st.session_state.abunderrs[1] = 1
                    with elecol3:
                        texth2o = st.text_input("H20", value=st.session_state.abundances[2])
                        try:
                            st.session_state.abundances[2] = float(texth2o)
                            st.session_state.abunderrs[2] = 0
                        except:
                            if len(texth2o) > 0:
                                st.markdown(f":red[MUST INPUT NUMBER]")
                                st.session_state.abunderrs[2] = 1
                    with elecol4:
                        textco = st.text_input("CO", value=st.session_state.abundances[3])
                        try:
                            st.session_state.abundances[3] = float(textco)
                            st.session_state.abunderrs[3] = 0
                        except:
                            if len(textco) > 0:
                                st.markdown(f":red[MUST INPUT NUMBER]")
                                st.session_state.abunderrs[3] = 1
                    with elecol5:
                        textch4 = st.text_input("CH4", value=st.session_state.abundances[4])
                        try:
                            st.session_state.abundances[4] = float(textch4)
                            st.session_state.abunderrs[4] = 0
                        except:
                            if len(textch4) > 0:
                                st.markdown(f":red[MUST INPUT NUMBER]")
                                st.session_state.abunderrs[4] = 1
                    with elecol6:
                        textn2o = st.text_input("N2O", value=st.session_state.abundances[5])
                        try:
                            st.session_state.abundances[5] = float(textn2o)
                            st.session_state.abunderrs[5] = 0
                        except:
                            if len(textn2o) > 0:
                                st.markdown(f":red[MUST INPUT NUMBER]")
                                st.session_state.abunderrs[5] = 1
                    with elecol7:
                        texto3 = st.text_input("O3", value=st.session_state.abundances[6])
                        try:
                            st.session_state.abundances[6] = float(texto3)
                            st.session_state.abunderrs[6] = 0
                        except:
                            if len(texto3) > 0:
                                st.markdown(f":red[MUST INPUT NUMBER]")
                                st.session_state.abunderrs[6] = 1
    
                if st.session_state.itmode:
                    if len(st.session_state.tried_abundances) > 0:
                        titlecol1, titlecol2 = st.columns([1,1.4])
                        with titlecol2:
                                st.markdown("Tried Abundances with Results:")
                        for i in range(len(st.session_state.tried_abundances)):
                            listcol1, listcol2, listcol3 = st.columns([1, 2, 1])
                            with listcol2:
                                abundancestring = f"\[O2: {st.session_state.tried_abundances[i][0]}, CO2: {st.session_state.tried_abundances[i][1]}, H2O: {st.session_state.tried_abundances[i][2]}, CO: {st.session_state.tried_abundances[i][3]}, CH4: {st.session_state.tried_abundances[i][4]}, N2O: {st.session_state.tried_abundances[i][5]}, O3: {st.session_state.tried_abundances[i][6]}\]"
                                logabund = f"O2: {st.session_state.tried_abundances[i][0]}, CO2: {st.session_state.tried_abundances[i][1]}, H2O: {st.session_state.tried_abundances[i][2]}, CO: {st.session_state.tried_abundances[i][3]}, CH4: {st.session_state.tried_abundances[i][4]}, N2O: {st.session_state.tried_abundances[i][5]}, O3: {st.session_state.tried_abundances[i][6]}"
                                if isinstance(st.session_state.chis[i], list):
                                    chistring = f",&emsp;X^2: {[float(z) for z in st.session_state.chis[i]]}"
                                    logchi = f"X^2: {[float(z) for z in st.session_state.chis[i]]}"
                                elif isinstance(st.session_state.chis[i], str):
                                    chistring = f",&emsp;X^2: {float(st.session_state.chis[i])}"
                                    logchi = f"X^2: {float(st.session_state.chis[i])}"
                                #if isinstance(st.session_state.dip_ratios[i], list):
                                    #dipstring = f",&emsp;Approximate Depth Correction: {[float(z) for z in st.session_state.dip_ratios[i]]}"
                                #elif isinstance(st.session_state.dip_ratios[i], str):
                                    #dipstring = f",&emsp;Approximate Depth Correction: {float(st.session_state.dip_ratios[i])}"
                                if st.session_state.correctselect:
                                    itchi = [float(z) for z in st.session_state.chis[i]]
                                    means = []
                                    for chichi in st.session_state.chis:
                                        ichi = [float(z) for z in chichi]
                                        means.append(np.mean(ichi))
                                    st.session_state.loglines.append(f"{logabund} | {logchi}\n")
                                    if float(np.mean(itchi)) == min(means):
                                        st.markdown(f":green[{abundancestring}{chistring}]")
                                    else:
                                        st.markdown(f"{abundancestring} {chistring}")
                                else:
                                    chichi = [float(z) for z in st.session_state.chis]
                                    st.session_state.loglines.append(f"{logabund} | {logchi}\n")
                                    if float(st.session_state.chis[i]) == min(chichi):
                                        st.markdown(f":green[{abundancestring}{chistring}]")
                                    else:
                                        st.markdown(f"{abundancestring}{chistring}")
                            with listcol3:
                                if st.button(ACTIVE_BUTTON_LABEL[i == st.session_state.appliedit], type = ACTIVE_BUTTON_COLOR[i == st.session_state.appliedit], key = f"rea{i}"):
                                    if i != st.session_state.appliedit:
                                        st.session_state.appliedit = i
                                        for z in range(5):
                                            st.session_state.abundances[z] = st.session_state.tried_abundances[i][z]
                                        st.session_state.tempfluxes = st.session_state.pastflux[i]
                                        st.session_state.temperrors = st.session_state.pasterr[i]
                                        st.rerun()
    
    
                        

if __name__ == "__main__":
    streamlit_molecfit()

