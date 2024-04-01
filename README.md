# Simple_MOLECFIT_GUI
Simple GUI for quickly and easily running ESO's MOLECFIT recipes. 

# Requirements
All MOLECFIT esorex recipes must be installed: https://www.eso.org/sci/software/cpl/esorex.html  
astropy==5.3.4  
matplotlib==3.6.2  
numpy==1.26.4  
streamlit==1.28.1  
streamlit_ace==0.1.1  

# File Structure
Data can be loaded in from any directory, so long as it is accurately pointed to in the code input section.
All outputs will be stored in the directory containing MOLGUI.py. 
The final output data will appear directly in the main directory, with the name "{sciencefile}\_CORRECTED.fits"
Other outputs (such as those from MOLECFIT directly) will appear in the outputs folder, typically somehow labelled with the date and time they were created.


# Tutorial
To run the code, run the command "streamlit run MOLGUI.py" in the directory containing MOLGUI.py, the "css" folder, and the "configs" folder containing esorex .rc files named as they are in this repository. This will open a browser window for the GUI.

To import your data, you will write a python script which populates specific variables.

Required variables:
"wavelengths": Data Wavelengths, can be a single list or a list of lists, segmented however you wish.
"fluxes": Data Fluxes, can be a single list or a list of lists, segmented however you wish.
"errors": Data Errors, can be a single list or a list of lists, segmented however you wish.
"wavunit": Units of Wavelength data, can be spelled out ("micron") or abbreviated ("nm"). If entry doesn't work, try different name
"sciencefile": Name of file containing original ESO data header. This header contains necessary information for running MOLECFIT

You can also define optional variables:
"telluricwavelengths": Wavelengths of Telluric Transmission Spectrum, can be a single list or a list of lists, ultimately makes no difference.
"tellurictransmission": Transmission Values of Transmission Spectrum, can be a single list or a list of lists, ultimately makes no difference

You will know you have successfully imported variables by a green confirmation message on the right of the screen, showing some of the imported information for reference.

After importing your data, navigate to the "MOLECFIT" tab at the top left.

Initially, you will see your first (or only) imported spectrum in the central plot. On the left, you have option to select between the segments as you imported them, in addition to being able to zoom in on your spectrum by changing the plotted lower and upper wavelength limits. NOTE: These options do not affect fitting, they only affect visualization. On the right of the plot, you have the option to toggle the telluric spectrum on and off.

When you have chosen the segment you would like to correct, click the "CORRECT THIS SEGMENT" button below the plot. More options will now be displayed.

The options to the left of the plot remain identical, although you can no longer change segment until you complete or cancel the correction. 
On the right, you will see "Run MOLECFIT," which runs MOLECFIT for the values that will be selected below the plot. 
Below this are 2 radio buttons. 
The first, "Correct entire segment (fit only selected region" is default MOLECFIT behavior, in which you can select wavelength regions to be fit, but the calculated correction will be universally applied to the selected segment. 
"Correct only selected region" will run MOLECFIT and apply the correction to ONLY the selected region. 
"Show Tellurics" is the same option to simply toggle telluric visualization. 
"Verbose Terminal" will have MOLECFIT print its usual output as it runs in the terminal. Toggling this off will display limited progress information.

Below the plot, you can select the Wavelength Region to be fit (or optionally corrected). This is done by specifying windows, each with the minimum wavelength value on the left, and the maximum on the right. 
Clicking "Add Window" will allow you to have multiple, seperate windows to be fit/corrected.
To delete a window, click the "X" button that will appear when you have more than 1 window.
Below wavelength specification, you will see the options for species abundances used for the correction.

Once wavelength windows and abundances have been specified, press the red "Run MOLECFIT" button the right to proceed. The code will now run MOLECFIT, which will take some time.

After MOLECFIT completes, the central plot will change. Now, instead of showing the base data with tellurics, it will now compare the original data (blue) to the corrected data (orange). By default, the orange spectrum will be shifted upwards for easier shape comparison. This can be toggled off with the "Overlap Old/New" option to the right of the plot.

Once you have begun applying correction, a previous correction log will appear below the abundance specification area. This will show the input abundances and the X^2 of the fit. The red "Current" button indicates the currently applied correction. After running MOLECFIT more than once, lines in the log will instead have a black "Reapply" buttons next to them. Clicking this button will again display the fit associated with those abundances.

When you have found your favored fit, make sure it has the red "Current" button next to it in the log, and click the black "Finalize Correction" button to the right of the plot. This will completely apply the correction and create a new output file named "{sciencefile}_CORRECTED.fits"

You can press the "Cancel Correction" button to the right of the plot at any time to not output any file and save no correction. 

This concludes the tutorial. The core premise of the GUI is to allow you to repeatedly iterate over several segments quickly and easily, without having to rerun MOLECFIT for every abundance you've tried for comparison. 
