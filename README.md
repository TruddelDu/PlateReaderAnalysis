# PlateReaderAnalysis

## General information
I wrote this program during my PhD to analyze my experimental data. I tried to keep it as general as possible but of course it still ended up to be best suited for analyzes of inhibitory concentrations of antibiotics under different environmental conditions and for different strains as this is what I used it for.

Unfortunately, documentation is rather lacking as I only learned the importance and how-to of documenting during the end of the development phase.

## Technologies
Project created with:
- python 3.8
- xlrd 1.2.0
- pandas 1.1.1
- numpy 1.19.1
- math
- sys
- numbers
- scipy 1.5.2
- matplotlib.pyplot 3.3.1
- seaborn 0.10.1
- sklearn

## Usage
The programm can prepare and plot data of growth experiments measured in a plate reader. Using metadata and user defined analysis-specifications, the raw data is automatically:
- read (staring from various input formats)
- background corrected
- luminescence data is normalized
- smoothing is applied
- sets of measurements that need to be analyzed together are determined (on both the level of dilution series as well as replicates)
- various analysis can be performed (both for single replicates as well as the average of all replicates)
  - inhibitory concentrations (IC)
  - interaction of two inducers on the IC
  - dose-response plots (on IC, growth rate and luminescence)
  - doubling time
  - timedependent OD & luminescence variations 

![AVGtimedependentOD_WT_Laspartomycin C](https://user-images.githubusercontent.com/68091502/139456853-2aae50b0-8e99-4302-9308-08de3450a8fc.png)

All files of this script must be contained in the same folder on your computer. All other files regarding your experiment can be in another file. The script only takes one input to read the data - the path of the AnalysisInformation.xlsx file, in which you specify which and how the data should be imported.

So far ouput of the following readers can be analyzed: victor, biotek, SPECTROstar (BMG), CLARIOstar (BMG), FLUOstar Omega (BMG).
The required format depends on the used reader. The respective functions to read the data describe the required output in more detail. The script recognizes automatically which device was used for the experiment.

For analysis each experiment needs a metainfo table placed in the same folder and with the same name but with an added '\_metainfo'. If luminescence data was measured and is required to be in another file (only CLARIOstar at this point) it also has to be in the same folder as the OD-data and with the same name but '\_lum' added.

