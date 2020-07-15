# PlateReaderAnalysis
This script plots data of growth experiments measured in a plate reader. 

All files of this script must be contained in the same folder on your computer. Python 3.X (best 3.7) is needed for execution. All other files regarding your experiment can be in another file. The script only takes one input to read the data - the path of the AnalysisInformation.xlsx file, in which you specify which and how the data should be imported.

So far ouput of the following readers can be analyzed: victor, biotek, SPECTROstar (BMG), CLARIOstar (BMG), FLUOstar Omega (BMG).
The required format depends on the used reader. The respective functions to read the data describe the required output in more detail. The script recognizes automatically which device was used for the experiment.

For analysis each experiment needs a metainfo table placed in the same folder and with the same name but with an added '\_metainfo'. If luminescence data was measured and is required to be in another file (only CLARIOstar at this point) it also has to be in the same folder as the OD-data and with the same name but '\_lum' added.

