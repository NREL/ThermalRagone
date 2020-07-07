# ThermalRagone

The numerical model used to generate the data in the manuscript Rate Capability and Ragone Plots for Thermal Energy Storage is provided in the GitHub repository (NREL/ThermalRagone). The code was developed in MATLAB R2019b. The MATLAB software is required to run the code, although no non-standard hardware is needed. Descriptions of the files are provided below:

•	NumericalModel_NatEng.m is the main code that sets up the numerical model and calculates the local temperatures, enthalpies and important global parameters such as the specific energy and power. More information about the numerical model is provided in the Methods section of the main paper.
•	ductflow.m calculates the friction factor and Nusselt number of the glycol-water mixture flowing through the flat channel (duct). The correlation used depends on the Reynolds number of the fluid.
•	The GlycolProperties folder includes MATLAB codes that calculate the specific heat, conductivity, viscosity, and density of a propylene-glycol water mixture using an artificial neural network.
•	The PCMProperties folder includes enthalpy and temperature data for different potential PCMs. For more information about the PCM properties, please see the Supplemental Information.

The steps to run the code are described below. The current code outputs sample data for discharging a thermal storage device with a PCM of tetradecane embedded in a graphite matrix (see Supplemental Information for material properties) at a C-rate of 1. The storage device has the baseline geometry used in the paper.

Instructions for Use
1.	Download the files from GitHub and place them all in the same folder.
2.	Open the NumericalModel_NatEng.m file in MATLAB (ideally in version R2019b or newer)
3.	Add the GlycolProperties and PCMProperties folder to the MATLAB path. One way to do this is to select the two folders in the “Current Folders” pane and right click on the files. Then select Add to Path – Selected Folders.
4.	Run the code

The code will output all calculated variables (see comments) and the rate capability and Ragone plots. The expected runtime for the demo case should be between 10 and 20 seconds. Many properties and operating conditions can be changed in the code, including the C-rate, geometry, PCM transport properties, h-T relationships, and the cutoff temperature. Please refer to the comments for more information. Three additional approximate h-T relationships are included in the files, which were generated based on the PCM properties described in the Supplemental Information.
