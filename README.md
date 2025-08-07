# air2water+ice

```
      .__       ________                  __                            .__              
_____  |__|______\_____  \__  _  _______ _/  |_  ___________     .__     |__| ____  ____  
\__  \ |  \_  __ \/  ____/\ \/ \/ /\__  \\   __\/ __ \_  __ \  __|  |___ |  |/ ___\/ __ \
 / __ \|  ||  | \/       \ \     /  / __ \|  | \  ___/|  | \/ /__    __/ |  \  \__\  ___/
(____  /__||__|  \_______ \ \/\_/  (____  /__|  \___  >__|       |__|    |__|\___  >___  >
     \/                  \/             \/          \/                           \/    \/
```

A model to predict Lake Surface Water Temperature (LSWT) and ice phenology, thickness, and composition using minimal inputs.\
**Version 1.0.0 - August 2025**

**Contact**: [s.piccolroaz@unitn.it](mailto\:s.piccolroaz@unitn.it)

## How to Cite

Fregona M., Lepparanta M., Jansen J., Ala-Konni J., Mammarella I., and Piccolroaz S. (2025), *Simulating lake ice phenology, thickness, and composition using minimal inputs: an extension of the air2water model*, Under Review for JAMES

---

## Preamble

This is an extension of the original hybrid physics-based/statistical air2water model by Piccolroaz et al. (2013), integrating an ice module. The extended model simulates ice cover dynamics, distinguishing between black ice (formed by direct lake water freezing) and white ice (from accumulated snow submerged and frozen over black ice). The extension preserves the original framework’s key advantages: minimal input requirements (air temperature and precipitation), low parametrization, and physically based governing equations.

**References**:

- Piccolroaz et al. (2013), DOI: [https://doi.org/10.5194/hess-17-3323-2013](https://doi.org/10.5194/hess-17-3323-2013)
- Toffolon et al. (2014), DOI: [https://doi.org/10.4319/lo.2014.59.6.2185](https://doi.org/10.4319/lo.2014.59.6.2185)
- Piccolroaz (2016), DOI: [https://doi.org/10.4081/aiol.2016.5791](https://doi.org/10.4081/aiol.2016.5791)

For operational details about how to run the model, process, and analyze results, refer to the GitHub repository: [https://github.com/spiccolroaz/air2water](https://github.com/spiccolroaz/air2water)

---

## Content of the Folder

- Precompiled executable files for Linux (`air2water+ice_v1.0.out`) and Windows (`air2water+ice_v1.0.exe`)
- Source files of the original Fortran 90 scripts (`src` folder)
- Example test cases for Lake Kilpisjärvi, Lake Kallavesi, and Lake Pyhäjärvi in Finland (see manuscript for data sources)

---

## Structure of Input and Output Files

**Important**: Do not change the encoding of input files. They must be ANSI text files. UTF-8 encoded files will not be read correctly.

### 1. Input Files

#### 1.1 `input.txt`

Must be in the same folder as the executable. This file contains the input information:


	! Input file
	Kallavesi	 					! LakeName: name of the lake
	air_ERA5						! IDair: name/ID of the air temperature station
	2014-2023zeros_no2020Final		! IDwater: name/ID of the water temperature station
	c								! type of series: c=continuous, m=mean year
	1d 			    				! time resolution: 1d=daily, nw=n weeks (n=1,2,...), 1m=monthly
	ICES							! model version: 1=a2w 4 par; 2=a2w 6 par, 3=a2w 8 par, ICES=6 par w/ ice module
	0								! Threshold temperature for ice formation
	NSE								! objective function: KGE, NSE, RMS. Use NSE if ICES is selected
	1								! beta (NSE weight for multiobjective calibration; 1=only LSWT; 0=only ice; 0.5: both equal weight)
	0								! calibration flag for the ice module (0=tot ice is ued for calibration; 1=black and white ice are used for calibration)
	0 								! snowfall flag (0=snowfall calculated from total precipitaion using a9 and a12; 1=snowfall from observations)
	8.9		    					! mean depth required by the ice module (meters) 19.5=Kilpisjarvi, 5.4=Pyhajarvi, 8.9=Kallavesi
	CRN								! mod_num : numerical method used to solve the LSWT model equation: EUL (Euler Explicit), RK2 (Runge-Kutta 2), RK4 (Runge-Kutta 4), CRN (Crank-Nicolson)
	PSO        						! run mode: PSO (Particle Swarm Optimization), LATHYP (Latin Hypercube), RANSAM (Random Sampling), or FORWARD
	0.01							! minimum percentage of data: 0...1. E.g., when using 1m time resolution, the monthly average is set to NaN if available data during one month are less than this percentage (1% in this case) 
	2000							! n_run: number of iterations
	-999							! minimum efficiency index (i.e., RMS, NSE, KGE). All model realization with efficiency index greater than this threshold are saved in file 0_...
	1 								! log_flag: 1--> parameters 2 and 3 are searched in a logarithmic space; 0--> parameters 2 and 3 are searched in a linear space
	0.9 							! Courant number 


#### 1.2 `PSO.txt`

Contains parameters for the PSO algorithm.

	! PSO parameters
	2000			! number of particles
	2 2	   			! c1, c2: constants known as cognitive and social learning factors (used in the equation of motion of particles)
	0.9  0.4    	! inertia max and min (used in the equation of motion of particles)
 
#### 1.3 Folder: `LakeName` (e.g., `Kallavesi`)

- Calibration input: `IDair_IDwater_typeofseriesc.txt` (e.g., air_ERA5_2014-2023zeros_no2020Final_cc.txt)
  The file has 11 columns (year, month, day, air temperature, water temperature, total ice thickness, black ice thickness, white ice thickness, snow thickness, total precipitation, snowfall):
  
			2014	1	1	 2.3604	-999	-999	-999	-999	-999	0.0036751628	0.00023376942
			2014	1	2	-0.0592	-999	-999	-999	-999	-999	0.00064960122	0.00036939979
			2014	1	3	-1.9642	-999	-999	-999	-999	-999	0.00025397539	0.00025603175
			2014	1	4	-0.5558	-999	-999	-999	-999	-999	0.00071430206	0.00071232021
			2014	1	5	 0.1198	-999	-999	-999	-999	-999	0.00061652064	0.00054576993
			2014	1	6	-0.0127	-999	-999	-999	-999	-999	0.00079214573	0.0004537776	
			...
    NOTES
  	- The series of observed air temperature and total precipitation must be complete. They cannot have gaps or no data.
    - The series of snowfall is used only if snowfall flag = 1. In this case, it must be complete. Otherwise, snowfall is calculated from total precipitation using the internal algorithm based on parameters a9 and a12.
    - The series of observed water temperature and ice thickness can contain no data (-999).
    - The series of data must start on the 1st of January. If data are available after that date, air temperature should be reconstructed and the value -999 assigned to water temperature.
    - Time series are always at daily time scale, as the equation of the model is solved with daily time step. The model automatically evaluates weekly, multi-weekly, or monthly averages (of water temperature) when using different time scales for model calibration.
    - Temperatures are in °C, ice ticknesses in m and precipitation and snowfall in m/day
      
- Validation input: `IDair_IDwater_typeofseriesv.txt`  (e.g., air_ERA5_2014-2023zeros_no2020Final_cv.txt)
  The file has the same structure as the calibration input file.
- Parameter range: `parameters.txt`
  The file contains the range of variation of each parameter of the model. The range of variation of the parameters can be defined on the basis of physical considerations and is strictly lake dependent, thus its reasonable a priori estimation is crucial to obtain reliable results. The range of parameters a1-a8 can be estimated using the equations presented in the Supplementary Material of Piccolroaz (2016) or, visually, from Figure 3 of the same manuscript. A pre-processing script in matlab is also available to evaluate the a priori physical range of model's parameter (see  [https://github.com/spiccolroaz/air2water](https://github.com/spiccolroaz/air2water)). The a priori range of variation of the ice model parameters a9-a12, can be reasonably fixed as reported below. The structure of the file 'parameters.txt' is as follows:

			-0.9972	0.0117	0.0199	1.00	0.00	0.00	0.00	0.0	-10	0	0.8	0
			 2.0000	0.2930	1.1595	46.52	17.07	1.00	26.52	0.5	 10	35	1.2	10

  The first line contains the minimum value of each parameters, while the second contains the maximum values. There are 12 columns, as the number of parameters. Parameters a9-a12 are not used in the original air2water model, while parameters a7 and a8 are not used when running the ice module.
- Forward run parameters: `parameters_forward.txt`
  The file contains a set of parameters to be used to run the model in forward mode, see run mode in the input file input.txt. The structure of the file is as follows:

			   0.153230   0.044833   0.049815  13.130569   0.301900   0.598710   0.000000   0.000000   2.011250   4.409121   1.015691   4.886216

### 2. Output Files

Located in `output_modelversion` (e.g., `output_ICES`):

- `0_PSO_*.out`: Binary file with parameters and efficiency metrics. The file contains a matrix with 17 columns. Each row refers to the set of parameters (columns 1-12) and the associated efficiency indexes (columns 13-17; total eff. index, eff. index LSWT, eff. index total ice, eff. index black ice, eff. index white ice) of each iteration performed by the optimization algorithm. Only the parameter sets that provide an efficiency larger than the minimum efficiency index defined in the file 'input.txt' are saved. Values are saved in double precision. This file allows: a) drawing the dotty plots of parameters to evaluate whether the a priori range of variations of 	parameters have been reasonably defined (i.e., not too narrow, not too large), b) evaluating whether the optimization (searching) algorithm converged towards the best solution, c) perform uncertainty analyses (when using LATHYP).
- `1_PSO_*.out`: ASCII file that contains the best set of parameters (1st line), the value of the total efficiency index during the calibration period (2nd line), and the value of the total efficiency index during the validation period (if any, 3rd line)
- `2_PSO_*.out`: Calibration results. Row represent time and columns are:
  
   year
 
   month

   day

   observed air temperature

   observed water temperature
	
   simulated water temperature
	
   observed total ice thickness  
	
   simulated total ice thickness  
	
   observed black ice thickness  
	
   simulated black ice thickness  	
	
   observed white ice thickness  
	
   simulated white ice thickness  		
	
   observed snow thickness  
	
   simulated snow thickness  	
 
   observed total precipitation
  
   observed snowfall
  
   simulated snowfall

  Note that the first year is replicated and is used as warm up year to mitigate the influence of initial conditions. 
  During the warm up year, year, month and day are set to -999.
- `3_PSO_*.out`: Validation results, with the same structure of file '2_PSO... .out'
- `4_PSO_*.out`: Delta (dimensionless mixed layer thickness) for calibration. Refer to Toffolon et al. (2014) and Piccolroaz et al., (2015) to see how to convert delta into the dimensional thickness of the well-mixed surface layer.
- `5_PSO_*.out`: Delta for validation

---

**Trento, 7 August 2025**\
Sebastiano Piccolroaz\
📧 [s.piccolroaz@unitn.it](mailto\:s.piccolroaz@unitn.it)

