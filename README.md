Source code provided by NOAA into the public domain for estimating geomagnetic fields based on geolocation.
This repo provides an easy way to build and test this code for educational and research purposes. 


# World Magnetic Model Software and support documents


### Sublibrary Files

- `GeomagnetismLibrary.c`              WMM Subroutine library, C functions
- `GeomagnetismHeader.h`               WMM Subroutine library, C header file 

### Main Programs

- `wmm_point.c`                 Command prompt version for single point computation
- `wmm_grid.c`                  Grid, profile and time series computation, C main function
- `wmm_file.c`                  C program which takes a coordinate file as input with WMM ISO formatted coefficients as input

### Data Files

- `WMM.COF`                 WMM Coefficients file


### Supporting Documents

- `Geomagnetism_Library_Software_Manual.pdf`      Technical documents of the WMM Model and softwares


### Excecutables

- `wmm_point`               Command Prompt program for single point 
- `wmm_grid`                Grid, time series or profile
- `wmm_file`                File processing 


### Test Files

- `WMM2015v2testvalues.pdf`         Test values for WMM2015
- `sample_input_file.txt`           Sample input file for program wmm_file.exe 
- `sample_output_file.txt`          Sample output file for program  wmmfile.exe


## Compiling with gcc

There is a Makefile in the source directory for building with gcc
`gcc -lm inputfile [dependencies] -o outputfile`
For example, the `wmm_cmd.c` can be compiled as
`gcc -lm wmm_cmd.c GeomagnetismLibrary.c -o wmm_cmd.exe`
Note: The library files `GeomagnetismLibrary.c` and `GeomagnetismHeader.h` should reside in the same directory for compiling.

## Compiling with cmake

- An experimental cmakefile is provided by this repo. 

## Executing the file processing program

The file processing program accepts this syntax:
`wmm_file.exe f INPUT_FILE.txt OUTPUT_FILE.txt`

Additionally to have uncertainty values appended to the output file add the "e" flag or "Errors".  For example:
`wmm_file.exe f e INPUT_FILE.txt OUTPUT_FILE.txt`

or

`wmm_file.exe f --Errors INPUT_FILE.txt OUTPUT_FILE.txt`


# Model Software Support


*  National Centers for Environmental Information (NCEI)
*  E/NE42 325 Broadway
*  Boulder, CO 80303 USA
*  Attn: Manoj Nair or Arnaud Chulliat
*  Phone:  (303) 497-4642 or -6522
*  Email:  geomag.models@noaa.gov
For more details about the World Magnetic Model visit 
http://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml
 




       
