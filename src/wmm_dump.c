#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>


#include "GeomagnetismHeader.h"
#include "EGM9615.h"

#define LAT_DEGREES_STEP 5.0
#define LON_DEGREES_STEP 2.5



int MAG_Grid_gen(
    MAGtype_CoordGeodetic minimum,
    MAGtype_CoordGeodetic maximum,
    double lat_step_size,
    double lon_step_size,
    double altitude_step_size,
    double time_step,
    MAGtype_MagneticModel *MagneticModel,
    MAGtype_Geoid *Geoid,
    MAGtype_Ellipsoid Ellip,
    MAGtype_Date StartDate,
    MAGtype_Date EndDate
    )
{
    int NumTerms;
    double a, b, c, d, PrintElement, ErrorElement = 0;

    MAGtype_MagneticModel *TimedMagneticModel;
    MAGtype_CoordSpherical CoordSpherical;
    MAGtype_MagneticResults MagneticResultsSph, MagneticResultsGeo, MagneticResultsSphVar, MagneticResultsGeoVar;
    MAGtype_SphericalHarmonicVariables *SphVariables;
    MAGtype_GeoMagneticElements GeoMagneticElements, Errors;
    MAGtype_LegendreFunction *LegendreFunction;
    MAGtype_Gradient Gradient;

    if (fabs(altitude_step_size) < 1.0e-10) altitude_step_size = 99999.0;
    if (fabs(time_step) < 1.0e-10) time_step = 99999.0;


    NumTerms = ((MagneticModel->nMax + 1) * (MagneticModel->nMax + 2) / 2);
    TimedMagneticModel = MAG_AllocateModelMemory(NumTerms);
    LegendreFunction = MAG_AllocateLegendreFunctionMemory(NumTerms); /* For storing the ALF functions */
    SphVariables = MAG_AllocateSphVarMemory(MagneticModel->nMax);
    // sets the loop initialization values
    a = minimum.HeightAboveGeoid;
    b = minimum.phi;
    c = minimum.lambda; //longitude
    d = StartDate.DecimalYear;

    printf("# min alt: %f lat: %f lon: %f year: %f \n",
        a,b,c,d
        );
    printf("# max alt: %f lat: %f lon: %f year: %f \n",
           maximum.HeightAboveGeoid,
           maximum.phi,
           maximum.lambda,
           EndDate.DecimalYear);

  printf("# step alt: %f lat: %f lon: %f year: %f \n",
         altitude_step_size,lat_step_size, lon_step_size, time_step);

  printf("# ===============\n");
  printf("# lat_deg, lon_deg, mag_x, mag_y, mag_z \n");

  unsigned total_points = 0;

    // Altitude loop
    for (minimum.HeightAboveGeoid = a;
      minimum.HeightAboveGeoid <= maximum.HeightAboveGeoid;
      minimum.HeightAboveGeoid += altitude_step_size)   {
        // Latitude loop
        for (minimum.phi = b; minimum.phi <= maximum.phi; minimum.phi += lat_step_size) {
            //Longitude loop
            for (minimum.lambda = c;
                 minimum.lambda <= maximum.lambda; minimum.lambda += lon_step_size) {
                if (Geoid->UseGeoid == 1)
                    MAG_ConvertGeoidToEllipsoidHeight(&minimum, Geoid); /* This converts the height above mean sea level to height above the WGS-84 ellipsoid */
                else
                    minimum.HeightAboveEllipsoid = minimum.HeightAboveGeoid;

                MAG_GeodeticToSpherical(Ellip, minimum, &CoordSpherical);
                MAG_ComputeSphericalHarmonicVariables(Ellip, CoordSpherical, MagneticModel->nMax,  SphVariables); /* Compute Spherical Harmonic variables  */
                MAG_AssociatedLegendreFunction(CoordSpherical, MagneticModel->nMax, LegendreFunction); /* Compute ALF  Equations 5-6, WMM Technical report*/

                // Year loop
                for (StartDate.DecimalYear = d;
                     StartDate.DecimalYear <= EndDate.DecimalYear; StartDate.DecimalYear += time_step)  {

                    MAG_TimelyModifyMagneticModel(StartDate, MagneticModel,
                                                  TimedMagneticModel); /*This modifies the Magnetic coefficients to the correct date. */
                    MAG_Summation(LegendreFunction, TimedMagneticModel, *SphVariables, CoordSpherical,
                                  &MagneticResultsSph); /* Accumulate the spherical harmonic coefficients Equations 10:12 , WMM Technical report*/
                    MAG_SecVarSummation(LegendreFunction, TimedMagneticModel, *SphVariables, CoordSpherical,
                                        &MagneticResultsSphVar); /*Sum the Secular Variation Coefficients, Equations 13:15 , WMM Technical report  */
                    MAG_RotateMagneticVector(CoordSpherical, minimum, MagneticResultsSph,
                                             &MagneticResultsGeo); /* Map the computed Magnetic fields to Geodetic coordinates Equation 16 , WMM Technical report */
                    MAG_RotateMagneticVector(CoordSpherical, minimum, MagneticResultsSphVar,
                                             &MagneticResultsGeoVar); /* Map the secular variation field components to Geodetic coordinates, Equation 17 , WMM Technical report*/
                    MAG_CalculateGeoMagneticElements(&MagneticResultsGeo,
                                                     &GeoMagneticElements); /* Calculate the Geomagnetic elements, Equation 18 , WMM Technical report */
                    MAG_CalculateGridVariation(minimum, &GeoMagneticElements);
                    MAG_CalculateSecularVariationElements(MagneticResultsGeoVar,
                                                          &GeoMagneticElements); /*Calculate the secular variation of each of the Geomagnetic elements, Equation 19, WMM Technical report*/
                    MAG_WMMErrorCalc(GeoMagneticElements.H, &Errors);


                    total_points++;
                    printf("%5.2f,%6.2f,\t%10.3f,%10.3f,%10.3f \n",
                        minimum.phi, minimum.lambda,
                        GeoMagneticElements.X, GeoMagneticElements.Y, GeoMagneticElements.Z);

                } // year loop
            } // Longitude Loop
        } // Latitude Loop
    } // Altitude Loop

    printf("# ========= Total Points: %u ====== \n", total_points);

    MAG_FreeMagneticModelMemory(TimedMagneticModel);
    MAG_FreeLegendreMemory(LegendreFunction);
    MAG_FreeSphVarMemory(SphVariables);

    return TRUE;
}

int main()
{
    MAGtype_MagneticModel * MagneticModels[1];
    MAGtype_Ellipsoid Ellip = {};
    MAGtype_CoordGeodetic minimum = {};
    MAGtype_CoordGeodetic maximum = {};
    MAGtype_Geoid Geoid = {};
    MAGtype_Date startdate = {};
    MAGtype_Date enddate =  {};
    int  i, epochs = 1;
    double lat_step_size = LAT_DEGREES_STEP;
    double lon_step_size = LON_DEGREES_STEP;
    double altitude_step_size = 500.0; //only one step
    double time_step_size = 1.0;
    char filename[] = "../bin/WMM.COF";
    char VersionDate[12];
    char ans[20];

    startdate.DecimalYear = 2019.1;
    enddate.DecimalYear = 2019.1;

    //longitude
  minimum.lambda = -180.0;
  maximum.lambda = 180.0;
  // geodetic latitude
  minimum.phi = -90.0;
  maximum.phi = 90.0;
  // height about geoid  in...meters?
  minimum.HeightAboveGeoid = 1.0;
  maximum.HeightAboveGeoid = 100.0;


    if(!MAG_robustReadMagModels(filename, &MagneticModels, 1)) {
        printf("\n '%s' not found. \n ", filename);
        return 1;
    }
    strncpy(VersionDate, VERSIONDATE_LARGE + 39, 11);
    VersionDate[11] = '\0';

    MAG_SetDefaults(&Ellip, &Geoid);
    /* Set EGM96 Geoid parameters */
    Geoid.GeoidHeightBuffer = GeoidHeightBuffer;
    Geoid.Geoid_Initialized = 1;


    MAG_Grid_gen(minimum, maximum,
        lat_step_size, lon_step_size,
        altitude_step_size, time_step_size,
        MagneticModels[0], &Geoid, Ellip, startdate, enddate);

    for(i = 0; i < epochs; i++) {
        MAG_FreeMagneticModelMemory(MagneticModels[i]);
    }

    return 0;
}
