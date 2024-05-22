#ifndef _PACER_GEOMETRY_H__
#define _PACER_GEOMETRY_H__

#include <math.h>
#include "pacer_imager_defs.h"
#include "pacer_imager.h"
// #include "slalib.h"
#include <star/pal.h>

// #define SPEED_OF_LIGHT 299792458.0        // speed of light in m/s

/**
 * 
   Sime functions come from corr2uvfits written by Randall Wayth.
   However, I am taking incremental approach here. Adding only what is needed and verifying each function / equations with my own calculations
    - at least that's the plan which will be realised as much as possible  
*/

class CObsMetadata;

class PacerGeometry
{
public:
        /**
         * Convert coordinator from local topocentric East, North, Height units to
         * 'local' XYZ units. Local means Z point north, X points through the equator
         * from the geocenter along the local meridian and Y is East.
         * Latitude is geodetic, in radian.
         */
        static void ENH2XYZ_local(double east, double north, double height, double lattitude_rad, double& x, double& y, double& z) {
           double latSin = sin(lattitude_rad),
           latCos = cos(lattitude_rad);

           x = -north * latSin + height * latCos;
           y = east;
           z = north * latCos + height * latSin;
           
           // COMMENT : these coordinate are nearly ok to calculate UVW, but 
           // WARNING : they still need to be "precessed" to J2000 or Current epoch !
        }
        
        static double GetMJD(int year, int month, double day, int refHour, int refMinute, double refSecond)
        {
                int res;
                double jdBase;
                palCaldj(year, month, day, &jdBase, &res); // get MJD for calendar day of obs
                double ret =  jdBase +
                        //2400000.5 +  // convert MJD to JD
                        refHour/24.0 + refMinute/1440.0 + refSecond/86400.0; // add intra-day offset
                        
                printf("DEBUG_MJD2 : jdBase = %.8f (from %d/%d/%.8f) + %d/24.0 + %d/1440.00 + %.4f/86400.00 -> mjd = %.8f\n",jdBase,year, month, day, refHour,refMinute,refSecond,ret);
                return ret;        
        }

        
        // WARNING : angle values in radians here :
        struct UVWTimestepInfo
        {
                double rmatpr[3][3];
                double ha2000; // HA angle of the direction to the observed source (or phase centre / MWA beam centre), corrected for precesion etc
                double lmst;   // local mean sidreal time 
                double dec_aber; // DEC angle of the direction to the observed source (or phase centre / MWA beam centre), corrected for precesion etc 
                double lmst2000; // local mean sidreal time J2000 (?)
        };
        
        // TODO : possibly add UVWTimestepInfo& obsInfo to class CObsMetadata
        static void PrepareTimestepUVW( UVWTimestepInfo& uvwInfo, CObsMetadata& obsMetaData )
        {
           double xprec, yprec, zprec;
           double raHrs = obsMetaData.raHrs;
           double decDeg = obsMetaData.decDegs;
           double mjd =  obsMetaData.dateRequestedMJD; // hardcoding for testing was 59019.68752315;
           double arrayLongitudeRad = obsMetaData.geo_long*(M_PI/180.00);
           double arrayLattitudeRad = obsMetaData.geo_lat*(M_PI/180.00);
           double lmst = palDranrm(palGmst(mjd) + arrayLongitudeRad); 
           printf("DEBUG lmst : mjd = %.8f , arrayLongitudeRad = %.8f [rad] -> lmst = %.8f\n",mjd,arrayLongitudeRad,lmst);
                                 
           double ra_app, dec_app;                 
           palMap(raHrs*(M_PI/12.0), decDeg*(M_PI/180.0), 0.0, 0.0, 0.0, 0.0, 2000.0, mjd, &ra_app, &dec_app);                
           
           double ra_aber = (raHrs*15.00)*(M_PI/180.00);
           double dec_aber = decDeg*(M_PI/180.00);
           aber_radec_rad(2000.0, mjd, raHrs*(M_PI/12.0), decDeg*(M_PI/180.0), &ra_aber,&dec_aber);

                 
           double rmattr[3][3];
           palPrenut(2000.0, mjd, rmattr);
                 
           mat_transpose(rmattr, uvwInfo.rmatpr);
                 
           double newarrlat, ha2000, lmst2000;
           ha_dec_j2000(uvwInfo.rmatpr,lmst,arrayLattitudeRad,ra_aber,dec_aber,&ha2000,&newarrlat,&lmst2000);
           printf("DEBUG_NEW : HA2000 = %.8f , LMST2000 = %.8f , mjd = %.8f , ra_app = %.8f dec_app = %.8f , raHrs = %.8f decDeg = %.8f , lmst = %.8f , lmst_2000 = %.8f\n",ha2000,lmst2000,mjd,ra_app,dec_app,raHrs,decDeg,lmst,lmst2000);
           
           uvwInfo.dec_aber = dec_aber;
           uvwInfo.ha2000 = ha2000;
           uvwInfo.lmst = lmst;
           uvwInfo.lmst2000 = lmst2000;
        }
        
        static void CalcUVW_NEW(UVWTimestepInfo &uvwInfo, double x, double y, double z, double &u, double &v, double &w)
        {
                static bool printed = false;
                double ha2000 = uvwInfo.ha2000;
                double lmst = uvwInfo.lmst;
                double dec_aber = uvwInfo.dec_aber;
                double lmst2000 = uvwInfo.lmst2000;
                double xprec, yprec, zprec;
                
                //---------------------------------------------------------------------------------------------------------
                // hardcode test :
                // DEBUG : (x,y,z) = (392.6971,-55.5600,-57.8339), ha2000 = 6.1717 , lmst = 4.75669697 , dec_aber = -0.0775, lmst2000 = 4.75121533 -> (u,v,w) = (-98.63851479,-27.86078831,387.47507206)
                // lmst = 4.75669697;
                // 2023-09-11 Commented these 2 lines :
                // ha2000 = 6.1717;
                // lmst2000 = 4.75121533;
                //------------------------------------------------------------------------------------------------------
                
                precXYZ(uvwInfo.rmatpr, x, y, z, lmst, &xprec, &yprec, &zprec, lmst2000);
                
                PRINTF_DEBUG("DEBUG : HA2000 = %.8f (hardcoded = 6.1717), LMST2000 = %.8f , LMST = %.8f\n",ha2000,lmst2000,lmst);

                CalcUVW(ha2000, dec_aber, xprec, yprec, zprec, u, v, w);
                PRINTF_DEBUG("DEBUG : (x,y,z) = (%.4f,%.4f,%.4f), ha2000 = %.4f , lmst = %.8f , dec_aber = %.4f, lmst2000 = %.8f -> (u,v,w) = (%.8f,%.8f,%.8f)\n",x,y,z,ha2000,lmst,dec_aber,lmst2000,u,v,w);

                if(!printed){
                    for(int i=0;i<3;i++){
                       for(int j=0;j<3;j++){
                          PRINTF_DEBUG("uvwInfo.rmatpr[%d][%d] = %.8f\n",i,j,uvwInfo.rmatpr[i][j]);
                       }
                    }
                    printed=true;
                }
        }

         
        // tested and working but hardcoded for obsID = 1103645160 , /media/msok/625ace3c-584e-49a0-b773-856d6fb8526f/pacer/imager/MWA/HydraA/1103645160/my/20220723/PACER/COTTER
        // see testing in /home/msok/Desktop/PAWSEY/PaCER/logbook/20220924_fixing_uvw_calculation.odt
        static void CalcUVW( double ha2000, double lmst, double dec_aber, double lmst2000, double x, double y, double z, double &u, double &v, double &w )
        {
                double xprec, yprec, zprec;
                double arrayLattitudeRad = -26.703319*(M_PI/180.00); // WARNING : sumulation is not for the MRO ???!!
//                double arrayLattitudeRad = -23.02*(M_PI/180.00); // alsma test
                double rmatpr[3][3];
                
                 
                 
                 // CASA simulation data:
                 double raHrs = 9.30158333;
                 double decDeg = -12.093469;
                 double mjd = 56798.40171;
                 
                 // obsID = 1103645160 : /media/msok/625ace3c-584e-49a0-b773-856d6fb8526f/pacer/imager/MWA/HydraA/1103645160/my/20220723/PACER/COTTER
                 // see also hard-coding / TODO in observation_metadata.cpp line 278
                 raHrs =  9.30160000;
                 decDeg = -12.0956;
                 mjd = 57017.67064814815;
//                 ha2000 = -20.89187500;                 
                 
                 double ra_aber = (raHrs*15.00)*(M_PI/180.00);
                 dec_aber = decDeg*(M_PI/180.00);
//                 double mjd = 56798.40171296296296296296;
                 double ra_app, dec_app;
                 
                 palMap(raHrs*(M_PI/12.0), decDeg*(M_PI/180.0), 0.0, 0.0, 0.0, 0.0, 2000.0, mjd, &ra_app, &dec_app);

                 
                 aber_radec_rad(2000.0, mjd, raHrs*(M_PI/12.0), decDeg*(M_PI/180.0), &ra_aber,&dec_aber);

                 
                 double rmattr[3][3];
                 palPrenut(2000.0, mjd, rmattr);
                 
                 mat_transpose(rmattr, rmatpr);
                 
                 double newarrlat; // ha2000, lmst2000;
                 ha_dec_j2000( rmatpr,lmst,arrayLattitudeRad,ra_aber,dec_aber,&ha2000,&newarrlat,&lmst2000);
                 printf("DEBUG : HA2000 = %.8f , LMST2000 = %.8f , LMST = %.8f , raHrs = %.8f , decDeg = %.8f , ra_app = %.8f , dec_app = %.8f, mjd = %.8f\n",ha2000,lmst2000,lmst,raHrs,decDeg,ra_app,dec_app,mjd);
                

                // TODO : This rotation matrix needs a proper fix !!!
                // palPrenut(2000.0, mjd, rmattr); 
                // mat_transpose(rmattr, uvwInfo.rmatpr);
                // see cotter/geometry.h PrepareTimestepUVW in /home/msok/mwa_software/Andre/anoko/cotter

/*                rmatpr[0][0] = 0.99999324;
                rmatpr[0][1] = 0.00337311;
                rmatpr[0][2] = 0.00146534;
                rmatpr[1][0] = -0.00337304;
                rmatpr[1][1] = 0.99999431;
                rmatpr[1][2] = -0.00004846;
                rmatpr[2][0] = -0.00146549;
                rmatpr[2][1] = 0.00004352;
                rmatpr[2][2] = 0.99999893;*/
                
// -----------------------------------------------------------------------------------------------------
                precXYZ(rmatpr, x, y, z, lmst, &xprec, &yprec, &zprec, lmst2000);
// temprorarily :
//                xprec = x;
//                yprec = y;
//                zprec = z;
// -----------------------------------------------------------------------------------------------------                
                CalcUVW(ha2000, dec_aber, xprec, yprec, zprec, u, v, w);
                
                PRINTF_DEBUG("DEBUG : (x,y,z)=(%.4f,%.4f,%.4f) , ha2000 = %.4f , lmst = %.8f , dec_aber = %.4f , lmst2000 = %.8f -> (u,v,w) = (%.4f,%.4f,%.4f)\n",xprec,yprec,zprec,ha2000,lmst,dec_aber,lmst2000,u,v,w);
        }

        static void CalcUVW(UVWTimestepInfo &uvwInfo, double x, double y, double z, double &u, double &v, double &w)
        {
                double ha2000 = uvwInfo.ha2000;
                double lmst = uvwInfo.lmst;
                double dec_aber = uvwInfo.dec_aber;
                double lmst2000 = uvwInfo.lmst2000; // local mean sidreal time - to be calculated 
                double xprec, yprec, zprec;

// -----------------------------------------------------------------------------------------------------
// TODO :                 precXYZ(uvwInfo.rmatpr, x, y, z, lmst, &xprec, &yprec, &zprec, lmst2000);
// temprorarily :
                xprec = x;
                yprec = y;
                zprec = z;
// -----------------------------------------------------------------------------------------------------                
               
                CalcUVW(ha2000, dec_aber, xprec, yprec, zprec, u, v, w);
        }

        // as in TMS chapter 4.1 and Figure 4.2 (see also 3.1 and 3.2 for explanations why this is done this way) :
        static void CalcUVW(double ha, double dec, double x, double y, double z, double &u, double &v, double &w)
        {
           const double sh = sin(ha), sd = sin(dec);
           const double ch = cos(ha), cd = cos(dec);
           u  = sh*x + ch*y;
           v  = -sd*ch*x + sd*sh*y + cd*z;
           w  = cd*ch*x  - cd*sh*y + sd*z;
        }
        
        static void mat_transpose(double rmat1[3][3], double rmat2[3][3])
        {
                int i, j;

                for(i=0;i<3;++i) {
                        for(j=0;j<3;++j) {
                                rmat2[j][i] = rmat1[i][j];
                        }
                }
        }

        static void aber_radec_rad(double eq, double mjd, double ra1, double dec1, double *ra2, double *dec2)
        {
                double v1[3], v2[3];

                palDcs2c(ra1, dec1, v1);
                stelaber(eq, mjd, v1, v2);
                palDcc2s(v2, ra2, dec2);
                *ra2 = palDranrm(*ra2);
        }

        static void stelaber(double eq, double mjd, double v1[3], double v2[3])
        {
                double amprms[21], v1n[3], v2un[3], w, ab1, abv[3], p1dv;
                int i;

                palMappa(eq,mjd,amprms);

                /* code from mapqk.c (w/ a few names changed): */

                /* Unpack scalar and vector parameters */
                ab1 = amprms[11];
                for ( i = 0; i < 3; i++ )
                {
                                abv[i] = amprms[i+8];
                }

                palDvn ( v1, v1n, &w );

                /* Aberration (normalization omitted) */
                p1dv = palDvdv ( v1n, abv );
                w = 1.0 + p1dv / ( ab1 + 1.0 );
                for ( i = 0; i < 3; i++ ) {
                                v2un[i] = ab1 * v1n[i] + w * abv[i];
                }
                /* normalize  (not in mapqk.c */
                palDvn ( v2un, v2, &w );
        }


        static void ha_dec_j2000(double rmat[3][3], double lmst, double lat_rad, double ra2000,
                  double dec2000, double *newha, double *newlat, double *newlmst)
        {
                double nwlmst, nwlat;

                rotate_radec(rmat,lmst,lat_rad,&nwlmst,&nwlat);
                *newlmst = nwlmst;
                *newha = palDranrm(nwlmst - ra2000);
                *newlat = nwlat;
        }
        
        /* rmat = 3x3 rotation matrix for going from one to another epoch                                      
        * ra1, dec1, ra2, dec2 are in radians                                                                 
        */
        static void rotate_radec(double rmat[3][3], double ra1, double dec1,
                                                double *ra2, double *dec2)
        {
                double v1[3], v2[3];

                palDcs2c(ra1,dec1,v1);
                palDmxv(rmat,v1,v2);
                palDcc2s(v2,ra2,dec2);
                *ra2 = palDranrm(*ra2);
        }



        /**
        * lmst, lmst2000 are the local mean sidereal times in radians                                         
        * for the obs. and J2000 epochs.                                                                      
        */
        static void precXYZ(double rmat[3][3], double x, double y, double z, double lmst,
                                        double *xp, double *yp, double *zp, double lmst2000)
        {
                double sep = sin(lmst);
                double cep = cos(lmst);
                double s2000 = sin(lmst2000);
                double c2000 = cos(lmst2000);

                /* rotate to frame with x axis at zero RA */
                double xpr = cep*x - sep*y;
                double ypr = sep*x + cep*y;
                double zpr = z;

                double xpr2 = (rmat[0][0])*xpr + (rmat[0][1])*ypr + (rmat[0][2])*zpr;
                double ypr2 = (rmat[1][0])*xpr + (rmat[1][1])*ypr + (rmat[1][2])*zpr;
                double zpr2 = (rmat[2][0])*xpr + (rmat[2][1])*ypr + (rmat[2][2])*zpr;

                /* rotate back to frame with xp pointing out at lmst2000 */
                *xp = c2000*xpr2 + s2000*ypr2;
                *yp = -s2000*xpr2 + c2000*ypr2;
                *zp = zpr2;
        }

        /* constants for WGS84 Geoid model for the Earth */
        #define EARTH_RAD_WGS84 6378137.0 /* meters */
        #define E_SQUARED 6.69437999014e-3

        /**
                * Convert Geodetic lat/lon/height to XYZ coords
                * uses constants from WGS84 system
                */
                
        // see also : /home/msok/bighorns/software/analysis/scripts/python$ ant_local2earth.py
        //  https://geodesy.noaa.gov/cgi-bin/xyz_getxyz.prl
        //  /home/msok/Desktop/PAWSEY/PaCER/logbook/20221126_image_simulation_FINAL.odt                
        static void Geodetic2XYZ(double lat_rad, double lon_rad, double height_meters, double &X, double &Y, double &Z) {
                double s_lat,s_lon,c_lat,c_lon;
                double chi;

                s_lat = sin(lat_rad); c_lat = cos(lat_rad);
                s_lon = sin(lon_rad); c_lon = cos(lon_rad);
                chi = sqrt(1.0 - E_SQUARED*s_lat*s_lat);

                X = (EARTH_RAD_WGS84/chi + height_meters)*c_lat*c_lon;
                Y = (EARTH_RAD_WGS84/chi + height_meters)*c_lat*s_lon;
                Z = (EARTH_RAD_WGS84*(1.0-E_SQUARED)/chi + height_meters)*s_lat;
        }



};

#endif
