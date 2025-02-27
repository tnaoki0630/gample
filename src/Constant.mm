#import "Constant.h"
// physical constants
double ec = 4.803204197e-10;     // elementary charge [statC = esu]
double c = 2.99792458e10;        // speed of light [cm/s]
double kb = 1.3806503e-16;       // Boltzmann's constant [erg/K]
// unit transform
double evtoerg = 1.6021764262e-12;       // from eV     to  erg
double ergtoev = 1/1.6021764262e-12;     // from erg    to  eV
double ergtoj = 1e-7;                    // from erg    to  Joule
double jtoerg = 1/1.0e-7;                // from Joule  to  erg
double ktoev = 1/11604;                  // from Kelvin to  eV
double evtok = 11604;                    // from eV     to  Kelvin
double patoatm = 1/101325;               // from Pa     to  atm
double atmtopa = 101325;                 // from atm     to  Pa
double VtoG = 1e6/c;                     // V/m to cgs-gauss [V/m = 1/(299*100) statV/cm]
double GtoV = 1/(1e6/c);                 // convert cgs-gauss to V/m [1 statV/cm = 299*100 V/m]
double VtosV = 1e8/c;                    // convert V to cgs-gauss [V = 1/(299) statV]
double sVtoV = 1/(1e8/c);                // convert cgs-gauss to V [1 statV = 299 V]
double AtosA = 1/(10/c);                 // convert A to cgs-gauss [1 A = 1/(10/c) statA]
double sAtoA = 10/c;                     // convert cgs-gauss to A [1 statA = 10/c A]
double JtosJ = 1/(10/c*1e4);             // convert A/m^2 to cgs-gauss (esu/s/cm^2) [1 A = 1/(10/c) statA]
double sJtoJ = 10/c*1e4;                 // convert cgs-gauss (esu/s/cm^2) to A/m^2 [1 statA = 10/c A]
double TtoG = 1e4;                       // convert T to to cgs-gauss G [1 G = 10^-4 T]
double GtoT = 1e-4;                      // convert T to to cgs-gauss G [1 G = 10^-4 T]