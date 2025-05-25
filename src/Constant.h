#import <Foundation/Foundation.h>
// physical constants
const double ec = 4.803204197e-10;     // elementary charge [statC = esu]
const double c = 2.99792458e10;        // speed of light [cm/s]
const double kb = 1.3806503e-16;       // Boltzmann's constant [erg/K]
// unit transform
const double evtoerg = 1.6021764262e-12;       // from eV     to  erg
const double ergtoev = 1/1.6021764262e-12;     // from erg    to  eV
const double ergtoj = 1e-7;                    // from erg    to  Joule
const double jtoerg = 1/1.0e-7;                // from Joule  to  erg
const double ktoev = 1/11604;                  // from Kelvin to  eV
const double evtok = 11604;                    // from eV     to  Kelvin
const double patoatm = 1/101325;               // from Pa     to  atm
const double atmtopa = 101325;                 // from atm     to  Pa
const double VtoG = 1e6/c;                     // V/m to cgs-gauss [V/m = 1/(299*100) statV/cm]
const double GtoV = 1/(1e6/c);                 // convert cgs-gauss to V/m [1 statV/cm = 299*100 V/m]
const double VtosV = 1e8/c;                    // convert V to cgs-gauss [V = 1/(299) statV]
const double sVtoV = 1/(1e8/c);                // convert cgs-gauss to V [1 statV = 299 V]
const double AtosA = 1/(10/c);                 // convert A to cgs-gauss [1 A = 1/(10/c) statA]
const double sAtoA = 10/c;                     // convert cgs-gauss to A [1 statA = 10/c A]
const double JtosJ = 1/(10/c*1e4);             // convert A/m^2 to cgs-gauss (esu/s/cm^2) [1 A = 1/(10/c) statA]
const double sJtoJ = 10/c*1e4;                 // convert cgs-gauss (esu/s/cm^2) to A/m^2 [1 statA = 10/c A]
const double TtoG = 1e4;                       // convert T to to cgs-gauss G [1 G = 10^-4 T]
const double GtoT = 1e-4;                      // convert T to to cgs-gauss G [1 G = 10^-4 T]