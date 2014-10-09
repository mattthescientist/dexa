// PathFinder.h
// Copyright (C) 2007-2009 Matthew Ruffoni
//
// This file is part of DEXA.
//
// DEXA is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DEXA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with DEXA.  If not, see <http://www.gnu.org/licenses/>.
//
#ifndef PF_PATHFINDER_H
#define PF_PATHFINDER_H

#define PF_NO_ERROR   0
#define PF_ERROR      1
#define PF_NOT_READY  2

#define PF_MODE_NONE        0
#define PF_MODE_POTENTIALS  1
#define PF_MODE_ATOMS       2
#define PF_MODE_END         3

#define PF_LENGTH_TOL       1e-4  // +/- error in total path length
#define PF_ANGLE_TOL        1e-4  // +/- error in any given scattering angle

#define PF_MED_VERBOSE 2

#ifndef PI
  #define PI (4.0 * atan (1.0))
#endif

#include <vector>
#include <string>
#include <cmath>
#include <iostream>

//------------------------------------------------------------------------------
// TypeDefs

typedef struct PF_Scatterer {
  double x;
  double y;
  double z;
  double Beta;              // scattering angle
  int Potential;            // index to correlate feff.inp and paths.dat atoms
} ScatterCoord;

typedef struct PF_Leg {
  int Atom;
  std::vector <PF_Leg*> *Next;
  PF_Leg () { Next = NULL; Atom = 0; }
} Leg;

typedef struct PF_Shell {
  std::vector <ScatterCoord> *Scatterers;
  
  int ChipIndex;            // eg. 1 corresponds to chip0001.dat, 2 to 0002 etc.
  int Legs;                 // nleg   Number of legs in the scattering path
  int Degeneracy;           // degen  Number of equivalent paths in the shell
  double Radius;            // r      Path Radius (1/2 path length for nleg > 2)
  PF_Shell () {
	  Scatterers = NULL;
  	  ChipIndex = 0;
  	  Legs = 0;
  	  Degeneracy = 0;
  	  Radius = 0.0;
  }
} Shell;

//------------------------------------------------------------------------------
// PathFinder Class

class PathFinder {

public:
  PathFinder () { MaxLegs = 0; MaxR = 0.0; };    // Empty default constructor
  ~PathFinder () {
    for (unsigned int i = 0; i < Shells.size(); i ++) { 
      if (Shells[i].Scatterers != NULL) delete Shells[i].Scatterers;
      Shells.clear ();
    }
    for (unsigned int i = 0; i < Paths.size(); i ++) { 
      deleteAllLegs(Paths[i]);
      Paths.clear ();
    }
  }
  
  // Public Commands
  int feffInp (const char *FeffInpName, int Verbose);
  int pathsDat (const char *PathsDatName, int Verbose);
    
  // Public GET Routines
  std::vector <Shell> getShells () { return Shells; }
  std::vector <ScatterCoord> getAtoms () { return Atoms; }
  
private:
  int definePaths (int Verbose);
  std::vector <Leg*>* getLegs 
    (int i, int j, int LegsLeft, double R, double *Angle, int *Potential);
  int countDegeneracy (std::vector <Leg*> *Leg);
  double vectorAngle (int i, int j, int k);
  void deleteAllLegs (std::vector <Leg*> *Leg);
  
  int MaxLegs;
  double MaxR;

protected:
  double vectorLength (double x, double y, double z) 
    { return pow (x * x + y * y + z * z, 0.5); }
  double abs (double a) { return a < 0 ? (-1.0 * a) : (a); }

  std::vector <std::string> Potentials;
  std::vector < std::vector <Leg*> *> Paths;
  std::vector <Shell> Shells;       // Information from FEFF's paths.dat file
  std::vector <ScatterCoord> Atoms; // Atom type and coords from FEFF's feff.inp
  std::vector <ScatterCoord> OriginalAtoms;
};

#endif // PF_MGSTSPECTRUM_H
