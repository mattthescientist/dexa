// PathFinder.cpp
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
#include "PathFinder.h"
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

// The basic premise of its operation is as follows. The user firstly, in some 
// higher level class or function, must call the 'feffInp' and 'pathsDat' class
// functions here. In each case, the user specifies the path to FEFF's feff.inp
// and paths.dat files respectively. Data pertaining to the position of atoms
// within the studied system and to each scattering path are then extracted and
// stored in the instance of this class.
//
// It doesn't matter which of these two functions is called first, but once both
// have been called, 'definePaths' will activate the path finding algorithm
// found in 'getLegs'. This will then search for and store every scattering path
// in the crystal that matches one of the forms listed in the paths.dat file,
// thus breaking the degeneracy of the scattering paths for any given shell so
// that every individual scattering path is know explictly.
//
// Take care: 'getLegs' is a RECURSIVE function. It works by exploiting the fact
// that any given scattering path may be successively broken down into a set of
// smaller paths and scattering legs. For example,
//   A tripe scattering path  = A double scattering path plus one leg
//   A double scattering path = A single scattering path plus one leg
//   A single scattering path = One scattering leg plus one leg
//
// With this in mind, 'definePaths' will call 'getLegs' and specify the number
// of legs, the path length, and the scattering angles and potentials for the
// sought paths.It will also state that the source and destination atoms for the
// current leg are zero, thus telling 'getLegs' to start at the absorbing atom. 
//
// 'getLegs' will then cycle through all the other atoms specified in the 
// feff.inp and see if the following conditions are met
//   1. Adding it to the path will not cause the number of legs to be exceeded.
//   2. Adding it to the path will not cause the path length to be exceeded.
//   3. The angle between the previous and proposed legs matches that listed
//      from paths.dat
//   4. The proposed scattering potential matches that listed from paths.dat
//
// If all four conditions are met, the leg formed between this and the source
// atom is added to a stack. 'getLegs' is then called again with the following
// changes
//   1. The number of legs left to find in the path is decremented by one.
//   2. The length of this new leg is subtracted from the total path length. 
//   3. The source atom set to what is currently the destination atom i.e. the
//      one at the end of this new leg.
//
// the whole process is thus repeated recursively until one of the four 
// conditions is not met. If at any time, condition 3 or 4 are violated, or if
// adding the proposed leg to the path exceeds the total path length, the 
// leg is automatically invalid and so discarded. 
// 
// Additionally, if all target atoms are searched, but none found to be 
// satisfactory, the source atom itself must be invalid (since the path cannot
// be completed) and so it too is deleted.
//
// If, however, the number of legs remaining equals zero AND the path length is
// also zero (within some tolerance), the proposed path looks good. A final
// check is performed to make sure that the last leg returns to the absorbing
// atom. If this is true, the proposed path is one of the ones sought and so 
// kept in the stack, which, at the end of the day, becomes a linked list.

// There is no functionality to read out the paths once they have been found.
// The idea is for this class to instead be inherited by some larger class that
// needs to use the path finder algorithm. That class will then have direct
// access to the paths and thus be able to manipulate them as required.

using namespace::std;

template <class T> bool fromString(T& t, const std::string& s)
{
  std::istringstream iss(s); return !(iss >> t).fail();
}


//------------------------------------------------------------------------------
// feffInp (char *) : Reads FEFF's feff.inp file and extracts the POTENTIALS and
// ATOMS sections, which are stored in class vectors Potentials and Atoms 
// respectively.
//
int PathFinder::feffInp (const char *Filename, int Verbose) {
  string theString, Temp; istringstream iss;
  string NewPotential; int PotentialIndex;
  ScatterCoord NewAtom;

  // Open the feff.inp specified in Filename. Abort if it fails to open.
  ifstream FeffInp (Filename, ios::in);  
  if (! FeffInp.is_open ()) {
    return PF_ERROR;
  }

  // Find the POTENTIALS and ATOMS sections of the feff.inp. These flags must be
  // in the first field of a line to be valid, so search only first fields.
  int Mode = PF_MODE_NONE;
  while (!FeffInp.eof()) { 
    getline (FeffInp, theString);
    iss.str (theString);
    iss >> theString;

    // Check to see if we've reached a new section. If a new section is found,
    // set theString to "*" so that the following section will treat it as a
    // comment and so avoid trying to extract data from the section header.
    if (theString == "POTENTIALS") { Mode = PF_MODE_POTENTIALS; theString ="*";}
    if (theString == "ATOMS") { Mode = PF_MODE_ATOMS; theString = "*"; }
    if (theString == "END") { Mode = PF_MODE_END; theString = "*"; }
    
    // Get data, avoiding comments, blank lines, and unwanted info respectively
    if (theString != "*" && !iss.eof () && Mode != PF_MODE_NONE) {
      switch (Mode) {
        case PF_MODE_POTENTIALS:
          iss >> Temp >> NewPotential; Potentials.push_back (NewPotential);
          break;
        case PF_MODE_ATOMS:
          fromString (NewAtom.x, theString);
          iss >> NewAtom.y >> NewAtom.z >> PotentialIndex;
          NewAtom.Potential = PotentialIndex;
          Atoms.push_back (NewAtom);
          OriginalAtoms.push_back (NewAtom);
          break;
      }
    }
    iss.clear ();
  }
  FeffInp.close ();

  // Using this list, and data from paths.dat, find necessary scattering paths
  definePaths (Verbose);

  return PF_NO_ERROR;
}


//------------------------------------------------------------------------------
// pathsDat (char *) : Reads FEFF's paths.dat file and extracts all necessary
// information concerning the scattering paths FEFF deemed to be significant.
//
// Each scattering path is read into a 'Shell'. Remember however that the 
// paths.dat only contains one varient of each type of scattering path. The 
// other, degenerate paths must be found maunally and added to the shell later.
//
// The following path properties are stored from each path's header:
//   ChipIndex  : the number needed to link the path to its chipnnnn.dat file.
//   Legs       : the number of scattering legs in the current path.
//   Radius     : Path radius. For MS paths, this is half the total path length.
//   Degeneracy : Number of scattering paths that are similar to the one listed.
// 
// Within each path, the following info is stored from each atom:
//   x, y, z    : Position of the atom with respect to the beamline frame.
//   Potential  : The atomic species off which the photoelectron is scattering.
//   Beta       : The scattering angle between the current leg and the next.
//
int PathFinder::pathsDat (const char *Filename, int Verbose) {
  string theString, Temp1, Temp2, Temp3, Temp4; 
  istringstream iss;
  Shell NewShell;
  ScatterCoord NewAtom;
  vector <ScatterCoord> *NewPath;
  int i;
  float FDegeneracy;

  // Open the paths.dat file. Abort with PF_ERROR if it fails to open.
  ifstream PathsDat (Filename, ios::in);  
  if (! PathsDat.is_open ()) {
    return PF_ERROR;
  }

  // Ditch the paths.dat header. Since different versions of FEFF have different
  // header sizes, it's not possible to simply discard x lines. Therefore,
  // search for the '-----' line at the end of the header.
  do {
    getline (PathsDat, theString);
  } while (theString[1] != '-');

  // Read each scattering path in turn.
  while (!PathsDat.eof ()) {
    getline (PathsDat, theString);
    iss.str (theString);
    iss >> NewShell.ChipIndex;
    if (!iss.eof ()) {    // Ensures we haven't read a blank line i.e. near EOF
    
      // Get the path header information
      iss >> NewShell.Legs >> FDegeneracy >> Temp1
          >> Temp2 >> Temp3 >> Temp4 >> NewShell.Radius;
      if (NewShell.Radius > MaxR) MaxR = NewShell.Radius;
      if (NewShell.Legs > MaxLegs) MaxLegs = NewShell.Legs;
      NewShell.Degeneracy = int(FDegeneracy);
      getline (PathsDat, theString);    // Path description line

      // Get data for each atom in the scattering sequence
      NewPath = new vector <ScatterCoord>;
      for (i = 0; i < NewShell.Legs; i ++) {
        getline (PathsDat, theString);
        iss.clear (); iss.str (theString);
        iss >> NewAtom.x >> NewAtom.y >> NewAtom.z >> NewAtom.Potential >> Temp1
            >> Temp2 >> Temp3 >> NewAtom.Beta;
        NewPath -> push_back (NewAtom);
      }

      // Store this shell and path variant in the Shells vector
      NewShell.Scatterers = NewPath;
      Shells.push_back (NewShell);
    }
  }

  // Now search for all varients of the scattering paths with aid of feff.inp
  definePaths (Verbose);
  
  return PF_NO_ERROR;
}

//------------------------------------------------------------------------------
// getLegs (int, int, int, double, double, int) : This routine basically IS the
// path finding algorithm. It is a RECURSIVE routine, so be careful when trying
// to read through it. Read the notes on its operation at the top of this file.
//
vector <Leg*>* PathFinder::getLegs (int i, int j, int LegsLeft, double S, double *Angle, int *Potential){
  vector <Leg*> *NewLegs = new vector <Leg*>;
  Leg *NewLeg;  

  // First update the variables that are tracking our current position within
  // the scattering path as a whole. Since we have moved to a new leg, we must:
  //   1. Decrement the counter indicating how many legs are left to find.
  //   2. Subtract the length of the leg from the total length of the path.
  //   3. Move to the next scattering angle in the path
  //   4. Move to the next scattering potential in the path
  LegsLeft --;
  Angle ++;
  Potential ++;
  S -= vectorLength (Atoms[j].x - Atoms[i].x, Atoms[j].y - Atoms[i].y,
                     Atoms[j].z - Atoms[i].z);

  // If we still have some path legs remaining AND S is still greater than zero
  if (LegsLeft > 0 && S > PF_LENGTH_TOL) {
  
    // Split the Atoms search into two 'for' loops to avoid legs of zero length
    // i.e. the physically meaningless leg between any atom and itself.
    for (int k = 0; k < j; k ++) {

      // Only accept the proposed leg if the scattering angle and potential are
      // those that are required.
      if (abs (*Angle - vectorAngle (i, j, k)) < PF_ANGLE_TOL 
        && Atoms[k].Potential == *Potential) {

        // The leg looks good, so remember its index, find all its child legs
        // and then, if it has child legs, store it in the NewLegs linked list.
        NewLeg = new Leg;
        NewLeg -> Atom = j;
        if ((NewLeg -> Next = 
          getLegs (j, k, LegsLeft, S, Angle, Potential)) != NULL) {
          NewLegs -> push_back (NewLeg);
        } else { delete NewLeg; }
      }
    }
    for (unsigned int k = j + 1; k < Atoms.size(); k ++) {
      if (abs (*Angle - vectorAngle (i, j, k)) < PF_ANGLE_TOL 
        && Atoms[k].Potential == *Potential) {
        NewLeg = new Leg;
        NewLeg -> Atom = j;
        if ((NewLeg -> Next = 
          getLegs (j, k, LegsLeft, S, Angle, Potential)) != NULL) {
          NewLegs -> push_back (NewLeg);
        } else { delete NewLeg; }
      }
    }

    // See if we have found some child legs. If so, return NewLegs. If not we're
    // at the end of the path. Return NULL so that the Next pointer of the
    // previous atom becomes the end of the path.
    if (NewLegs -> size () == 0) {
      delete NewLegs;
      return NULL;
    }
    return NewLegs;

  // LegsLeft = 0 and S is zero (or near to it), so the path seems good.
  // Just make sure that the last atom in the path is the absorbing atom.
  } else if (LegsLeft == 0 && abs (S) <= PF_LENGTH_TOL && j == 0) {
    NewLeg = new Leg;
    NewLeg -> Atom = j;
    NewLeg -> Next = NULL;
    NewLegs -> push_back (NewLeg);
    return NewLegs;
  
  // Finally, either LegsLeft = 0 or S is zero (or near to it) but not both of 
  // them, or the last atom in the proposed path is not the absorbing atom. 
  // Whatever the case, the path is invalid so delete it.
  } else {
    delete NewLegs;
    return NULL;
  }
}


//------------------------------------------------------------------------------
// vectorAngle (int, int, int) : Returns the angle between atoms i, j, and k,
// which are stored in the Atoms array. The fudge-factor of 1e-12 is present
// here to avoid a curious error caused when the machine rounds off a number 
// ever so slightly greater than 1 in the argument of acos; a condition that
// will occur often on a backscattering event where beta = 180deg. This
// cannot be caught with an 'if (x > 1.0)' type statement since the 'if' sees
// x as being exactly 1 where as acos sees x as slightly greater than 1. Strange
//
double PathFinder::vectorAngle (int i, int j, int k) {
  double x1 = Atoms[j].x - Atoms[i].x;
  double y1 = Atoms[j].y - Atoms[i].y;
  double z1 = Atoms[j].z - Atoms[i].z;
  double x2 = Atoms[k].x - Atoms[j].x;
  double y2 = Atoms[k].y - Atoms[j].y;
  double z2 = Atoms[k].z - Atoms[j].z;
  double l1 = vectorLength (x1, y1, z1);
  double l2 = vectorLength (x2, y2, z2);

  if (l1 > 0.0 && l2 > 0.0) {
    return acos ((x1 * x2 + y1 * y2 + z1 * z2) / ((l1 * l2)+1e-12)) / PI *180.0;
  } else {
    return 0;
  }
}

//------------------------------------------------------------------------------
// countDegeneracy (vector <Leg*> *) : Walks through the specifed list of 
// scattering legs and counts the path degeneracy.
//
int PathFinder::countDegeneracy (vector <Leg*> *Leg) {
  int PathCount = 0;
  if (Leg == NULL) {
    return 1;
  } else {
    for (unsigned int i = 0; i < Leg -> size (); i ++) {
      PathCount += countDegeneracy (Leg -> at(i) -> Next);
    }
    return PathCount;
  }
}

//------------------------------------------------------------------------------
// deleteAllLegs (vector <Leg*> *) : Walks through the specifed scattering path 
// and deletes all its legs. This is called by the class destructor to clean up
// the memory before it too is deleted.
//
void PathFinder::deleteAllLegs (vector <Leg*> *Leg) {
  for (unsigned int i = 0; i < Leg -> size (); i ++) {
    if (Leg -> at(i) -> Next != NULL) {
      deleteAllLegs (Leg -> at(i) -> Next);
    }
    delete Leg -> at(i);
  }
  delete Leg;
}

//------------------------------------------------------------------------------
// definePaths () : Goes through each path type found in paths.dat and generates
// arrays containing the scattering angles and potentials before calling 
// 'getLegs'. If the path finder in 'getLegs' returns a different number of
// of paths to that specified by the degeneracy in paths.dat, a warning is
// printed on stdout, and the degeneracy changed to that found by 'getLegs'.
//
int PathFinder::definePaths (int Verbose) {

  // First make sure that both the feff.inp and paths.dat have been read. If
  // not return PF_NOT_READY.
  if (!Shells.size () || !Atoms.size ()) { return PF_NOT_READY; }
  double *Angles, *AnglesReverse;  
  int *thePotentials, *PotentialsReverse;
  int Degeneracy;

  // Go through each scattering path type based on the information in Shells
  // that was obtained from the paths.dat file.
  for (unsigned int i = 0; i < Shells.size (); i ++) {

    // Create a list from the scattering angles for the current path.
    // Angles[0] and Angles[1] MUST be zero to handle the special condition
    // generated by the very first call to 'getLegs'.
    Angles = new double [Shells[i].Scatterers -> size() + 1];
    AnglesReverse = new double [Shells[i].Scatterers -> size() + 1];
    Angles [0] = 0.0; AnglesReverse [0] = 0.0;
    Angles [1] = 0.0; AnglesReverse [1] = 0.0;
    for (unsigned int j = 0; j < Shells[i].Scatterers -> size() - 1; j ++) {
      Angles [j + 2] = Shells[i].Scatterers -> at(j).Beta;
      AnglesReverse [Shells[i].Scatterers -> size() - j] = Angles [j + 2];
    }

    // Create a list from the scattering potentials for the current path.
    // Potentials[0] MUST be zero to handle the special condition
    // generated by the very first call to 'getLegs'.
    thePotentials = new int [Shells[i].Scatterers -> size() + 1];
    PotentialsReverse = new int [Shells[i].Scatterers -> size() + 1];
    thePotentials [0] = 0; PotentialsReverse [Shells[i].Scatterers->size()] = 0;
    for (unsigned int j = 0; j < Shells[i].Scatterers -> size(); j ++) {
      thePotentials [j + 1] = Shells[i].Scatterers -> at(j).Potential;
      PotentialsReverse [Shells[i].Scatterers -> size() - 1 - j] 
        = thePotentials[j + 1];
    }

    // Print some useful information to stdout and then call 'getLegs'. Once
    // 'getLegs' has finished, push the list of Legs into the Paths vector. 
    if (Verbose >= PF_MED_VERBOSE) {
      cout << "Finding legs for path " << i + 1 << " (";
      for (unsigned int j = 0; j < Shells[i].Scatterers -> size (); j ++) {
        cout << Potentials[Shells[i].Scatterers -> at (j).Potential];
        if (j < Shells[i].Scatterers -> size () - 1) { cout << "->"; }
      }
      cout << ") ... " << flush;
    }
    
    vector <Leg*> *thePaths = getLegs 
      (0, 0, Shells[i].Legs + 1, Shells[i].Radius * 2.0, Angles, thePotentials);
    Paths.push_back (thePaths);
    Degeneracy = countDegeneracy (Paths[Paths.size() - 1]);    

    // For double and higher order scattering paths, we must go around the path
    // in both directions since each yields unique information. The only 
    // exception are co-linear paths, which in some cases are identical in both
    // directions.
    if (Shells[i].Degeneracy != Degeneracy && Shells[i].Scatterers->size() > 2){
      if (Shells[i].Scatterers -> size () % 2 == 0) {
        for (unsigned int j = 0; j < Shells[i].Scatterers -> size() + 1; j++) {
          if (Angles[j] != 0.0 && Angles[j] != 180.0) {
            vector <Leg*> *thePathsReverse = getLegs (0, 0, Shells[i].Legs + 1, 
              Shells[i].Radius * 2.0, AnglesReverse, PotentialsReverse);
            for (unsigned int k = 0; k < thePathsReverse -> size (); k ++) {
              Paths[Paths.size () - 1] -> push_back (thePathsReverse -> at (k));
            }
            break;
          }
        }
      } else {
        vector <Leg*> *thePathsReverse = getLegs (0, 0, Shells[i].Legs + 1, 
          Shells[i].Radius * 2.0, AnglesReverse, PotentialsReverse);
        for (unsigned int k = 0; k < thePathsReverse -> size (); k ++) {
          Paths[Paths.size () - 1] -> push_back (thePathsReverse -> at (k));
        }
      }
    }

    // Finally, check that the number of paths found by 'getLegs' matches the
    // degneracy specified in the paths.dat file. If the don't print a warning
    // to stdout and update the degeneracy to match that from 'getLegs'.
    Degeneracy = countDegeneracy (Paths[Paths.size() - 1]);
    if (Verbose >= PF_MED_VERBOSE) 
      cout << "done" << endl;
    if (Degeneracy != Shells[i].Degeneracy) {
      cout << "   Warning: degeneracy is " << Degeneracy 
           << " but paths.dat says it should be " << Shells[i].Degeneracy <<
            endl;
      Shells[i].Degeneracy = Degeneracy;
    }
    delete [] Angles;
    delete [] thePotentials;
    delete [] AnglesReverse;
    delete [] PotentialsReverse;
  }
  return PF_NO_ERROR;
}
