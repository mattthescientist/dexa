// ScriptEventsGeneral.cpp
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
//------------------------------------------------------------------------------
// ScriptEventsGeneral.cpp : Contains all the handler functions for DEXA script
// commands that are common to both the Minuit and GSL versions of the program.
// Some of these functions are not required when the tensor calculations are 
// removed for analysing amorphous samples. In this case, ScriptReader.h defines
// MGST_AMORPHOUS. Compiler preprocessor directives are then used to detect this
// and omit the unneeded routines.
//

#include "ScriptReader.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

// File name sufficies
#define K_SUFFIX "_KSpace"
#define R_SUFFIX "_RSpace"
#define Q_SUFFIX "_QSpace"

//------------------------------------------------------------------------------
// GLOBAL Command List
// These commands affect the DEXA fitting environment as a whole.
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Command : create <spectrum name>
// Creates a magnetostriction spectrum object with the specified <spectrum name>
// This is then pushed onto the end of the Spectra vector.
//
int create (stringstream &args, vector<struct spectrum> &Spectra) {
  struct spectrum NewSpectrum;
  int Verbose;
  args >> NewSpectrum.Name >> Verbose;
  if (args.eof ()) {
    NewSpectrum.SubStructure[0].Spectrum.setVerbose (Verbose);
    Spectra.push_back (NewSpectrum);
    if (Verbose >= MIN_VERBOSE)
      cout << "Adding " << NewSpectrum.Name << " to list of spectra..." << endl;
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}

//------------------------------------------------------------------------------
// Command : saveparams <file name>
// Saves the current values of all the DEXA fit parameters to the file specified
// at arg1.
int saveparams (stringstream &args, vector<struct spectrum> &Spectra) {
  string Filename;
  int ErrNo, Verbose;
  
  args >> Filename >> Verbose;
  if (args.eof ()) {
    ofstream ParamFile (Filename.c_str(), ios::out);
    if (! ParamFile.is_open()) {
      cout << "Unable to write parameters to " << Filename 
        << ". Parameters NOT saved" << endl;
      return SCRIPT_ERROR;
    }
    ParamFile.precision (SCRIPT_PRECISION);
    ErrNo = printParams (ParamFile, Spectra);
    ParamFile.close ();
    return ErrNo;
  }
  return SCRIPT_BAD_ARGS;
}

//------------------------------------------------------------------------------
// Command : preforientation <x> <y> <z>
// Sets the preferential orientation of the crystallites within a polycrystal
// sample. <x>, <y>, and <z> must combine to form a unit vector that points in a
// given direction with respect to the beamline frame i.e. +x is horizontally
// towards the machine, +y is up, and +z is along the direction of beam 
// propogation.
//
int preforientation (stringstream &args, vector<struct spectrum> &Spectra) {
  string ExpSpectrum, SpectrumName;
  double x, y, z, Degree;
  int err, Verbose;
  
  args >> x >> y >> z >> Degree >> Verbose;
  if (args.eof () && !args.fail ()) {
    for (unsigned int i = 0; i < Spectra.size (); i ++) {
      for (unsigned int j = 0; j < Spectra[i].SubStructure.size (); j ++) {
        err = Spectra[i].SubStructure[j].Spectrum.addPrefOrientation (x, y, z, Degree);
        if (err != DX_NO_ERROR) {
          cout << "Failed to add a preferential orientation in ("
            << x << ", " << y << ", " << z << ") with a degree of " 
            << Degree * 100.0 << "%" << endl;
          return SCRIPT_ERROR;
        }
      }
    }
    if (Verbose >= MIN_VERBOSE)
      cout << "Crystallite preferential orientation added in (" 
        << x << ", " << y << ", " << z << ") with a degree of " 
        << Degree * 100.0 << "%" << endl;
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}


//------------------------------------------------------------------------------
// Command : magnetisation <x1> <y1> <z1> <x2> <y2> <z2>
// Sets the two sample magnetisation vectors across which the DiffXAS 
// measurement was made. Each set, <x1>, <y1>, <z1>, and <x2>, <y2>, <z2> must
// combine to form a unit vector that points in a given direction with respect 
// to the beamline frame i.e. +x is horizontally towards the machine, +y is up, 
// and +z is along the direction of beam propogation.
//
int magnetisation (stringstream &args, vector<struct spectrum> &Spectra) {
  string ExpSpectrum, SpectrumName;
  double x1, y1, z1, x2 ,y2 ,z2;
  int err, Verbose;
  
  args >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> Verbose;
  if (args.eof ()) {
    for (unsigned int i = 0; i < Spectra.size (); i ++) {
      for (unsigned int j = 0; j < Spectra[i].SubStructure.size (); j ++) {
        err = Spectra[i].SubStructure[j].Spectrum.magnetisation1 (x1, y1, z1);
        if (!err) err = Spectra[i].SubStructure[j].Spectrum.magnetisation2 (x2, y2, z2);
        if (err != DX_NO_ERROR) {
          cout << "Failed to set the magnetisation to "
            << " M1 = (" << x1 << ", " << y1 << ", " << z1 << ")"
            << " and M2 = (" << x2 << ", " << y2 << ", " << z2 << ")";
          return SCRIPT_ERROR;
        }
      }
    }
    if (Verbose >= MIN_VERBOSE)
      cout << "Magnetisation set to"
        << " M1 = (" << x1 << ", " << y1 << ", " << z1 << ")"
        << " M2 = (" << x2 << ", " << y2 << ", " << z2 << ")" << endl;
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}


//------------------------------------------------------------------------------
// Command : polarisation <x> <y> <z>
// Sets the x-ray polarisation vector. This must be a linear polarisation.
// <x>, <y>, and <z> must combine to form a unit vector that points in a given 
// direction with respect to the beamline frame i.e. +x is horizontally towards 
// the machine, +y is up, and +z is along the direction of beam propogation.
//
int polarisation (stringstream &args, vector<struct spectrum> &Spectra) {
  double x, y, z;
  int err, Verbose;
  
  args >> x >> y >> z >> Verbose;
  if (args.eof ()) {
    for (unsigned int i = 0; i < Spectra.size (); i ++) {
      for (unsigned int j = 0; j < Spectra[i].SubStructure.size (); j ++) {
        err = Spectra[i].SubStructure[j].Spectrum.polarisation (x, y, z);
        if (err != DX_NO_ERROR) {
          cout << "Failed to set the polarisation orientation to ("
            << x << ", " << y << ", " << z << ")" << endl;
          return SCRIPT_ERROR;
        }
      }
    }
    if (Verbose >= MIN_VERBOSE)
      cout << "Polarisation set to (" << x <<", " << y <<", " << z <<")" <<endl;
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}


//------------------------------------------------------------------------------
// Command : numavesteps <steps>
// Sets the number of individual orientations to be used when averaging the
// magnetostriction about a defined preferential axis. A larger number will give
// a more accurate average but at the expense of a longer calculation time. This
// function is not needed in the tensor-free version of the code and is
// therefore omitted by the compiler preporocessor.
//
int numavesteps (stringstream &args, vector<struct spectrum> &Spectra) {
  int Steps, err, Verbose;
  
  args >> Steps >> Verbose;
  if (args.eof ()) {
    if (Steps < 1) {
      cout << "ERROR: The number of steps must be a positive, non-zero integer"
        << endl;
      return SCRIPT_ERROR;
    }
    for (unsigned int i = 0; i < Spectra.size (); i ++) {
      for (unsigned int j = 0; j < Spectra[i].SubStructure.size (); j ++) {
        err = Spectra[i].SubStructure[j].Spectrum.numSteps (Steps);
        if (err != DX_NO_ERROR) {
          cout << "Failed to set the number of averaging steps to "
            << Steps << endl;
          return SCRIPT_ERROR;
        }
      }
    }
    if (Verbose >= MIN_VERBOSE)
      cout << "Number of averaging steps set to " << Steps << endl;
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}





//------------------------------------------------------------------------------
// SPECTRUM Command List
// These commands affect only a specified spectrum object. The readScript()
// function in ScriptReader.cpp will pass the name of the spectrum to each of
// these handler functions at the end of the list of arguments.
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Command : addstructure <name>
// Adds a new sub-structure to the specified spectrum with the name <name>.
//
int addstructure (stringstream &args, vector<structure> &Structures) {
  structure NewStructure;
  string SpectrumName;
  int Verbose;
  args >> NewStructure.Name >> SpectrumName >> Verbose;
  if (args.eof ()) {
    NewStructure.Spectrum.setVerbose (Verbose);
    if (Structures[0].Name == "") { Structures[0].Name = NewStructure.Name; }
    else { Structures.push_back (NewStructure); }
    if (Verbose >= MIN_VERBOSE)
      cout << "Adding " << NewStructure.Name << " structure to " << SpectrumName 
      << "..." << endl;
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}

//------------------------------------------------------------------------------
// Command : spectrum.feffinp <feff inp>
// Loads the contents of the FEFF feff.inp file <feff inp> into the specified
// spectrum object.
//
int feffinp (stringstream &args, vector<structure> &Structures) {
  string FeffInp, SpectrumName;
  int ErrNo, Verbose;
  
  args >> FeffInp >> SpectrumName >> Verbose;
  if (args.eof ()) {
    for (unsigned int i = 0; i < Structures.size(); i ++) {
      ErrNo = Structures[i].Spectrum.feffInp (FeffInp.c_str(), Verbose);
      if (ErrNo == PF_ERROR) {
        cout << "Unable to read FEFF data from \"" << FeffInp << "\"." << endl;
        return SCRIPT_ERROR;
      }
    }
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}


//------------------------------------------------------------------------------
// Command : spectrum.pathsdat <paths dat>
// Loads the contents of the FEFF paths.dat file <paths dat> into the specified
// spectrum object.
//
int pathsdat (stringstream &args, vector<structure> &Structures) {
  string PathsDat, SpectrumName;
  int ErrNo, Verbose;
  
  args >> PathsDat >> SpectrumName >> Verbose;
  if (args.eof ()) {
    for (unsigned int i = 0; i < Structures.size(); i ++) {
      ErrNo = Structures[i].Spectrum.pathsDat (PathsDat.c_str(), Verbose);
      if (ErrNo == PF_ERROR) {
        cout << "Unable to read paths data from \"" << PathsDat << "\"." <<endl;
        return SCRIPT_ERROR;
      }
    }
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}


//------------------------------------------------------------------------------
// Command : spectrum.spectrum <filename> <header> <x col> <y col>
// Loads the experimental DiffXAS spectrum in <filename> into the specified
// spectrum object. The first <header> lines of the file are discarded, and then
// x-axis data (energy/k-space etc.) is loaded from the column <x col>, and the
// dChi data from <y col>. If the specified spectrum is divided into sub-
// structures (each of which will have their own MgstSpectrum object in the
// code), the data will be copied to all structures using the MgstSpectrum
// copySpectrum () function.
//
int spectrum (stringstream &args, vector<structure> &Structures) {
  string ExpSpectrum, SpectrumName;
  int Header, x, y, err, Verbose;
  
  args >> ExpSpectrum >> Header >> x >> y >> SpectrumName >> Verbose;
  if (args.eof ()) {
    cout << "Loading " << ExpSpectrum << " into " << SpectrumName << endl;
    err = Structures[0].Spectrum.loadSpectrum 
      (ExpSpectrum.c_str(), Header, x, y);
    if (err != DX_NO_ERROR) {
      cout << "ERROR: Failed to open " << ExpSpectrum 
        << " when trying to add it to " << SpectrumName << endl;
      return SCRIPT_ERROR;
    }
    for (unsigned int i = 1; i < Structures.size(); i ++) {
      err = Structures[i].Spectrum.copySpectrum (Structures[0].Spectrum);
    }
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}


//------------------------------------------------------------------------------
// Command : spectrum.savespectrum <Filename>
// Saves the final theory spectrum generated by DEXA from the FEFF input files.
// If more than one sub-structure has been used, the individual sub-structure
// components will also be saved to their own files, suffixed with the structure
// name, and within those files, the individual path contributions. Also, if 
// Fourier filtering has been used in the analysis, this function will save the
// corresponding intermediate Fourier transform diagnostic spectra, which are
// calculated within DiffXasSpectrum::fourierFilter() and saved internally.
// Thus, there are potentially three more output files, which are only saved if
// the 'filter' command is used. They have <Filename> as their root name, but
// are suffixed with K_SUFFIX, R_SUFFIX, and Q_SUFFIX respectively for data in
// k-, R-, and q-space. 
//
// In all cases, the files are structured such that the x ordinate is saved in
// column 1, followed by a tab delimiter, and then the y ordinate in column 2.
// For the k- and R- space spectra, a third column contains the y ordinate data
// after the application of the Hann filter window, a fourth column with the
// real FT component, and a fifth column with the complex FT component. None of
// the files have a header.
//
// Since a DEXA theory spectrum may potentially be constructed from many sub-
// structures, each of which is described by its own MgstSpectrum object, the
// final, sum-total spectrum must be calculated here, at the level of the script
// engine. Note also that the total FT for all the structures summed together is
// calculated here from the sums of the Real and Complex components. This is
// used rather than the .getTheoryRData() function of DiffXasSpectrum as the
// later, which acts on each sub-structure before summation, introduces rounding
// errors.
//
int savespectrum (stringstream &args, vector<structure> &Structures) {
  unsigned int i, j;
  int Verbose;
  string Filename, SpectrumName;
  ostringstream oss;
  vector<Coordinate> StructureDiffXas, TotalDiffXas;
  vector<Coordinate> KData, RData, QData, KWindowed, RWindowed, RReal, RComplex;
  vector<Coordinate> TotalKData, TotalRData, TotalQData;
  vector<Coordinate> TotalKWindowed, TotalRWindowed, TotalRReal, TotalRComplex;
  
  args >> Filename >> SpectrumName >> Verbose;
  if (args.eof ()) {
  
    // Cycle through each of the sub-structures used in the calculation, and sum
    // all the components together to make the final theory spectrum.
    for (i = 0; i < Structures.size (); i ++) {

      // First save this individual sub-structure's contribution. In the case of
      // there being just one structure in the calculation, this component will 
      // also effectively be the sum over all contributions.
      oss.str (""); oss.clear ();
      oss << Filename;
      if (Structures[i].Name != "") { oss << "_" << Structures[i].Name; }
      StructureDiffXas.clear ();
      Structures[i].Spectrum.saveTheory (oss.str().c_str(), StructureDiffXas);

      // Add this structure's contribution to the running total over all 
      // structures
      if (i == 0) {
        TotalDiffXas = StructureDiffXas;      
      } else {
        for (j = 0; j < TotalDiffXas.size (); j ++) {
          TotalDiffXas[j].y += StructureDiffXas[j].y;
        }
      }
      
      // Check to see if the Fourier filter is also being used. If so, obtain
      // all the intermediate Fourier filter spectra and sum them over all
      // sub-structures as well.
      if (Structures[0].Spectrum.usingFourierFilter ()) {
        KData = Structures[i].Spectrum.getTheoryKData ();
        QData = Structures[i].Spectrum.getTheoryQData ();
        RReal = Structures[i].Spectrum.getTheoryRDataReal ();
        RComplex = Structures[i].Spectrum.getTheoryRDataComplex ();
        KWindowed = Structures[i].Spectrum.getTheoryKWindowed ();
        RWindowed = Structures[i].Spectrum.getTheoryRWindowed ();
        if (i == 0) {
          TotalKData = KData;
          TotalQData = QData;
          TotalRReal = RReal;
          TotalRComplex = RComplex;
          TotalKWindowed = KWindowed;
          TotalRWindowed = RWindowed;
        } else {
          for (j = 0; j < TotalKData.size (); j ++) {
            TotalKData[j].y += KData[j].y;
            TotalKWindowed[j].y += KWindowed[j].y;
          }
          for (j = 0; j < TotalRReal.size (); j ++) {
            TotalRReal[j].y += RReal[j].y;
            TotalRComplex[j].y += RComplex[j].y;
          }
          for (j = 0; j < TotalRWindowed.size (); j ++) {
            TotalRWindowed[j].y += RWindowed[j].y;
          }
          for (j = 0; j < TotalQData.size (); j ++) {
            TotalQData[j].y += QData[j].y;
          }
        }
      }
    }

    // Save the full, optimised theory DiffXAS (sum over paths and structures)
    // This is only needed if there is more than one sub-structure as otherwise
    // the full spectrum will just be the first structure component saved above.
    ofstream OutputFile (Filename.c_str(), ios::out);
    for (i = 0; i < TotalDiffXas.size (); i ++) {
      OutputFile << showpos << scientific << TotalDiffXas[i].x 
        << '\t' << TotalDiffXas[i].y << endl;
    }
    OutputFile.close ();

    // If the Fourier filter is being used, save the intermediate filter files.
    if (Structures[0].Spectrum.usingFourierFilter ()) {

      // Save the K-space data
      if (Verbose >= MED_VERBOSE) {
        oss.clear(); oss.str("");
        oss << Filename << K_SUFFIX;
        ofstream KFile (oss.str().c_str(), ios::out);
        for (i = 0; i < TotalKData.size (); i ++) {
          KFile << showpos << scientific << TotalKData[i].x << '\t' 
            << TotalKData[i].y << '\t' << TotalKWindowed[i].y << endl;
        };
        KFile.close ();
      }

      // Save the R-space data
      if (Verbose >= MIN_VERBOSE) {
        oss.clear(); oss.str("");
        oss << Filename << R_SUFFIX;
        ofstream RFile (oss.str().c_str(), ios::out);
        for (i = 0; i < TotalRReal.size (); i ++) {
          RFile << showpos << scientific << TotalRReal[i].x << '\t' 
            << sqrt (pow (TotalRReal[i].y, 2) + pow (TotalRComplex[i].y, 2))
            << '\t' << sqrt (pow (TotalRWindowed[i].y, 2) + pow (TotalRWindowed[i + TotalRWindowed.size() / 2].y, 2)) << '\t' << TotalRReal[i].y << '\t' 
            << TotalRComplex[i].y << endl;
        }
        RFile.close ();
      }
    }
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}


//------------------------------------------------------------------------------
// Command : filter <kmin> <kmax> <dk> <rmin> <rmax> <dr>
// Activates Fourier filtering for the data in the present analysis. Since the
// filter acts on the DiffXAS spectrum as a whole, all sub-structures, as well
// as the experimental spectrum itself, will be filtered between <kmin> <= k <=
// <kmax> in k-space, and <rmin> <= R <= <rmax> in r-space, using a Hann window.
// The additional two parameters <dk> and <dR> define the width of the 
// transition regions at each end of the Hann window. For example, if <kmin> is
// 1 and <dk> is 0.1, the window amplitude will be zero at 0.95, and progress to
// unity at 1.05.
//
// If 'filter' is used, all fitting performed by DEXA will then take place in 
// back-transformed q-space.
//
int filter (stringstream &args, vector<structure> &Structures) {
  string SpectrumName;
  double kMin, kMax, dK, rMin, rMax, dR;
  int Verbose;

  args >> kMin >> kMax >> dK >> rMin >> rMax >> dR >> SpectrumName >> Verbose;
  if (args.eof () && !args.fail ()) {
    if (Verbose >= MIN_VERBOSE)
    cout << "Fourier filtering data to " << kMin << " <= k <= " << kMax 
      << ", " << rMin << " <= R <= " << rMax << "\n  with a Hann window, where"
      << " dK = " << dK << " and dR = " << dR << endl;
    for (unsigned int i = 0; i < Structures.size(); i ++) {
      Structures[i].Spectrum.setFilter 
        (kMin, kMax, dK, rMin, rMax, dR);
    }
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}

int fitrange (stringstream &args, vector<structure> &Structures) {
  string SpectrumName;
  double qMin, qMax;
  int Verbose;

  args >> qMin >> qMax >> SpectrumName >> Verbose;
  if (args.eof () && !args.fail ()) {
    if (Verbose >= MIN_VERBOSE)
      cout << "Setting fit range to " << qMin << " <= q <= " << qMax << endl;
    for (unsigned int i = 0; i < Structures.size(); i ++) {
      Structures[i].Spectrum.setFitRange (qMin, qMax);
    }
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}


//------------------------------------------------------------------------------
// Command : spectrum.saveexperiment <Filename>
// Saves the processed experiment spectrum. This may just be the same as that
// loaded by 'spectrum', but if any processing has been done, such as Fourier
// filtering, the spectrum saved here will reflect those modifications.
//
// If sub-structures have been used, it is the case that the experimental 
// spectrum will be identical in all MgstSpectrum objects, so here, it is only 
// necessary to consider the experimental spectrum in the first sub-structure.
// This also covers the case where no sub-structures are used since the total
// structure will be in the form of just one, single 'sub-structure'.
//
// If Fourier filtering has been used in the analysis, this function will also
// save the corresponding intermediate Fourier transform diagnostic spectra,
// which are calculated within DiffXasSpectrum::fourierFilter() and saved 
// internally. Thus, there are potentially three extra output files, which are
// only saved if the 'filter' command is used. They have <Filename> as their 
// root name, but are suffixed with K_SUFFIX, R_SUFFIX, and Q_SUFFIX
// respectively for data in k-, R-, and q-space. 

// In all cases, the files are structured such that the x ordinate is saved in
// column 1, followed by a tab delimiter, and then the y ordinate in column 2.
// For the k- and R- space spectra, an additional, third column contains the y
// ordinate data after the application of the Hann filter window. None of the
// files have a header.
//
// Since a DEXA theory spectrum may potentially be constructed from many sub-
// structures, each of which is described by its own MgstSpectrum object, the
// final, sum-total spectrum must be calculated here, at the level of the script
// engine. 
//
//
int saveexperiment (stringstream &args, vector<structure> &Structures) {
  string Filename, SpectrumName;
  int Verbose;
  
  ostringstream oss;
  args >> Filename >> SpectrumName >> Verbose;
  if (args.eof ()) {
  
    // Save the processed experimental spectrum
    Structures[0].Spectrum.saveSpectrum (Filename.c_str());
    
    // If the Fourier filter has been used, the rest of the function from
    // here down will save the intermediate Fourier filter data.
    if (Structures[0].Spectrum.usingFourierFilter ()) {
      vector <Coordinate> KData = Structures[0].Spectrum.getExptKData ();
      vector <Coordinate> RData = Structures[0].Spectrum.getExptRData ();
      vector <Coordinate> QData = Structures[0].Spectrum.getExptQData ();
      vector <Coordinate> RReal = Structures[0].Spectrum.getExptRDataReal ();
      vector <Coordinate> RComplex = Structures[0].Spectrum.getExptRDataComplex ();
      vector <Coordinate> KWindowed = Structures[0].Spectrum.getExptKWindowed ();
      vector <Coordinate> RWindowed = Structures[0].Spectrum.getExptRWindowed ();

      // Save the K-space data
      if (Verbose >= MED_VERBOSE) {
        oss.clear(); oss.str("");
        oss << Filename << K_SUFFIX; ofstream KFile (oss.str().c_str(), ios::out);
        for (unsigned int i = 0; i < KData.size (); i ++) {
          KFile << showpos << scientific << KData[i].x << '\t' << KData[i].y 
            << '\t' << KWindowed[i].y << endl;
        };
        KFile.close ();
      }

      // Save the R-space data
      if (Verbose >= MIN_VERBOSE) {
        oss.clear(); oss.str("");
        oss << Filename << R_SUFFIX;
        ofstream RFile (oss.str().c_str(), ios::out);
        for (unsigned int i = 0; i < RData.size (); i ++) {
          RFile << showpos << scientific << RData[i].x << '\t' << RData[i].y 
            << '\t' << RWindowed[i].y << '\t' << RReal[i].y << '\t' << RComplex[i].y << endl;
        }
        RFile.close ();
      }
    }    
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}


//------------------------------------------------------------------------------
// Command : spectrum.loadsigmasqr <name> <header> <xcol> <ycol>
// Since the statistical component of the noise in a DiffXAS spectrum - which
// can be obtained from the high frequency components of the signal - is often
// more than an order of magnitude smaller than the external systematic noise,
// this function allows a custom noise spectrum to be loaded from <name>. This
// will then be used to weight the least-squares regression. This noise spectrum
// can, for example, be calculated from MPR's swerror code. Basically, the file
// should contain data in columns delimited by C white-space (space, tab, etc.).
// The x ordinate is then loaded from column <xcol> and the y ordinate from 
// <ycol>. If the file contains it header, it may be discarded by specifying how
// many lines in contains in <header>.
//
int loadsigmasqr (stringstream &args, vector<structure> &Structures) {
  string Filename, SpectrumName;
  int Header, x, y, err, Verbose;
  
  args >> Filename >> Header >> x >> y >> SpectrumName >> Verbose;
  if (args.eof () && !args.fail ()) {
    if (Verbose >= MIN_VERBOSE)
      cout << "Loading sigma^2 from " << Filename << " into " << SpectrumName << endl;
    for (unsigned int i = 0; i < Structures.size(); i ++) {
      err = Structures[i].Spectrum.loadSigmaSqr 
        (Filename.c_str(), Header, x, y);
      if (err != DX_NO_ERROR) {
        cout << "ERROR: Failed to open " << Filename 
          << " when trying to add sigma^2 to " << SpectrumName << endl;
        return SCRIPT_ERROR;
      }
    }
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}


//------------------------------------------------------------------------------
// createParameterArrays : Creates the two vector<double> parameter arrays
// required by Minuit - one for the starting parameter values and one for the
// errors associated with them. A third array for parameter names is created so
// as to identify parameters that have the same name in the input script and
// which hence need to be linked to the same internal Minuit parameter.
//
// For each parameter added to Minuit's UserParams array, the index of the array
// element is stored alongside the fit parameter in Spectra. This allows the
// optimised values to be read out and applied to the spectra correctly from 
// both the Fcn function and from applyOptimisedParameters below.
//
double createParameterArrays (vector<double> &UserParams, vector<double> 
  &MinusLimits, vector<double> &PlusLimits, vector<double> &ParamErrors, 
  vector<string> &ParamNames, vector<struct spectrum> &Spectra, int Verbose) {
  
  double NumIndepPoints = 0.0;
  for (unsigned int i = 0; i < Spectra.size (); i ++) {
    NumIndepPoints +=Spectra[i].SubStructure[0].Spectrum.numIndependentPoints();
    for (unsigned int l = 0; l < Spectra[i].SubStructure.size (); l ++) {
      // Deal with parameters that apply to the crystal as a whole first
      for (unsigned int j = 0; j < Spectra[i].SubStructure[l].Crystal.size (); j ++) {
        for (unsigned int k = 0; k < ParamNames.size (); k ++) {

          //If a parameter of the same name already exists in the Minuit arrays,
          //link the current parameter to it by assigning the same index. No new
          // array element is created.
          if (Spectra[i].SubStructure[l].Crystal[j].Name == ParamNames [k]) {
            Spectra[i].SubStructure[l].Crystal[j].Index = k;
            break;
          }
        }

        // The current parameter doesn't already exist in the Minuit arrays so
        // create a new element for it.
        if (Spectra[i].SubStructure[l].Crystal[j].Index == -1) {
          Spectra[i].SubStructure[l].Crystal[j].Index = UserParams.size ();
          if (Verbose >= MAX_VERBOSE) {
            cout << "Adding " << Spectra[i].SubStructure[l].Crystal[j].Name 
              << ": " << "Spectrum=" << i << " Struct=" << l << " Crystal=" 
              << j << " Index=" << Spectra[i].SubStructure[l].Crystal[j].Index 
              << endl;
          }
          ParamNames.push_back (Spectra[i].SubStructure[l].Crystal[j].Name);
          MinusLimits.push_back (Spectra[i].SubStructure[l].Crystal[j].MinusLimit);
          PlusLimits.push_back (Spectra[i].SubStructure[l].Crystal[j].PlusLimit);
          UserParams.push_back 
            ((Spectra[i].SubStructure[l].Spectrum.*Spectra[i].SubStructure[l].Crystal[j].ParamGet)());
          ParamErrors.push_back (PARAM_ERRORS);
        }
      }

      // Now handle the parameters that depend upon specific scattering paths in
      // the same way as those above.
      for (unsigned int j = 0; j < Spectra[i].SubStructure[l].Path.size (); j ++) {
        for (unsigned int k = 0; k < ParamNames.size (); k ++) {
          if (Spectra[i].SubStructure[l].Path[j].Name == ParamNames [k]) {
            Spectra[i].SubStructure[l].Path[j].Index = k;
            break;
          }
        }
        if (Spectra[i].SubStructure[l].Path[j].Index == -1) {
          Spectra[i].SubStructure[l].Path[j].Index = UserParams.size ();
          ParamNames.push_back (Spectra[i].SubStructure[l].Path[j].Name);  
          MinusLimits.push_back (Spectra[i].SubStructure[l].Path[j].MinusLimit);
          PlusLimits.push_back (Spectra[i].SubStructure[l].Path[j].PlusLimit);
          UserParams.push_back 
            ( (Spectra[i].SubStructure[l].Spectrum.*Spectra[i].SubStructure[l].Path[j].ParamGet)
              (Spectra[i].SubStructure[l].Path[j].PathNum) );
          ParamErrors.push_back (PARAM_ERRORS);
        }
      }
    }
  }
  return NumIndepPoints;
}


//------------------------------------------------------------------------------
// STRUCTURE Command List
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Command : spectrum.structure.addpath <chipfile> <header>
// Loads the scattering path found in <chipfile>, one of FEFF's chipnnnn.dat
// files. Since the header size of these files varies across calculations, the
// script must specify the header size in <header>. These lines are discarded
// before the chip data is read.
//
int addpath (stringstream &args, structure &Structure) {
  string Chipfile, SpectrumName;
  int Header, Err, Verbose;

  args >> Chipfile >> Header >> SpectrumName >> Verbose;
  if (args.eof ()) {
    Err = Structure.Spectrum.addPath (Chipfile.c_str(), Header);
    if (Err == MS_NO_ERROR) {
      if (Verbose >= MIN_VERBOSE)
        cout << "Added " << Chipfile << " to " << SpectrumName << endl;
      return SCRIPT_NO_ERROR;
    } else {
      cout << "Failed to add "<< Chipfile << " to " << SpectrumName << endl;
      return SCRIPT_ERROR;
    }
  }
  return SCRIPT_BAD_ARGS;
}


//------------------------------------------------------------------------------
// Command : spectrum.updatepath <chipfile> <parameter> <value>
//
int updatepath (stringstream &args, structure &Structure) {
  string Chipfile, Property, SpectrumName;
  double Value;
  int Err, Verbose;

  args >> Chipfile >> Property >> Value >> SpectrumName >> Verbose;
  if (args.eof ()) {
    Err = Structure.Spectrum.updatePath 
      (Chipfile.c_str(), Property.c_str(), Value);
    if (Err == MS_NO_ERROR) {
      if (Verbose >= MIN_VERBOSE)
        cout << SpectrumName << "." 
        << Chipfile << "." << Property << " = " << Value << endl;
      return SCRIPT_NO_ERROR;
    } else {
      cout << "Failed to set " << Property << " of " 
        << SpectrumName << "." << Chipfile << " to " << Value << endl;
      return SCRIPT_ERROR;
    }
  }
  return SCRIPT_BAD_ARGS;
}


//------------------------------------------------------------------------------
// Command : spectrum.setparameter <parameter name> [<path>]
//
int setparameter (stringstream &args, structure &Structure) {
  string Parameter, SpectrumName;
  double Value;
  int Err, Path = -1, Verbose;
  bool CommandExecuted = false;
  
  // There are two acceptable forms for this command: 
  //   string >> int >> int >> string   AND   string >> int >> string
  // Try reading an int into the third parameter first, and if that fails, we
  // must be looking at the second condition, so reset stream and read a string
  args >> Parameter >> Value;
  if (!(args >> Path)) { args.unget(); args.clear(); }
  args >> SpectrumName >> Verbose;  
  
  if (args.eof ()) {
    // Identify the command on the current line and execute it
    for (unsigned int i = 0; i < sizeof(Parameters)/sizeof(*Parameters);i ++){
      if (Parameter.compare (Parameters[i].Name) == 0) {
        if (Path >= 0) {
          Err = (Structure.Spectrum.*Parameters[i].FnPathSet) 
            (Path, Value);
          CommandExecuted = true;
        } else {
          Err = (Structure.Spectrum.*Parameters[i].FnSet) (Value);
          CommandExecuted = true;
        }
        break;
      }
    }
    if (!CommandExecuted) {
      cout << "Cannot set " << Parameter 
        << ". It is not a valid parameter." << endl;
      return SCRIPT_PARAMETER_NOT_FOUND;
    }

    // Examine the parameter get function error code and return the 
    // appropriate scripte error code.
    if (Err == MS_NO_ERROR) {
      if (Verbose >= MIN_VERBOSE) {
        cout << SpectrumName << "." << Parameter;
        if (Path >= 0) cout << " (" << Path << ")";
        cout << " = " << Value << endl;
      }
      return SCRIPT_NO_ERROR;
    } else {
      cout << "Failed to set " << SpectrumName << "." << Parameter;
      if (Path >= 0) cout << " (" << Path << ")";
      cout << " to " << Value << " (error code " << Err << ")" << endl;
      return SCRIPT_ERROR;
    }
  }
  return SCRIPT_BAD_ARGS;
}


//------------------------------------------------------------------------------
// Command : spectrum.fitparameter <parameter name> [<path>]
// Adds the <parameter name> of spectrum to the list of parameters to be fitted.
// If it's a parameter that acts on a single path, the <path> index (base 0)
// must also be specified. 
//
// The halfVoigt parameter is handled somewhat explicitly since the "path"
// property here does in fact refer to a tensor element, and thus does not
// conform to the path number limits checking.
//
// The list of fitting parameters itself is just a collection of function 
// pointers. Each fit parameter requires a pointer to its 'get' function, and
// another to its 'set' function. A global list of all valid parameters is 
// specified in the Parameters array, defined in ScriptReader.h. This function
// takes the appropriate pointers from that array and stores them in Structure.
//
int fitparameter (stringstream &args, structure &Structure) {
  string Parameter, SpectrumName, ParamName;
  int Path = -1;
  int ArgPos = 0, Verbose;
  double MinusLimit = 0.0;
  double PlusLimit = 0.0;
  
  // There are two acceptable forms for this command: 
  //   string >> int >> string   AND   string >> string
  // Try reading an int into the second parameter first, and if that fails, we
  // must be looking at the second condition, so reset stream and read a string
  args >> ParamName >> Parameter;
  ArgPos = args.tellg ();
  args >> Path >> SpectrumName >> Verbose;
  if (args.fail() || !args.eof()) {   
    args.clear(); args.seekg (ArgPos, ios::beg); Path = -1;
    args >> MinusLimit >> PlusLimit >> SpectrumName >> Verbose;
    if (args.fail() || !args.eof()) {
      args.clear(); args.seekg (ArgPos, ios::beg);
      args >> Path >> MinusLimit >> PlusLimit >> SpectrumName >> Verbose;
      if (args.fail()) {       
        args.clear(); args.seekg (ArgPos, ios::beg);
        args >> SpectrumName >> Verbose;
      }
    }
  }
  
  if (args.eof () && !args.fail ()) {
    // Identify the parameter on the current line and add it to the fit params
    for (unsigned int i = 0; i < sizeof(Parameters)/sizeof(*Parameters);i ++){
      if (Parameter.compare (Parameters[i].Name) == 0) {

        // First make sure that the specified path is in bounds. Abort if not.
        if (Path >= Structure.Spectrum.numPaths () && Parameter != "halfVoigt"){
          cout << "Invalid path index: " << Path << endl;
          return SCRIPT_BAD_ARGS;
        }

        // Set the Get and Set function pointers according to whether or not
        // they act on individual paths (first case) or the spectrum as a
        // whole (second case).
        if (Path >= 0) {
          pathParam NewPathParam;
          NewPathParam.ParamGet = Parameters[i].FnPathGet;
          NewPathParam.ParamSet = Parameters[i].FnPathSet;
          NewPathParam.PathNum = Path;
          NewPathParam.PlusLimit = PlusLimit;
          NewPathParam.MinusLimit = MinusLimit;
          NewPathParam.Name = ParamName;
          NewPathParam.Index = -1;
          Structure.Path.push_back (NewPathParam); 
        } else {
          crystalParam NewCrystalParam;
          NewCrystalParam.ParamGet = Parameters[i].FnGet;
          NewCrystalParam.ParamSet = Parameters[i].FnSet;
          NewCrystalParam.Name = ParamName;
          NewCrystalParam.PlusLimit = PlusLimit;
          NewCrystalParam.MinusLimit = MinusLimit;
          NewCrystalParam.Index = -1;
          Structure.Crystal.push_back (NewCrystalParam);
        }
        cout << "Fitting: " << SpectrumName << "." << Parameter;
        if (Parameter == "halfVoigt") cout << " element " << Path;
        else if (Path >= 0) cout << " (" <<
          Structure.Spectrum.getPathj(Path).ChipName << ")";
        cout << endl;
        return SCRIPT_NO_ERROR;
      }
    }

    // To reach here, some sort of error must have occurred. Identify which
    // kind of error and return the appropriate error code.
    cout << "Invalid parameter " << SpectrumName << "." << Parameter;
    if (Path >= 0 && Parameter != "halfVoigt") cout << " (" <<
      Structure.Spectrum.getPathj(Path).ChipName << ")";
    cout << endl;
    return SCRIPT_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}





int fitfulltensor (stringstream &args, structure &Structure) {
  string SpectrumName;
  int Verbose;
  
  args >> SpectrumName >> Verbose;
  ostringstream oss;
  if (args.eof () && !args.fail ()) {
    for (unsigned int i = 0; i < MGST_TENSOR_SIZE; i ++) {
      pathParam NewPathParam;
      NewPathParam.ParamGet = &MgstSpectrum::tensorElement;
      NewPathParam.ParamSet = &MgstSpectrum::tensorElement;
      NewPathParam.PathNum = i;
      oss.str (""); oss << "Tensor" << i;
      NewPathParam.Name = oss.str();
      NewPathParam.Index = -1;
      Structure.Path.push_back (NewPathParam); 
    }
    cout << "Fitting: ALL " << MGST_TENSOR_SIZE 
      << " Magnetostriction Tensor Elements" << endl;
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}

int fitfullvoigt (stringstream &args, structure &Structure) {
  string SpectrumName;
  int Verbose;
  
  args >> SpectrumName >> Verbose;
  ostringstream oss;
  if (args.eof () && !args.fail ()) {
    for (unsigned int i = 0; i < VOIGT_SIZE; i ++) {
      pathParam NewPathParam;
      NewPathParam.ParamGet = &MgstSpectrum::voigtElement;
      NewPathParam.ParamSet = &MgstSpectrum::voigtElement;
      NewPathParam.PathNum = i;
      oss.str (""); oss << "Voigt" << i;
      NewPathParam.Name = oss.str();
      NewPathParam.Index = -1;
      Structure.Path.push_back (NewPathParam); 
    }
    cout << "Fitting: ALL " << VOIGT_SIZE << " Voigt Matrix Elements" << endl;
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}


/*int fithalfvoigt (stringstream &args, structure &Structure) {
  string SpectrumName;
  int Element;
  int FittedElements = -1;
  streampos FirstArgument;
  
  FirstArgument = args.tellg ();
  while (!args.eof ()) { FittedElements ++; args >> SpectrumName; }
  args.seekg (FirstArgument, ios::beg);
  args.clear ();
  ostringstream oss;
  if (!args.fail ()) {
    int SpectrumRef = findSpectrum (SpectrumName, Structure);
    if (SpectrumRef != SCRIPT_SPECTRUM_NOT_FOUND) {
      if (FittedElements > 0) {
        cout << "Fitting: HALF Voigt Matrix (elements " << flush;
        for (int i = 0; i < FittedElements; i ++) {
          args >> Element;
          if (Element < 0 || Element >= HALF_VOIGT_SIZE) {
            cout << "ERROR: Invalid Voigt element. Index must be between 0 and "
              << HALF_VOIGT_SIZE << endl;
            return SCRIPT_ERROR;
          }
          cout << Element << ", " << flush;
          pathParam NewPathParam;
          NewPathParam.ParamGet = &MgstSpectrum::halfVoigtElement;
          NewPathParam.ParamSet = &MgstSpectrum::halfVoigtElement;
          NewPathParam.PathNum = Element;
          oss.str (""); oss << "Voigt" << Element;
          NewPathParam.Name = oss.str();
          NewPathParam.Index = -1;
          Structure[SpectrumRef].Path.push_back (NewPathParam); 
        }
        cout << "\b\b)" << endl;
      } else {
        cout << "Fitting: HALF Voigt Matrix (all " 
          << HALF_VOIGT_SIZE << " elements)" << endl;
        for (unsigned int i = 0; i < HALF_VOIGT_SIZE; i ++) {
          pathParam NewPathParam;
          NewPathParam.ParamGet = &MgstSpectrum::halfVoigtElement;
          NewPathParam.ParamSet = &MgstSpectrum::halfVoigtElement;
          NewPathParam.PathNum = i;
          oss.str (""); oss << "Voigt" << i;
          NewPathParam.Name = oss.str();
          NewPathParam.Index = -1;
          Structure[SpectrumRef].Path.push_back (NewPathParam); 
        }
        FittedElements = HALF_VOIGT_SIZE;
      }
      return SCRIPT_NO_ERROR;
    }
    return SCRIPT_SPECTRUM_NOT_FOUND;
  }
  return SCRIPT_BAD_ARGS;
}*/



int symmetry (stringstream &args, structure &Structure) {
  string Symmetry, SpectrumName;
  int Err, Verbose;
  
  args >> Symmetry >> SpectrumName >> Verbose;
  if (args.eof () && !args.fail ()) {
    if (Symmetry == "cubic") { Err = 
      Structure.Spectrum.symmetry (MT_TENSOR_CUBIC); 
    }
    else if (Symmetry == "isotropic") { Err =
      Structure.Spectrum.symmetry (MT_TENSOR_ISOTROPIC);
    }
    else if (Symmetry == "cylindrical") { Err = 
      Structure.Spectrum.symmetry (MT_TENSOR_CYLINDRICAL);
    }
    else if (Symmetry == "tetragonal") { Err = 
      Structure.Spectrum.symmetry (MT_TENSOR_TETRAGONAL_T1); 
    }
    else if (Symmetry == "tetragonal2") { Err = 
      Structure.Spectrum.symmetry (MT_TENSOR_TETRAGONAL_T2); 
    }
    else {
      cout << "Unrecognised crystal symmetry: " << Symmetry << endl;
      return SCRIPT_ERROR;
    }
    if (Verbose >= MIN_VERBOSE)
      cout << SpectrumName << ".symmetry = " << Symmetry << endl;
    return SCRIPT_NO_ERROR;
  }
  return SCRIPT_BAD_ARGS;
}
