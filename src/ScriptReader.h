// ScriptReader.h
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
#ifndef SCRIPT_READER_H
#define SCRIPT_READER_H

#include <string>
#include "ThermalSpectrum.h"

using namespace::std;

#define PARAM_ERRORS 1.0e-3

// Script error codes
#define SCRIPT_NO_ERROR             0
#define SCRIPT_LOAD_ERROR           1
#define SCRIPT_SPECTRUM_NOT_FOUND  -1
#define SCRIPT_STRUCTURE_NOT_FOUND -2
#define SCRIPT_SYNTAX_ERROR         2
#define SCRIPT_ERROR                3
#define SCRIPT_BAD_ARGS             4
#define SCRIPT_COMMAND_NOT_FOUND    5
#define SCRIPT_PARAMETER_NOT_FOUND  6

// General script definitions
#define SCRIPT_PRECISION           6  // Decimal places for output of parameters
#define COMMENT_DELIMITER '*'         // Ignore lines where * is the first field
#define END_OF_FILE       '�'         // Used to prevent reading past EOF
#define END_OF_LINE       '\0'        // Used to prevent reading past EOL

#define MGST_TENSOR_SIZE  81
#define VOIGT_SIZE        36
#define HALF_VOIGT_SIZE   21

// Verbose output options
#define NO_VERBOSE  0
#define MIN_VERBOSE 1
#define MED_VERBOSE 2
#define MAX_VERBOSE 3
#define ERR_VERBOSE -1

// The following structs, pathParam and crystalParam, store the information
// required to manipulate the fit parameters in a ThermalSpectrum. Both contain two
// function pointers, which are refer to the parameter's get and set functions, 
// and a name for the parameter. When a parameter with a matching name is 
// identified in the input script, the function pointers are copied from the
// "const struct parameter Parameters []" array below. All parameters must be in
// this array or they will not be identified.
//
// Both structs also contain an int Index, which stores the array index used by
// the parameter within the solver's (GSL or Minuit) own list of fit parameters.
// Normally, a new parameter is just tagged onto the end of the solver's 
// parameter list and so this is just assigned the next sequential index. 
// However, if two or more parameters are LINKED, only the first will be given a
// new element in the solver's parameter list. The rest will just have Index set
// to match that of the first, so that they extract their value from the same
// array element.
//
// The two structs differ only in that pathParam contains an additional int (int
// PathNum), which stores the ScatteringPaths index (from DiffXasSpectrum.h) of 
// the path upon which the parameter acts.
//
struct pathParam {
  double (ThermalSpectrum::*ParamGet) (int);
  int (ThermalSpectrum::*ParamSet) (int, double);
  string Name;
  int PathNum;
  int Index;
  double PlusLimit, MinusLimit;
  
  ~pathParam () { ParamGet = NULL; ParamSet = NULL; Name.clear (); }
};

struct crystalParam {
  double (ThermalSpectrum::*ParamGet) ();
  int (ThermalSpectrum::*ParamSet) (double);
  string Name;
  int Index;
  double PlusLimit, MinusLimit;

  ~crystalParam () { ParamGet = NULL; ParamSet = NULL; Name.clear (); }
};

// When a new spectrum is created in the input script, a struct Spectrum object
// is created and stored in the "vector <struct spectrum> Spectra" vector that
// is passed to all the command handlers defined below. The struct creates a
// ThermalSpectrum object (i.e. the actual spectrum itself), and two vectors for
// the associated fit parameters: Crystal for parameters derived from 
// crystalParams, which act on the crystal as a whole, and Path for parameters
// derived from pathParam, that act on single scattering paths only.
//
typedef struct td_structure {
  ThermalSpectrum Spectrum;
  vector <struct crystalParam> Crystal;
  vector <struct pathParam> Path;
  string Name;
  
  ~td_structure () { Crystal.clear (); Path.clear (); }
} structure;

struct spectrum {
  vector<structure> SubStructure;
  string Name;
  
  spectrum () { Name = ""; SubStructure.resize (1); }
  ~spectrum () { SubStructure.clear (); }
};


struct parameter {
  string Name;
  double (ThermalSpectrum::*FnGet) ();
  int (ThermalSpectrum::*FnSet) (double);
  double (ThermalSpectrum::*FnPathGet) (int);
  int (ThermalSpectrum::*FnPathSet) (int, double);
  
  static parameter add (string NewName, 
    double (ThermalSpectrum::*NewFnGet) (),
    int (ThermalSpectrum::*NewFnSet) (double)) {
    parameter NewParameter;
    NewParameter.Name = NewName;
    NewParameter.FnGet = NewFnGet; NewParameter.FnPathGet = NULL;
    NewParameter.FnSet = NewFnSet; NewParameter.FnPathSet = NULL;
    return NewParameter;
  }
  
  static parameter addPath (string NewName, 
    double (ThermalSpectrum::*NewFnPathGet) (int),
    int (ThermalSpectrum::*NewFnPathSet) (int, double)) {
    parameter NewParameter;
    NewParameter.Name = NewName;
    NewParameter.FnPathGet = NewFnPathGet; NewParameter.FnGet = NULL;
    NewParameter.FnPathSet = NewFnPathSet; NewParameter.FnSet = NULL;
    return NewParameter;
  }
};

// In order for the readScript function to be able to identify commands, a look
// up list must created. This list, "const struct command Commands []" below, 
// is an array of struct command objects, which contain the command name and a 
// function pointer to its associated handler. The struct also contains the
// static command add (...) function to define a new instance of itself when
// populating the Commands array.
//
struct command_global {
  string Name;
  int (*Function)(stringstream &, vector<struct spectrum> &);

  static command_global add (const char* NewName,
    int (*Fn)(stringstream &, vector<struct spectrum> &)) {
    command_global NewCommand;
    NewCommand.Name = NewName;
    NewCommand.Function = Fn;
    return NewCommand;
  }
};

struct command_spectrum {
  string Name;
  int (*Function)(stringstream &, vector<structure> &);

  static command_spectrum add (const char* NewName,
    int (*Fn)(stringstream &, vector<structure> &)) {
    command_spectrum NewCommand;
    NewCommand.Name = NewName;
    NewCommand.Function = Fn;
    return NewCommand;
  }
};

struct command_structure {
  string Name;
  int (*Function)(stringstream &, structure &);

  static command_structure add (const char* NewName,
    int (*Fn)(stringstream &, structure &)) {
    command_structure NewCommand;
    NewCommand.Name = NewName;
    NewCommand.Function = Fn;
    return NewCommand;
  }
};


int printParams (ostream &Stream, vector<struct spectrum> &Spectra);

// The following functions are those which actually execute the operations
// required of each script command. They must follow the function prototype
// given in struct command above. Just to make life easier, the functions
// have the same names as their associated commands.
int create (stringstream &args, vector<struct spectrum> &Spectra);
int startfit (stringstream &args, vector<struct spectrum> &Spectra);
int saveparams (stringstream &args, vector<struct spectrum> &Spectra);
int preforientation (stringstream &args, vector<struct spectrum> &Spectra);
int magnetisation (stringstream &args, vector<struct spectrum> &Spectra);
int polarisation (stringstream &args, vector<struct spectrum> &Spectra);
int numavesteps (stringstream &args, vector<struct spectrum> &Spectra);

int addstructure (stringstream &args, vector<structure> &Structures);
int feffinp (stringstream &args, vector<structure> &Structures);
int pathsdat (stringstream &args, vector<structure> &Structures);
int spectrum (stringstream &args, vector<structure> &Structures);
int savespectrum (stringstream &args, vector<structure> &Structures);
int filter (stringstream &args, vector<structure> &Structures);
int fitrange (stringstream &args, vector<structure> &Structures);
int saveexperiment (stringstream &args, vector<structure> &Structures);
int loadsigmasqr (stringstream &args, vector<structure> &Structures);

int addpath (stringstream &args, structure &Structure);
int updatepath (stringstream &args, structure &Structure);
int setparameter (stringstream &args, structure &Structure);
int fitparameter (stringstream &args, structure &Structure);
int symmetry (stringstream &args, structure &Structure);
int fitfulltensor (stringstream &args, structure &Structure);
int fitfullvoigt (stringstream &args, structure &Structure);


// Define the list of valid commands for the script. The Commands array is used
// in readScript to identify each command and then call its executing function.
// The commands are case sensitive, and should be entered below exactly how you
// wish them to be written in the script. The add function of 'struct command' 
// is used to set all the properties of each command. To define a command, add
// it to the following list, declare its executing function just above, and 
// finally define that function (in ScriptReader.cpp). The syntax below is:
//   command::add (<cmd name>, <num args>, &<executing function>);
const struct command_global GlobalCommands [] = {
  command_global::add ("create",          &create),
  command_global::add ("startfit",        &startfit),
  command_global::add ("saveparams",      &saveparams)/*,
  command_global::add ("preforientation", &preforientation),
  command_global::add ("magnetisation",   &magnetisation),
  command_global::add ("polarisation",    &polarisation),
  command_global::add ("numavesteps",     &numavesteps)*/
};

const struct command_spectrum SpectrumCommands [] = {
  command_spectrum::add ("addstructure",   &addstructure),
  command_spectrum::add ("feffinp",        &feffinp),
  command_spectrum::add ("pathsdat",       &pathsdat),
  command_spectrum::add ("spectrum",       &spectrum),
  command_spectrum::add ("savespectrum",   &savespectrum),
  command_spectrum::add ("filter",         &filter),
  command_spectrum::add ("fitrange",       &fitrange),
  command_spectrum::add ("saveexperiment", &saveexperiment),
  command_spectrum::add ("loadsigmasqr",   &loadsigmasqr)
};

const struct command_structure StructureCommands [] = {
  command_structure::add ("addpath",         &addpath),
  command_structure::add ("updatepath",      &updatepath),
  command_structure::add ("setparameter",    &setparameter),
  command_structure::add ("fitparameter",    &fitparameter)/*,
  command_structure::add ("symmetry",        &symmetry),
  command_structure::add ("fitfulltensor",   &fitfulltensor),
  command_structure::add ("fitfullvoigt",    &fitfullvoigt)*/
};


// Define the list of spectrum parameters that may either be set in the script
// or added to the list of fit parameters. There are two types of function
// pointer to the executing code. The first takes just one double as an argument
// and acts on the spectrum as a whole. The second, takes and int and a double,
// where the int defines a path index. The double then only acts on that path.
// Only one or the other may be used. For the first type, use parameter:add, and
// for the second use parameter::addPath. setparameter(...) will work out which
// of the two pointers to use by looking at the arguments to setparameter in the
// input script i.e. does it have a path number as well as a parameter value, or
// just a parameter value?
const struct parameter Parameters [] = {
  parameter::add ("dE",       &ThermalSpectrum::dE,       &ThermalSpectrum::dE),
  parameter::addPath ("dSig2j", &ThermalSpectrum::dSig2,  &ThermalSpectrum::dSig2),
  parameter::addPath ("dS", &ThermalSpectrum::dS,      &ThermalSpectrum::dS)
};


int readScript (const char* ScriptName, vector <struct spectrum> &Spectra, int Verbose);
int findSpectrum (string Name, vector<struct spectrum> &Spectra);
void errorHandler (int Error, string LineOfData, int LineNumber);
double createParameterArrays (vector<double> &UserParams, vector<double> 
  &MinusLimits, vector<double> &PlusLimits, vector<double> &ParamErrors, 
  vector<string> &ParamNames, vector<struct spectrum> &Spectra, int Verbose);
void catchSigInt (int sig);

#endif // SCRIPT_READER_H
