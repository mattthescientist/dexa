// ScriptReader.cpp
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
#include "ScriptReader.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

// This file contains a general script reading engine for use with DiffXAS
// analysis, but which could easily be transferred to other applications by
// removing/modifying the printParams and findSpectrum functions. The engine
// here relies heavily on function pointers. The ScriptReader.h file contains
// a list of valid script commands, and with each, a pointer to its associated
// handler function. The readScript function below scans the input script for
// these commands, and, when it finds one, calls its associated handler. This
// handler should then perform all the operations required for the associated
// command. 
//
// Note that no argument checking is performed by readScript (). The arguments
// should be checked in the event handler when they are parsed. In the event of
// a failure, the appropriate error code should be returned. For a list of error
// codes, see ScriptReader.h. Some of these codes generate error messages on
// stdout when errorHandler() is called. For others, the message should be
// generated in the event handler if deemed necessary.
//
// General DiffXAS event handlers are found in ScriptEventsGeneral.cpp,
// those specific to the GSL backend in ScriptEventsGsl.cpp, and those specific
// to the Minuit backend in ScriptEventsMinuit.cpp.

//------------------------------------------------------------------------------
// printParams (ostream &, vector <struct spectrum> &) : Generates parameter
// output information and sends it to the output stream at arg1. This would
// generally either be cout or a file stream that has been initialised elsewhere
//
int printParams (ostream &Stream, vector<struct spectrum> &Spectra) {
  Stream << "The parameter values are:" << endl;
  for (unsigned int i = 0; i < Spectra.size (); i ++) {
    Stream << endl << "For " << Spectra[i].Name << ":" << endl << endl;
    for (unsigned int k = 0; k < Spectra[i].SubStructure.size (); k ++) {

      for (unsigned int j = 0; j < Spectra[i].SubStructure[k].Crystal.size (); j ++) {
        Stream << "  "  << Spectra[i].SubStructure[k].Crystal[j].Name << " = " 
          << (Spectra[i].SubStructure[k].Spectrum.*Spectra[i].SubStructure[k].Crystal[j].ParamGet) () << endl;
      }
      for (unsigned int j = 0; j < Spectra[i].SubStructure[k].Path.size (); j ++) {
        Stream << "  " << Spectra[i].SubStructure[k].Path[j].Name
          << " (" << /*Spectra[i].SubStructure[k].Spectrum.getPathj(*/Spectra[i].SubStructure[k].Path[j].PathNum/*).ChipName*/ << ") = "
          << (Spectra[i].SubStructure[k].Spectrum.*Spectra[i].SubStructure[k].Path[j].ParamGet) 
            (Spectra[i].SubStructure[k].Path[j].PathNum) << endl;
      }
    }
  }
  return SCRIPT_NO_ERROR;
}


//------------------------------------------------------------------------------
// findSpectrum (string, vector<struct spectrum> &) : Given a vector of spectrum
// objects, this function will compare the spectrum names against the input Name
// and return the vector index of the corresponding match,if found. If Name does
// not exist in the vector, SCRIPT_SPECTRUM_NOT_FOUND is returned instead.
//
int findSpectrum (string Name, vector<struct spectrum> &Spectra) {
  for (unsigned int i = 0; i < Spectra.size (); i ++) {
    if (!Name.compare (Spectra[i].Name)) {
      return i;
    }
  }
  cout << "There is no spectrum named " << Name << endl;
  return SCRIPT_SPECTRUM_NOT_FOUND;
}

//------------------------------------------------------------------------------
// findStructure (string, vector<struct spectrum> &) : Given a vector of
// sub-structures, this function will compare the structure names against the
// input Name. If this name is found, the index of that structure in the vector
// is returned. Failing that, SCRIPT_STRUCTURE_NOT_FOUND is returned.
//
int findStructure (string Name, vector<structure> &Structures) {
  for (unsigned int i = 0; i < Structures.size (); i ++) {
    if (!Name.compare (Structures[i].Name)) {
      return i;
    }
  }
  cout << "There is no structure named " << Name << endl;
  return SCRIPT_STRUCTURE_NOT_FOUND;
}

//------------------------------------------------------------------------------
// readScript (const char*, vector <struct spectrum>) : This function starts the
// script reading process. The vector Spectra, passed in by reference, should be
// empty. It will be filled with the spectra loaded from the input script given
// by Filename.
//
// The list of valid script commands is obtained from 'Commands', which can be
// found in ScriptReader.h. For each command, this array holds a function 
// pointer. For every command that is read, this associated function is called,
// with the aim of actually executing the operations required by the command.
//
// Lines which are blank, or where the first non-whitespace character is a
// COMMENT_DELIMITER, are ignored. The first field on all other lines must be
// a valid command, which is then followed by its arguments. If no valid command
// is found on a given line, the script is aborted with SCRIPT_SYNTAX_ERROR and 
// a message to stdout stating the line on which the error occurred.The validity
// of the arguments to a command are checked in the command's own function.
//
int readScript (const char* Filename, vector <struct spectrum> &Spectra, int Verbose) {
  
  string LineOfData;
  stringstream ScriptStream;
  stringstream iss;
  string Command;
  string SpectrumName;
  string StructureName;
  string TempFile;
  size_t DelimiterPos;
  int Error;
  int LineNumber = 0;

  // Open the specified script file and verify it opened correctly
  ifstream InputScript (Filename, ios::in);
  if (! InputScript.is_open()) {
    cout << "ERROR: Unable to open " << Filename << 
      ". Check the file exists and is readable." << endl;
    return SCRIPT_LOAD_ERROR;
  }
  ScriptStream << InputScript.rdbuf();
  InputScript.close();
  ScriptStream.seekg (ios_base::beg);
  cout <<
 
"================================================================================"
<< endl;
  cout << "Loading parameters from " << Filename << "..." << endl << endl;
  
  // Now scan through each line of the input script. Find the command name, and
  // (where applicable) extract the spectrum to which the it applies.
  while (!ScriptStream.eof()) {
  
    // Get the next line of data and load it into an istringstream for formatted
    // output later on.
    getline (ScriptStream, LineOfData);
    if (Verbose >= MAX_VERBOSE) {
    	cout << LineOfData << endl;
    }
    LineNumber ++;
    iss.clear ();
    iss.str (LineOfData);
    
    // First catch any blank or comment lines and ignore them
    if (LineOfData[0] != END_OF_LINE) {
      iss >> Command;
      if (Command[0] != COMMENT_DELIMITER) {
      
        // If the current command is to be applied to a given spectrum, separate
        // the command name and spectrum names here. Put the spectrum name onto
        // the end of the list of arguments so that the command can find it.
        Error = SCRIPT_COMMAND_NOT_FOUND;
        if ((DelimiterPos = Command.find ('.')) != string::npos) {
          iss.clear ();
          iss.seekp (LineOfData.length ());
          iss << " " << Command.substr (0, DelimiterPos);
          SpectrumName = Command.substr (0, DelimiterPos);
          Command = Command.substr (DelimiterPos + 1, Command.length () - DelimiterPos - 1);
          
          int SpectrumRef = findSpectrum (SpectrumName, Spectra);
          if (SpectrumRef != SCRIPT_SPECTRUM_NOT_FOUND) {  
            // Now find the structure name, if given
            if ((DelimiterPos = Command.find ('.')) != string::npos) {
              iss << "." << Command.substr (0, DelimiterPos);
              StructureName = Command.substr (0, DelimiterPos);
              Command = Command.substr (DelimiterPos + 1, Command.length () - DelimiterPos - 1);
              
              int StructureRef = findStructure (StructureName, Spectra[SpectrumRef].SubStructure);
              if (StructureRef != SCRIPT_STRUCTURE_NOT_FOUND) {
              
                // There is both a spectrum and a structure name so the command
                // must be a structure specific command. Handle it with the
                // StructureCommands list.
                for (unsigned int i = 0; i < sizeof (StructureCommands)
                  / sizeof (*StructureCommands); i ++) {
                  if (Command.compare (StructureCommands[i].Name) == 0) {
                    iss << " " << Verbose;
                    Error = (StructureCommands[i].Function) 
                      (iss, Spectra[SpectrumRef].SubStructure[StructureRef]);
                    break;
                  }
                }
              } else { Error = SCRIPT_STRUCTURE_NOT_FOUND; }
            } else {
            
              // There is no structure name so the command must either be a
              // spectrum specific command, or a structure specific command used
              // in the context of a single-structure spectrum. Handle it with 
              // both the SpectrumCommands list and the StructureCommands list.
              // In the latter case, the StructureRef will be zero.
              for (unsigned int i = 0; 
                i < sizeof(SpectrumCommands) / sizeof(*SpectrumCommands); i ++){
                if (Command.compare (SpectrumCommands[i].Name) == 0) {
                  iss << " " << Verbose;
                  Error = (SpectrumCommands[i].Function) 
                    (iss, Spectra[SpectrumRef].SubStructure);
                  break;
                }
              }
              for (unsigned int i = 0; 
                i < sizeof(StructureCommands)/sizeof(*StructureCommands); i ++){
                if (Command.compare (StructureCommands[i].Name) == 0) {
                  iss << " " << Verbose;
                  Error = (StructureCommands[i].Function) 
                    (iss, Spectra[SpectrumRef].SubStructure[0]);
                  break;
                }
              }
            }
          } else { Error = SCRIPT_SPECTRUM_NOT_FOUND; }
        } else { 

          // There is no spectrum name so the command must be a global
          // command. Handle it with the GlobalCommands list.
          iss.clear ();
          iss.seekp (LineOfData.length ());
          iss << " " << Verbose;
          for (unsigned int i = 0; 
            i < sizeof(GlobalCommands) / sizeof(*GlobalCommands); i ++){
            if (Command.compare (GlobalCommands[i].Name) == 0) {
              Error = (GlobalCommands[i].Function) (iss, Spectra);
              break;
            }
          } 
        }
        if (Error) {
          errorHandler (Error, LineOfData, LineNumber);
          cout << "Aborting!" << endl;
          return Error;
        }
      }
    }// else { cout << endl; }
  }
  cout << "Successfully finished reading input script" << endl;
  cout <<
 
"================================================================================"
<< endl;
  return SCRIPT_NO_ERROR;
}


//------------------------------------------------------------------------------
// errorHandler (int, string, int) : Prints error information to standard output
// The error codes may be found at the top of ScriptReader.h. Not all codes are
// handled explicitly here, since some types of error generate their own output
// in the command handlers above. In those cases, it is not necessary to provide
// the user with any more information here other than the script line which
// caused the error.
void errorHandler (int Error, string LineOfData, int LineNumber) {
  switch (Error) {
    case SCRIPT_BAD_ARGS:
      cout << "Bad command arguments";
      break;
    case SCRIPT_COMMAND_NOT_FOUND:
      cout << "Unrecognised command found"; 
      break;
    case SCRIPT_SPECTRUM_NOT_FOUND:
      cout << "Bad spectrum specified";
      break;
    default:
      cout << "An error occurred (code " << Error << ")";
      break;
  }
  cout << " on line " << LineNumber << ": \"" << LineOfData << "\"" << endl;
}
