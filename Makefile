# Makefile
# Copyright (C) 2007-2009 Matthew Ruffoni
#
# This file is part of DEXA.
#
# DEXA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DEXA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with DEXA.  If not, see <http:#www.gnu.org/licenses/>.
#
# ------------------------------------------------------------------------------
# 
# By default this file will attempt to compile DEXA using g++, and will look for
# GSL and Minuit in /usr/local/gsl and /usr/local/minuit. If these differ on
# your system, change the CC variable just below and the paths in the 'Flags'
# section.
#
# Valid commands:
#
#   make        : Makes just the Minuit version of DEXA
#   make minuit : Makes just the Minuit version of DEXA
#   make all    : Makes both GSL and Minuit versions of DEXA
#   make gsl    : Makes just the GSL version (usful if Minuit is not available)
#   make clean  : Removes source object files after compilation
#
# After compilation, the executables DexaMinuit and DexaGsl will appear. These
# can then be copied to a system location for all users to access if desired.

CC = g++
SRC_DIR := src

# Output files
BIN_GSL := dexa_gsl
BIN_MIN := dexa

# Source files: COM common, GSL Gsl only, MIN Minuit only
_OBJ_COM := Fourier.o DiffXasSpectrum.o PathFinder.o ThermalSpectrum.o ScriptReader.o ScriptEventsGeneral.o
_OBJ_GSL := ScriptEventsGsl.o MgstFcnGsl.o MgstFit.o
_OBJ_MIN := ScriptEventsMinuit.o MgstFcnMinuit.o MgstFit.o

OBJ_COM := $(patsubst %,$(SRC_DIR)/%,$(_OBJ_COM))
OBJ_GSL := $(patsubst %,$(SRC_DIR)/%,$(_OBJ_GSL))
OBJ_MIN := $(patsubst %,$(SRC_DIR)/%,$(_OBJ_MIN))

# Flags
C_FLAGS := -I/usr/local/gsl/include -I/usr/local/minuit/include -Wall -g
FLAGS_GSL := -L/usr/local/gsl/lib -lgsl -lgslcblas -o $(BIN_GSL)
FLAGS_MIN := -L/usr/local/gsl/lib -lgsl -lgslcblas -L/usr/local/minuit/lib -static -llcg_Minuit -o $(BIN_MIN)

# General object dependencies
%.o: %.cpp %.h
	$(CC) -c -o $@ $< $(C_FLAGS)

# Rules for building FitMgst in its various guises
.PHONY: all gsl minuit clean

minuit: $(OBJ_COM) $(OBJ_MIN)
	$(CC) $(OBJ_COM) $(OBJ_MIN) $(FLAGS_MIN)

all:    gsl minuit

gsl:    $(OBJ_COM) $(OBJ_GSL)
	$(CC) $(OBJ_COM) $(OBJ_GSL) $(FLAGS_GSL)

clean:
	rm -f $(OBJ_COM) $(OBJ_GSL) $(OBJ_MIN)

# Explicit declariation of dependencies for src objects that are not satisfied
# by the general declaration (%.o:...) above. i.e. classes that inherit others
# and source files that include headers with different root names.
$(SRC_DIR)/MgstFcnGsl.o: $(SRC_DIR)/MgstFcnGsl.cpp $(SRC_DIR)/MgstFcnGsl.h \
   $(SRC_DIR)/ScriptReader.h $(SRC_DIR)/Deriv.cpp
	$(CC) -c -o $(SRC_DIR)/MgstFcnGsl.o $(SRC_DIR)/MgstFcnGsl.cpp $(C_FLAGS)
  
$(SRC_DIR)/DiffXasSpectrum.o: $(SRC_DIR)/DiffXasSpectrum.cpp \
  $(SRC_DIR)/DiffXasSpectrum.h $(SRC_DIR)/Resample.cpp $(SRC_DIR)/Fourier.cpp
	$(CC) -c -o $(SRC_DIR)/DiffXasSpectrum.o $(SRC_DIR)/DiffXasSpectrum.cpp $(C_FLAGS)
  
$(SRC_DIR)/ThermalSpectrum.o: $(SRC_DIR)/ThermalSpectrum.cpp $(SRC_DIR)/ThermalSpectrum.h \
  $(SRC_DIR)/DiffXasSpectrum.o $(SRC_DIR)/PathFinder.o
	$(CC) -c -o $@ $< $(C_FLAGS)                

$(SRC_DIR)/ScriptEventsGeneral.o: $(SRC_DIR)/ScriptEventsGeneral.cpp $(SRC_DIR)/ScriptReader.h
	$(CC) -c -o $@ $< $(C_FLAGS)

$(SRC_DIR)/ScriptEventsMinuit.o: $(SRC_DIR)/ScriptEventsMinuit.cpp $(SRC_DIR)/ScriptReader.h
	$(CC) -c -o $@ $< $(C_FLAGS)

$(SRC_DIR)/ScriptEventsGsl.o: $(SRC_DIR)/ScriptEventsGsl.cpp $(SRC_DIR)/ScriptReader.h
	$(CC) -c -o $@ $< $(C_FLAGS)

#$(SRC_DIR)/MgstFcnGsl.o: $(SRC_DIR)/MgstFcnGsl.cpp $(SRC_DIR)/MgstFcnGsl.h $(SRC_DIR)/ScriptReader.h
#	$(CC) -c -o $@ $< $(C_FLAGS)
