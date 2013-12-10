# THIRD-PARTY LIBRARIES <<-- CHANGE AS APPROPRIATE -->>

PATH_PROFIL  = /home/bchachua/Programs/ThirdParty/Profil-2.0.8
LIB_PROFIL = -L$(PATH_PROFIL)/lib -lProfilPackages -lProfil -lBias -llr
INC_PROFIL = -I$(PATH_PROFIL)/include

PATH_FILIB = /opt/filib++
LIB_FILIB = -L$(PATH_FILIB)/lib -lprim
INC_FILIB = -I$(PATH_FILIB)/include/
FLAGS_FILIB = -frounding-math -ffloat-store

LIB_LAPACK = -llapack
LIB_CPPUNIT = -lcppunit

# COMPILATION <<-- CHANGE AS APPROPRIATE -->>

# PROF = -pg
OPTIM = -O3
DEBUG = -g
WARN  = -Wall

CC = gcc
CPP = g++
#CPP = icpc
FLAGS_CPP = $(DEBUG) $(OPTIM) $(WARN) $(FLAGS_FILIB)

LINK = $(CPP)
FLAGS_LINK = 
