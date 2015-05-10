########################################
# MAKEFILE TO COMPILE SPHPPM
###################################
KNOWN_SYSTEM=no

SYSTEM=jostalkingcomputer

ifeq ($(SYSTEM), ice)
-include Makefile.ice
KNOWN_SYSTEM=yes
endif
ifeq ($(SYSTEM), mpp)
-include Makefile.mpp
KNOWN_SYSTEM=yes
endif
ifeq ($(SYSTEM), miyuko)
-include Makefile.miyuko
KNOWN_SYSTEM=yes
endif
ifeq ($(SYSTEM), supermuc)
-include Makefile.supermuc
KNOWN_SYSTEM=yes
endif
ifeq ($(SYSTEM), laptop)
-include Makefile.laptop
KNOWN_SYSTEM=yes
endif
ifeq ($(SYSTEM), mitsuko)
-include Makefile.mitsuko
KNOWN_SYSTEM=yes
endif
ifeq ($(SYSTEM), jostalkingcomputer)
-include Makefile.jostalkingcomputer
KNOWN_SYSTEM=yes
endif


#----------------------------------------------------------------------
# Folders
#----------------------------------------------------------------------
SRC_DIR = src
OBJ_DIR = obj
MOD_DIR = mod

#----------------------------------------------------------------------
#  Makefile for sph_ppm client 
#----------------------------------------------------------------------
TARGET      			= mcf$(TAG)

#----------------------------------------------------------------------
#  Compiler settings (Linux)
#----------------------------------------------------------------------
DEFINE       	 = -D__MPI 
DIM_FLAG    	 = -D__DIM=$(DIM)
KERNEL_FLAG  	 = -D__KERNEL=$(KERNEL)
OTHER            = -D__COLLOID_SEPARATE_FILE -D__COLLOID_NOACCE
COMPILERSETTINGS = $(DEFINE) $(DIM_FLAG) $(KERNEL_FLAG) $(OTHER)

#----------------------------------------------------------------------
#  Compiler settings
#----------------------------------------------------------------------
INC = -I$(PPM_BASE)/include/ppm -I$(MOD_DIR) -I. -Isrc/boundary -Isrc/colloid -Isrc/control -Isrc/debug -Isrc/io -Isrc/kernel -Isrc/marching -Isrc/particles -Isrc/physics -Isrc/random -Isrc/rhs -Isrc/stateEquation -Isrc/statistic -Isrc/technique -Isrc/tool -Isrc/tracers
FFLAGS1	= $(INC) $(COMPILERSETTINGS) $(COMPILEROPTION) 
#IFORT FLAG:
#FFLAGS2	= -module $(MOD_DIR)
FFLAGS2	= -J $(MOD_DIR)
#GFORTRAN FLAG:
FFLAGS	= $(FFLAGS1) $(FFLAGS2)
LIBS    = -lppm -lmetis -lmpich
LDFLAGS	= -L$(PPM_BASE)/lib -L$(METIS_BASE)/lib -L$(MPI_BASE)/lib 

#----------------------------------------------------------------------
#  Object files / RULES                                                        
#----------------------------------------------------------------------
-include Make.depend

#--------------------------------------------------------------------
#  PATTERN RULES
#-------------------------------------------------------------------
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.F90
	$(FC) -c $(FCFLAGS) $(FFLAGS) $< -o $@

#----------------------------------------------------------------------
#  All
#----------------------------------------------------------------------
all: $(TARGET)

#--------------------------------------------------------------------
#  Make dependency has to be done if Make.depend does not exist
#--------------------------------------------------------------------
depend: 
	$(MAKEDEPF90_BASE)/bin/makedepf90 -I src/ -o $(TARGET) -b $(OBJ_DIR) $(SRC_DIR)/*.F90 > Make.depend;\
	mkdir -p $(OBJ_DIR);\
	mkdir -p $(MOD_DIR)

#----------------------------------------------------------------------
#  Clean                                                               
#----------------------------------------------------------------------
clean:
	rm -rf $(OBJ_DIR)/*;\
	rm -rf $(MOD_DIR)/*

err:
	$(error aborting);
