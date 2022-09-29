#
# die .c - Dateien werden NICHT kompiliert und in die
# Bibliothek eingebunden !!
#

#include ../../../Makefile.sys

AR          =ar r
RM          =rm -rf
CP          =cp


AKUSTIK_INCLUDES += -I./include
SRCDIR    = ./src
OBJDIR    = ./src

vpath %.cpp $(SRCDIR)
vpath %.c   $(SRCDIR)
vpath %.o   $(SRCDIR)

# ----------------------------------------------------------------------------
#   compilation rules
%.h: %.h

%.o: %.cpp
	@echo compiling '( $(BUILDINFO) ) '$< ...
	@$(CXX) $(PCC_FLAGS) $(AKUSTIK_DEFS) $(AKUSTIK_INCLUDES) $(COPTFLAGS) $(CFLAGS) $(CCPPFLAGS) -c $< -o $(OBJDIR)/$(@F)

%.d: %.cpp
	@echo creating dependencies '( $(BUILDINFO) )' $< ...
	@$(CXX) $(DEPFLAGS) $(AKUSTIK_DEFS)  $(PCC_FLAGS) $(AKUSTIK_INCLUDES) $(COPTFLAGS) $(CFLAGS) $(CCPPFLAGS) $< > $(DEPDIR)/$(@F)


DEST_DIR    = ./

LIB       =libxmlio.a

HDR_DEST_DIR=./
LIB_DEST_DIR=./

CSRCS= xmlfile_c.c xmlinput_c.c 
CXXSRCS= xmlfile.cpp xmlinput.cpp xmloutput.cpp
OBJLIB=  $(CXXSRCS:.cpp=.o) $(CSRCS:.c=.o)
HDRS= $(CXXSRCS:.cpp=.h) xmlconfig.h xmlinputp.h xmlstream.h
LINKLIB = $(addprefix $(OBJDIR)/,$(OBJLIB))

all: $(OBJLIB)
	@echo create $(LIB)
	@$(AR) $(LIB) $(LINKLIB)
	@echo $(LIB) is now up-to-date

clean:
	@echo cleaning xmlio-0.92
	@$(RM) $(LINKLIB)

wipe: clean
	@echo delete $(LIB)
	@$(RM) $(LIB)

xmlfile_c.o: xmlfile_c.c 
	@echo compiling '( $(PETSC_ARCH) )' src/xmlfile_c.c ...
	@$(CC) $(CC_FLAGS) $(AKUSTIK_INCLUDES) $(COPTFLAGS) $(CFLAGS) $(CCPPFLAGS)  -c src/xmlfile_c.c -o $(OBJDIR)/xmlfile_c.o

xmlinput_c.o: xmlinput_c.c 
	@echo compiling '( $(PETSC_ARCH) )' src/xmlinput_c.c ...
	@$(CC) $(CC_FLAGS) $(AKUSTIK_INCLUDES) $(COPTFLAGS) $(CFLAGS) $(CCPPFLAGS) -c src/xmlinput_c.c -o $(OBJDIR)/xmlinput_c.o

