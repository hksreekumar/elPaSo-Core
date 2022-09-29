
# =============================================================================
#   Makefile used to build gzstream
# =============================================================================

AR          =ar r
RM          =rm -rf

SRCDIR  = src
vpath %.cpp $(SRCDIR)
vpath %.c   $(SRCDIR)
SRCS    = gzstream.cpp
OBJDIR  = ./src
OBJLIB  = $(SRCS:.cpp=.o)
HDRS    = gzstream.h$
LINKLIB = $(addprefix $(OBJDIR),/$(OBJLIB))

LIB = libgzstream.a

# ------------------------------------------------------------
# -- rules
%.h: %.h

%.o: %.cpp
	@echo compiling $< ...
	@$(CXX) $(PCC_FLAGS) -I./include $(COPTFLAGS) $(CFLAGS) $(CCPPFLAGS) -c $< -o $(OBJDIR)/$(@F)

# ------------------------------------------------------------
# -- targets
all: $(OBJLIB)
	@echo create $(LIB)
	@$(AR) $(AR_FLAGS) ./lib/$(LIB) $(LINKLIB)
	@echo $(LIB) is now up-to-date
	@echo build of gzstream finished

clean:
	@echo cleaning gzstream ...
	@$(RM) ./lib/$(LIB)
	@$(RM) ./$(OBJDIR)/*.o

new: wipe all

wipe: distclean

distclean: clean

