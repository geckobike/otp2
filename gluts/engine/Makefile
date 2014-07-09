#!/bin/make
#------------------------------------------------------------
#
#  Generic make file for the tools
#
#  usage:   clean - remove objects and executable
#	    rebuild - clean and then make
#
#  Alex Weighell
#
#------------------------------------------------------------

#ifneq ($(DEBUG),TRUE)
#DEBUGCFLAG = -O2
#endif

TARGET	   = engine.a
SOURCES    = engine.cpp

DEBUGCFLAG  =
ADDRCFLAGS  =
ADDRINCS    =-I/usr/include/GL

ifndef SOURCES
.PHONY: srcerror
srcerror:
	@echo "-----------------------"
	@echo "! No Sources Specified"
	@echo "-----------------------"
endif

ifndef TARGET
.PHONY: targeterror
targeterror:
	@echo "----------------------"
	@echo "! No Target Specified"
	@echo "----------------------"
endif

OBJDIR	    = obj

# Include directories

INCS	=   -I. -I/usr/include $(ADDRINCS)

# Flags for compiler (add -pg for profile)
CFLAGS += -g -c $(DEBUGCFLAG) -MMD -Wall -DDSO -DLINUX -DFUNCPROTO -DSTYLEAPP -DAL_STYLE -D_REENTRANT \
	    -D_GNU_SOURCE  -DLINUX_IA32 -DTOOLS_COMPILE -DPLATFORM_LINUX=1 \
	    -DNDEBUG \
	    -Wno-comment \
	    -Wno-unused \
	    $(ADDRCFLAGS)
	    
LDFLAGS = -lglut

GCC     := gcc
GLD     := g++
RM      := /bin/rm -f

# Derived variables

CFILES	=   $(filter %.c, $(SOURCES))
CPPFILES =  $(filter %.cpp, $(SOURCES))
COBJS	=   $(foreach file,$(CFILES),$(OBJDIR)/$(subst .c,.o,$(notdir $(file))))
CPPOBJS	=   $(foreach file,$(CPPFILES),$(OBJDIR)/$(subst .cpp,.o,$(notdir $(file))))
OBJECTS =   $(COBJS) $(CPPOBJS)
CDEPS	=   $(COBJS:%.o=%.d)
CPPDEPS	=   $(CPPOBJS:%.o=%.d)
 
# Rules

.PHONY: all notify

all: $(TARGET) notify

notify:
	@echo "-- building engine --"


$(COBJS) $(CDEPS):	EXT = c
$(CPPOBJS) $(CPPDEPS):	EXT = cpp
DEP = $(basename $(notdir $@)).d
OBJ = $(basename $@).o


$(TARGET): $(SOURCES) $(OBJECTS)
	@echo -- linking $@ --
#	@$(GLD) -o $@ $(OBJECTS) $(LDFLAGS)
	@$(RM) $(TARGET)
	ar rcs $(TARGET) $(OBJECTS) $(LIBS)

%.d:
	@mkdir -p $(@D)
	@$(RM) $*.o
	@echo "# empty" > $@

%.h:
	@

$(COBJS) $(CPPOBJS): SRCBASE = $(basename $(notdir $@)).$(EXT)
$(COBJS) $(CPPOBJS): SRC = $(filter %/$(SRCBASE) $(SRCBASE), $(SOURCES))
$(COBJS): _CFLAGS = $(CFLAGS) -Werror-implicit-function-declaration
$(CPPOBJS): _CFLAGS = $(CFLAGS)

$(COBJS) $(CPPOBJS):
	@echo -- compiling $(SRC) --
	$(GCC) $(INCS) $(_CFLAGS) -c $(SRC) -o $(OBJ)
#	@if [ -e $(DEP) ]; then \
#	sed -e 's/.*:/$$(OBJDIR)\/$(DEP):/g' $(DEP) > $(OBJDIR)/$(DEP); \
#	rm -f $(DEP); \
#	fi

clean:
	@$(RM) obj/*
	@$(RM) $(TARGET)

# Silently include dependancies
sinclude $(CDEPS) $(CPPDEPS)
