#!/bin/make
#------------------------------------------------------------
#
#  Generic make file for the tools
#
#  usage:   clean - remove objects and executable
#	    rebuild - clean and then make
#
#------------------------------------------------------------

#ifneq ($(DEBUG),TRUE)
#DEBUGCFLAG = -O2
#endif

#TARGET	    = crap
#SOURCES	    = sleevetest.cpp
ENGINEROOT  = /d/programming/gluts/engine
ENGINELIB   = $(ENGINEROOT)/engine.a
DEBUGCFLAG  =
ADDRCFLAGS  =
ADDRINCS    =-I/d/programming/gluts/freeglut/include -I$(ENGINEROOT)/../

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
CFLAGS := -Wall -Werror

CFLAGS += -g -c $(DEBUGCFLAG) -MMD -Wall -DDSO -DLINUX -DFUNCPROTO -DSTYLEAPP -DAL_STYLE -D_REENTRANT \
	    -D_GNU_SOURCE  -DLINUX_IA32 -DTOOLS_COMPILE -DPLATFORM_LINUX=1 \
	    -DNDEBUG \
	    -Wno-comment \
	    -Wno-unused \
	    -Wno-missing-braces \
	    $(ADDRCFLAGS)

LDFLAGS = -L/d/programming//gluts/freeglut/lib/
LDFLAGS += -lgdi32 -lopengl32 -lglu32 -lfreeglut
# -lGLU -lX11 -lXi -lXmu -lfreeglut
#-lGL

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

.PHONY: all $(ENGINELIB) notify

all: $(ENGINELIB) notify $(TARGET)

notify:
	@echo -- building $(TARGET) --

$(ENGINELIB):
	@make -s -C $(ENGINEROOT)


$(COBJS) $(CDEPS):	EXT = c
$(CPPOBJS) $(CPPDEPS):	EXT = cpp
DEP = $(basename $(notdir $@)).d
OBJ = $(basename $@).o

$(TARGET): $(SOURCES) $(OBJECTS)
	@echo -- linking $@ --
	@$(GLD) -o $@ $(OBJECTS) $(ENGINELIB) $(LDFLAGS)
	@cp -f /d/programming/gluts/freeglut/bin/freeglut.dll .

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
	@$(GCC) $(INCS) $(_CFLAGS) -c $(SRC) -o $(OBJ)
#	@if [ -e $(DEP) ]; then \
#	sed -e 's/.*:/$$(OBJDIR)\/$(DEP):/g' $(DEP) > $(OBJDIR)/$(DEP); \
#	rm -f $(DEP); \
#	fi

clean:
	@$(RM) obj/*
	@$(RM) $(TARGET)

# Silently include dependancies
sinclude $(CDEPS) $(CPPDEPS)
