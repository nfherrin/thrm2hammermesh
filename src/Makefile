# Defining variables
FC:=gfortran -Wall -Wextra -Wconversion -fcheck=all
EXE:=thrm2hammermesh
EXT:=.exe

all:			TYPE :=
opt:			TYPE :=
debug:    TYPE := _Debug
coverage: TYPE := _Coverage

all:      OPT := -O3
debug:    OPT := -O0 -g -ffpe-summary=all
coverage: OPT := -O0 -g --coverage

#-----------------------------------------------------------------------------------
#Add source files as necessary
#If the files are not comiled using the genieric rules, add commands under their
#rule declaration. Add these items to FILTER
#-----------------------------------------------------------------------------------

SRC := thrm2hammer.f95

OBJ := $(SRC:.f95=.o)
MOD := $(OBJ:.o=.mod)

FILTER := \

OBJ_FILTER := $(FILTER:.f95=.o)
MOD_FILTER := $(FILTER:.f95=.mod)

#-----------------------------------------------------------------------------------
#Phony targets for cleaning and building
#-----------------------------------------------------------------------------------
.PHONY: all clean reset

all: $(EXE)

debug: $(EXE)

coverage: $(EXE)

#Intended to clean up compilation artifacts but leave executable & coverage
clean:
	rm -f $(OBJ) $(CONT_OBJ) $(LIB_OBJ)
	rm -f $(MOD) $(CONT_MOD) $(LIB_MOD)
	rm -f *$(EXT)
	rm -f *.mod
	rm -f *.o
	rm -f $(COMP_DEP)

#Intended to reset directory to fresh state with no exe or artifacts
reset: clean
	rm -f *.gcno *.gcda
	rm -f -r $(EXE)*.dSYM
	rm -f ../$(EXE)*
	rm -f *$(EXT)
	rm -f *.mod
	rm -f *.o

#-----------------------------------------------------------------------------------
#Generics for source files
#-----------------------------------------------------------------------------------
$(filter-out $(OBJ_FILTER), $(OBJ)): %.o:	%.f95
	$(FC) -c $(OPT) $<

$(filter-out $(MOD_FILTER), $(MOD)):	%.mod:	%.f95
	$(FC) -c $(OPT) $<

$(EXE): $(OBJ) $(CONT_OBJ) $(LIB_OBJ)
	$(FC) -o $@$(TYPE)$(EXT) $(OPT) $(OBJ) $(CONT_OBJ) $(LIB_OBJ) -llapack
	mv ./$(EXE)$(TYPE)$(EXT) ../

#-----------------------------------------------------------------------------------
#Dependency List
#Use [gfortran -M -cpp *.f95] repeatedly until clean compile to update rules below
#-----------------------------------------------------------------------------------
thrm2hammer.o: thrm2hammer.f95
