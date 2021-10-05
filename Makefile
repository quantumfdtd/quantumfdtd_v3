DIR_MAIN       = ./
DIR_BUILD      = $(DIR_MAIN)build/
DIR_H          = $(DIR_BUILD)include/
DIR_OBJ        = $(DIR_BUILD)obj/
DIR_SRC        = $(DIR_BUILD)src/

DEBUG =
OPTIMIZATION = -O2 -funroll-loops -finline-functions
FLOWTRACE =
CFLAGS = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE)
COMPILER = mpicxx
LIBS = -lfftw3_mpi -lfftw3 -lgsl -lgslcblas
INCLUDES = -I build/include/

CPP := $(wildcard $(DIR_SRC)*.cpp)
OBJ  = $(CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)

EXE =\
mpisolve

$(EXE): $(OBJ)
	echo "Linking:   $@ ($(COMPILER))"
	$(COMPILER) -o $@ $^ $(LIBS)

$(DIR_OBJ)%.o: $(DIR_SRC)%.cpp
	@[ -d $(DIR_OBJ) ] || mkdir -p $(DIR_OBJ)
	@echo "Compiling: $< ($(COMPILER))"
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

run:
	mpirun -np 8 --use-hwthread-cpus mpisolve

# clean up misc files
clean:
	@echo "Object files and executable deleted"
	if [ -d "$(DIR_OBJ)" ]; then rm -f $(EXE) $(DIR_OBJ)/*; rmdir $(DIR_OBJ); fi

cleandata:
	@echo "Output and log files deleted"
	./scripts/cleandatafiles.sh

.SILENT :
