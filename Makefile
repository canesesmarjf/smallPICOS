COMPILER=g++
BIN=bin
SRC=src
OBJ=obj

INCL=includes/
INCL_ARMA=lib_ARMA/include/
INCL_HDF5=lib_HDF5/include/
LIBS_HDF5=lib_HDF5/lib/
LIBS_ARMA=lib_ARMA/lib/

CPPFLAGS= -O3 -I$(INCL) -I$(INCL_ARMA) -I$(INCL_HDF5)

BINARIES = $(addprefix $(OBJ)/,\
main.o	       	\
initialize.o		\
types.o         \
)

# NOTE:
# If I add -std=c++11 or ++17 to the linker step, then I really need to link the dynamic library to main.o, otherwise the compilation process will fail.
# If we succesfully compile the program using the linking process when we include -std=c++XX, we will still get and error when we attempt to run the binary because the loader will not find the location of the dylib correctly even if we have specfied it on the compilation process.
# To solve this in MacOS we need to use the following command:
# install_name_tool -change @rpath/libarmadillo.9.dylib arma_libs/lib/libarmadillo.9.dylib "bin/"$OUTPUT
# https://www.fullstaq.com/knowledge-hub/blogs/an-alternative-to-macos-dyld-library-path
# To solve this in LINUX, you need to add the search path on LD_LIBRARY_PATH

smallPICOS: $(BIN)/smallPICOS++

$(BIN)/smallPICOS++: $(BINARIES)
	$(COMPILER) -o $@ $^ -L$(LIBS_ARMA) -larmadillo -std=c++17 -DARMA_USE_HDF5 -L$(LIBS_HDF5) -lhdf5

	# For MacOS use the following:
	install_name_tool -add_rpath @executable_path/../$(LIBS_ARMA) ./bin/smallPICOS++
	install_name_tool -add_rpath @executable_path/../$(LIBS_HDF5) ./bin/smallPICOS++

	# The above command adds an additional search path to the binary so it can find the dylib located where we have it in this project.

$(OBJ)/%.o: $(SRC)/%.cpp
	$(COMPILER) -c $< -o $@ $(CPPFLAGS) -std=c++17

clean: cleanBin cleanObj

cleanBin:
	-rm -r $(BIN)/PICOS++

cleanObj:
	-rm -r $(OBJ)/*.o
