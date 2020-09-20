########################################################################
# - Compile ProteinEvolverABC -
########################################################################
# Makefile ProtABC
# 	Compile: make all
# 	Clean: make clean
#
PROGRAM_NAME = ProteinEvolverABC


help:
	@echo ""
	@echo "To use the $(PROGRAM_NAME) Makefile, type"
	@echo "     make all           to compile the whole package but not install it"
	@echo "                          or remove the object files"
	@echo "     make clean         to remove all object files and executables from the"
	@echo "                          current directory"
	@echo " "


all: 
	@echo "Compiling Phi .."
	$(MAKE) -C source/PhiPack/src clean
	$(MAKE) -C source/PhiPack/src Phi
	@cp source/PhiPack/Phi bin/
	@echo "Done!"
	@echo ""
	@echo "Compiling ProteinEvolverProtABC .."
	$(MAKE) -C source/ProteinEvolverProtABC clean
	$(MAKE) -C source/ProteinEvolverProtABC all
	@cp source/ProteinEvolverProtABC/ProteinEvolverProtABC1.2.0 bin/
	@echo "Done!"
	@echo ""
	@echo "Compiling ProteinEvolverABC_GUI .."
	$(MAKE) -C source/ProteinEvolverABC_GUI clean
	$(MAKE) -C source/ProteinEvolverABC_GUI all
	@cp source/ProteinEvolverABC_GUI/ProteinEvolverABC_GUI.jar bin/
	@echo "Done!"
	@echo ""
	@echo "Note that R will require the following libraries: lattice, MCMCpack, ape, graphics, abc"
	@echo "These libraries can be installed from R by typing:"
	@echo "install.packages("lattice")"
	@echo "install.packages("MCMCpack")"
	@echo "install.packages("ape")"
	@echo "install.packages("graphics")"
	@echo "install.packages("abc")"
	@echo "See the documentation for additional details about $(PROGRAM_NAME)"
	@echo "Compilation completed!"
	@echo ""
	

clean:
	@echo "Removing executables .."
	$(MAKE) -C source/PhiPack/src clean
	@rm -f bin/Phi
	@rm -f source/PhiPack/Phi
	$(MAKE) -C source/ProteinEvolverProtABC clean
	@rm -f bin/ProteinEvolverProtABC1.2.0
	$(MAKE) -C source/ProteinEvolverABC_GUI clean
	@rm -f bin/ProteinEvolverABC_GUI.jar
	@echo "Done!"

