########################################################################
# - Compile ProteinEvolverABC -
########################################################################
# Makefile ProteinEvolverABC
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
	@cp source/ProteinEvolverProtABC/ProteinEvolverProtABC2.0.0 bin/
	@echo "Done!"

	@echo "Compiling SSCPE .."
	@unzip source/SSCPE-main.zip -d ./source/
	@unzip -o source/SSCPE-main/SSCPE.zip -d ./source/SSCPE-main/
	@cp -r source/SSCPE-main/* .
	@chmod u+x script_install_SSCPE.sh
	@./script_install_SSCPE.sh
	@cp Prot_evol bin/
	@cp tnm bin/
	@rm SSCPE.zip
	@mkdir SSCPE
	@mv *.in SSCPE
	@mv DIR_PROT_EVOL SSCPE
	@mv DIR_TNM SSCPE
	@mv jtt.txt SSCPE
	@mv lg.txt SSCPE
	@mv k2_mat.pl SSCPE
	@mv Prot_evol SSCPE
	@mv qsubmit.pl SSCPE
	@mv README_SSCPE SSCPE
	@mv script_* SSCPE
	@mv SSCPE.pl SSCPE
	@mv tnm SSCPE
	@mv wag.txt SSCPE
	@mv *.zip SSCPE
	@mv README.md SSCPE
	@mv raxml-ng SSCPE
	@mv ALI SSCPE
	@mv PDB SSCPE
	@mv SSCPE bin
	@rm -r source/SSCPE-main
	@echo "Done!"

	@echo "Compiling Prot_evol3 .."
	@unzip source/SSCPE-main.zip -d ./source/
	@unzip -o source/SSCPE-main/SSCPE.zip -d ./source/SSCPE-main/
	@mkdir source/SSCPE-main/Prot_evol3
	@unzip -o source/SSCPE-main/Prot_evol.zip -d ./source/SSCPE-main/Prot_evol3/
	@cp ./source/PathsProt_evol.pl ./source/SSCPE-main/Prot_evol3/
	@perl ./source/SSCPE-main/Prot_evol3/PathsProt_evol.pl ./source/SSCPE-main/Prot_evol3/Prot_evol.c
	$(MAKE) -C source/SSCPE-main/Prot_evol3
	@cp source/SSCPE-main/Prot_evol3/Prot_evol source/SSCPE-main/Prot_evol3/Prot_evol3
	@cp source/SSCPE-main/Prot_evol3/Prot_evol3 bin/
	@rm -r source/SSCPE-main
	@echo "Done!"

	@echo "Compiling Prot_evol_ProtASR2 .."
	$(MAKE) -C source/Prot_evol_ProtASR2_src
	@cp source/Prot_evol_ProtASR2_src/Prot_evol source/Prot_evol_ProtASR2_src/Prot_evol_ProtASR2_exe
	@cp source/Prot_evol_ProtASR2_src/Prot_evol_ProtASR2_exe bin/
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
	@rm -f bin/ProteinEvolverProtABC2.0.0
	@rm -f bin/tnm
	@rm -f bin/Prot_evol
	@rm -f -r bin/SSCPE
	@rm -f bin/Prot_evol3
	$(MAKE) -C source/Prot_evol_ProtASR2_src clean
	@rm -f bin/Prot_evol_ProtASR2_exe
	@rm -f source/Prot_evol_ProtASR2_src/Prot_evol_ProtASR2_exe
	
	@echo "Done!"

