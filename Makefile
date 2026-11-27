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
	@cp -r source/SSCPE/SOURCE/* .
	@chmod u+x script_install_SSCPE.sh
	@./script_install_SSCPE.sh
	@cp Prot_evol bin/
	@cp tnm bin/
	@cp raxml-ng bin/

	@mkdir -p SSCPE
	@mv Input_Prot_evol.in SSCPE
	@mv Input_TNM.in SSCPE
	@mv Mutation_para.in SSCPE
	@mv structures.in SSCPE
	@mv DIR_PROT_EVOL SSCPE
	@mv DIR_TNM SSCPE
	@mv TNM_DATA SSCPE
	@mv ALI SSCPE
	@mv PDB SSCPE
	@mv jtt.txt SSCPE
	@mv lg.txt SSCPE
	@mv wag.txt SSCPE
	@mv k2_mat.pl SSCPE
	@mv script_dir_Prot_evol.pl SSCPE
	@mv script_run_models.pl SSCPE       
	@mv script_write_dir.pl SSCPE
	@mv qsubmit.pl SSCPE  
	@mv script_read_tree.pl SSCPE     
	@mv script_run_SSCPE_models.pl SSCPE
	@mv SSCPE.pl SSCPE
	@mv script_install_SSCPE.sh SSCPE
	@mv Prot_evol.zip SSCPE
	@mv tnm.zip SSCPE
	@mv Prot_evol SSCPE
	@mv tnm SSCPE
	@mv raxml-ng SSCPE
	@mv README_SSCPE SSCPE
	@mv SSCPE bin
	@echo "Done!"

	@echo "Compiling Prot_evol3 .."
	@mkdir -p source/SSCPE/Prot_evol3
	@unzip -o source/SSCPE/SOURCE/Prot_evol.zip -d ./source/SSCPE/Prot_evol3/
	@cp ./source/PathsProt_evol.pl ./source/SSCPE/Prot_evol3/
	@perl ./source/SSCPE/Prot_evol3/PathsProt_evol.pl ./source/SSCPE/Prot_evol3/Prot_evol.c
	$(MAKE) -C source/SSCPE/Prot_evol3
	@cp source/SSCPE/Prot_evol3/Prot_evol source/SSCPE/Prot_evol3/Prot_evol3
	@cp source/SSCPE/Prot_evol3/Prot_evol3 bin/
	@echo "Done!"

	@echo "Compiling Prot_evol_ProtASR2 .."
	$(MAKE) -C source/Prot_evol_ProtASR2_src
	@cp source/Prot_evol_ProtASR2_src/Prot_evol source/Prot_evol_ProtASR2_src/Prot_evol_ProtASR2_exe
	@cp source/Prot_evol_ProtASR2_src/Prot_evol_ProtASR2_exe bin/
	@echo "Done!"

	@echo ""
	@echo "Note that R will require the following libraries: lattice, MCMCpack, ape, graphics, abc"
	@echo "These libraries can be installed from R by typing:"
	@echo 'install.packages(c("lattice", "MCMCpack", "ape", "abc"))'
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
	@rm -rf source/SSCPE/Prot_evol3
	@rm -f bin/raxml-ng
	
	@echo "Done!"

