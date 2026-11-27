#!/bin/sh

#unzip *.zip

pwd
MAIN=$PWD

RAXML="raxml-ng"
chmod u+x $RAXML
echo "Directory where $RAXML is installed: $MAIN"

PROG="SSCPE.pl"
echo "Directory where $PROG is installed: $MAIN"

# Compile the TNM program
DIR="DIR_TNM"
echo "Compiling tnm in folder $DIR and moving to $MAIN"
rm -rf $DIR
mkdir $DIR
cp -f tnm.zip $DIR
cd $DIR
unzip tnm.zip
make
mv -f tnm ../
mv -f Mutation_para.in ../
mv Input_TNM.in  ../
rm -rf *.o
cd ..

# Compile the Prot_evol program
DIR="DIR_PROT_EVOL"
echo "Compiling Prot_evol in folder $DIR and moving to current folder"
rm -rf $DIR
mkdir $DIR
cp Prot_evol.zip $DIR
cp script_dir_Prot_evol.pl $DIR
cd $DIR
chmod u+x script_dir_Prot_evol.pl
unzip Prot_evol.zip
echo "Running script_dir_Prot_evol.pl"
./script_dir_Prot_evol.pl
make
mv -f Prot_evol ../
mv structures.in ../
mv Input_Prot_evol.in  ../
rm -f *.o
cd ..

echo "Assigning dependencies of $PROG"
chmod u+x $PROG
chmod u+x script_write_dir.pl
echo "Running the program script_write_dir.pl"
./script_write_dir.pl ./

echo "Executable file $PROG is in folder: "
pwd

