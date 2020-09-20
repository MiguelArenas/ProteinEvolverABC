# Makefile ProteinEvolverABC_GUI.java
# 	Compile: make all
# 	Clean: make clean
#

JAVAC=javac
sources = $(wildcard *.java)
classes = $(sources:.java=.class)

all: $(classes)


clean :
	rm -f *.class
	rm -f *.jar

%.class : %.java
	$(JAVAC) $<
	jar cmf ProteinEvolverABC_GUI.txt ProteinEvolverABC_GUI.jar *.class
	rm -f *.class