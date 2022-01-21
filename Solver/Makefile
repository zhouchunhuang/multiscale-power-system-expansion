EXDIR = /home/mchen/zhouchun/CEP/Model/Model
GRBPATH = /apps/gurobi/gurobi-9.1.2
EXINC         = $(EXDIR)/inc
EXSRCCPP      = $(EXDIR)/src
EXDATA = $(EXDIR)/../../runfolder
CPP_EX = Model

# ---------------------------------------------------------------------
# Compiler options
# ---------------------------------------------------------------------
CCC = mpicxx
CCOPT = -m64 -O -fPIC -fexceptions -D_NDEBUG
CCFLAGS = $(CCOPT) -I$(GRBPATH)/include -I$(EXINC)

# ---------------------------------------------------------------------
# Link options and libraries
# ---------------------------------------------------------------------
CPPLIB   = -L$(GRBPATH)/lib -lgurobi_g++5.2 -lgurobi91 -lm

#------------------------------------------------------------
#  make all      : to compile the code.
#  make execute  : to compile and execute the code.
#------------------------------------------------------------
all:
	make all_cpp

execute: all
	make execute_cpp

all_cpp: $(CPP_EX)

execute_cpp: $(CPP_EX)
	 ./Model $(EXDATA) -d
 # ------------------------------------------------------------
clean :
	/bin/rm -rf *.o *~ *.class
	/bin/rm $(CPP_EX)
	/bin/rm -rf *.mps *.ord *.sos *.lp *.sav *.net *.msg *.log *.clp

# ------------------------------------------------------------
Model: main.o model.o parameters.o read.o Benders.o ColGen.o direct.o NestColGen.o test.o
	$(CCC) $(CCFLAGS) main.o model.o parameters.o read.o Benders.o ColGen.o direct.o NestColGen.o test.o -o $(CPP_EX) $(CPPLIB)
main.o: $(EXSRCCPP)/main.cpp $(EXINC)/model.h
	$(CCC) $(CCFLAGS) -c $(EXSRCCPP)/main.cpp $(CPPLIB)
parameters.o: $(EXSRCCPP)/parameters.cpp $(EXINC)/parameter.h
	$(CCC) $(CCFLAGS) -c $(EXSRCCPP)/parameters.cpp $(CPPLIB)
read.o: $(EXSRCCPP)/read.cpp $(EXINC)/model.h
	$(CCC) $(CCFLAGS) -c $(EXSRCCPP)/read.cpp $(CPPLIB)
model.o: $(EXSRCCPP)/model.cpp $(EXINC)/model.h
	$(CCC) $(CCFLAGS) -c $(EXSRCCPP)/model.cpp $(CPPLIB)
Benders.o: $(EXSRCCPP)/Benders.cpp $(EXINC)/model.h
	$(CCC) $(CCFLAGS) -c $(EXSRCCPP)/Benders.cpp $(CPPLIB)
ColGen.o: $(EXSRCCPP)/ColGen.cpp $(EXINC)/model.h
	$(CCC) $(CCFLAGS) -c $(EXSRCCPP)/ColGen.cpp $(CPPLIB)
direct.o: $(EXSRCCPP)/direct.cpp $(EXINC)/model.h
	$(CCC) $(CCFLAGS) -c $(EXSRCCPP)/direct.cpp $(CPPLIB)
NestColGen.o: $(EXSRCCPP)/NestColGen.cpp $(EXINC)/model.h
	$(CCC) $(CCFLAGS) -c $(EXSRCCPP)/NestColGen.cpp $(CPPLIB)
test.o: $(EXSRCCPP)/test.cpp $(EXINC)/model.h
	$(CCC) $(CCFLAGS) -c $(EXSRCCPP)/test.cpp $(CPPLIB)

# mode: makefile
# End:
