OBJ=main.o util_funcs.o SuffixArray.o BranchPointGroups.o Reads.o GenomeMapper.o string.o
EXE=ICSMuFin
CXX=g++
COMPFLAGS=-Wall -ggdb -MMD -std=c++11 -pthread 
OBJDIR=./objects/

$(EXE):$(OBJ)
	$(CXX) $(COMPFLAGS) $(OBJ) -o $(EXE) -lz -lboost_regex
	mv *.o *.d ./obj

%.o: %.cpp
	$(CXX) $(COMPFLAGS) -c $<
-include $(OBJ:.o=.d)	

.PHONY: clean

clean:
	rm ./*.o
	rm ./*.d

cleaner:
	rm ./$(EXE)

