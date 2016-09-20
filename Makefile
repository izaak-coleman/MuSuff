OBJ=main.o helper_functions.o SuffixArray.o BranchPointGroups.o Reads.o GenomeMapper.o
EXE=ICSMuFin
CXX=g++
COMPFLAGS=-Wall -ggdb -MMD -std=c++11 -pthread
OBJDIR=./objects/

$(EXE):$(OBJ)
	$(CXX) $(COMPFLAGS) $(OBJ) -o $(EXE) -lz

%.o: %.cpp
	$(CXX) $(COMPFLAGS) -c $<

-include $(OBJ:.o=.d)	

.PHONY: clean

clean:
	rm ./*.o
	rm ./*.d

cleaner:
	rm ./$(EXE)

