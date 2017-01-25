OBJ=main.o string.o build_suffix_array.o parse_coordinate_data.o
EXE=extract
CXX=g++
COMPFLAGS=-Wall -ggdb -MMD -std=c++11 -pthread 
OBJDIR=./objects/

$(EXE):$(OBJ)
	$(CXX) $(COMPFLAGS) $(OBJ) -o $(EXE) -lz -lboost_regex

%.o: %.cpp
	$(CXX) $(COMPFLAGS) -c $<
-include $(OBJ:.o=.d)	

.PHONY: clean

clean:
	rm ./*.o
	rm ./*.d

cleaner:
	rm ./$(EXE)

