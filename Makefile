CXX = g++
CXXFLAGS = -Wall -std=c++11 $(DEBUG)

burnman: burnman.o models.o extended_margules.o slb.o aa.o birch_murnaghan.o electronic.o einstein.o debye.o global.o
	$(CXX) $(CXXFLAGS) -g -o burnman models.o extended_margules.o burnman.o slb.o aa.o birch_murnaghan.o electronic.o einstein.o debye.o global.o

burnman.o: burnman.cpp burnman.hpp
	$(CXX) $(CXXFLAGS) -g -c burnman.cpp 

models.o: models.cpp models.hpp
	$(CXX) $(CXXFLAGS) -g -c models.cpp

extended_margules.o: solution_models/extended_margules.cpp solution_models/extended_margules.hpp
	$(CXX) $(CXXFLAGS) -g -c solution_models/extended_margules.cpp

electronic.o: eos/electronic.cpp eos/electronic.hpp
	$(CXX) $(CXXFLAGS) -g -c eos/electronic.cpp 

einstein.o: eos/einstein.cpp eos/einstein.hpp
	$(CXX) $(CXXFLAGS) -g -c eos/einstein.cpp

debye.o: eos/debye.cpp eos/debye.hpp
	$(CXX) $(CXXFLAGS) -g -c eos/debye.cpp

slb.o: eos/slb.cpp eos/slb.hpp
	$(CXX) $(CXXFLAGS) -g -c eos/slb.cpp

aa.o: eos/aa.cpp eos/aa.hpp
	$(CXX) $(CXXFLAGS) -g -c eos/aa.cpp

birch_murnaghan.o: eos/birch_murnaghan.cpp eos/birch_murnaghan.hpp
	$(CXX) $(CXXFLAGS) -g -c eos/birch_murnaghan.cpp

global.o: global.cpp global.hpp
	$(CXX) $(CXXFLAGS) -g -c global.cpp

clean: 
	$(RM) count *.o *~ 

tidy:
	rm *.o
