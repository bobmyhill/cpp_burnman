CXX = g++
CXXFLAGS = -Wall -std=c++11 $(DEBUG)

btst: btst.o fsolve.o
	$(CXX) $(CXXFLAGS) -g -o btst btst.o fsolve.o

btst.o: btst.cpp
	$(CXX) $(CXXFLAGS) -g -c btst.cpp

fsolve.o: fsolve.cpp
	$(CXX) $(CXXFLAGS) -g -c fsolve.cpp 

clean: 
	$(RM) count *.o *~ 

tidy:
	rm *.o
