LIBS= 
COMPILE=g++
FLAGS    = -O3 -lm

fdm: main.o solver.o myData.o  
	$(COMPILE) $(FLAGS) -o fdm main.o solver.o myData.o $(LIBS)

main.o: main.cpp myData.h solver.h
	$(COMPILE) $(FLAGS) -c main.cpp 

solver.o: solver.cpp solver.h 
	$(COMPILE) $(FLAGS) -c solver.cpp

myData.o: myData.cpp myData.h 
	$(COMPILE) $(FLAGS) -c myData.cpp

clean:
	rm -f fdm *.o
 # -fopenmp -lm
