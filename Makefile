#use DPRINT to run in debug mode with lost of messages
#PARAM=-O3 -DPRINT -g 
PARAM=-O3

uniform_sampling:matrix_base.o StringTokenizer.o random.o uniform_freq.cpp database.h pattern.h adj_matrix.h pattern_factory.h subgraph_iso.h graph_iso_check.h uniform_freq_random_walk.h
	g++ $(PARAM) -o uniform_sampling.out matrix_base.o StringTokenizer.o random.o uniform_freq.cpp


matrix_base.o: matrix_base.cpp matrix_base.h
	g++ $(PARAM) -c matrix_base.cpp matrix_base.h

StringTokenizer.o : StringTokenizer.cpp StringTokenizer.h
	g++ $(PARAM) -c StringTokenizer.cpp StringTokenizer.h

random.o: random.h random.cpp
	g++ $(PARAM) -c random.cpp random.h

clean:
	rm -f *.gch
	rm -f *.o
	rm -f *.out
