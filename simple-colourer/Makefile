CPPFLAGS = -Wall -Wno-unused-function -pthread -O3 -march=native

graph_colour: graph_colour.cpp ../graph_colour_solver.h ../graph_colour_solver.cpp
	g++ $(CPPFLAGS) -std=c++14 -Wall -Wno-unused-function -o graph_colour graph_colour.cpp ../graph_colour_solver.cpp -lm

clean:
	rm -f graph_colour
