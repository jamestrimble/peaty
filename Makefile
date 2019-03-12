GC2_HEADERS = vc_solver/graph_colour_solver.h
GC2_C_FILES = simple-colourer/graph_colour.cpp \
              vc_solver/graph_colour_solver.cpp

CPPFLAGS = -Wall -Wno-unused-function -pthread
VC_HEADERS = \
             vc_solver/graph.h \
             vc_solver/util.h \
             vc_solver/bitset.h \
             vc_solver/int_stack_without_dups.h \
             vc_solver/sparse_graph.h \
             vc_solver/sequential_solver.h \
             vc_solver/colourer.h \
             vc_solver/root_node_processing.h \
             vc_solver/params.h \
             vc_solver/graph_colour_solver.h
VC_C_FILES = \
             vc_solver/graph.cpp \
             vc_solver/util.cpp \
             vc_solver/sparse_graph.cpp \
             vc_solver/solve_mwc.cpp \
             vc_solver/sequential_solver.cpp \
             vc_solver/colourer.cpp \
             vc_solver/root_node_processing.cpp \
             vc_solver/graph_colour_solver.cpp

all: simple-colourer/graph_colour vc_solver/solve_vc

vc_solver/solve_vc: $(VC_C_FILES) $(VC_HEADERS)
	g++ -O3 -march=native $(CPPFLAGS) -std=c++14 -o vc_solver/solve_vc -DNDEBUG $(VC_C_FILES)

simple-colourer/graph_colour: $(GC2_C_FILES) $(GC2_HEADERS)
	g++ -O3 -march=native $(CPPFLAGS) -std=c++14 -Wall -Wno-unused-function -o simple-colourer/graph_colour $(GC2_C_FILES) -lm

test:
	./test-script
clean:
	rm -f simple-colourer/graph_colour vc_solver/solve_vc
