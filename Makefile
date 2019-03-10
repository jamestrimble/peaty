GC_HEADERS = simple-graph-colourer/graph_colour_solver.h \
             simple-graph-colourer/c_program_timing.h
GC_C_FILES = simple-graph-colourer/graph_colour.c \
             simple-graph-colourer/graph_colour_solver.c \
             simple-graph-colourer/c_program_timing.c

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
             vc_solver/params.h
VC_C_FILES = \
             vc_solver/graph.c \
             vc_solver/util.c \
             vc_solver/sparse_graph.c \
             vc_solver/solve_mwc.c \
             vc_solver/sequential_solver.c \
             vc_solver/colourer.c \
             vc_solver/root_node_processing.c

all: simple-graph-colourer/graph_colour vc_solver/solve_vc

vc_solver/solve_vc: $(VC_C_FILES) $(VC_HEADERS)
	g++ -O3 -march=native $(CPPFLAGS) -std=c++14 -o vc_solver/solve_vc -DNDEBUG $(VC_C_FILES)

simple-graph-colourer/graph_colour: $(GC_C_FILES) $(GC_HEADERS)
	gcc -O3 -march=native -Wall -std=c11 -Wno-unused-function -o simple-graph-colourer/graph_colour $(GC_C_FILES) -lm

test:
	./test-script
clean:
	rm -f simple-graph-colourer/graph_colour vc_solver/solve_vc
