import sys
import subprocess

def read_instance(lines):
    edges = []
    for line in lines:
        tokens = line.strip().split()
        if tokens:
            if tokens[0] == "p":
                n = int(tokens[2])
            elif tokens[0] != "c":
                v, w = int(tokens[0]) - 1, int(tokens[1]) - 1
                if v > w:
                    v, w = w, v
                edges.append((v, w))
    return n, edges

def write_graph_to_proc(proc, n, edges):
    proc.stdin.write("p edge {} {}\n".format(n, len(edges)))
    for v, w in edges:
        proc.stdin.write("e {} {}\n".format(v+1, w+1))
    proc.stdin.close()

def get_colouring_bound(vv, new_edges, fractional_level=1, seed=1, time_limit=999999):
    proc = subprocess.Popen(
            ['simple-graph-colourer/graph_colour', '-l'+str(time_limit), '-f'+str(fractional_level), '-r'+str(seed), 'DUMMY_FILENAME'],
            stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    write_graph_to_proc(proc, len(vv), new_edges)
    result = proc.stdout.read().strip().split("\n")
    timeout = False
    for line in result:
        tokens = line.strip().split()
        if len(tokens)==4 and tokens[3] == "TIMEOUT":
            timeout = True
        if len(tokens)==3 and tokens[2] == "SATISFIABLE":
            bound = int(tokens[0]) / fractional_level
    proc.wait()

    return 0 if timeout else bound

def run_vc_solver(vv, new_edges, bound, timelimit):
    proc = subprocess.Popen(
            ['vc_solver/solve_vc', '-a0', '-c3', '-m-1', '-l'+str(timelimit), '-i'+str(bound), 'DUMMY_FILENAME'],
            stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    write_graph_to_proc(proc, len(vv), new_edges)
    result = proc.stdout.read().strip().split("\n")
    timeout = False
    for line in result:
        tokens = line.strip().split()
        if tokens and tokens[0] == "TIMEOUT":
            timeout = True
        if tokens and tokens[0] == "VertexCover":
            vc = [int(token) for token in tokens[1:]]
    proc.wait()

    return None if timeout else [vv[v] for v in vc]

def solve_subgraph(n, adj_lists, vv):
    vv = sorted(vv)
    new_adj_lists = [[] for v in vv]
    old_to_new_vtx = [-1] * n
    for i, v in enumerate(vv):
        old_to_new_vtx[v] = i
    new_edges = []
    for new_v, old_v in enumerate(vv):
        for old_w in adj_lists[old_v]:
            if old_w > old_v:   # so that we do each edge only once
                new_w = old_to_new_vtx[old_w]
                if new_w != -1:
                    new_edges.append((new_v, new_w))

    if len(vv) > 1000:
        return run_vc_solver(vv, new_edges, 0, 999999)
    else:
        result = run_vc_solver(vv, new_edges, 0, 5)
        if result is not None:
            return result
        for seed in range(40):
            bound = get_colouring_bound(vv, new_edges, 2, seed, 1 if seed<35 else 3)
            if bound == 0:
                print "c colouring bound timed out with seed", seed
            else:
                break
        return run_vc_solver(vv, new_edges, bound, 999999)

def vv_are_clique(vv, adj_lists):
    for i, v in enumerate(vv):
        for j in range(i+1, len(vv)):
            w = vv[j]
            if w not in adj_lists[v]:
                return False
    return True

def isolated_vertex_removal(n, adj_lists, in_cover, deleted):
    made_a_change = False
    for v in range(n):
        if not deleted[v]:
            neighbours = [w for w in adj_lists[v] if not deleted[w]]
            if vv_are_clique(neighbours, adj_lists):
                deleted[v] = True
                made_a_change = True
                for w in neighbours:
                    deleted[w] = True
                    in_cover[w] = True
                    for u in adj_lists[w]:
                        adj_lists[u].remove(w)
                    adj_lists[w] = []

    return made_a_change

def vertex_folding(n, adj_lists, in_cover, deleted, deg_2_reductions):
    made_a_change = False

    for i in range(n):
        for v in adj_lists[i]:
            if i not in adj_lists[v]:
                print "c ???", v, i

    for v in range(n):
        if len(adj_lists[v]) == 2:
            w = adj_lists[v][0]
            x = adj_lists[v][1]

            # for this reduction, w and x must not be adjacent
            if x not in adj_lists[w]:
                adj_lists[w].remove(v)
                adj_lists[x].remove(v)
                for u in adj_lists[x]:
                    adj_lists[u].remove(x)
                    if u not in adj_lists[w]:
                        adj_lists[w].append(u)
                        adj_lists[u].append(w)
                adj_lists[v] = []
                adj_lists[x] = []
                deleted[v] = True
                deleted[x] = True
                deg_2_reductions.append((v, w, x))
                made_a_change = True
    return made_a_change


def make_list_of_components(adj_lists):
    vertex_used = [len(adj_lists[i])==0 for i in range(n)]
    component_sizes = []
    components = []
    for i in range(n):
        if not vertex_used[i]:
            component = [i]
            to_explore = [i]
            vertex_used[i] = True
            while to_explore:
                v = to_explore[-1]
                del to_explore[-1]
                for w in adj_lists[v]:
                    if not vertex_used[w]:
                        component.append(w)
                        to_explore.append(w)
                        vertex_used[w] = True
#            print "adjlistlen", sorted([len(adj_lists[v]) for v in component])

            components.append(component)
            component_sizes.append(len(component))
    print "c component sizes", component_sizes
    print "c sum of component sizes", sum(component_sizes)
    return components


def presolve(n, edge_list):
    deleted = [False] * n   # has vertex i been deleted?
    in_cover = [False] * n   # has vertex i been added to the vertex cover?
    for v, w in edge_list:
        if v == w:
            in_cover[v] = True

    edge_list = [(v, w) for (v, w) in edge_list if not in_cover[v] and not in_cover[w]]

    adj_lists = [[] for i in range(n)]
    for v, w in edge_list:
        adj_lists[v].append(w)
        adj_lists[w].append(v)

    # deg_2_reductions is a list of tuples (v, w, x), where v is the vertex of deg 2,
    # w is the kept neighbour and x is the removed neighbour
    deg_2_reductions = []

    while True:
        a = isolated_vertex_removal(n, adj_lists, in_cover, deleted)
        b = vertex_folding(n, adj_lists, in_cover, deleted, deg_2_reductions)
        if not a and not b:
            break

    components = make_list_of_components(adj_lists)

    for component in components:
        vc = solve_subgraph(n, adj_lists, component)
        if vc is None:
            print "c TIMEOUT"
            exit()
        for v in vc:
            in_cover[v] = True

    for v, w, x in reversed(deg_2_reductions):
        if in_cover[w]:
            in_cover[x] = True
        else:
            in_cover[v] = True

    print "c VC SIZE", sum(in_cover)
    print "s vc {} {}".format(n, sum(in_cover))
    for v in range(n):
        if in_cover[v]:
            print v + 1
#    print " ".join(str(i+1) for i in range(n) if in_cover[i])

if __name__ == "__main__":
    with open(sys.argv[1], "r") as f:
        n, edge_list = read_instance([line for line in f.readlines()])

    edge_set = set(edge_list)
    edge_list = list(edge_set)
    edge_list.sort()

    presolve(n, edge_list)
