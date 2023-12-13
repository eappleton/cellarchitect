import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import itertools
import meshio
import networkx as nx
import random, pdb
import operator
import pickle


def unique_rows_counts(arr):
    """
    Count efficiently the number of occurence of a np.array in a 2D np.array
    :param arr: 2d numpy array
    :return:
    :arr[unq_idx]: unique arrays in arr
    :counts: counts of the unique arrays
    """
    # from stackoverflow

    # Calculate linear indices using rows from arr
    lidx = np.ravel_multi_index(arr.T,arr.max(0)+1 )

    # Get the unique indices and their counts
    _, unq_idx, counts = np.unique(lidx, return_index = True, return_counts=True)

    return arr[unq_idx], counts

def distance_3d(a,b):
    """
    Calculate the euclidean distance between two points in 3d
    :param a,b: points
    :return: euclidean distance
    """
    return np.linalg.norm(a-b)

def d3_to_d2(arr):
    """
    Convert a 3d np array to a 2d array by squeezing the 2nd dim
    :param arr: array to convert
    :return: converted array
    """
    return arr.reshape(-1, arr.shape[-1])

def unique_rows(arr):
    """
    Find the unique rows in an array
    :param arr: array
    :return: unique rows
    """
    # from stackoverflow
    
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def plot_formula(formula, x_range, x_label=None, y_label=None): 
    """
    Plot a forumula in str format
    :param formula: formula to plot, e.g. '2*x+1'
    :param x_range: x values to be evaluated
    :params x_label, y_label: label to add on x and/or y axis
    """
    x = np.array(x_range)  
    y = eval(formula)
    plt.plot(x, y)
    
    ax = plt.gca()
    if x_label is not None:
        ax.set_xlabel(x_label)
    if y_label is not None:
        ax.set_ylabel(y_label)
        
    plt.show()
    
def plot_histo_distribution_edges_length(path_to_file, mesh_constraints):
    """
    Plot the distribtion of edges length from a .msh file
    :param path_to_file: path to the .msh
    :param mesh_constraints: mesh size constraints used to create the mesh (int)
    """    
    points, cells, _, _, _ = meshio.read(path_to_file)
    
    tetra = cells['tetra']
    tetra_all = np.empty([tetra.shape[0], tetra.shape[1], 3])
    
    for k in range(tetra.shape[0]):
        for i,val in enumerate(tetra[k]):
            tetra_all[k][i] = points[val]
            
    all_length = []
    for i in range(tetra.shape[0]):
        for a, b in itertools.combinations(tetra_all[i], 2):
            all_length.append(distance_3d(a,b))
      
    plt.axvline(x=mesh_constraints, color='r')
    plt.hist(all_length)
    plt.xlabel('Length of edges')
    plt.ylabel('Number of edges')

    plt.show()
    
def print_edges_and_vertices(graph):
    
    print('Number of nodes {} and number of edges {}'.format(graph.number_of_nodes(), graph.number_of_edges()))
    
    
    
####################################
#                                  #
#          Connectors algo         #
#                                  #
####################################

def principal_vertex(elements, repeats):
    """ Not sure we will need this one anymore """
       
    index = repeats.argmax()
    elem = elements[index]
    
    print("Max elem is " + str(elem) + " with " + str(repeats[index]) + " repetitions.")

    return elem

def convert_mesh_to_graph(path_to_file):
    """
    Convert a thetra based mesh (.msh) file into a graph representation where
    vertices are building blocks (bb) and edges represent neighbor bb
    (i.e bb sharing at least one common vertex)
    :param path_to_file: path to the .msh file
    :return: the graph
    """
    
    points, cells, _, _, _ = meshio.read(path_to_file)
    
    # Get back the tetra information
    tetra = cells['tetra']
    tetra_all = np.empty([tetra.shape[0], tetra.shape[1], 3])
    
    for k in range(tetra.shape[0]):
        for i,val in enumerate(tetra[k]):
            tetra_all[k][i] = points[val]
          
    # create the graph
    G=nx.Graph()
    
    # add all nodes to the graph
    for idx,coordinate in enumerate(tetra_all):
        G.add_node(idx, coord=coordinate)  
    # note: G has now a node for each tetra. And each node has an attribute coord
    # which corresponds to the 2d matrix of the vertices of the tetra.    
    
    # compare all tet which each other
    comp_list = itertools.combinations(enumerate(tetra_all), 2) 
    #comparison_list is a list with the index of the tet compared (as in tetra_all) and the vertices of the two
    #compared arrays in two different arrays.
    comp_lit_with_vertices = [ [[idx1, idx2], [arr1,arr2]] for ((idx1, arr1),(idx2,arr2)) in comp_list ]
    
    # add the edges to the graph by looking into tet that shared common vertex
    for comparison in comp_lit_with_vertices:
        C = np.array([x for x in set(tuple(x) for x in comparison[1][0]) & set(tuple(x) for x in comparison[1][1])])
        if C.size != 0:
            G.add_edge(comparison[0][0],comparison[0][1])
            
    print_edges_and_vertices(G)
            
    return G


def k_connected_from_none(graph, start, k):
      
    # We use MC sampling to find a very close to optimial solution
    # to the optimization probltem
        
    n_vertices_under_k = sum( x < k for x in graph.degree().values() )
    if n_vertices_under_k != 0: 
        print('{} nodes have less thank K connections in the original graph!'.format(n_vertices_under_k))
    
    
    n_connections = []
    best_graph = None
    min_connection = 100000
    
    # MC sampling to find a close-to-optimal minimum
    for i in range(200):
        
        visited, stack = set(), [start]
        G = nx.Graph()
        G.add_nodes_from(graph.nodes())
        
        while stack:
            vertex = stack.pop()
            if vertex not in visited:
                visited.add(vertex)
                stack.extend(set(graph[vertex]) - visited)
                vertex_to_connect = []
            
                diff = k-len(G[vertex])
                if diff <= k and diff > 0:
                    available_nodes = [x for x in graph[vertex] if len(G[x])<k]                  
                    L = len(available_nodes)
                    if L >= diff:
                        vertex_to_connect = random.sample(available_nodes,diff)
                    else:
                        vertex_to_connect = random.sample(available_nodes,L) # add connection from a node under k connected
                
                to_connect = []
                for element in vertex_to_connect:
                    to_connect.append((vertex,element))
                G.add_edges_from(to_connect)
                          
        n_edges = G.number_of_edges()
        n_connections.append(n_edges)
        
        if n_edges < min_connection and nx.is_connected(G):
            min_connection = n_edges
            best_graph = G
            
    # some vertices might still have less than k connections, so
    # we break the degree() constraint at k to make k connections
    # everywhere
    
    nodes_under_k = [x for x in best_graph if len(best_graph[x])<k]
    for node in nodes_under_k:
        last_connections = []
        diff = k-len(best_graph[node])
        available = list(set(graph[node]) - set(best_graph[node]))
        
        if diff < 0:
            diff = 0
        elif len(available) < diff: 
            diff = len(available)
        
        last_connections = random.sample(available,diff)
        to_connect_under = []
        for element in last_connections:
            to_connect_under.append((node,element))
        best_graph.add_edges_from(to_connect_under)    
    
         
    print_edges_and_vertices(best_graph)
   
    return best_graph,n_vertices_under_k
    
    
def k_connected_from_mst(mst_graph, original_graph, k):
    """
    Create a k_connected graph from a minimum spanning tree graph and the original graph
    note: the original graph is not modified, a new one is created
    :param mst_graph: minimum spanning tree graph
    :param original_graph: fully connected graph
    :return: k connected from mst graph
    """
        
    # TODO: this is a preliminary version; sufficient to present the results
    # Need to find if there is a pattern in the random choice of vertices
    # such that it minimizes the number of edges (~15% variation now with MC sampling)
    # No success to find an algo yet -> look in graph theory paper    
    
    
    n_vertices_under_k = sum( x < k for x in original_graph.degree().values() )
    if n_vertices_under_k != 0: 
        print('{} nodes have less thank K connections in the original graph!'.format(n_vertices_under_k))
    
    
    n_connections = []
    best_graph = None
    min_connection = 100000
    
    # MC sampling to find a close-to-optimal minimum
    for i in range(200):
        
        T = mst_graph.copy()
        
        for vertex in T:            
            diff = k-len(T[vertex])
            if diff>0:
                vertex_to_connects = [x for x in list(original_graph[vertex]) if x not in list(T[vertex])][:diff]                
                for vertex_to_connect in vertex_to_connects:
                    T.add_edge(vertex,vertex_to_connect)
        
        n_edges = T.number_of_edges()
        n_connections.append(n_edges)
        if n_edges < min_connection:
            min_connection = n_edges
            best_graph = T
            
    print_edges_and_vertices(best_graph)        
            
    return best_graph, n_vertices_under_k


def moving_part(mst_graph, original_graph, L):
    
    T = mst_graph.copy()
        
    # Find all nodes with less than or equal than 2 connections
    nodes = [n for n in T if len(mst_graph[n]) <= 2]
    
    # Create new graph with nodes found above
    G = nx.Graph()
    G.add_nodes_from(nodes)
    
    # Add the edges that connected the node from the mst graph
    comp_list = itertools.combinations(nodes, 2)
    for comp in comp_list:
        if mst_graph.has_edge(comp[0],comp[1]):
            G.add_edge(comp[0],comp[1])
    pdb.set_trace()
    # Find all connected subgraph from G
    graphs_connected = list(nx.connected_component_subgraphs(G))
    
    # find subgraph having more than L nodes
    graphs_to_modify = [g for g in graphs_connected if len(g)>L]
    
    # add connections to satify the L constraint:
    for graph in graphs_to_modify:
        # find a dead end:
        node_start = [node for node in graph if len(graph[node])==1][0]
        # find to which nodes a connection will be added
        nodes_to_add_connections = find_nodes_to_add_connections(graph,node_start,L)
        # add connections
        connect_node_to_new_node(mst_graph,original_graph,nodes_to_add_connections,T)
        
    print_edges_and_vertices(T)
        
    return T
    
def connect_node_to_new_node(mst_graph,original_graph,nodes,new_graph):
    """
    Algo to add connectors according to the moving_part() algo
    :param mst_graph: original mst graph
    :param original_graph: original graph with all connections
    :param nodes: list of nodes which need the additional connections
    :param new_graph: the moving_part graph
    """
        
    for node in nodes:
        while True:
            vertex_to_connect = random.choice(list(original_graph[node]))
            if  vertex_to_connect not in list(mst_graph[node]):
                break
        new_graph.add_edge(node,vertex_to_connect)
    
        
def find_nodes_to_add_connections(graph, start, L):
    """
    Algo to find the nodes which need a additional connection
    Implementation modify from Depth-First search
    :param graph: sub graph with nodes having maximum 2 edges to which we
    want to add connectors
    :param start: starting node; a dead end. ie a node with only one edge
    :param L: moving part parameters
    """
    
    counter = 0
    nodes_to_add_connection = []
    
    visited, stack = set(), [start]
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            stack.extend(list(set(graph[vertex]) - visited))
            counter +=1
        
        if counter == L:
            counter = 0
            nodes_to_add_connection.append(vertex)
            
    return nodes_to_add_connection

def check_L_moving_part(graph,L):
    
    nodes = [n for n in graph if len(graph[n]) <= 2]
    G = nx.Graph()
    G.add_nodes_from(nodes)
    
    # Find all connected subgraph from graph
    graphs_connected = list(nx.connected_component_subgraphs(G))
    
    # find subgraph having more than L nodes
    graphs_to_modify = [g for g in graphs_connected if len(g)>L]
    
    if graphs_to_modify:
        print('Moving part did NOT worked')
    else:
        print('Moving part worked!')
        
def bool_L_moving_part(graph,L):
    
    nodes = [n for n in graph if len(graph[n]) <= 2]
    G = nx.Graph()
    G.add_nodes_from(nodes)
    
    # Find all connected subgraph from graph
    graphs_connected = list(nx.connected_component_subgraphs(G))
    
    # find subgraph having more than L nodes
    graphs_to_modify = [g for g in graphs_connected if len(g)>L]
    
    if graphs_to_modify:
        return False
    else:
        return True
    
    
def check_k_connections(graph,k):
    """
    Check how many vertices have less than k connections in a graph
    :param graph:
    :param k:
    """
    
    n_connection_under_k = 0
    for i in range(len(graph)):
        if len(graph[i]) < k: n_connection_under_k+=1
    print('n connections under k: ', n_connection_under_k)
    
def n_connections_under_k(graph,k):
    """
    return how many vertices have less than k connections in a graph
    :param graph:
    :param k:
    """
    
    n_connection_under_k = 0
    for i in range(len(graph)):
        if len(graph[i]) < k: n_connection_under_k+=1
    
    return n_connection_under_k
    
    
def plot_analysis_algo(algo_names, n_edges, n_edges_percent_of_origin):
    """
    Create a bar plot to display how many edges each graph created by
    algorithms
    :param algo_names: a list with the algo names (string)
    :param n_edges: a list with the number of edges (int)
    :param n_edges_percent_of_origin: a list with the % the new number
    of edges represent from the orginal fully connected graph (int or float)
    note: the three list must be in the same order
    """
    
    pos = np.arange(len(algo_names))
    width = 1.0
    
    fig, ax = plt.subplots(figsize=(4, 5), dpi=100)    
    ax = plt.axes()
    
    ax.set_xticks(pos - (width / 2.5))
    ax.set_xticklabels(algo_names, rotation=45)
    ax.set_ylabel('N of connections')
        
    plt.bar(pos, n_edges, width * 0.5, color='c')
    
    rects = ax.patches
    labels = [ '%.2f' % elem for elem in n_edges_percent_of_origin ]
    
    for rect, label in zip(rects, labels):
        height = rect.get_height()
        #ax.text(rect.get_x() + rect.get_width()/2, height + 5, str(label)+' %', ha='center', va='bottom')
        ax.text(rect.get_x() + rect.get_width()/2, height, str(label)+' %', ha='center', va='bottom')
    
    
    plt.show()
            
        
def plot_absolute_and_relative_n_connectors(graph_pruned, algo_names, original_graph):
    
    n_edges = []
    n_edges_percent_of_origin = []
    n_edges_origin = original_graph.number_of_edges()
    
    for graph in graph_pruned:
        n_edges.append(graph.number_of_edges())
        n_edges_percent_of_origin.append( (graph.number_of_edges()/n_edges_origin)*100 )
        
    plot_analysis_algo(algo_names,n_edges,n_edges_percent_of_origin)

    

        
    
def write_output_file(graph, directory, filename):
    """
    Read a graph and write it into an external txt file. Header:
    tet_id_1, tet_id_2, connector_id, vertex_coordinate (belong to tet_id_1)
    :param graph: 
    :param filename: name of the output file (string), without extension. eg: 'test'
    """    
    # note: here connector id will represent a face id, ie three cells of a given face will have this id.
    # we start the id at 1. id 0 is kept for the leaf of the mesh structure; they all have the same
    # external connector.
    connectors_id = 1
    output_file = []
    already_in = []
    
    for node in graph.nodes():
        voisins = list(graph[node].keys())
        current_vertices = graph.node[node]['coord']
        
        # for every neighbor
        for voisin in voisins:
            increase = False
            for vertex in graph.node[voisin]['coord']:
                if any(np.equal(current_vertices,vertex).all(1)):
                    list_a_to_b = [node, voisin]
                    list_a_to_b.extend(list(vertex))
                    list_a_to_b.append(connectors_id)
                    list_b_to_a = [voisin, node]
                    list_b_to_a.extend(list(vertex))
                    list_b_to_a.append(connectors_id)
                                    
                    if list_a_to_b[:3] not in already_in:
                    #if list_a_to_b[:3] not in already_in and list_b_to_a[:3] not in already_in:
                        
                        already_in.append(list_a_to_b[:3])
                        already_in.append(list_b_to_a[:3])
                        
                        # check if the vertex is a leaf
                        if graph.degree(voisin)==1 or graph.degree(node)==1:
                            # in this case, we use the unique connector
                            # because cells are part of a leaf block
                            output_file.append(list_a_to_b[:-1] + [0])
                            output_file.append(list_b_to_a[:-1] + [0])
                        else:
                            output_file.append(list_a_to_b)
                            output_file.append(list_b_to_a)
                            increase = True
                        
            if increase:
                connectors_id +=1
                    
    save_pickle(directory, filename, output_file)
    nx.write_gpickle(graph, directory + 'Graph_' + filename)
            
            
####################################
#                                  #
#               other              #
#                                  #
####################################            
            
            
            

def save_pickle(directory, name, list_):
    with open(directory + name + '.txt', 'wb') as fp:
        pickle.dump(list_, fp)

def load_pickle(path):
    with open(path, 'rb') as f:
        list_ = pickle.load(f)
    return list_            
            
                
            
####################################
#                                  #
#            dev program           #
#                                  #
####################################

def create_actions_lists(n_d, a_dict, dif):
    
    main_a_list = []
    rest_a_list = []

    main_a_list.extend([a_dict['divide'], a_dict['increase_count']]*n_d)
    
    if dif > 0:
        rest_a_list.extend(main_a_list)
        rest_a_list.extend([a_dict['divide'], a_dict['increase_count']])
    
        main_a_list.append(a_dict['increase_count'])
        
    # create blocks
    create_blocks = [a_dict['divide_and_stuck'], a_dict['increase_count']]*3
    main_a_list.extend(create_blocks)
    rest_a_list.extend(create_blocks)

    print('main_actions_list: ', main_a_list) 
    print('rest_actions_list: ',  rest_a_list)
    
    return main_a_list, rest_a_list
    
    
def express_pbp(tets, m_action_list, r_action_list, n_div, diff, power2):
    
    # create a list of lists; each list represents the programm of one tet, shared
    # until this point by all cells
    
    # this code is more elegant than the two for loops, but can't
    # get a way to make a copy of the list in it; it refers to the 
    # original list and messed up everything
    #tets_program = list(itertools.repeat(m_action_list, power2[n_div]))
    #tets_program.extend(list(itertools.repeat(r_action_list, diff)))
    
    tets_program = []
    for i in range(power2[n_div]):
        tets_program.append(m_action_list[:])
        
    for i in range(diff):
        tets_program.append(r_action_list[:])
    
    # now we add to the tet program the program to express a specific binders
    # format: last 4 entries are the binder each cell has to express
    for idx,val in enumerate(tets_program):
        val.extend(tets[idx])
    
    return tets_program
    
    
    
    
    
    
        
    