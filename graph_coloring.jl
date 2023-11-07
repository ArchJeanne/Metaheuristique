### Import packages
using Random
using BenchmarkTools

### Load and transform graph structure if necessary
Base.@kwdef struct Graph_color
    nbr_vertices::Int
    adj::Matrix{Int}
    color::Vector{Int}
end

function read_file(path::String, permute::Bool) 
    """Returns the graph corresponding to the instance"""
    file = open(path, "r")
    # Determine the number of vertices
    n = 0 
    for line in eachline(file)
        if line[1]=='p' 
            n = parse(Int, match(r"\d+", line).match) #number of vertices
        end
    end
    close(file)
    # Initialisation of the adj_matrix and color vector
    adj_matrix = zeros(Int, n, n)
    color_vector = ones(Int, n)
    #Update of the adj matrix
    file = open(path, "r")
    for line in eachline(file)  
        if line[1]=='e'
            vertices = split(line)
            vertice1 = parse(Int, vertices[2])  # First vertex
            vertice2 = parse(Int, vertices[3])  # Second vertex
            adj_matrix[vertice1, vertice2] = 1
            adj_matrix[vertice2, vertice1] = 1
        end
    end
    close(file)
    G = Graph_color(;nbr_vertices = n, adj = adj_matrix, color = color_vector)
    if permute
        G = permute_vertices(G)
    end
    return(G)
end

function permute_vertices(G::Graph_color) 
    """Changes randomly the vertices order considered for the coloring"""
    adj_matrix = G.adj
    colors = G.color
    perm = randperm(size(adj_matrix, 1)) #permutation of the order of vertices

    adj_matrix_permuted = adj_matrix[perm, perm]
    colors_permuted = colors[perm]

    G_permuted = Graph_color(;nbr_vertices = G.nbr_vertices, adj = adj_matrix_permuted, color = colors_permuted)
    return G_permuted
end

### General graph functions

function evaluate_single_instance(G::Graph_color) 
    """Evaluates the number of conflicts in the coloring of a single instance"""
    conflicts = 0
    for i in 2:G.nbr_vertices
            conflicts += sum(G.adj[k,i]*(G.color[i]==G.color[k]) for k=1:(i-1))
    end
    return conflicts
end

function neighbors(G::Graph_color, vertex::Int)
    """Returns neighbors of vertex"""
    neighbor_indices = []
    for v in 1:G.nbr_vertices
        if G.adj[vertex, v] == 1
            push!(neighbor_indices, v)
        end
    end
    return neighbor_indices
end


### Random assignements heuristics

function random_assignment(G::Graph_color, nbr_colors::Int) 
    """Assigns to each vertices a random color"""
    for i=1:G.nbr_vertices
        G.color[i] = rand(1:nbr_colors)
    end
end

function greedy_random_assignment(G::Graph_color, nbr_colors::Int)
    """
    Assigns to each vertex a random color that is not its neighboors' if possible, and a random color otherwise
    Returns corresponding graph
    """
    for i=1:G.nbr_vertices
        adj_colors = [G.color[v] for v=1:(i-1) if G.adj[v,i]==1] #colors of neighboors
        no_conflicts_color = [k for k=1:nbr_colors if !(k in(adj_colors))] #colors of neighboors that are not in conflict with current color
        if length(no_conflicts_color) > 0
            G.color[i] = no_conflicts_color[rand(1:length(no_conflicts_color))] #pick one color out of the colors that don't create a conflict
        else
            G.color[i] = rand(1:nbr_colors) #otherwise pick a random color
        end
    end
    return G
end

### Tabu search : 1-vertix neighborhood structure

function diff_evaluate_one_vertice(G::Graph_color, v::Int, old_color::Int, new_color::Int)
    diff_conflicts = 0
    for i=1:G.nbr_vertices
        if G.adj[i,v] == 1
            if old_color == G.color[i]
                diff_conflicts = diff_conflicts - 1
            end
            if new_color == G.color[i]
                diff_conflicts = diff_conflicts + 1
            end
        end
    end
    return diff_conflicts
end


function tabu_search_1_vertice(G::Graph_color, nbr_colors::Int, nbr_max_iter::Int, tabu_memory_iter::Int)
    best_color = copy(G.color)
    best_value = evaluate_single_instance(G)
    current_value = best_value
    iter = 0
    tabu_list = zeros(Int, G.nbr_vertices, nbr_colors)
    current_vertex = rand(1:G.nbr_vertices)
    while iter <= nbr_max_iter
        local_best_value = 10000000
        local_best_color = G.color[current_vertex]
        tabu_list[current_vertex, G.color[current_vertex]] += round(tabu_memory_iter * (iter//nbr_max_iter))
        for c=1:nbr_colors
            if iter > tabu_list[current_vertex, c]
                neighboor_value = current_value + diff_evaluate_one_vertice(G, current_vertex, G.color[current_vertex], c)
                if neighboor_value < local_best_value
                    local_best_color = c
                    local_best_value = neighboor_value
                end
            end
        end
        if local_best_value < 10000000
            current_value = local_best_value
        end
        G.color[current_vertex] = local_best_color
        if local_best_value < best_value
            for v in 1:G.nbr_vertices
                best_color[v] = G.color[v]
            end
            best_value = local_best_value
        end
        iter += 1
        current_vertex = rand(1:G.nbr_vertices)
    end
    for v in 1:G.nbr_vertices
        G.color[v] = best_color[v]
    end
end

### Tabu search : 2-vertix-swap neighborhood structure

function diff_evaluate_two_vertices(G::Graph_color, v1::Int, v2::Int, color1::Int, color2::Int)
    diff_conflicts = 0
    old_color_1 = color1 # old color of v1
    new_color_1 = color2 # new color of v1
    old_color_2 = color2 # old color of v2
    new_color_2 = color1 # new color of v2
    
    # conflicts for v1
    for i in neighbors(G,v1) 
        if i!=v2    
            if old_color_1 == G.color[i]
                diff_conflicts = diff_conflicts - 1
            end
            if new_color_1 == G.color[i]
                diff_conflicts = diff_conflicts + 1
            end
        end
    end
    # conflicts for v2
    for i in neighbors(G,v2)  
        if i!=v1      
            if old_color_2 == G.color[i]
                diff_conflicts = diff_conflicts - 1
            end
            if new_color_2 == G.color[i]
                diff_conflicts = diff_conflicts + 1
            end
        end
    end
    return diff_conflicts
end


function swap_colors_in_neighborhood(G::Graph_color, tabu_list,  iter::Int, tabu_memory_iter, nbr_max_iter, current_value)
    """Swaps a random vertice with all its neighbors and evaluates if the current_value is improved"""
    # Two random vertices
    v1 = rand(1:G.nbr_vertices)
    neighbors_v1 = neighbors(G, v1)
    color_v1 = G.color[v1] 
    best_color_v1 = color_v1
    best_v2 = false # by default, if no color swap improves the objective, we don't swap vertices
    best_color_v2 = false # by default, if no color swap improves the objective, we don't swap vertices

    for v2 in neighbors_v1
        color_v2 = G.color[v2]

        if iter > tabu_list[v1,v2] #couple v1 v2 allowed   
            # Change node colors
            G.color[v1] = color_v2
            G.color[v2] = color_v1
            # Value after exchange
            new_value = current_value + diff_evaluate_two_vertices(G,v1,v2,color_v1,color_v2)
            if new_value < current_value
                current_value = new_value
                best_color_v1 = color_v2
                best_color_v2 = color_v1
                best_v2 = v2
            end
            # get back to previous graph
            G.color[v1] = color_v1
            G.color[v2] = color_v2
        end    
    end
    if best_v2!=false    # we found a vertice best_v2 that strictly improves the obj function
        G.color[v1] = best_color_v1
        G.color[best_v2] = best_color_v2
        tabu_list[v1, best_v2] += round(tabu_memory_iter * (iter//nbr_max_iter))
    end
    return G, current_value, tabu_list
end


function tabu_search_2_vertices(G::Graph_color, nbr_colors::Int, nbr_max_iter::Int, tabu_memory_iter::Int)
    # Initialisation
    G = greedy_random_assignment(G,nbr_colors)
    iter = 0
    tabu_list = zeros(Int, G.nbr_vertices, G.nbr_vertices)
    best_value = evaluate_single_instance(G)
    best_G = G
    current_value = evaluate_single_instance(G)

    while iter <= nbr_max_iter
        G, current_value, tabu_list = swap_colors_in_neighborhood(G, tabu_list, iter, tabu_memory_iter, nbr_max_iter, current_value)
        
        if current_value < best_value
            best_G = G
            best_value = current_value
        end
        iter+=1
    end
    return best_G, best_value
end

### Evaluation functions

function evaluate_file(file_txt::String, nbr_iters::Int, nbr_colors::Int, tabu::Bool, permute::Bool)
    """Evaluates for a given file the number of conflicts nbr_iters times
    Returns the min number of conflicts and corresponding execution time"""
    conflicts = zeros(Int,nbr_iters) #list of conflicts for each iteration
    execution_times = zeros(Float64,nbr_iters) #list of execution times for each iteration
    G_0 = read_file(file_txt, permute)
    for i=1:nbr_iters
        start_execution_time = time()
        if tabu
            if neighborhood_structure=="1_vertice"
                tabu_search_1_vertice(G_0, nbr_colors, nbr_max_iter_tabu, iter_tabu_memory) #Graph_color, nbr_colors, nbr_max_iter, tabu_memory_iter
            end
            if neighborhood_structure=="2_vertice"
                tabu_search_2_vertices(G_0, nbr_colors, nbr_max_iter_tabu, iter_tabu_memory) #Graph_color, nbr_colors, nbr_max_iter, tabu_memory_iter
            end
        else
            greedy_random_assignment(G_0, nbr_colors)
        end
        nbr_conflicts = evaluate_single_instance(G_0)
        conflicts[i] = nbr_conflicts
        end_execution_time = time()
        execution_times[i] = end_execution_time - start_execution_time
    end
    index_of_min = argmin(conflicts) #index of min number of conflicts in list of conflicts
    println(minimum(conflicts), " : min nb of conflicts")
    println(round(mean(conflicts),digits = 4), " : mean nb of conflicts")
    println(maximum(conflicts), " : max nb of conflicts")
    println("Best solution execution time : ", round(execution_times[index_of_min], digits = 4), " secondes")
end


function evaluate_all_files()
    """Evaluates all files, each file is evaluated nbr_iters times"""
    for i in eachindex(list_paths)
        time_start = time()
        path = list_paths[i]
        nbr_colors = list_nb_colors[i]

        println("------------------------------------")
        println(path, " - ", nbr_colors, " colors",if tabu " - Tabu search" else " - Greedy heuristic" end, if permute " - Permuted" else " - Not permuted" end)
        evaluate_file(path,nbr_iters,nbr_colors,tabu,permute)

        time_end = time()
        execution_time = round(time_end - time_start,digits = 4)
        println("Total execution time : $execution_time secondes")
        if tabu == false
            nb_solutions_evaluated_per_sec = round(nbr_iters/(time_end - time_start), digits=2)
            println("Number of solutions evaluated per second : $nb_solutions_evaluated_per_sec")
        end
    end
end
