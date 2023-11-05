using Random
using BenchmarkTools

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
            vertice1 = parse(Int, vertices[2])  # First vertice
            vertice2 = parse(Int, vertices[3])  # Second vertice
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

function evaluate_single_instance(G::Graph_color) 
    """Evaluates the number of conflicts in the coloring of a single instance"""
    conflicts = 0
    for i in 2:G.nbr_vertices
            conflicts += sum(G.adj[k,i]*(G.color[i]==G.color[k]) for k=1:(i-1))
    end
    return conflicts
end

function permute_vertices(G::Graph_color) 
    """Changes randomly the vertixes order considered for the coloring"""
    adj_matrix = G.adj
    colors = G.color
    perm = randperm(size(adj_matrix, 1)) #permutation of the order of vertices

    adj_matrix_permuted = adj_matrix[perm, perm]
    colors_permuted = colors[perm]

    G_permuted = Graph_color(;nbr_vertices = G.nbr_vertices, adj = adj_matrix_permuted, color = colors_permuted)
    return G_permuted
end

function random_assignment(G::Graph_color, nbr_colors::Int) 
    """Assigns to each vertixes a random color"""
    for i=1:G.nbr_vertices
        G.color[i] = rand(1:nbr_colors)
    end
end

function greedy_random_assignment(G::Graph_color, nbr_colors::Int)
    """
    Assigns to each vertix a random color that is not its neighboors' if possible, and a random color otherwise
    """
    for i=1:G.nbr_vertices
        adj_colors = [G.color[v] for v=1:(i-1) if G.adj[v,i]==1]
        no_conflicts_color = [k for k=1:nbr_colors if !(k in(adj_colors))]
        if length(no_conflicts_color) > 0
            G.color[i] = no_conflicts_color[rand(1:length(no_conflicts_color))]
        else
            G.color[i] = rand(1:nbr_colors)
        end
    end
end

function tabu_search(G::Graph_color, nbr_colors::Int, nbr_max_iter::Int, tabu_memory_iter::Int)
    best_color = copy(G.color)
    best_value = evaluate_single_instance(G)
    iter = 0
    tabu_list = zeros(Int, G.nbr_vertices, nbr_colors)

    current_vertice = rand(1:G.nbr_vertices)
    local_best_color = rand(1:nbr_colors)
    local_best_value = 100000
    while iter <= nbr_max_iter
        local_best_value = 100000
        local_best_color = G.color[current_vertice]
        tabu_list[current_vertice, G.color[current_vertice]] = iter
        for c=1:nbr_colors
            if iter > tabu_list[current_vertice, c]
                G.color[current_vertice] = c
                neighboor_value = evaluate_single_instance(G)
                if neighboor_value < local_best_value
                    local_best_color = c
                    local_best_value = neighboor_value
                end
            end
        end
        G.color[current_vertice] = local_best_color
        tabu_list[current_vertice, local_best_color] += tabu_memory_iter
        if local_best_value < best_value
            for v in 1:G.nbr_vertices
                best_color[v] = G.color[v]
            end
            best_value = evaluate_single_instance(G)
        end
        iter += 1
        current_vertice = rand(1:G.nbr_vertices)
    end
    for v in 1:G.nbr_vertices
        G.color[v] = best_color[v]
    end
end


function evaluate_file(file_txt::String, nbr_iters::Int, nbr_colors::Int, tabu::Bool, permute::Bool)
    """Evaluates for a given file the number of conflicts nbr_iters times
    Returns the min number of conflicts and corresponding execution time"""
    conflicts = zeros(Int,nbr_iters)
    execution_times = zeros(Float64,nbr_iters)
    G_0 = read_file(file_txt, permute)
    for i=1:nbr_iters
        start_execution_time = time()
        if tabu
            # println("Executing Tabu Search ...") 
            tabu_search(G_0, nbr_colors, 50000, 500) #Graph_color, nbr_colors, nbr_max_iter, tabu_memory_iter

        else
            # println("Executing Greedy Heuristic ...")
            greedy_random_assignment(G_0, nbr_colors)
        end
        nbr_conflicts = evaluate_single_instance(G_0)
        conflicts[i] = nbr_conflicts
        end_execution_time = time()
        execution_times[i] = end_execution_time - start_execution_time
    end
    index_of_min = argmin(conflicts)
    println(minimum(conflicts), " : min nb of conflicts")
    println(round(mean(conflicts),digits = 4), " : mean nb of conflicts")
    println(maximum(conflicts), " : max nb of conflicts")
    # println(execution_times)
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
        nb_solutions_evaluated_per_sec = round(nbr_iters/(time_end - time_start), digits=2)
        # mean_execution_time = round(execution_time/nbr_iters,digits = 4)
        println("Total execution time : $execution_time secondes")
        # println("Mean execution time : $mean_execution_time secondes")
        if tabu == false
            println("Number of solutions evaluated per second : $nb_solutions_evaluated_per_sec")
        end
    end
end


#--------------------------------------------------------------------------------------

# Parameters
list_paths = ["Fichiers/dsjc125.1.col.txt","Fichiers/dsjc125.9.col.txt","Fichiers/dsjc250.1.col.txt","Fichiers/dsjc250.5.col.txt","Fichiers/dsjc250.9.col.txt","Fichiers/dsjc1000.5.col.txt","Fichiers/dsjc1000.5.col.txt","Fichiers/dsjc1000.5.col.txt","Fichiers/flat300_26_0.col.txt","Fichiers/le450_15c.col.txt"]
list_nb_colors = [5, 44, 8, 28, 72, 86, 85, 84, 26, 15]
nbr_iters = 10 #nbr of iterations for each file (we compute the nb of conflicts nbr_iters time and compute the min nbr of conflicts)
tabu = false
permute = false

# Evaluation
evaluate_all_files()








# #Tests
# println(typeof(G))
# println(typeof(G.nbr_vertices))
# println(G.nbr_vertices)
# println(typeof(G.adj))
# println(size(G.adj))
# println(typeof(G.color))
# println(size(G.color))
# println("number of edges : ", sum(G.adj[i,j] for i=1:G.nbr_vertices, j=1:G.nbr_vertices)/2)

# println(evaluate_single_instance(G))
# random_assignment(G, 5)
# println(evaluate_single_instance(G))
# greedy_random_assignment(G, 5)
# println(evaluate_single_instance(G))
