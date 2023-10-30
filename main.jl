Base.@kwdef struct Graph_color
    nbr_vertices::Int
    adj::Matrix{Int}
    color::Vector{Int}
end

function read_file(path::String)
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
    return(G)
end

function evaluate(G::Graph_color)
    conflicts = 0
    for i in 2:G.nbr_vertices
            conflicts += sum(G.adj[k,i]*(G.color[i]==G.color[k]) for k=1:(i-1))
    end
    return conflicts
end


function random_assignment(G::Graph_color, nbr_colors::Int)
    for i=1:G.nbr_vertices
        G.color[i] = rand(1:nbr_colors)
    end
end

function greedy_random_assignment(G::Graph_color, nbr_colors::Int)
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

#--------------------------------------------------------------------------------------

fichier_texte = "Fichiers/dsjc125.1.col.txt"

println("------------------")
G = read_file(fichier_texte)

#Tests
println(typeof(G))
println(typeof(G.nbr_vertices))
println(G.nbr_vertices)
println(typeof(G.adj))
println(size(G.adj))
println(typeof(G.color))
println(size(G.color))
println("number of edges : ", sum(G.adj[i,j] for i=1:G.nbr_vertices, j=1:G.nbr_vertices)/2)

println(evaluate(G))
random_assignment(G, 5)
println(evaluate(G))
greedy_random_assignment(G, 5)
println(evaluate(G))

