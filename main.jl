Base.@kwdef struct Graph_color
    nbr_vertices::Int
    adj::Matrix{Int}
    color::Vector{Int}
end

function read_file(path::String, nbr_color::Int)
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
    color_vector = zeros(Int, nbr_color)
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

#--------------------------------------------------------------------------------------

fichier_texte = "Fichiers/dsjc125.1.col.txt"

println("------------------")
G = read_file(fichier_texte, 21)
println(typeof(G))
println(typeof(G.nbr_vertices))
println(G.nbr_vertices)
println(typeof(G.adj))
println(size(G.adj))
println(typeof(G.color))
println(size(G.color))

println("number of edges : ", sum(G.adj[i,j] for i=1:G.nbr_vertices, j=1:G.nbr_vertices)/2)