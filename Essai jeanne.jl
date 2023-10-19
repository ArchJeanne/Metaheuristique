function read_file(path::String)

    # Dictionnary of vertices
    # Key = name of vertice ; Values = vertexes adjacents to the vertex
    Dict_vertices = Dict() 
    vertices_int = 0

    # Open the file
    file = open(path, "r")
    # Read file, fill vertices dictionnary
    for line in eachline(file)

        # get number of vertices
        if line[1]=='p' 
            vertices_int = parse(Int, match(r"\d+", line).match) #number of vertices
        end

        # # fill dictionnary keys
        # for i in 1:vertices_int
        #     Dict_vertices[i]=[] 
        # end

        # fill dictionnary values
        if line[1]=='e'
            vertices = split(line)
            vertice1 = parse(Int, vertices[2])  # First vertice
            vertice2 = parse(Int, vertices[3])  # Second vertice
            # push!(Dict_vertices[vertice1], vertice2)    
            Dict_vertices[vertice1] = vcat(Dict_vertices[vertice1],vertice2)
            Dict_vertices[vertice2]= vcat(Dict_vertices[vertice2],vertice1)    
   
            # push!(Dict_vertices[vertice2], vertice1)  
            # println(vertice1," ", vertice2)
        end
    end
    print("ok")
    print(Dict_vertices)
   
    close(file)
end

fichier_texte = "Fichiers/dsjc125.1.col.txt"

read_file(fichier_texte)