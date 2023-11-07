module Config
    export list_paths, list_nb_colors, nbr_iters, tabu, neighborhood_structure, permute, nbr_max_iter_tabu, iter_tabu_memory

    # Paramètres
    const list_paths = ["Fichiers/dsjc125.1.col.txt","Fichiers/dsjc125.9.col.txt","Fichiers/dsjc250.1.col.txt","Fichiers/dsjc250.5.col.txt","Fichiers/dsjc250.9.col.txt","Fichiers/dsjc1000.5.col.txt","Fichiers/dsjc1000.5.col.txt","Fichiers/dsjc1000.5.col.txt","Fichiers/flat300_26_0.col.txt","Fichiers/le450_15c.col.txt"]
    const list_nb_colors = [5, 44, 8, 28, 72, 86, 85, 84, 26, 15]
    const nbr_iters = 10
    const tabu = true
    const neighborhood_structure = "2_vertice" #["1_vertice","2_vertice"]
    const permute = false
    const nbr_max_iter_tabu = 500
    const iter_tabu_memory = 50
end