include("config.jl")
using .Config  # Importez le module Config du fichier config.jl en utilisant le chemin absolu

include("graph_coloring.jl")

# Exécution du code
evaluate_all_files()
