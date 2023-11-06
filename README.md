# Graph Coloring Solver

Ce programme résout le problème de la coloration de graphes à l'aide de l'algorithme de recherche tabou et d'une heuristique gloutonne. Il peut lire des fichiers d'instances de graphes et évaluer la qualité de la coloration pour plusieurs itérations.
Il est codé en Julia.

## Utilisation

1. Mdifiez les paramètres d'exécution en éditant le fichier `config.jl` pour définir les valeurs souhaitées, notamment pour les chemins d'accès aux instances à étudier.

2. Exécutez le fichier `main.jl` pour lancer le programme.

## Structure du Code

Le code est organisé en deux fichiers principaux :
- `graph_coloring.jl` : Ce fichier contient toutes les fonctions liées à la résolution du problème de coloration de graphes.
- `config.jl` : Ce fichier est utilisé pour définir les paramètres nécessaires à l'exécution du code.
- `main.jl` : Ce fichier est utilisé pour exécuter le code en important les fonctions du fichier `graph_coloring.jl`.

## Paramètres

Les paramètres d'exécution du programme sont définis dans le fichier `config.jl`. Vous pouvez modifier les valeurs des paramètres en éditant ce fichier. Voici les paramètres disponibles :

- `list_paths` : Une liste des chemins des fichiers d'instances de graphes à évaluer.
- `list_nb_colors` : Une liste du nombre de couleurs à utiliser pour chaque fichier d'instance.
- `nbr_iters` : Le nombre d'itérations à exécuter pour chaque fichier d'instance.
- `tabu` : Booléen pour spécifier l'utilisation de la recherche tabou.
- `permute` : Booléen pour spécifier la permutation des sommets du graphe.
- `nbr_max_iter_tabu` : Le nombre maximum d'itérations pour la recherche tabou.
- `iter_tabu_memory` : La mémoire de la liste tabou pour la recherche tabou.

Jeanne ARCHAMBAULT & Cyril BERTRAND