# Projet_ALG_Laval_Lesbats_M2_BI

## Description
Le projet consiste à implémenter et analyser une structure d’indexation de génomes appelée graphe de de Bruijn coloré (cDBG).  
L’objectif est de stocker les k-mers d’une collection de génomes tout en conservant l’information de leur appartenance aux différents génomes (couleurs), puis de calculer des mesures de similarité alignment-free entre des séquences requêtes et cette collection.

Dans le projet nous aurons propose deux implémentations :
- une version naïve
- une version avancée, exploitant la redondance des couleurs au sein des unitigs
- 

## Auteurs
- Adrien Laval
- Cédric Lesbats


## Organisation du dépôt

src/
├── naive/
│   ├── build.py
│   ├── query.py
│   └── dbg_indexer.py
└── advanced/
    ├── build.py
    ├── query.py
    └── dbg_indexer.py



## Principe général

### Construction de l’index
À partir d’un fichier listant les chemins vers plusieurs fichiers FASTA (un génome par ligne), le programme construit un graphe de de Bruijn coloré :
- les k-mers sont extraits des génomes,
- chaque k-mer est associé à une ou plusieurs couleurs, correspondant aux génomes dans lesquels il apparaît,
- l’ordre des fichiers dans la liste définit l’ordre des couleurs.

Le graphe est ensuite sérialisé à l’aide de la librairie pickle.



### Requête
À partir d’un graphe sérialisé et d’un fichier FASTA de séquences requêtes, le programme calcule, pour chaque requête, un ratio de k-mers partagés avec chacun des génomes de la collection.



## Versions implémentées

### Version naïve
- Structure : dictionnaire {k-mer → liste de couleurs}
- Avantages :
  - implémentation simple
  - requêtes rapides
- Inconvénients :
  - redondance importante de l’information de couleur
  - structure plus volumineuse



### Version avancée
- Structure : dictionnaire {unitig → liste de couleurs}
- Hypothèse exploitée :
  - tous les k-mers d’un même unitig partagent les mêmes couleurs
- Avantages :
  - structure plus compacte en mémoire et sur disque
- Inconvénients :
  - requêtes plus lentes (recherche d’un k-mer par parcours des unitigs)


## Utilisation

### Construction de l’index

Version naïve :

Version avancée :


### Requête

Version naïve :

Version avancée :


## Entrées et sorties

### Fichier d’entrée des génomes
Le fichier passé avec l’option `-i` doit contenir un chemin vers un fichier FASTA par ligne.  
L’ordre des lignes définit l’ordre des couleurs dans les résultats.

### Sortie des requêtes
Pour chaque séquence requête, le programme produit une ligne contenant :
- l’identifiant de la séquence
- une liste de similarités (ratios de k-mers partagés), une valeur par génome


## Environnement et dépendances
- Python 3
- Bibliothèques :
  - bibliothèque standard Python
  - BioPython (module `SeqIO`)

