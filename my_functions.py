# import basic libraries
import numpy as np 
import pandas as pd
import gudhi as gd
import math
import matplotlib.pyplot as plt

# Extract list of genes in a particular pathway based on pathway ID
def get_pathway_genes(kegg_map, pathway_id):
    gene_list = []
    for i in range(len(kegg_map)):
        kg = kegg_map.iloc[i, 2]
        if kegg_map.iloc[i, 1] == pathway_id:
            gene_list.append(kg)
 
    return gene_list
    
# Extract a sub-data (for a particular pathway) from the experession data based on a pathway ID
def extract_pathway_data(data, kegg_map, pathway_id): 
    gene_list = []
    for i in range(len(kegg_map)):
        kg = kegg_map.iloc[i, 2]
        if kegg_map.iloc[i, 1] == pathway_id:
            gene_list.append(kg)
    
    result = data[data.index.isin(gene_list)]
    return result

# Get Pearson correlation matrix
def corr_matrix(data, corr_type="Similarity", by="Sample"):
    
    corr_data=[]
    
    if by == "Sample":
        if corr_type == "Similarity":
            corr_data = data.corr()
        elif corr_type == "Dissimilarity":
            corr_data = 1 - data.corr()
        elif corr_type == "Absolute Dissimilarity":
            corr_data = 1 - abs(data.corr())
        else: None
    elif by == "Gene":
        data = data.T
        if corr_type == "Similarity":
            corr_data = data.corr()
        elif corr_type == "Dissimilarity":
            corr_data = 1 - data.corr()
        elif corr_type == "Absolute Dissimilarity":
            corr_data = 1 - abs(data.corr())
        else: None
    else: None
    
    return corr_data

# Plot the correlation matrix
def plot_corr_matrix(correlation_matrix, by="Sample", size=(8,7)):
          
    if by == "Sample":
        plt.figure(figsize=size)
        sns.heatmap(correlation_matrix)
        plt.title(f"Correlation matrix")
        plt.xlabel("Samples")
        plt.ylabel("Samples")
        plt.show()
    elif by == "Gene":
        plt.figure(figsize=size)
        sns.heatmap(correlation_matrix)
        #plt.title(f"Correlation matrix for {corr_type}")
        plt.title(f"Correlation matrix")
        plt.xlabel("Genes")
        plt.ylabel("Genes")
        plt.show()
    else: None

# Get persistence diagram from a given data
def get_persistence_diagram(data, max_dim=2, edge_length=10):
    # Convert data to numpy array
    data = np.array(data)
    
    # Calculate persistence diagrams
    rips_complex = gd.RipsComplex(points=data, max_edge_length=edge_length)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=max_dim)
    persistence = simplex_tree.persistence()
    
    return persistence

# Plot the persistence diagram
def plot_PD(persistence, dim=-1):
    grouped_data = {}
    for tuple in persistence:
        dimension = tuple[0]
        if dimension not in grouped_data:
            grouped_data[dimension] = []
        grouped_data[dimension].append(tuple)
        
    colormap = None
 
    if dim == 0:
        persistence = grouped_data[0]
        gd.plot_persistence_diagram(persistence, colormap = colormap)  
    elif dim == 1:
        persistence = grouped_data[1]
        gd.plot_persistence_diagram(persistence, colormap = colormap)
    elif dim == -1:
        gd.plot_persistence_diagram(persistence, colormap = colormap)
    else:
        None
        
    plt.grid(False)  #

# Plot the persistence barcode
def plot_PB(persistence, dim=-1):
    grouped_data = {}
    for tuple in persistence:
        dimension = tuple[0]
        if dimension not in grouped_data:
            grouped_data[dimension] = []
        grouped_data[dimension].append(tuple)
        
    # colors = [custom_cmap(dim) for dim, _ in persistence]
    colormap = None
    
    if dim == 0:
        persistence = grouped_data[0]
        gd.plot_persistence_barcode(persistence, colormap = colormap)  
    elif dim == 1:
        persistence = grouped_data[1]
        gd.plot_persistence_barcode(persistence, colormap = colormap)
    elif dim == -1:
        gd.plot_persistence_barcode(persistence, colormap = colormap)
    else:
        None
        
    plt.grid(False) #

# Get the Betti numbers of a simplicial complex from a given data
def get_betti_numbers(data, max_dim=2, edge_length=10, filtration_values = np.arange(0, 1.1, 0.1)):
    # Convert data to numpy array
    data = np.array(data)
    
    # Calculate persistence diagrams
    rips_complex = gd.RipsComplex(points=data, max_edge_length=edge_length)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=max_dim)
    persistence = simplex_tree.persistence()
    
    betti_numbers_list=[]
    
    # Compute Betti numbers at each filtration value
    for filtration in filtration_values:
        simplex_tree_copy = gd.SimplexTree()
        simplex_tree_copy.set_dimension(simplex_tree.dimension())
    
        for s, f in simplex_tree.get_filtration():
            if f <= filtration:
                simplex_tree_copy.insert(s, filtration)
    
        persistence = simplex_tree_copy.persistence()
        betti_numbers = simplex_tree_copy.betti_numbers()
        betti_numbers_list.append(betti_numbers)
        #print(f"Betti numbers at filtration {filtration}: {betti_numbers}")
    betti_numbers_list = pd.DataFrame(betti_numbers_list)
    
    return betti_numbers_list

# Get sum of death times of a PD
def sum_death(persistence, pers_dim=0):
    sums = [[] for _ in range(5)]
    dims = [[] for _ in range(5)]
    
    for i in range(len(persistence)):
        dim = persistence[i][0]
        death = persistence[i][1][1]
        
        if dim < 5:
            dims[dim].append(death)
            sums[dim].append(death)
            
    sums = [np.array(sum_list) for sum_list in sums]
    
    if pers_dim in [0, 1, 2, 3, 4]:
        death_sum = np.sum(sums[pers_dim][~np.isinf(sums[pers_dim])])
        return death_sum
    else:
        return None
    
# Get sum of birth times of a PD
def sum_birth(persistence, pers_dim=0):
    sums = [[] for _ in range(5)]
    dims = [[] for _ in range(5)]
    
    for i in range(len(persistence)):
        dim = persistence[i][0]
        birth = persistence[i][1][0]
        
        if dim < 5:
            dims[dim].append(birth)
            sums[dim].append(birth)
            
    sums = [np.array(sum_list) for sum_list in sums]
    
    if pers_dim in [0, 1, 2, 3, 4]:
        birth_sum = np.sum(sums[pers_dim][~np.isinf(sums[pers_dim])])
        return birth_sum
    else:
        return None

# Get persistence (lifespan) of a certain dim features of a PD
def length_persistence(persistence, pers_dim=0):
    lengths = [[] for _ in range(5)]
    
    for i in range(len(persistence)):
        dim = persistence[i][0]
        length = persistence[i][1][1] - persistence[i][1][0]
        
        if dim < 5:
            lengths[dim].append(length)
            
    if pers_dim in [0, 1, 2, 3, 4]:
        return lengths[pers_dim]
    else:
        return None
    
# Get sum of persistence (lifespan) of a certain dim features of a PD
def sum_persistence(persistence, pers_dim=0):
    sum_lengths = [0, 0, 0, 0, 0]
    
    for i in range(len(persistence)):
        dim = persistence[i][0]
        length = persistence[i][1][1] - persistence[i][1][0]
        
        if dim < 5 and not math.isinf(length):
            sum_lengths[dim] += length
            
    if pers_dim in [0, 1, 2, 3, 4]:
        return sum_lengths[pers_dim]
    else:
        return None
    
# Get maximum persistence (lifespan) amongst features in a certain dim
def max_persistence(persistence, pers_dim=0):
    lengths = []
    dimensions = [[] for _ in range(5)]
    
    for i in range(len(persistence)):
        dim = persistence[i][0]
        length = persistence[i][1][1] - persistence[i][1][0]
        
        if dim < 5 and not math.isinf(length):
            dimensions[dim].append(length)
            lengths.append(length)
            
    if pers_dim in [0, 1, 2, 3, 4]:
        if dimensions[pers_dim]:
            return np.max(dimensions[pers_dim])
        else:
            return 0  # or return None if desired
    else:
        return None
    
# Get minimum persistence (lifespan) amongst features in a certain dim
def min_persistence(persistence, pers_dim=0):
    lengths = []
    dimensions = [[] for _ in range(5)]
    
    for i in range(len(persistence)):
        dim = persistence[i][0]
        length = persistence[i][1][1] - persistence[i][1][0]
        
        if dim < 5 and not math.isinf(length):
            dimensions[dim].append(length)
            lengths.append(length)
            
    if pers_dim in [0, 1, 2, 3, 4]:
        if dimensions[pers_dim]:
            return np.min(dimensions[pers_dim])
        else:
            return 0  # or return None if desired
    else:
        return None
    
# Get number of features based on dimension using persistence homology
def number_features(persistence, dim=0):
    dims = [[] for _ in range(5)]
    
    for i in range(len(persistence)):
        d = persistence[i][0]
        pt = persistence[i][1]
        
        if d < 5:
            dims[d].append(pt)
            
    if dim in [0, 1, 2, 3, 4]:
        num_features = len(dims[dim])
        return num_features
    else:
        return None
    
# Get number of features based on a threshold and dimension using persistence homology
def count_features(persistence, dim=0, threshold=17):
    count = 0
    
    for item in persistence:
        if item[0] == dim and item[1][-1] <= threshold:
            count += 1
    
    return count

# Get the classical Euler characteristics
def euler_characteristic(persistence):
    dims = [[] for _ in range(5)]
    
    for i in range(len(persistence)):
        d = persistence[i][0]
        bd = persistence[i][1]
        
        if d < 5:
            dims[d].append(bd)
    
    b0 = len(dims[0])
    b1 = len(dims[1])
    b2 = len(dims[2])
    b3 = len(dims[3])
    b4 = len(dims[4])

    euler = b0 - b1 + b2 - b3 + b4
    
    return euler

# Get topological descriptors with significant values
def significants(data, adjp):
    if len(adjp) == len(data):
        data['Both_Adj_P_Values'] = adjp
    else:
        # Adjust the length of the list to match the DataFrame
        adjp = adjp[:len(data)]
        data['Both_Adj_P_Values'] = adjp
    
    significant_desc = data[data['Both_Adj_P_Values'] < 0.05]
    
    return print(significant_desc.index.to_list())

