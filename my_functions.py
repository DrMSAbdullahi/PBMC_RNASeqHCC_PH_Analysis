import numpy as np
import pandas as pd

# Extract gene set given a specific pathway
def get_pathway_genes(kegg_map, pathway_id):
    gene_list = []
    for i in range(len(kegg_map)):
        kg = kegg_map.iloc[i, 2]
        if kegg_map.iloc[i, 1] == pathway_id:
            gene_list.append(kg)
    #len(pathway) 
    return gene_list

# Extract the genes from the experession data based on a pathway given
# HCC hsa05225 is set as default for the function
def extract_pathway_data(data, kegg_map, pathway_id='hsa05225'): 
    gene_list = []
    for i in range(len(kegg_map)):
        kg = kegg_map.iloc[i, 2]
        if kegg_map.iloc[i, 1] == pathway_id:
            gene_list.append(kg)
    
    result = data[data.index.isin(gene_list)]
    return result

# Get the persistence diagram from a given data
def get_persistence_diagram(data, max_dim=2, edge_length=10):
    # Convert data to numpy array
    data = np.array(data)
    
    # Calculate persistence diagrams
    rips_complex = gd.RipsComplex(points=data, max_edge_length=edge_length)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension=max_dim)
    persistence = simplex_tree.persistence()
    
    return persistence

# Get the persistence diagram of a large data through the lazy witness based on landmark points 
def get_witness_PD(data, max_dim=2, max_alpha=20, num_points=50):
    # Convert data to numpy array
    data = np.array(data)
    
    #Get the landmark points
    landmarks= gd.pick_n_random_points(points=data, nb_points=num_points)
    
    # Calculate witnessess and get persistence diagrams
    wit_complex = gd.EuclideanWitnessComplex(witnesses=data, landmarks=landmarks)
    wit_tree = wit_complex.create_simplex_tree(max_alpha_square=max_alpha,limit_dimension=max_dim)
    wit_persistence = wit_tree.persistence()
    
    return wit_persistence

def witness_PD(data, max_dim=2, max_alpha=20, num_points=50):
    n_rows = data.shape[0]
    n_cols = data.shape[1]
    metric = np.zeros([n_cols, n_cols])

    i = 0
    while i < n_cols:
        metric[i, i] = 1 / np.std(data.iloc[:, i])
        i += 1

    DataPrime = []
    j = 0
    while j < n_rows:
        DataPrime.append(metric.dot(data.iloc[j, :]))
        j += 1
    
    # Convert data to numpy array
    data = np.array(DataPrime)
    
    # Get the landmark points
    landmarks = gd.pick_n_random_points(points=data, nb_points=num_points)
    
    # Calculate witnesses and get persistence diagrams
    wit_complex = gd.EuclideanWitnessComplex(witnesses=data, landmarks=landmarks)
    wit_tree = wit_complex.create_simplex_tree(max_alpha_square=max_alpha, limit_dimension=max_dim)
    wit_persistence = wit_tree.persistence()
    
    return wit_persistence

# Plot the persistence diagram
def plot_PD(persistence, dim=-1):
    grouped_data = {}
    for tuple in persistence:
        dimension = tuple[0]
        if dimension not in grouped_data:
            grouped_data[dimension] = []
        grouped_data[dimension].append(tuple)
 
    if dim == 0:
        persistence = grouped_data[0]
        gd.plot_persistence_diagram(persistence)  
    elif dim == 1:
        persistence = grouped_data[1]
        gd.plot_persistence_diagram(persistence)
    elif dim == -1:
        gd.plot_persistence_diagram(persistence)
    else:
        None

# Plot the persistence barcode
def plot_PB(persistence, dim=-1):
    grouped_data = {}
    for tuple in persistence:
        dimension = tuple[0]
        if dimension not in grouped_data:
            grouped_data[dimension] = []
        grouped_data[dimension].append(tuple)
 
    if dim == 0:
        persistence = grouped_data[0]
        gd.plot_persistence_barcode(persistence)  
    elif dim == 1:
        persistence = grouped_data[1]
        gd.plot_persistence_barcode(persistence)
    elif dim == -1:
        gd.plot_persistence_barcode(persistence)
    else:
        None

# Get the betti numbers of a simplicial complex from a given data
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

# Get the sum of death times of a PD
def sum_death(persistence, pers_dim=0):
    sums = [[] for _ in range(3)]
    dims = [[] for _ in range(3)]
    
    for i in range(len(persistence)):
        dim = persistence[i][0]
        death = persistence[i][1][1]
        
        if dim < 3:
            dims[dim].append(death)
            sums[dim].append(death)
            
    sums = [np.array(sum_list) for sum_list in sums]
    
    if pers_dim in [0, 1, 2]:
        death_sum = np.sum(sums[pers_dim][~np.isinf(sums[pers_dim])])
        return death_sum
    else:
        return None
    
# Get the sum of birth times of a PD
def sum_birth(persistence, pers_dim=0):
    sums = [[] for _ in range(3)]
    dims = [[] for _ in range(3)]
    
    for i in range(len(persistence)):
        dim = persistence[i][0]
        birth = persistence[i][1][0]
        
        if dim < 3:
            dims[dim].append(birth)
            sums[dim].append(birth)
            
    sums = [np.array(sum_list) for sum_list in sums]
    
    if pers_dim in [0, 1, 2]:
        birth_sum = np.sum(sums[pers_dim][~np.isinf(sums[pers_dim])])
        return birth_sum
    else:
        return None

# Get the persistence of a certain dim features of a PD
def length_persistence(persistence, pers_dim=0):
    lengths = [[] for _ in range(3)]
    
    for i in range(len(persistence)):
        dim = persistence[i][0]
        length = persistence[i][1][1] - persistence[i][1][0]
        
        if dim < 3:
            lengths[dim].append(length)
            
    if pers_dim in [0, 1, 2]:
        return lengths[pers_dim]
    else:
        return None
    
# Get the sum of persistence of a certain dim features of a PD
def sum_persistence(persistence, pers_dim=0):
    sum_lengths = [0, 0, 0]
    
    for i in range(len(persistence)):
        dim = persistence[i][0]
        length = persistence[i][1][1] - persistence[i][1][0]
        
        if dim < 3 and not math.isinf(length):
            sum_lengths[dim] += length
            
    if pers_dim in [0, 1, 2]:
        return sum_lengths[pers_dim]
    else:
        return None
    
# Get the maximum persistence amongst features in a certain dim
def max_persistence(persistence, pers_dim=0):
    lengths = []
    dimensions = [[] for _ in range(3)]
    
    for i in range(len(persistence)):
        dim = persistence[i][0]
        length = persistence[i][1][1] - persistence[i][1][0]
        
        if dim < 3 and not math.isinf(length):
            dimensions[dim].append(length)
            lengths.append(length)
            
    if pers_dim in [0, 1, 2]:
        if dimensions[pers_dim]:
            return np.max(dimensions[pers_dim])
        else:
            return 0  # or return None if desired
    else:
        return None
    
# Get the number of traced by the persistence homology based on dimension of the features
def number_features(persistence, dim=0):
    dims = [[] for _ in range(3)]
    
    for i in range(len(persistence)):
        d = persistence[i][0]
        pt = persistence[i][1]
        
        if d < 3:
            dims[d].append(pt)
            
    if dim in [0, 1, 2]:
        num_features = len(dims[dim])
        return num_features
    else:
        return None
    
# Get the number of traced by the persistence homology at a threshold based on dimension of the features
def count_features(persistence, dim=0, threshold=17):
    count = 0
    
    for item in persistence:
        if item[0] == dim and item[1][-1] <= threshold:
            count += 1
    
    return count

#Get the Euler characteristics of a persistence
def euler_characteristic(persistence):
    dims = [[] for _ in range(3)]
    
    for i in range(len(persistence)):
        d = persistence[i][0]
        bd = persistence[i][1]
        
        if d < 3:
            dims[d].append(bd)
    
    vertices = len(dims[0])
    edges = len(dims[1])
    faces = len(dims[2])

    euler = vertices - edges + faces
    return euler

#Get the correlation matrix
def corr_matrix(data, corr_type="Similarity", by="Sample"):
    
    corr_data=[]
    
    if by == "Sample":
        if corr_type == "Similarity":
            corr_data = data.corr()
        elif corr_type == "Dissimilarity":
            corr_data = 1 - data.corr()
        else: None
    elif by == "Gene":
        data = data.T
        if corr_type == "Similarity":
            corr_data = data.corr()
        elif corr_type == "Dissimilarity":
            corr_data = 1 - data.corr()
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
