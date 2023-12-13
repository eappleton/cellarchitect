import numpy as np
from scipy.spatial.distance import euclidean
from scipy.spatial import ConvexHull
import pylab as plt
import pandas as pd
import networkx as nx
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#from skimage import measure
#import ShapeStatistics
import math

"""
Eigenvector Similarity

calculate the Laplacian eigenvalues for the adjacency matrices of each of the graphs.
For each graph, find the smallest k such that the sum of the k largest eigenvalues
constitutes at least 90% of the sum of all of the eigenvalues.
If the values of k are different between the two graphs, then use the smaller one.
The similarity metric is the sum of the squared differences between the largest k eigenvalues
between the graphs. This will produce a similarity metric in the range [0, inf),
where values closer to zero are more similar.

Further analysis must be done before using this method, e.g. is the number of nodes/edges

A combination of grapgh similarity metrics, number of the connected components,
sphericity(Isoperimetric Quotient=(36*np.pi*V**2)/(S**3)) can be used for general case. 
To Do: check a case with negative eigen values 
"""
def read_points(filename):
    """
    Function: read csv file and extract coordinates as a 2D list.
    Returns points from the largest connected components. 
    
    Steps:
        1. Find the largest connected componet of cells.
        2. Translate and scale to unit sphere. issue with edge effects

    Returns
    ----------
    network - networkx grapgh for using the eigen-value similarity
    #For spherical Harmonics
    points_list - list of lists. [[x,y,z],[x,y,z],[x,y,z]]
    points_norm - points scaled and translated to unit coordinates [-1,1]
    """
    points_list = []
    cell_radius = 6.0
    
    #read CSV file
    dataframe = pd.read_csv(filename,header=0)
    
    ### Find largest connected compoent. We want largest single object.
    # 1A. create cell neighborhood network
    network = nx.Graph()
    neighbor_threshold = 20.0
    for m, row_m in dataframe.iterrows():
        for n, row_n in dataframe.iterrows():
            coord_m = [row_m["x"],row_m["y"],row_m["z"]]
            coord_n = [row_n["x"],row_n["y"],row_n["z"]]
            if euclidean(coord_m,coord_n) <= neighbor_threshold:
                network.add_edge(m,n)
            else:
                network.add_nodes_from([m,n])
    
    # 1B. Find largest connected component
    all_connected_componets = sorted(nx.connected_components(network), key = len, reverse=True)
    max_connected_component = all_connected_componets[0]
    #plt.figure()
    #nx.draw_networkx(network)    
    #print(str(len(all_connected_componets ))+"connected components")      
    if len(max_connected_component)*2<len(network):
        return None, None, None, None
    # 1C. Extract 3D points from dataframe
    for idx, row in dataframe.iterrows():
        if idx in max_connected_component:
            points_list.append([row["x"],row["y"],row["z"]])
    
    # 2. translate and scale input to unit sphere
    points_norm = np.array(points_list)
    center = np.mean(points_norm, axis=0)
    points_norm = points_norm - center
    maxd = np.max(np.sqrt(np.sum(points_norm**2, axis=1)))
    maxd += 2.0*cell_radius
    points_norm /=  maxd# normalized to unit circle [-R,R]^3
    points_norm += 1 #recenter in positive euclidean space of size [0,2R]^3
    
    scale_factor = neighbor_threshold/maxd
    
    #num_connected_componets=len(all_connected_componets )
    
    return network,points_list,points_norm,scale_factor#,num_connected_componets

def select_k(spectrum, minimum_energy = 0.9):
    running_total = 0.0
    total = sum(spectrum)
    if total == 0.0:
        return len(spectrum)
    for i in range(len(spectrum)):
        running_total += spectrum[i]
        if running_total / total >= minimum_energy:
            return i + 1
    return len(spectrum)

def asSpherical(points):
   
    r_theta_phi = np.zeros(np.array(points).shape)
    for idx,xyz in enumerate(points):
        #takes list xyz (single coord)
        x = xyz[0]
        y = xyz[1]
        z = xyz[2]
        r =  np.sqrt(x*x + y*y + z*z) #radial [0,r]
        phi = math.atan2(y,x) #azimuthal. [0,2pi]
        theta = math.acos(z/r) #polar/zenith [0,pi]
        r_theta_phi[idx]=[r,theta,phi]
    return r_theta_phi
def SphericalHarmonicFeatures(m,l,n,theta,phi):
    f=0 #To Do: find how to compute f
    fn=[f[np.random.randint(len(theta))][np.random.randint(len(phi))] for i in range(3)]
    SHF=[[np.sum(fn) for mm in range(-ll+1,ll)] for ll in range(l)]

    return SHF  

def Sphericity(points): 
    
    hull = ConvexHull(np.array(points))
    V=hull.volume
    S=hull.area
    IQ=(36*np.pi*V**2)/(S**3)
    return IQ
    
    
if __name__ == "__main__":
    filename1 = "/Users/noushinm/Downloads/2018_10_26_Sim_Runs/Lim_Green_Red_2_2018-10-26_12-54-50/cell_positions.csv.160"#"/Users/noushinm/Downloads/Shape metric/cell_positions.csv.70"
    filename2 = "/Users/noushinm/Downloads/2018_10_26_Sim_Runs/Lim_Green_Red_2_2018-10-26_12-42-44/cell_positions.csv.160"#"/Users/noushinm/Downloads/Shape metric/cell_positions.csv.70"
    
    network1,points_list_1,points_norm_1,scaled_diameter_1 = read_points(filename1)
    network2,points_list_2,points_norm_2,scaled_diameter_2 = read_points(filename2)


    laplacian1 = nx.spectrum.laplacian_spectrum(network1)
    #laplacian2 = nx.spectrum.laplacian_spectrum(network1) #self similarity is 0
    laplacian2 = nx.spectrum.laplacian_spectrum(network2)
    
    k1 = select_k(laplacian1)
    k2 = select_k(laplacian2)
    k = min(k1, k2)
    # We can set a threshold for dropping graphs with larger similarity, e.g. similarity > 50
    similarity = sum((laplacian1[:k] - laplacian2[:k])**2)
    print("similarity is " +str(similarity))
    # Clustering based on Sphericity 
    IQ1=Sphericity(points_list_1)
    IQ2=Sphericity(points_list_2)
    print("Sphericity of grapghs are %s and %s" %(IQ1,IQ2))