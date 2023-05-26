import numpy as np
from numpy import sqrt, exp, dot, conj, min,max, abs
import matplotlib.pyplot as plt
import lsquant as qtr

"""
Let us first define a hopping function, which will set a hopping
for two sites depending on their type and position. For graphene, it would
be t for nearest neighbors, which we define as sites where the module of 
the differences |R_i - R_j| is less equal d_nn = 2.46  
"""

import numpy as np
from scipy.spatial import KDTree


def get_neighbor_list(tree, points, cutoff_radius):
    # Get all neighbors for a given cutoff the tree for points within the cutoff distance
    neighbors_list = tree.query_ball_point(points, cutoff_radius)

    # Determine maximum number of neighbors 
    max_neighbors_number = max([len(n) for n in neighbors_list])

    # Create a structured numpy array to load the neighbor list 
    # if the number of neighbors for a ith-entry is less than the maximum
    #fill with -1, else use correct values
    result = np.full((len(neighbors_list), max_neighbors_number), -1)
    for i, neighbors in enumerate(neighbors_list):
        result[i, :len(neighbors)] = neighbors
        
    return result


def get_minimal_supercell( lat_vec, radius, pbc=(False, False, False) ):
    scdim =[1,1,1]
    for i,(ai, periodic) in enumerate(zip(lat_vec.T, pbc)):
        if periodic:
            scdim[i] = np.ceil(radius/np.linalg.norm(ai)).astype(int)
    return scdim


def extend_cell_to_radius( lat_vec, orb_pos, radius, pbc=(False, False, False) ):
    scdim =[1,1,1]
    for i,(ai, periodic) in enumerate(zip(lat_vec.T, pbc)):
        if periodic:
            scdim[i] = np.ceil(radius/np.linalg.norm(ai)).astype(int)
    return create_supercell( scdim, lat_vec, orb_pos)

def create_supercell( scdim, lat_vec, orb_pos):
    cell_pos = lat_vec@([ xi.flatten() for xi in np.mgrid[0:scdim[0],0:scdim[1], 0:scdim[2]]])
    positions = cell_pos.T + orb_pos[0]
    for ro in orb_pos[1:]:
        positions = np.concatenate((positions, cell_pos.T + ro))
    return positions


ncell    = (4,4,1)
orb_pos  = np.array([[0,0,0], [-1,0,0] ])
lat_vec  = np.matrix([ [ 3/2, sqrt(3)/2,0 ], [ 3/2,-sqrt(3)/2,0 ], [ 0,0,1 ] ]).T

max_hop_dist = 4;
pbc=(True, True, False)
min_scdim = get_minimal_supercell( lat_vec, max_hop_dist, pbc=pbc ) 
ext_scdim = [ n+1 if periodic else n for n, periodic in zip(min_scdim, pbc) ]

ext_positions = create_supercell(ext_scdim,lat_vec, orb_pos) 

center = np.mean(ext_positions, axis=0);

# Construct a KDTree
tree = KDTree(ext_positions)
neighbors_list = get_neighbor_list(tree, ext_positions, cutoff_radius= max_hop_dist) 
print( neighbors_list)








def hopping_function(site_i, site_j, model_params):
    pass
    #lat_vec = [ [ 1/2, sqrt(3)/2,0 ], [ 1/2,-sqrt(3)/2,0 ], [ 0,0,1 ] ]
#def hamiltonian(k):
#    a_0, a_1, a2  = lat_vec;
#    hop = 2.8;
#    f_k = hop*( 1 + exp( -1j*dot(k,a_0)) + exp( -1j*dot(k,a_1)) );
#    return [ [ 0        , f_k],
#             [ conj(f_k),  0 ]
#            ];

#The band class requires lattice vectors and the hamiltonian function
#graphene = k.bandstructure(lat_vec, hamiltonian );

#To plot a desire band-path you pass it to the class as a list of tuples
#npts= 100;
#bandpath = [ ("K", (1/3,2/3,0), 111 ), ("G", (0,0,0), 35) , ("M",(1/2,1/2,0),55), ("K'",(2/3,1/3,0),1) ];
#graphene.set_bandpath(bandpath);

#The computing the band structure is as simmple as
#sigma_x,sigma_y = [ [[0,1],[1,0]], [[0,-1j],[1j,0]] ];
#bandstructure = graphene.compute_bands( proj_ops=[ sigma_x,sigma_y] );

#fig= plt.gcf();
#ax = plt.gca();
#xaxis = graphene.Xaxis();
#for proj_band in bandstructure:
#    band,sigma_x,sigma_y = proj_band;
#    z   = sigma_x;
#    s   =  20*abs(z);
#    plt.plot(xaxis,band, c="k");
#    im  = ax.scatter(xaxis,band,s=s, c=z,cmap="coolwarm",vmin=-1, vmax=1);
#fig.colorbar(im, ax=ax);

#xlabels= graphene.XLabels();
#ax.set_xticks(xlabels[0])
#ax.set_xticklabels(xlabels[1])
#plt.savefig('proj_band_sigma_x.pdf');
#plt.show();


