#include "cell_list.hpp"
#include <iostream>
int main() {
    //The box will be described by the lattice vectors
    const int nx=100;
    const int ny=100;
    const int nz=10;

    std::cout<<"Doing a search in "<<nx*ny*nz<<" sites"<<std::endl;
    std::vector<Site> site_list = std::vector<Site>(nx*ny*nz);

    const LatticeCell  lattice_vectors( {1.0*nx,0.0,0.0},
                                        {0.0,1.0*ny,0.0},
                                        {0.0,0.0,1.0*nz});

    //We will put some points
    int n=0;
    for( auto iz=0; iz < nz; iz++)
    for( auto iy=0; iy < ny; iy++)
    for( auto ix=0; ix < nx; ix++)
    {
        site_list[n]= Site(ix,iy,iz);
        n=n+1;      
    }
    double cutoff_radius = 1.0;
    CellList cell_list(lattice_vectors, cutoff_radius);
    cell_list.SetSites(site_list);
    cell_list.GetNeighborList();
    
    return 0;
}   

