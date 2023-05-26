#include "neighbors.hpp"



void CellList::GetNeighborList()
{
    std::vector< std::vector<int> > neighborStruct( this->GetNumberOfSites() );
    auto cellmap = this->GetGrid();

    for (const auto& current_cell_pair  :cellmap ) 
    {
        //Iterate over all sites in each cell 
        auto current_cell = current_cell_pair.second;
        for( auto& current_site : current_cell)
        {
            const auto current_cell_triple = current_site.GetHostCell();
            const auto current_index = current_site.GetIndex();
 
             //Iterate over all sites in neighbor cells
            for(int di0=-1; di0<=1; di0++ )
            for(int di1=-1; di1<=1; di1++ )
            for(int di2=-1; di2<=1; di2++ )
            {
                auto neighbor_cell_triple= current_cell_triple+ Triplet(di0,di1,di2);
                auto neighboh_cell_key = TripletToString(current_cell_triple);
                auto neighbor_cell = cellmap[neighboh_cell_key];
                for (auto& neighbor_site : neighbor_cell ) 
                {
                    const auto neighbor_index = current_site.GetIndex();
                    Real distance = current_site.Distance(neighbor_site);
                    if ( distance < this->GetCutoffRadius() )
                        neighborStruct[current_index].emplace_back(neighbor_index);
                }
            }
        }
    }
}


const int spatial_dim =3;
void CellList::SetCellVectors(const LatticeCell lattice_vectors, const Real cutoff_radius){
    //cutoff guard condition ->    if( cutoff_radius<= 0 ) { exit(EXIT_FAILURE);   }
    //Lattice vector guard condition -> 
    this->cutoff_radius_ = cutoff_radius;
    this->lattice_vectors_= lattice_vectors;
    cell_vectors_ = lattice_vectors;

    //To get the minium cell that fix the cutoff sphere
    // I will extend the region to a cube
    // of vertex  (2*cutoff,0,0), (0,2*cutoff,0), (0,2*cutoff,0)
    // to get the cell, we can just express these vectors 
    // in fractionall coordinates of the lattice vectors
    // and the the maximum n
    
    const auto D = 2*cutoff_radius;
    const LatticeCell  cutoff_cube( {D,0.0,0.0},
                                    {0.0,D,0.0},
                                    {0.0,0.0,D});
    for( int dir=0; dir < 3; dir++)
    {
        const auto fracV = lattice_vectors.FractionalCoords(cutoff_cube.GetVector(dir));
        std::cout<<fracV[0]<<" "<<fracV[1]<<" "<<fracV[2]<<std::endl; 
    }
}

Triplet CellList::CellTriplet(const Site& site)
{
    const auto pos = site.GetPosition();
    Triplet cell_triplet;
    
    auto fraccoord = this->GetCellVectors().FractionalCoords(pos);
    for(auto i=0; i< spatial_dim; i++){
        cell_triplet[i] = (size_t) std::floor(fraccoord[i]);
    }
    return cell_triplet;
}

//Checke
void CellList::SetSites(std::vector<Site> & sites){
    // Guard condition-> Here
    for(auto & site: sites)
        Insert(site);
    return;
}

void CellList::Insert(Site & site){
    num_sites_+=1; // increase number of sites

    //Make the class determine the cell triplet
    // and pass it to the site as host cell 
    auto cell_triplet= this->CellTriplet(site);
    site.SetHostCell(cell_triplet);

    //Convert the triplet into a tag and use it as key
    //in the map
    auto cell_tag    = TripletToString(cell_triplet);
    if(this->GetGrid().count(cell_tag) == 0)
        this->GetGrid()[cell_tag]= std::vector<Site>();
    this->GetGrid()[cell_tag].emplace_back(site);
    return;
}

int main() {
    //The box will be described by the lattice vectors
    const int nx=10;
    const int ny=10;
    const int nz=10;

    std::vector<Site> site_list;
    const LatticeCell  lattice_vectors({10.0,0.0,0.0},{0.0,10.0,0.0},{0.0,0.0,10.0});

    //We will put some points
    for( auto ix=0; ix < nx; ix++)
    for( auto iy=0; iy < ny; iy++)
    for( auto iz=0; iz < nz; iz++)
      site_list.push_back( Site(ix,iy,iz));
      
      double cutoff_radius = 1.0;
      CellList cell_list(lattice_vectors, cutoff_radius);
      cell_list.SetSites(site_list);
      cell_list.GetNeighborList();
    
    return 0;
}   

