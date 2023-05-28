#include "neighbors.hpp"

const int spatial_dim =3;
//default tag ? 
std::string TripletToString(const Triplet& arr) {
    std::ostringstream oss;
    for (size_t i = 0; i < 3; ++i) 
        oss << arr[i]<<" ";
    return oss.str();
}
//    for(int i=0; i<spatial_dim; i++)
//    {
 //       const auto xi = cell.GetVector(i);
  //      std::cout<<xi[0]<<" "<<xi[1]<<" "<<xi[2]<<std::endl;
   // }


/////****************BEGIN: SET SITES FUNCTIONS********************///

/**
 * @brief Sets the sites for neighbor searching and inserts them into the cell grid.
 * @param sites The collection of sites to be considered.
 *
 * This function sets the sites for neighbor searching and inserts them into the cell grid. It iterates over
 * the given sites and inserts each site into the appropriate cell in the cell grid using the Insert function.
 */
void CellList::SetSites(std::vector<Site> & sites){ // Function verified
    // Guard condition-> Here
    for(auto & site: sites)
        Insert(site);
    return;
}

/**
 * @brief Inserts a site into the cell grid.
 * @param site The site to be inserted.
 *
 * This function inserts a site into the cell grid. It increases the number of sites in the
 * cell list by one and determines the site's host cell using the CellTriplet function. The host cell information
 * is then assigned to the site. The site is inserted into the appropriate cell in the cell grid using a cell tag
 * derived from the host cell's triplet.
 */
void CellList::Insert(Site & site){ // Function verified
    num_sites_+=1; // increase number of sites
    this->SetSiteHostCell(site); //Set the host cell index
    const auto& cell_triplet = site.GetHostCell();

    //Convert the triplet into a tag and use it as key
    //in the map
    const auto& cell_key = TripletToString(cell_triplet);
    auto& cell_map = this->GetCellMap();
    
    // gf the cell does not exists 
    // create an empty vector
    if(cell_map.count(cell_key) == 0) 
        cell_map[cell_key]= std::vector<Site>();
        
    //Add a site to the cell
    cell_map[cell_key].emplace_back(site);
    return;
}

/**
 * @brief Computes the cell triplet for a given site.
 * @param site The site for which to compute the cell triplet.
 * @return The cell triplet.
 *
 * This function computes the cell triplet for a given site. It extracts the position of the site using GetPosition,
 * and then calculates the fractional coordinates of the site using the GetCellVectors function. The fractional coordinates
 * are converted to the nearest integer values using floor, and the resulting values are stored in the cell triplet.
 * The computed cell triplet represents the cell in the cell grid to which the site belongs.
 */
void CellList::SetSiteHostCell(Site& site) // Function verified
{
    const auto& cell = this->GetCell();  
    const auto& fpos = cell.FractionalCoords( site.GetPosition() );

    Triplet triplet;
    for(auto i=0; i< spatial_dim; i++)
        triplet[i] = (size_t) std::floor(fpos[i]);
    site.SetHostCell(triplet);
}

/////****************END: SET SITES FUNCTIONS********************///


void CellList::GetNeighborList()
{
    std::vector< std::vector<int> > neighborStruct( this->GetNumberOfSites() );
    
    const auto& cell_map = this->GetCellMap();
    for (const auto& current_cell_pair  :cell_map ) 
    {
        //Iterate over all sites in each cell 
        const auto current_cell = current_cell_pair.second;
        for(const auto& current_site : current_cell)
        {
            const auto& current_hostcell = current_site.GetHostCell();
            const auto& current_index = current_site.GetIndex();

             //Iterate over all sites in neighbor cells
            for(int di0=-1; di0<=1; di0++ )
            for(int di1=-1; di1<=1; di1++ )
            for(int di2=-1; di2<=1; di2++ )
            {
                //PERIODICITY SHOULD BE ENFORCED HERE
                const auto neighbor_hostcell = current_hostcell+ Triplet(di0,di1,di2);
                const auto neighboh_cellkey = TripletToString(neighbor_hostcell);

                //We exclude all cells that does not exists.
                //this is equivalent to assume open boundary conditions
                if( cell_map.count(neighboh_cellkey) !=0 )
                {
                    const auto& neighbor_cell = cell_map.at(neighboh_cellkey);
                    for (auto& neighbor_site : neighbor_cell ) 
                    {
                        const auto& neighbor_index = current_site.GetIndex();
                        const auto& distance = current_site.Distance(neighbor_site);
                        if ( distance < this->GetCutoffRadius() )
                            neighborStruct[current_index].emplace_back(neighbor_index);
                    }
                }
            }
        }
    }
}

void CellList::SetCellVectors(const LatticeCell lattice, const Real cutoff_radius){
    //cutoff guard condition ->    if( cutoff_radius<= 0 ) { exit(EXIT_FAILURE);   }
    //Lattice vector guard condition -> 
    this->cutoff_radius_ = cutoff_radius;
    this->lattice_= lattice;

    //To get the minium cell that fix the cutoff sphere
    // I will extend the region to a cube
    // of vertex  (2*cutoff,0,0), (0,2*cutoff,0), (0,2*cutoff,0)
    // to get the cell, we can just express these vectors 
    // in fractionall coordinates of the lattice vectors
    // Then, we can just choose the maximum length
    // in fraccional space to know a potential 
    // box that fit the sphere 
    
    const auto D = 2*cutoff_radius;
    const LatticeCell  cutoff_cube( {D,0.0,0.0},
                                    {0.0,D,0.0},
                                    {0.0,0.0,D});
    Vector3D maxLen(0,0,0);
    //printf("Cutof Cube in lattice vectors\n");
    for( int dir=0; dir < 3; dir++)
    {
        const auto& fracV = lattice.FractionalCoords(cutoff_cube.GetVector(dir));        
        //std::cout<<fracV[0]<<" "<<fracV[1]<<" "<<fracV[2]<<std::endl;

        //Get the maximum length 
        for( int comp=0; comp < 3; comp ++)
        if( maxLen[comp]< std::abs(fracV[comp]) ) maxLen[comp]= std::abs(fracV[comp]);
    }
    //std::cout<<std::endl;

    //We know use this length to construct an adequate cell
    for( int i=0; i < 3; i++)
        cell_.SetVector(i, lattice.GetVector(i)*maxLen[i] );

    //printf("Cutof Cube in cell vectors\nShould be larger than lattice less or equal to one\n");
    for( int dir=0; dir < 3; dir++)
    {
        const auto& fracV = cell_.FractionalCoords(cutoff_cube.GetVector(dir));      
        //std::cout<<fracV[0]<<" "<<fracV[1]<<" "<<fracV[2]<<std::endl;
    }
    //std::cout<<" FIN OF THE SPHEREW"<<std::endl;

}

int main() {
    //The box will be described by the lattice vectors
    const int nx=1000;
    const int ny=1000;
    const int nz=1;

    std::vector<Site> site_list;
    const LatticeCell  lattice_vectors( {1.0*nx,0.0,0.0},
                                        {0.0,1.0*ny,0.0},
                                        {0.0,0.0,1.0*nz});

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

