#include "cell_list.hpp"

const int spatial_dim =3;


/////****************BEGIN: SET SITES FUNCTIONS********************///

/**
 * @brief Sets the sites for neighbor searching and inserts them into the cell grid.
 * @param sites The collection of sites to be considered.
 *
 * This function sets the sites for neighbor searching and inserts them into the cell grid. It iterates over
 * the given sites and inserts each site into the appropriate cell in the cell grid using the Insert function.
 */
void CellList::SetSites(const std::vector<Site> & sites){ // Function verified
    // Guard condition-> Here

    //Modifies the site and make it a cell site
    auto& supercell = this->GetSuperCell();
    num_sites_ = 0;
    for(const auto & site: sites)
    {
        auto cell_site = site;
        cell_site.SetIndex(num_sites_); 
        supercell.InsertSite(cell_site);
        num_sites_++;
    }
    return;
}


/////****************END: SET SITES FUNCTIONS********************///


std::vector< std::vector<int> > CellList::GetNeighborList()
{
    auto& supercell = this->GetSuperCell();
    std::vector< std::vector<int> > neighborStruct( this->GetNumberOfSites() );

    //Here we determine the maximum number of neighbors
    Integer max_num_neighbors = 0 ;
    for (const auto& cell: supercell.GetCells())
    if (max_num_neighbors < cell.size() )
        max_num_neighbors = cell.size();
    std::cout<<"max_num_neighbors"<<max_num_neighbors<<std::endl;        

    //We reserve the neighbor list
    for( auto& site_neighbors: neighborStruct )
        site_neighbors.reserve(max_num_neighbors);



    for (const auto& cell: supercell.GetCells())
        for(const auto& site: cell)
        {
            //std::cout<<"Querying site"<<std::endl;
            //std::cout<<"SiteID: "<<site.GetIndex()<<std::endl;
            //std::cout<<"SitePos: "<<site.GetPosition()[0]<<" "<<site.GetPosition()[1]<<" "<<site.GetPosition()[2]<<std::endl;
            //std::cout<<"SiteCell: "<<supercell.GetSiteHostCell(site)(0)<<" "<<supercell.GetSiteHostCell(site)(1)<<" "<<supercell.GetSiteHostCell(site)(2)<<std::endl;
            //std::cout<<"CellIndex: "<<supercell.GetCellIndex(site)<<std::endl<<std::endl;

            //Now we need to iterate over sites at neighboring cells
             //Iterate over all sites in neighbor cells
            for(int di0=-1; di0<=1; di0++ )
            for(int di1=-1; di1<=1; di1++ )
            for(int di2=-1; di2<=1; di2++ )
            {
                Triplet neighbor_triplet = supercell.GetSiteHostCell(site) + Triplet(di0,di1,di2);
                if( supercell.ValidCell(neighbor_triplet) )
                for(const auto& neighbor_site: supercell.GetCell(neighbor_triplet))
                    if ( site.Distance(neighbor_site) < this->GetCutoffRadius() )
                        neighborStruct[site.GetIndex()].emplace_back(neighbor_site.GetIndex());                    
            }
        }        
return neighborStruct;
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
    Vector3D minCellLen(0,0,0);
    //printf("Cutof Cube in lattice vectors\n");
    for( int dir=0; dir < 3; dir++)
    {
        const auto& fracV = lattice.FractionalCoords(cutoff_cube.GetVector(dir));        
        //Get the maximum length 
        for( int comp=0; comp < 3; comp ++)
        if( minCellLen[comp]< std::abs(fracV[comp]) ) 
            minCellLen[comp]= std::abs(fracV[comp]);
    }
    //Before, we determine the minimum cell length
    //Nevertheless, we must pointed out, that we need an integer number
    //of the lattice cell. Therefore, we will take the ceil of the
    //inverse and the invert it again. That way, it will always 
    //be an integer. The ceil of the inverse is nothing
    //but the number of cells
    Triplet numCells;
    for( int dir=0; dir < 3; dir++)
    {
        numCells[dir] = (Integer)std::ceil(1/minCellLen[dir]);
        minCellLen[dir] = 1.0/(Real)numCells[dir];
    }
    //std::cout<<std::endl;

    //We know use this length to construct an adequate cell
    for( int i=0; i < 3; i++)
        cell_.SetVector(i, lattice.GetVector(i)*minCellLen[i] );
    supercell_ = LatticeSupercell(numCells,cell_);

    //printf("Cutof Cube in cell vectors\nShould be larger than lattice less or equal to one\n");
    for( int dir=0; dir < 3; dir++)
    {
        const auto& fracV = supercell_.FractionalCoords(cutoff_cube.GetVector(dir));      
        //std::cout<<fracV[0]<<" "<<fracV[1]<<" "<<fracV[2]<<std::endl;
    }
    
    
    //std::cout<<" FIN OF THE SPHEREW"<<std::endl;



}
