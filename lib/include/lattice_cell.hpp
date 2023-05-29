#ifndef LATTICE_CELL_HPP
#define LATTICE_CELL_HPP
#include <list>
#include "typedef.hpp" //Define the types: Matrix3D, Vector3D, Triplet, Float, Real, Integer;
#include "site.hpp"

/**
 * @brief Class representing a lattice cell.
 */
class LatticeCell {
public:
    /**
     * @brief Constructs a lattice cell with the given lattice vectors.
     * @param lat0 The first lattice vector (default: [1, 0, 0]).
     * @param lat1 The second lattice vector (default: [0, 1, 0]).
     * @param lat2 The third lattice vector (default: [0, 0, 1]).
     */
    explicit 
    LatticeCell(const std::array<Real,3> lat0={1,0,0},
                const std::array<Real,3> lat1={0,1,0},
                const std::array<Real,3> lat2={0,0,1})
    {
        for(int i=0; i < 3; i++)
        {
            lattice_vectors_(i,0) = lat0[i];
            lattice_vectors_(i,1) = lat1[i];
            lattice_vectors_(i,2) = lat2[i];
        }
    }

    /**
     * @brief Returns the lattice vector at the specified indices.
     * @param i The row index.
     * @param j The column index.
     * @return The lattice vector.
     */
    Real& operator()(const Integer i, const Integer j) { return lattice_vectors_(i,j);}

    /**
     * @brief Returns the lattice vectors of the cell.
     * @return The lattice vectors.
     */
    const Matrix3D& GetVectors() const {
        return lattice_vectors_;
    }

    /**
     * @brief Returns the lattice vector at the specified index.
     * @param i The column index.
     * @return The lattice vector.
     */    
    const Vector3D GetVector(const Integer i) const {
        return lattice_vectors_.col(i);
    }

    /**
     * @brief Sets the lattice vector at the specified index.
     * @param i The column index.
     * @param x_i The new lattice vector.
     */
    void SetVector(const Integer i, const Vector3D& x_i){
        lattice_vectors_.col(i) = x_i;
    }

    /**
     * @brief Computes the norms of the lattice vectors.
     * @return The norms of the lattice vectors.
     */
    std::array<Real, 3> Norms() const 
    {
        std::array<Real, 3> norms;
        for (size_t i = 0; i < 3; ++i) {
            norms[i] = lattice_vectors_.col(i).norm();
        }
        return norms;
    }

    /**
     * @brief Converts Cartesian coordinates to fractional coordinates.
     * @param X The Cartesian coordinates.
     * @return The fractional coordinates.
     */
    Vector3D FractionalCoords(const Vector3D& X) const {
        return lattice_vectors_.inverse()*X;
    }

    /**
     * @brief Set Periodic Boundary conditions
     * @param i0 True if the system is periodic along the 0 direction
     * @param i1 True if the system is periodic along the 1 direction
     * @param i2 True if the system is periodic along the 2 direction
     */
    void SetPBC(const bool i0,const bool i1,const bool i2)
    {
        pbc_= {i0,i1,i2};
    }

    /**
     * @brief Get true if the system is perdiodic along the i direction
     * @param i0 Direction to check for periodicity
     */
    bool IsPeriodicAt(const int i) const
    {
        return pbc_[i];
    }
    
private:
    Matrix3D lattice_vectors_; ///< The lattice vectors of the cell.
    std::array<bool,3> pbc_={true,true,true} ;//array containing the periodic boundary conditions
};


/**
 * @brief Class representing a Lattice constructed as a super cell.
 */
class LatticeSupercell : public LatticeCell {
public:
    /**
     * @brief Constructs a supercell with the given lattice vectors and supercell dimensions.
     * @param lat0 The first lattice vector (default: [1, 0, 0]).
     * @param lat1 The second lattice vector (default: [0, 1, 0]).
     * @param lat2 The third lattice vector (default: [0, 0, 1]).
     * @param dim0 The supercell dimension along the first lattice vector (default: 1).
     * @param dim1 The supercell dimension along the second lattice vector (default: 1).
     * @param dim2 The supercell dimension along the third lattice vector (default: 1).
     */
    explicit LatticeSupercell(  const Integer dim0 = 1,
                                const std::array<Real, 3> lat0 = {1, 0, 0},
                                const Integer dim1= 1,
                                const std::array<Real, 3> lat1 = {0, 1, 0},
                                const Integer dim2=1,
                                const std::array<Real, 3> lat2 = {0, 0, 1})
        : LatticeCell(lat0, lat1, lat2), scdim_(dim0,dim1,dim2)
        {
            cell_sites_= std::vector< std::list<Site> >(dim0*dim1*dim2);
        }

    /**
     * @brief Constructs a supercell using a cell and supercell dimensions.
     * @param dim0 The supercell dimension along the first lattice vector (default: 1).
     * @param dim1 The supercell dimension along the second lattice vector (default: 1).
     * @param dim2 The supercell dimension along the third lattice vector (default: 1).
     * @param cell A instance of the latticeCell 
     */    
    explicit LatticeSupercell(const Triplet dim, const LatticeCell& lat_cell)
        : LatticeCell(lat_cell), scdim_(dim){
            cell_sites_= std::vector< std::list<Site> >(dim[0]*dim[1]*dim[2]);
        }

    
    /**
     * @brief Insert a site at a particular cell based on the site position
     * @param site Site to be inserted.
     */    
    void InsertSite(const Site& site)
    {
        auto index = this->GetCellIndex(site);
        cell_sites_[index].insert(cell_sites_[index].end(), site);


//        std::cout<<"SiteID: "<<site.GetIndex()<<std::endl;
//        std::cout<<"SitePos: "<<site.GetPosition()[0]<<" "<<site.GetPosition()[1]<<" "<<site.GetPosition()[2]<<std::endl;
//        std::cout<<"SiteCell: "<<this->GetSiteHostCell(site)(0)<<" "<<this->GetSiteHostCell(site)(1)<<" "<<this->GetSiteHostCell(site)(2)<<std::endl;
//        std::cout<<"CellIndex: "<<this->GetCellIndex(site)<<std::endl<<std::endl;

    }

    /**
     * @brief Get a site at a particular cell based on the site position
     * @return site Site to be inserted.
     */    
    const std::vector< std::list<Site> >& GetCells()
    {
        return cell_sites_;
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
    Triplet GetSiteHostCell(const Site& site) // Function verified
    {
        const auto& fpos = this->FractionalCoords( site.GetPosition() );
        Triplet triplet;
        for(auto i=0; i< 3; i++)
            triplet[i] = (size_t) std::floor(fpos[i]);
    return this->FoldToSuperCell(triplet);
    }

    /**
     * @brief Fold a particular cell triplet inside the supercell
     * @Triplet cell triplet to be folded.
     */    
    Triplet FoldToSuperCell(const Triplet& host_cell_triplet) const
     {
        auto triplet = host_cell_triplet;
        for(int i = 0 ; i <3; i++)
        if(LatticeCell::IsPeriodicAt(i) )
        {
            const auto& scdim = this->GetSuperCellNumber(i);
            triplet[i] = (triplet[i]+ scdim)%scdim;
        }
        return triplet;
    }

    /**
     * @brief Get the Cell index of a given site based on its host cell triplet
     * @Site Site to be checked
     */    
    size_t ValidCell(const Triplet& triplet) 
     {
        bool valid = true;
        for(int i = 0; i < 3; i++)
        if( ! LatticeCell::IsPeriodicAt(i) )
            valid *= (triplet[i]>=0 && triplet[i]< this->GetSuperCellNumber(i) );  
        return valid;
    }


    /**
     * @brief Get the Cell index of a given site based on its host cell triplet
     * @Site Site to be checked
     */    
    size_t GetCellIndex(const Triplet& triplet) 
     {
        auto fold_triplet = this->FoldToSuperCell(triplet);      
        return ( fold_triplet[2]*GetSuperCellNumber(1) + fold_triplet[1]) *GetSuperCellNumber(0) + fold_triplet[0];
    }


    /**
     * @brief Get the Cell index of a given site based on its host cell triplet
     * @Site Site to be checked
     */    
    size_t GetCellIndex(const Site& site) 
     {
        auto triplet = this->FoldToSuperCell(this->GetSiteHostCell(site));      
        return ( triplet[2]*GetSuperCellNumber(1) + triplet[1]) *GetSuperCellNumber(0) + triplet[0];
    }


    /**
     * @brief Get the Cell index of a given site based on its host cell triplet
     * @Site Site to be checked
     */    
    const std::list<Site>& GetCell(const Triplet& cell_index) 
    {
        const size_t index = this->GetCellIndex(cell_index);
        return cell_sites_[index];
    }


    /**
     * @brief Returns the supercell dimension along the first lattice vector.
     * @return The supercell dimension.
     */
    Integer GetSuperCellNumber(const Integer dir) const { return scdim_[dir]; }


    /**
     * @brief Sets the supercell dimension along the first lattice vector.
     * @param dim The supercell dimension.
     */
    void SetSuperCellNumber(const Integer dir,const Integer scdim) { scdim_[dir] = scdim; }

    /**
     * @brief Sets the supercell dimension along the second lattice vector.
     * @param dim The supercell dimension.
     */



private:
    Triplet scdim_; ///< The supercell dimension along the first lattice vector.
    std::vector< std::list<Site> > cell_sites_;
};



#endif