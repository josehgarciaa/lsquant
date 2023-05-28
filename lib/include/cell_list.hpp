#ifndef CELL_LIST_HPP
#define CELL_LIST_HPP

#include <iostream>   // for std::cout, std::endl, etc.
#include <vector>     // for std::vector
#include <cstddef>    // for size_t
#include <unordered_map> // for std::unordered_map


#include "typedef.hpp" //Define the types: Matrix3D, Vector3D, Triplet, Float, Real, Integer;
#include "site.hpp" //Defines the Sites object
#include "lattice_cell.hpp" //Defines the Sites object



typedef std::unordered_map<std::string,std::vector<Site> > Map;

/**
 * @brief Class for computing nearest neighbors using the CellList algorithm.
 */
class CellList {
 public:
  /**
   * @brief Constructs a CellList object.
   * @param lattice_vectors The lattice vectors defining the periodic boundaries.
   * @param cutoff_radius The cutoff radius for neighbor searching.
   */
  explicit CellList(const LatticeCell& lattice, 
                    const double cutoff_radius)
                    :cutoff_radius_(cutoff_radius) 
  {
    // We will cover the space with another lattice
    // where each cells contains all points within a 
    // sphere with radius R= cutoff_radius.
    // This guarantee that it is only necessary to look for
    // neighbors at the nearest neighbors of this cell
    this->SetCellVectors(lattice, cutoff_radius);
   }
 
  /**
   * @brief Sets the lattice vectors and updates the cell vectors.
   * @param lattice_vectors The lattice vectors defining the periodic boundaries.
   * @param cutoff_radius The cutoff radius for neighbor searching.
   */ 
 void SetCellVectors(const LatticeCell lattice_, const Real cutoff_radius_);

  /**
   * @brief Sets the sites for neighbor searching.
   * @details This function modifies the Site objects by setting their Host Cell variable.
   * @param sites The collection of sites to be considered.
   * 
   * \note In this current implementation, the object sites will be modified.
   */
 void SetSites(std::vector<Site>& sites);
 
  /**
   * @brief Computes the neighbor list using the CellList algorithm.
   */ 
 void GetNeighborList();
 
  /**
   * @brief Computes the triplet identifying a host site and set at site.
   * @param site The site for which to compute the cell triplet.
   */
 void SetSiteHostCell(Site& site);
  /**
   * @brief Inserts a site into the cell grid.
   * @param site The site to be inserted.
   */
void Insert(Site & site);

  /**
   * @brief Returns a reference to the cell grid.
   * @return The cell grid.
   */
inline
Map& GetCellMap(){ return cell_map_; };

  /**
   * @brief Returns the lattice vectors.
   * @return The lattice vectors.
   */
inline
const LatticeCell& GetLattice() const { return lattice_;}

  /**
   * @brief Returns the cell vectors.
   * @return The cell vectors.
   */
inline
const LatticeCell& GetCell() const { return cell_;}

  /**
   * @brief Returns the number of sites.
   * @return The number of sites.
   */
inline
const size_t& GetNumberOfSites() const { return num_sites_;}

  /**
   * @brief Returns the cutoff radius.
   * @return The cutoff radius.
   */
inline
const Real& GetCutoffRadius() const { return cutoff_radius_;}


 private:
  Real cutoff_radius_;  ///< The cutoff radius for neighbor searching.
  LatticeCell lattice_; ///< The lattice vectors defining the periodic boundaries.
  LatticeCell cell_; ///< The cell vectors used for neighbor searching.
  Map cell_map_; ///< The grid storing the sites in cells.
  size_t num_sites_ = 0; ///< The number of sites.
};

#endif