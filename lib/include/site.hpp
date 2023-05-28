#ifndef SITE_HPP
#define SITE_HPP

#include "typedef.hpp" //Define the types: Matrix3D, Vector3D, Triplet, Float, Real, Integer;
#include <iostream>   // for std::cout, std::endl, etc.
#include <cstddef>    // for size_t

/**
 * @brief Class representing a site in the system.
 */
class Site {
 public:

  /**
   * @brief Default constructor.
   * @param X The position vector of the site (default: zero vector).
   */
 explicit Site(Vector3D X = Vector3D()) : position_(X){};


  /**
   * @brief Initializes a site with the given coordinates.
   * @param x The x-coordinate of the site.
   * @param y The y-coordinate of the site.
   * @param z The z-coordinate of the site.
   */ 
 explicit Site( const Real x, 
                const Real y, 
                const Real z) : position_(Vector3D(x,y,z)){ tag_="X";};

  /**
   * @brief Default destructor.
   */
 ~Site() = default;


  /**
   * @brief Computes the distance between this site and another site.
   * @param A The other site.
   * @return The distance between the sites.
   */
inline Real Distance(const Site& A) const {
    const auto D = A.GetPosition() - this->GetPosition() ;
    return D.norm();
}


// Getters/Setters with CamelCase.
  /**
   * @brief Returns the position vector of the site.
   * @return The position vector .
   */
inline const Vector3D& GetPosition() const { return position_; }



  /**
   * @brief Sets the position vector of the site.
   * @param X The new position vector.
   */
inline void SetPosition(Vector3D X) { position_ = X; }


  /**
   * @brief Returns the tag of the site.
   * @return The tag.
   */
inline const std::string& GetTag() const { return tag_; }

  /**
   * @brief Sets the tag of the site.
   * @param Id The new tag.
   */
inline void SetTag(std::string Id) { tag_ = Id; }

  /**
   * @brief Sets the host cell of the site.
   * @param host_cell The host cell.
   */
inline void SetHostCell(const Triplet host_cell) { host_cell_ = host_cell; }

  /**
   * @brief Returns the host cell of the site.
   * @return The host cell.
   */
inline const Triplet& GetHostCell() const { return host_cell_; }

  /**
   * @brief Returns the index of the site.
   * @return The index.
   */
inline const size_t& GetIndex() const { return id_; }

private:
    std::string tag_;
    Triplet host_cell_; 
    Vector3D position_;
    size_t id_;
};

#endif