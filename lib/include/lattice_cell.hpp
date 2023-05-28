#ifndef LATTICE_CELL_HPP
#define LATTICE_CELL_HPP

#include "typedef.hpp" //Define the types: Matrix3D, Vector3D, Triplet, Float, Real, Integer;

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
    
private:
    Matrix3D lattice_vectors_; ///< The lattice vectors of the cell.
};

#endif