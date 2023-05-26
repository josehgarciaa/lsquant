
#ifndef SPATIAL
#define SPATIAL

#include <iostream>   // for std::cout, std::endl, etc.
#include <vector>     // for std::vector
#include <unordered_map> // for std::unordered_map
#include <cstddef>    // for size_t
#include <Eigen/Dense> // for Eigen::Vector3d and Eigen::Vector3i

typedef Eigen::Matrix3d Matrix3D;
typedef Eigen::Vector3d Vector3D; //.norm() method
typedef Eigen::Vector3d Position; //.norm() method
typedef Eigen::Vector3i Triplet;
typedef double Float;
typedef double Real;
typedef int Integer;

//default tag ? 
std::string TripletToString(const Triplet& arr) {
    std::ostringstream oss;
    for (size_t i = 0; i < 3; ++i) 
        oss << arr[i]<<" ";
    return oss.str();
}


class LatticeCell {
private:
    Matrix3D lattice_vectors_;

public:
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

    Real operator()(const int i, const int j) const { return lattice_vectors_(i,j);}


    Matrix3D GetVectors() const {
        return lattice_vectors_;
    }
    Vector3D GetVector(const int i) const {
        return lattice_vectors_.col(i);
    }


    std::array<Real, 3> Norms() const 
    {
        std::array<Real, 3> norms;
        for (size_t i = 0; i < 3; ++i) {
            norms[i] = lattice_vectors_.col(i).norm();
        }
        return norms;
    }

    Vector3D FractionalCoords(const Vector3D& X) const {
        return lattice_vectors_.inverse()*X;
    }
};


class Site {
 public:

 // Default constructor 
 explicit Site(Vector3D X = Position()) : position_(X){};

// Initialize by 
 explicit Site( const double x, 
                const double y, 
                const double z) : position_(Vector3D(x,y,z)){ tag_="X";};

 // Copy constructor
 explicit Site(const Site& A) : position_(A.position_), tag_(A.tag_){};
 
 // Move constructor
 Site(Site&& A) noexcept : position_(std::move(A.position_)), tag_(std::move(A.tag_)) {}

 // Default destructor
 ~Site() = default;


// Add const for function and arguments, follow naming conventions and use getters.
Real Distance(const Site& A) const {
    Position dx = A.GetPosition();
    dx -= GetPosition();
    return dx.norm();
}

// Add const for function and arguments, follow naming conventions and use getters.
Position Displacement(const Site& A) const {
    Position dx = A.GetPosition();
    dx -= GetPosition();
    return dx;
}

// Getters/Setters with CamelCase.
inline Vector3D& GetPosition()  { return position_; }
inline Vector3D GetPosition() const { return position_; }
inline void SetPosition(Vector3D X) { position_ = X; }
inline std::string GetTag() const { return tag_; }
inline void SetTag(std::string Id) { tag_ = Id; }
inline void SetHostCell(Triplet host_cell) { host_cell_ = host_cell; }
inline Triplet GetHostCell() const { return host_cell_; }
inline size_t GetIndex() const { return id_; }

private:
    std::string tag_;
    Triplet host_cell_; 
    Vector3D position_;
    size_t id_;
};



typedef std::unordered_map<size_t,std::vector<Site> > Map;
class CellList {
 public:
 
 //explicit 
 CellList(const LatticeCell lattice_vectors,
          const double cutoff_radius)
    {
    //We will create cells that fully contain 
    //all points within a given raiuds. 
    //Therefore, we will only need to look for
    //nearest cell neighbors.  
    this->SetCellVectors(lattice_vectors, cutoff_radius);
    }
 
 void SetCellVectors(const LatticeCell lattice_vectors, const Real cutoff_radius_);

 void SetSites(std::vector<Site>& sites);
 
 void GetNeighborList();
 
 Triplet CellTriplet(const Site& site);

void Insert(Site & site);

inline
std::unordered_map<std::string,std::vector<Site> >& GetGrid(){ return grid_; };

inline
LatticeCell GetLatticeVectors() const { return lattice_vectors_;}

inline
LatticeCell GetCellVectors() const { return cell_vectors_;}

inline
size_t GetNumberOfSites(){ return num_sites_;}

inline
Real GetCutoffRadius(){ return cutoff_radius_;}


 private:
  Real cutoff_radius_;
  LatticeCell lattice_vectors_;
  LatticeCell cell_vectors_;
  std::unordered_map<std::string, std::vector<Site> > grid_;
  size_t num_sites_ = 0;
};

#endif
