
#Armadillo backend
#g++ -shared -o lib/armadillo_matvec.so lib/sparse/backends/armadillo/matvec.cpp  -larmadillo -fPIC

#eigen3
#g++ -I/usr/include/eigen3 -I./lib/include/ -shared -o lib/eigen3_spalg.so lib/sparse/backends/eigen3/sparse_algebra.cpp -fPIC -lcblas

g++ -I/usr/include/eigen3 -I./lib/include lib/src/neighbors.cpp -fPIC -O3 -o neighbors
