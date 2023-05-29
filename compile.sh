
#Armadillo backend
#g++ -shared -o lib/armadillo_matvec.so lib/sparse/backends/armadillo/matvec.cpp  -larmadillo -fPIC -O3

#eigen3
#g++ -I/usr/include/eigen3 -I./lib/include/ -shared -o lib/eigen3_spalg.so lib/sparse/backends/eigen3/sparse_algebra.cpp -fPIC -lcblas -O3

g++ -pg -I/usr/include/eigen3 -Ilib/include -shared -o lib/libcell_list.so lib/src/cell_list.cpp -fPIC -O3
g++ -pg -I/usr/include/eigen3 -Ilib/include -Llib/ lib/src/neighbors.cpp -o neighbors -lcell_list -fPIC -O3 

export LD_LIBRARY_PATH=./lib:$LD_LIBRARY_PATH

