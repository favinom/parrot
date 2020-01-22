cd ~/projects
git clone --recurse-submodules https://bitbucket.org/zulianp/utopia.git

export LIBMESH_DIR=MOOSE_DIR/libmesh/installed/

cd utopia
mkdir build
export UTOPIA_DIR=~/projects/utopia/build

cd utopia
mkdir bin
cd bin
cmake .. -DCMAKE_INSTALL_PREFIX=$UTOPIA_DIR
make -j 4 && make install

cd ../../utopia_fe/

mkdir bin
cd bin
cmake .. -DUTOPIA_DIR=$UTOPIA_DIR -DLIBMESH_DIR=$LIBMESH_DIR -DCMAKE_INSTALL_PREFIX=$UTOPIA_DIR -DMOONOLITH_INSTALL_PREFIX=$UTOPIA_DIR

make -j 4 && make install
