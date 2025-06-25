# Get all of this working, comment it out, then prototype the commands that follow.

mkdir .build
cd .build
rm -rf xtb
git clone https://github.com/grimme-lab/xtb

cd xtb
cmake -B build \
 -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_C_COMPILER="gcc-15" \
 -DCMAKE_Fortran_COMPILER="gfortran-15" \
 -DWITH_TBLITE=OFF \
 -DWITH_CPCMX=OFF \
 -DBLA_VENDOR=Apple
make -C build -j8
