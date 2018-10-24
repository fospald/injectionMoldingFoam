
rm polyMesh
ln -s polyMesh.$1 polyMesh
cd ..
decomposePar -force
