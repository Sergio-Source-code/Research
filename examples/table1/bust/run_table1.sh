
if [ ! -f "LMO" ]; then
    mkdir LMO
fi
cd LMO
cp ../input.vtk .
# This is very important
export PATH=../../../../code/bin:$PATH

# Preprocessing
MeshScale input.vtk scaled.vtk
# Untangling
LocalMeshOpt  input=scaled.vtk alpha=100 beta=100 gamma=1.0 anisotropy=0.6 allowBigStep=true minScaledJacobian=0.0
cp BestLocalOpt.vtk untangled.vtk
# Inversion-free deformation
SlimHex input=untangled.vtk target=scaled.vtk result=slim.vtk iters=20
# MSJ Improving
LocalMeshOpt  input=slim.vtk alpha=1000 beta=1000 gamma=1.0 anisotropy=0.6 allowBigStep=true minScaledJacobian=0.18  cosangle=0.935925826  useSmallBlock=true blockSize=100 localIters=20 iters=20 targetSurface=scaled.vtk projectToTargetSurface=true
cp BestLocalOpt.vtk MinSJ_0.2.vtk
# Post Processing
SmoothMap orig=scaled.vtk input=MinSJ_0.2.vtk output=output.vtk iters=10  preserveQuality=true preserveSharpFeature=true treatSharpFeatureAsRegular=false cosangle=0.984807753 2>err.log
MeshScale output.vtk result.vtk orig=input.vtk
extract_surface result.vtk result.surface.off 1
extract_surface input.vtk input.surface.off 1
# Hausdorff Distance
metro input.surface.off result.surface.off
# Scaled Jacobian Measurement
MeshQuality result.vtk

