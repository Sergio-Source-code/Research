
if [ ! -f "LMO" ]; then
    mkdir LMO
fi
cd LMO
cp ../input.vtk .
# This is very important
export PATH=../../../../code/bin:$PATH

# Preprocessing
MeshScale input.vtk scaled.vtk
SmoothMap orig=scaled.vtk input=scaled.vtk output=smooth.vtk iters=10 preserveSharpFeature=true treatSharpFeatureAsRegular=true cosangle=0.965925826
# Untangling
LocalMeshOpt  input=smooth.vtk alpha=100 beta=100 gamma=1.0 anisotropy=0.4 allowBigStep=true minScaledJacobian=0.0 useSmallBlock=true blockSize=200 localIters=20 iters=20 targetSurface=smooth.vtk projectToTargetSurface=true
cp BestLocalOpt.vtk untangled.vtk
# Inversion-free deformation
SlimHex input=untangled.vtk target=smooth.vtk result=slim.vtk iters=50
# MSJ Improving
LocalMeshOpt input=slim.vtk alpha=500 beta=500 gamma=1.0 anisotropy=0.2 allowBigStep=false minScaledJacobian=0.16 useSmallBlock=true blockSize=100 localIters=20 iters=20 targetSurface=smooth.vtk projectToTargetSurface=true  cosangle=0.965925826
cp BestLocalOpt.vtk final.vtk
# Post Processing
SmoothMap orig=scaled.vtk input=final.vtk output=output.vtk iters=20 preserveQuality=true preserveSharpFeature=true treatSharpFeatureAsRegular=false 2>err.log
MeshScale output.vtk result.vtk orig=input.vtk
extract_surface result.vtk result.surface.off 1
extract_surface input.vtk input.surface.off 1
# Hausdorff Distance
metro input.surface.off result.surface.off
# Scaled Jacobian Measurement
MeshQuality result.vtk

