
if [ ! -f "LMO" ]; then
    mkdir LMO
fi
cd LMO
cp ../input.vtk .
# This is very important
export PATH=../../../../code/bin:$PATH

# Preprocessing
MeshScale input.vtk scaled.vtk
SmoothMap orig=scaled.vtk input=scaled.vtk output=out.vtk iters=10 preserveSharpFeature=true treatSharpFeatureAsRegular=false cosangle=0.965925826
# Untangling
LocalMeshOpt  input=out.vtk alpha=100 beta=100 gamma=1.0 anisotropy=0.6 allowBigStep=true minScaledJacobian=0.0 cosangle=0.965925826  useSmallBlock=true blockSize=100 localIters=20 iters=20  targetSurface=out.vtk projectToTargetSurface=true useSmallBlock=true blockSize=100
cp BestLocalOpt.vtk untangled.vtk
# Inversion-free deformation
SlimHex input=untangled.vtk target=out.vtk result=slim.vtk iters=20
# MSJ Improving
LocalMeshOpt  input=slim.vtk alpha=500 beta=500 gamma=1.0 anisotropy=0.4 allowBigStep=true minScaledJacobian=0.1 cosangle=0.965925826  useSmallBlock=true blockSize=100 localIters=20 iters=20 targetSurface=out.vtk projectToTargetSurface=true useSmallBlock=true blockSize=100
cp BestLocalOpt.vtk MinSJ_0.1.vtk
LocalMeshOpt  input=MinSJ_0.1.vtk alpha=500 beta=500 gamma=1.0 anisotropy=0.4 allowBigStep=true minScaledJacobian=0.2 cosangle=0.965925826  useSmallBlock=true blockSize=100 localIters=20 iters=20 targetSurface=out.vtk projectToTargetSurface=true useSmallBlock=true blockSize=100
cp BestLocalOpt.vtk MinSJ_0.2.vtk
LocalMeshOpt  input=MinSJ_0.2.vtk alpha=500 beta=500 gamma=1.0 anisotropy=0.4 allowBigStep=true minScaledJacobian=0.3 cosangle=0.965925826  useSmallBlock=true blockSize=100 localIters=20 iters=20 targetSurface=out.vtk projectToTargetSurface=true useSmallBlock=true blockSize=100
cp BestLocalOpt.vtk MinSJ_0.3.vtk
LocalMeshOpt  input=MinSJ_0.3.vtk alpha=500 beta=500 gamma=1.0 anisotropy=0.4 allowBigStep=true minScaledJacobian=0.4 cosangle=0.965925826  useSmallBlock=true blockSize=100 localIters=20 iters=20 targetSurface=out.vtk projectToTargetSurface=true useSmallBlock=true blockSize=100
cp BestLocalOpt.vtk MinSJ_0.4.vtk
LocalMeshOpt  input=MinSJ_0.4.vtk alpha=500 beta=500 gamma=1.0 anisotropy=0.4 allowBigStep=true minScaledJacobian=0.5 cosangle=0.965925826  useSmallBlock=true blockSize=100 localIters=20 iters=20 targetSurface=out.vtk projectToTargetSurface=true useSmallBlock=true blockSize=100
cp BestLocalOpt.vtk MinSJ_0.5.vtk
LocalMeshOpt  input=MinSJ_0.5.vtk alpha=500 beta=500 gamma=1.0 anisotropy=0.4 allowBigStep=true minScaledJacobian=0.53 cosangle=0.965925826  useSmallBlock=true blockSize=100 localIters=20 iters=20 targetSurface=out.vtk projectToTargetSurface=true useSmallBlock=true blockSize=100
cp BestLocalOpt.vtk MinSJ_0.53.vtk
LocalMeshOpt  input=MinSJ_0.53.vtk alpha=500 beta=500 gamma=1.0 anisotropy=0.4 allowBigStep=true minScaledJacobian=0.58 cosangle=0.965925826  useSmallBlock=true blockSize=100 localIters=20 iters=20 targetSurface=out.vtk projectToTargetSurface=true useSmallBlock=true blockSize=100
cp BestLocalOpt.vtk final.vtk
SmoothMap orig=scaled.vtk input=final.vtk output=output.vtk iters=10 preserveQuality=true preserveSharpFeature=true treatSharpFeatureAsRegular=false 2>err.log cosangle=0.965925826
MeshScale output.vtk result.vtk orig=input.vtk cosangle=0.965925826
extract_surface result.vtk result.surface.off 1
extract_surface input.vtk input.surface.off 1
metro input.surface.off result.surface.off
MeshQuality result.vtk

