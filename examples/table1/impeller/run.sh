
if [ ! -f "LMO" ]; then
    mkdir LMO
fi
cd LMO
cp ../input.vtk .
# This is very important
export PATH=../../../../code/bin:$PATH

# Preprocessing
MeshScale input.vtk scaled.vtk
# Global Untangling
MeshOpt  input=scaled.vtk alpha=2000 beta=2000 gamma=1.0 anisotropy=0.02 allowBigStep=true minScaledJacobian=0.66 cosangle=0.999  stepSize=0.999 useProjection=false iters=8
# Inversion-free deformation
SlimHex input=opt.vtk target=scaled.vtk result=slimopt.vtk iters=20
# MSJ Improving
LocalMeshOpt  input=slimopt.vtk alpha=500 beta=1250 gamma=1.0 anisotropy=0.1 allowBigStep=true minScaledJacobian=0.2 cosangle=0.9994807753   useSmallBlock=true blockSize=100 localIters=20 iters=20 targetSurface=scaled.vtk projectToTargetSurface=true
cp BestLocalOpt.vtk final1.vtk
# Inversion-free deformation
time SlimHex input=final1.vtk target=scaled.vtk result=final.vtk iters=5
# Post Processing
SmoothMap orig=scaled.vtk input=final.vtk output=smooth.vtk iters=20 preserveQuality=true preserveSharpFeature=true treatSharpFeatureAsRegular=false 2>err.log
MeshScale smooth.vtk result.vtk orig=input.vtk
extract_surface result.vtk result.surface.off 1
extract_surface input.vtk input.surface.off 1
# Hausdorff Distance
metro input.surface.off result.surface.off
# Scaled Jacobian Measurement
MeshQuality result.vtk

