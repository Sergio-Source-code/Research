
if [ ! -f "LMO" ]; then
    mkdir LMO
fi
cd LMO
cp ../input.vtk .
# This is very important
export PATH=../../../../code/bin:$PATH

# Preprocessing
MeshScale input.vtk scaled.vtk scale=1.0
# Global MSJ Improving
MeshOpt  input=scaled.vtk alpha=20 beta=20 gamma=1.0 anisotropy=1.2 allowBigStep=true minScaledJacobian=0.66 cosangle=0.999  stepSize=0.999 useProjection=false
cp opt.vtk final.vtk
SmoothMap orig=scaled.vtk input=final.vtk output=smooth.vtk iters=20 preserveQuality=true preserveSharpFeature=true treatSharpFeatureAsRegular=false 2>err.log
MeshScale smooth.vtk result.vtk orig=input.vtk scale=1.0
cp mapback.vtk result.vtk
extract_surface result.vtk result.surface.off 1
extract_surface input.vtk input.surface.off 1
# Hausdorff Distance
metro input.surface.off result.surface.off
# Scaled Jacobian Measurement
MeshQuality result.vtk

