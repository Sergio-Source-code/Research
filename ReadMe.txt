# 1. Install Docker

# 2. Go to this directory and run:
docker build . -t quadoptimization
to build the container image

# 3. To run the code with an input file, run:
docker run -v <directory_path_on_host>:<directory_path_on_container> quadoptimization /code/bin/optimize <input_file_path_in_container> <output_file_path_in_container>

For example, to run the quadoptimization container image with the file C:\Users\sergi\Research\examples\2d\simplified_patch7_2D.vtk, I would mount
the host directory: C:\Users\sergi\Research\examples to the container directory: /examples, and then pass the container file paths as
arguments to /code/bin/optimize. Running the command below, the output file will be in /examples/output.vtk in the container
which is C:\Users\sergi\Research\examples\output.vtk on the host.

docker run -v C:\Users\sergi\Research\examples:/examples quadoptimization /code/bin/optimize /examples/2d/simplified_patch7_2D.vtk /examples/output.vtk