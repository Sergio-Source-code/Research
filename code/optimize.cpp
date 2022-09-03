#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "Solver.h"
#include <iostream>

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage: optimize <input/file/path> <output/file/path>\n";
        return -1;
    }

    std::cout << "Running optimize.cpp\n";
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();

    Solver solver(mesh);
    solver.getHessian();
    solver.optimize();

    MeshFileWriter fileWriter(solver.mesh, argv[2]);
    fileWriter.WriteFile();
    std::cout << "Done\n";
    
    return 0;
}
