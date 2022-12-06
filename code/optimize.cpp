#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "Solver.h"
#include <iostream>
#include <ctime>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: optimize <input/file/path> <output/file/path>\n";
        return -1;
    }

    std::cout << "Running optimize.cpp\n";
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();

    if (mesh.m_cellType != QUAD) {
        std::cout << "Mesh must have QUAD cell type, is " << mesh.m_cellType << std::endl;
        return -1;
    } else if (mesh.V.size() == 0) {
        std::cout << "Error reading mesh file, 0 vertices found" << std::endl;
        return -1;
    }

    std::cout << "-------------------------" << std::endl;
    std::cout << "Initial quality" << std::endl;
    mesh.process();
    mesh.printQuality();
    Solver solver(mesh);

    // std::cout << "Running global optimization" << std::endl;
    clock_t globalTime = clock();
    solver.runGlobalOptimization();
    globalTime = clock() - globalTime;

    std::cout << "Intermediate quality" << std::endl;
    solver.mesh.printQuality();
    std::string fileName2 = std::string("/examples/globalResult.vtk");
    MeshFileWriter fw2(solver.mesh, fileName2.c_str());
    fw2.WriteFile();

    // std::cout << "Running local optimizations" << std::endl;
    clock_t localTime = clock();
    solver.runLocalOptimizations();
    localTime = clock() - localTime;

    std::cout << "Final quality" << std::endl;
    solver.mesh.printQuality();
    MeshFileWriter fileWriter(solver.mesh, argv[2]);
    fileWriter.WriteFile();

    std::cout << "Global time: " << (float)globalTime/CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "Local  time: " << (float)localTime/CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "Total  time: " << ((float)globalTime + localTime)/CLOCKS_PER_SEC  << " seconds" << std::endl;
    std::cout << "Done" << std::endl;
    std::cout << "-------------------------" << std::endl;
    return 0;
}
