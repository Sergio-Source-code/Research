#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "Solver.h"
#include <iostream>
#include <ctime>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "Usage: optimize <input/file/path> <output/file/path>\n";
        return -1;
    }

    std::cout << "Running optimize.cpp\n";
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();

    mesh.process();
    mesh.printQuality();
    Solver solver(mesh);

    clock_t globalTime = clock();
    solver.runGlobalOptimization();
    globalTime = clock() - globalTime;

    solver.mesh.printQuality();

    clock_t localTime = clock();
    solver.runLocalOptimizations();
    localTime = clock() - localTime;

    solver.mesh.printQuality();

    MeshFileWriter fileWriter(solver.mesh, argv[2]);
    fileWriter.WriteFile();

    std::cout << "Done\n";
    std::cout << "Global time: " << (float)globalTime/CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "Local  time: " << (float)localTime/CLOCKS_PER_SEC << " seconds" << std::endl;
    
    return 0;
}
