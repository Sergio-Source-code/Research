#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "Solver.h"
#include <iostream>
#include <ctime>

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage: optimize <input/file/path> <output/file/path>\n";
        return -1;
    }

    std::cout << "Running optimize.cpp\n";
    clock_t time_req = clock();
    MeshFileReader reader(argv[1]);
    Mesh& mesh = (Mesh&)reader.GetMesh();
    mesh.BuildAllConnectivities();
    mesh.ExtractBoundary();
    mesh.orderNeighboringVertices();
    mesh.orderVerticesInFaces();

    // Vertex& v1 = mesh.V.at(std::stoi(argv[3]));
    // Vertex& v2 = mesh.V.at(std::stoi(argv[4]));
    // for (int i = 0; i < mesh.E.size(); i++) {
    //     Edge& e = mesh.E.at(i);
    //     if ((e.Vids.at(0) == v1.id && e.Vids.at(1) == v2.id) || (e.Vids.at(0) == v2.id && e.Vids.at(1) == v1.id)) {
    //         std::cout << "Collapsing edge " << e.id << std::endl;
    //         mesh.collapseEdge(e.id);
    //         break;
    //     }
    // }
    // mesh.expandFace(std::stoi(argv[3]));


    mesh.printQuality();
    Solver solver(mesh);
    
    // std::vector<int> verticesToIgnoreDELETE = {496, 488, 487, 505, 507, 480, 481, 467, 468};
    // for (int i = 0; i < verticesToIgnoreDELETE.size(); i++) {
    //     Vertex& v = solver.mesh.V.at(verticesToIgnoreDELETE.at(i));
    //     v.isBoundary = true;
    // }
    // solver.iterativeOptimization();


    solver.getHessian();
    solver.optimize();

    // MeshFileWriter fileWriter1(solver.mesh, "/examples/globalResult.vtk");
    // fileWriter1.WriteFile();
    // solver.pastMeshes.push_back(solver.mesh);

    // std::vector<int> freeVertexIndices = solver.mesh.selectBoundaryForIterativeOptimization();
    // solver.iterativeOptimization();
    // MeshFileWriter fileWriter2(solver.mesh, "/examples/iterativeResult.vtk");
    // fileWriter2.WriteFile();
    // solver.mesh.ExtractBoundary();
    // // solver.mesh.fixInvertedFaces();
    // solver.getHessian();
    // solver.getInitialTerm(freeVertexIndices);
    // solver.optimize();

    solver.mesh.printQuality();

    // for (int i = 0; i < solver.mesh.F.size(); i++) {
    //     Face& f = solver.mesh.F.at(i);
    //     if (i == 558 || i == 376) {
    //         std::cout << "\n " << f.id << ": ";
    //         for (int j = 0; j < f.Vids.size(); j++) {
    //             std::cout << f.Vids.at(j) << ", ";
    //         }
    //     }
    // }

    // solver.iterativeOptimizeCornerVertices();
    clock_t time_req1 = clock();
    solver.iterativeOptimizeMinSJ();
    time_req1 = clock() - time_req1;

    MeshFileWriter fileWriter(solver.mesh, argv[2]);
    fileWriter.WriteFile();
    solver.mesh.printQuality();
    std::cout << "Done\n";
	time_req = clock() - time_req;
    std::cout << "Time: " << (float)time_req/CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << "Iterative time: " << (float)time_req1/CLOCKS_PER_SEC << " seconds" << std::endl;
    
    return 0;
}
