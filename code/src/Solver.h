#ifndef SOLVER_H_
#define SOLVER_H_

#include "Mesh.h"
#include </usr/include/eigen3/Eigen/Sparse>

class Solver {
public:
    Mesh mesh;
    double DELTA;

public:
    Solver() : mesh() {}
    Solver(const Solver& solver) : mesh(solver.mesh) {}
    Solver(const Mesh& rMesh) : mesh(rMesh) {}
    ~Solver() {}

    void runGlobalOptimization();
    void runLocalOptimizations();

private:
    void getGlobalHessianAndVector(Eigen::SparseMatrix<double>& hessian, Eigen::VectorXd& vector);
    void addCoefficientsForIndex(std::vector<Eigen::Triplet<double>>& hessianCoefficients, Eigen::VectorXd& vector, int index1, int index2, float value);
    Eigen::VectorXd solveSystem(Eigen::SparseMatrix<double>& hessian, Eigen::VectorXd& vector);

    void runLocalOptimization(std::vector<int> localFaceRegion, int& counter);
    void setDelta(std::vector<std::vector<int>> localFaceRegions);
    std::vector<double> getLocalGradient(std::vector<int> localFaceRegion);
    double getLocalStepSize(std::vector<int> localFaceRegion, std::vector<double> gradient, double maxStepSize);
    double getLocalStepSizeFromEdges(std::vector<int> localFaceRegion, std::vector<double> gradient);
    std::vector<double> estimateLocalPartialDerivativesForV(size_t vId, std::vector<int> localFaceRegion);

    void updateVertices(std::vector<Vertex> newV);
    void updateVertices(Eigen::VectorXd newV);
};

#endif /* SOLVER_H_ */