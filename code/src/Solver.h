#ifndef SOLVER_H_
#define SOLVER_H_

#include "Mesh.h"
#include </usr/include/eigen3/Eigen/Sparse>

class Solver {
public:
    Solver()
    : mesh()
    {}
    Solver(const Solver& solver)
    : mesh(solver.mesh)
    {}
    Solver(const Mesh& rMesh)
    : mesh(rMesh)
    {
        int size = 2 * mesh.V.size();
        standardHessianCoefficients.clear();
        standardHessian = Eigen::SparseMatrix<double>(size, size);
        standardVector = Eigen::VectorXd::Zero(size);
    }
    ~Solver()
    {}
    void getHessian();
    void getHessianForEdgeTerms();
    void getInitialTerm(std::vector<int> vIds);
    void getInitialTermWithConstantWeight(std::vector<int> vIds, float weight);
    void getInitialTermWithDepthWeight(std::vector<int> vIds);
    void getInitialTermWithLength(std::vector<int> vIds);
    void addInitialWeightForV(int vId, float weight);

    void addAverageTermForVertex(int vId);
    void addAverageTermForFace(int fId);
    void optimize();
    void iterativeOptimization();
    void iterativeOptimizeCornerVertices();
    void iterativeOptimizeMinSJ();
    std::vector<double> estimateGradientForFV(size_t fId, size_t vId, float initialSJ);
    std::vector<double> estimateMinSJGradient(size_t vId, std::vector<int> localFaceRegions);

private:
    void addCoefficientsForIndex(int index1, int index2, float value);
    std::vector<Vertex> iterativeSmoothInnerVertices();
    std::vector<Vertex> iterativeSmoothBoundaryVertices();
    void updateVertices(std::vector<Vertex> newV);

public:
    Mesh mesh;
    std::vector<Mesh> pastMeshes;
    std::vector<Eigen::Triplet<double>> standardHessianCoefficients;
    Eigen::SparseMatrix<double> standardHessian;
    Eigen::VectorXd standardVector;
};

#endif /* SOLVER_H_ */