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
    void optimize();

private:
    void addCoefficientsForIndex(int index1, int index2, float value);

public:
    Mesh mesh;
    std::vector<Eigen::Triplet<double>> standardHessianCoefficients;
    Eigen::SparseMatrix<double> standardHessian;
    Eigen::VectorXd standardVector;
};

#endif /* SOLVER_H_ */