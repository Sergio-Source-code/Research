#include "Solver.h"
#include <iostream>

using namespace std;

int getIndexOf(int num, vector<unsigned long int> vec);

void Solver::getHessian() {
    cout << "In getHessian()\n";
    standardHessianCoefficients.clear();
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex v = mesh.V.at(i);
        for (int k = 0; k < v.N_Fids.size(); k++) {
            Face& f = mesh.F.at(v.N_Fids.at(k));
            for (int j = 0; j < f.Vids.size(); j++) {
                Vertex u = mesh.V.at(f.Vids.at(j));
                int vIndex = getIndexOf(v.id, f.Vids);
                int uIndex = getIndexOf(u.id, f.Vids);
                int difference = (uIndex - vIndex + 4) % 4;
                float value = difference % 2 == 0 ? 1 : -1;
                if (difference == 0)
                    value = 3;
                else if (difference == 1)
                    value = -2;
                else if (difference == 2)
                    value = 1;
                else
                    value = -2;
                addCoefficientsForIndex(v.id, u.id, value);
            }
        }
    }
    standardHessian.setFromTriplets(standardHessianCoefficients.begin(), standardHessianCoefficients.end());
}

void Solver::addCoefficientsForIndex(int index1, int index2, float value) {
    Vertex v = mesh.V.at(index1);
    Vertex u = mesh.V.at(index2);
    if (v.isBoundary) {
        return;
    }
    if (u.isBoundary) {
        standardVector(2 * index1) += value * u.x;
        standardVector(2 * index1 + 1) += value * u.y;
    } else {
        standardHessianCoefficients.push_back(Eigen::Triplet<double>(2 * index1, 2 * index2, value));
        standardHessianCoefficients.push_back(Eigen::Triplet<double>(2 * index1 + 1, 2 * index2 + 1, value));
    }
}

void Solver::optimize() {
    cout << "In optimize()\n";

    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
    Eigen::SparseMatrix<double> A = standardHessian;
    Eigen::VectorXd b = standardVector;

    solver.analyzePattern(A);
    solver.factorize(A);

    Eigen::VectorXd x(b.size());
    x = solver.solve(-b);

    for (int i = 0; i < mesh.V.size(); i++) {
        if (mesh.V.at(i).isBoundary)
            continue;
        mesh.V.at(i).x = x(2 * i);
        mesh.V.at(i).y = x(2 * i + 1);
    }
}

int getIndexOf(int num, vector<unsigned long int> vec) {
    for (int i = 0; i < vec.size(); i++) {
        if (vec.at(i) == num)
            return i;
    }
    return -1;
}