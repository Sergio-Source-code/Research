#include "Solver.h"
#include "MeshFileWriter.h"
#include <iostream>

#include <iostream>
#include <iomanip>

using namespace std;

void Solver::runGlobalOptimization() {
    int n = 2 * mesh.V.size();
    Eigen::SparseMatrix<double> hessian = Eigen::SparseMatrix<double>(n, n);
    Eigen::VectorXd vector = Eigen::VectorXd::Zero(n);

    getGlobalHessianAndVector(hessian, vector);
    Eigen::VectorXd newVertexPositions = solveSystem(hessian, vector);
    updateVertices(newVertexPositions);
}

void Solver::getGlobalHessianAndVector(Eigen::SparseMatrix<double>& hessian, Eigen::VectorXd& vector) {
    std::vector<Eigen::Triplet<double>> hessianCoefficients;

    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex v = mesh.V.at(i);

        for (int k = 0; k < v.N_Fids.size(); k++) {
            Face& f = mesh.F.at(v.N_Fids.at(k));

            for (int j = 0; j < f.Vids.size(); j++) {
                Vertex u = mesh.V.at(f.Vids.at(j));

                int vIndex = mesh.getLocalIndexOf(v.id, f.Vids);
                int uIndex = mesh.getLocalIndexOf(u.id, f.Vids);
                int difference = (uIndex - vIndex + 4) % 4;
                float value = difference % 2 == 0 ? 1 : -1;
                
                // float a = 1;
                if (difference == 0)
                    value = 3;// * a * a - 2 * a + 1;
                else if (difference == 1)
                    value = -2;// * a * a;
                else if (difference == 2)
                    value = 1;// * a * a + 2 * a - 1;
                else
                    value = -2;// * a * a;

                addCoefficientsForIndex(hessianCoefficients, vector, v.id, u.id, value);
            }
        }
    }

    hessian.setFromTriplets(hessianCoefficients.begin(), hessianCoefficients.end());
}

void Solver::addCoefficientsForIndex(std::vector<Eigen::Triplet<double>>& hessianCoefficients, Eigen::VectorXd& vector, int index1, int index2, float value) {
    Vertex v = mesh.V.at(index1);
    Vertex u = mesh.V.at(index2);
    if (v.isBoundary) {
        return;
    }
    if (u.isBoundary) {
        vector(2 * index1) += value * u.x;
        vector(2 * index1 + 1) += value * u.y;
    } else {
        hessianCoefficients.push_back(Eigen::Triplet<double>(2 * index1, 2 * index2, value));
        hessianCoefficients.push_back(Eigen::Triplet<double>(2 * index1 + 1, 2 * index2 + 1, value));
    }
}

Eigen::VectorXd Solver::solveSystem(Eigen::SparseMatrix<double>& A, Eigen::VectorXd& b) {
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> systemSolver;
    systemSolver.analyzePattern(A);
    systemSolver.factorize(A);

    return systemSolver.solve(-b);
}

void Solver::updateVertices(std::vector<Vertex> newV) {
    if (mesh.V.size() != newV.size()) {
        std::cerr << "ERROR: mesh.V.size() != newV.size() in updateVertices()\n";
        std::exit(1);
    }
    for (int i = 0; i < mesh.V.size(); i++) {
        mesh.V.at(i) = newV.at(i);
    }
}

void Solver::updateVertices(Eigen::VectorXd newV) {
    if (newV.size() != 2 * mesh.V.size()) {
        std::cerr << "ERROR: newV.size() != 2 * mesh.V.size() in updateVertices()\n";
        std::exit(1);
    }
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        if (!v.isBoundary) {
            v.x = newV(2 * i);
            v.y = newV(2 * i + 1);
        }
    }
}

void Solver::runLocalOptimizations() {
    std::vector<std::vector<int>> localFaceRegions = mesh.getLocalFaceRegions();
    setDelta(localFaceRegions);
    int counter = 0;
    // std::cout << "Number of local face regions: " << localFaceRegions.size() << std::endl;
    for (int i = 0; i < localFaceRegions.size(); i++) {
        std::vector<int> localFaceRegion = localFaceRegions.at(i);
        // std::cout << "Optimizing local face region " << i << " with " << mesh.getNumOfInvertedElements(localFaceRegion) << " inversions" << std::endl;
        runLocalOptimization(localFaceRegion, counter);
        // std::cout << "done in " << counter << " iterations" << std::endl;
    }
    // std::cout << "All done in " << counter << " iterations\n";
    // std::cout << "-------------------------\n";
}

void Solver::setDelta(std::vector<std::vector<int>> localFaceRegions) {
    double minimumEdgeLength = 1000000000;
    for (int i = 0; i < localFaceRegions.size(); i++) {
        std::vector<int> localFaceRegion = localFaceRegions.at(i);

        for (int j = 0; j < localFaceRegion.size(); j++) {
            Face& f = mesh.F.at(localFaceRegion.at(j));

            for (int k = 0; k < f.Eids.size(); k++) {
                Edge& e = mesh.E.at(f.Eids.at(k));

                double currentEdgeLength = glm::length(mesh.V.at(e.Vids.at(1)) - mesh.V.at(e.Vids.at(0)));
                if (currentEdgeLength < minimumEdgeLength) {
                    minimumEdgeLength = currentEdgeLength;
                }
            }
        }
    }
    DELTA = minimumEdgeLength / 10;
}

void Solver::runLocalOptimization(std::vector<int> localFaceRegion, int& counter) {
    double stepSize;
    double maxStepSize;
    bool jumpedOnLast = false;
    double lastEnergy = 0.0;
    std::vector<int> innerVertices = mesh.getInnerVerticesOfFaces(localFaceRegion);

    for (int i = 0; i < 100; i++) {
        std::vector<double> gradient = getLocalGradient(localFaceRegion, innerVertices);

        double tempCurrentEnergy = mesh.getIterativeEnergyOfRegion(localFaceRegion);

        if (i == 0) {
            stepSize = getLocalStepSizeFromEdges(localFaceRegion, gradient);
            maxStepSize = stepSize;
        }
        stepSize = getLocalStepSize(localFaceRegion, gradient, maxStepSize);
        bool isReturnedStepSizeZero = stepSize == 0;
        if (stepSize == 0) {
            if (mesh.getMinimumScaledJacobian(localFaceRegion) <= 0) {
                // std::cout << "jump at iteration " << counter << std::endl;
                stepSize = maxStepSize;
            } else {
                // std::cout << "Iterative optimization done at iteration " << counter << " because stepSize == 0 and no inversions" << std::endl;
                break;
            }
        }

        for (int j = 0; j < mesh.V.size(); j++) {
            Vertex& v = mesh.V.at(j);
            v.x -= stepSize * gradient.at(2 * v.id);
            v.y -= stepSize * gradient.at(2 * v.id + 1);
        }

        double tempNextEnergy = mesh.getIterativeEnergyOfRegion(localFaceRegion);
        // if (jumpedOnLast && tempNextEnergy <= lastEnergy) {
        //     std::cout << "Iterative optimization done at iteration " << counter << " because jump start failed" << std::endl;
        //     break;
        // }
        // jumpedOnLast = isReturnedStepSizeZero;// (stepSize == 0) && mesh.getMinimumScaledJacobian(localFaceRegion) <= 0;
        // if (jumpedOnLast)
        //     lastEnergy = tempNextEnergy;

        std::string fileName = std::string("/examples/_iterativeMinSJ_") + std::to_string(counter) + ".vtk";
        MeshFileWriter fw(mesh, fileName.c_str());
        // fw.WriteFile();
        counter++;
    }
}

double Solver::getLocalStepSizeFromEdges(std::vector<int> localFaceRegion, std::vector<double> gradient) {
    // make sure that max( length(v - newV) / length(min(v.N_E)) ) <= c
    // where newV = v - stepSize * gradient
    // c = 0.2
    // means that vertex moving most relative to its minimum edge is moving at c length of minimum edge
    double stepSize = 100000000;
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex v = mesh.V.at(i);
        Vertex newV = mesh.V.at(i);
        double minEdgeLength = mesh.getMinimumEdgeLength(v);

        // select s such that length(v - newV) = sqrt( (s*g(x))^2 + (sg(y))^2 ) = target length
        // => s^2 = targetLength^2 / (g(x)^2 + g(y)^2)
        double targetLength = 0.2 * minEdgeLength;
        double s = std::sqrt( (targetLength * targetLength) / (std::pow(gradient.at(2 * v.id), 2) + std::pow(gradient.at(2 * v.id + 1), 2)) );
        if (s < stepSize) {
            // std::cout << "step size determined by vertex " << v.id << std::endl;
            stepSize = s;
        }
    }
    return stepSize;
}

std::vector<double> Solver::getLocalGradient(std::vector<int> localFaceRegion) {
    std::vector<double> gradient(2 * mesh.V.size());
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex& v = mesh.V.at(i);
        if (!v.isBoundary) {
            std::vector<double> estimatedPartialDerivatives = estimateLocalPartialDerivativesForV(v.id, localFaceRegion);
            gradient.at(2 * v.id) += estimatedPartialDerivatives.at(0);
            gradient.at(2 * v.id + 1) += estimatedPartialDerivatives.at(1);
        }
    }
    return gradient;
}

std::vector<double> Solver::getLocalGradient(std::vector<int> localFaceRegion, std::vector<int> innerVertices) {
    std::vector<double> gradient(2 * mesh.V.size());
    for (int i = 0; i < innerVertices.size(); i++) {
        Vertex& v = mesh.V.at(innerVertices.at(i));
        if (!v.isBoundary) {
            std::vector<double> estimatedPartialDerivatives = estimateLocalPartialDerivativesForV(v.id, localFaceRegion);
            gradient.at(2 * v.id) += estimatedPartialDerivatives.at(0);
            gradient.at(2 * v.id + 1) += estimatedPartialDerivatives.at(1);
        }
    }
    return gradient;
}

double Solver::getLocalStepSize(std::vector<int> localFaceRegion, std::vector<double> gradient, double maxStepSize) {
    double t = 0;
    for (int i = 0; i < gradient.size(); i++) {
        t += std::pow(gradient.at(i), 2);
    }
    if (t <= 0) {
        return 0;
    }
    t *= -1;
    double a = maxStepSize;
    double currentEnergy = mesh.getIterativeEnergyOfRegion(localFaceRegion);
    double nextEnergy;
    for (int j = 0; j < 10; j++) {
        for (int k = 0; k < mesh.V.size(); k++) {
            Vertex& v = mesh.V.at(k);
            v.x -= a * gradient.at(2 * k);
            v.y -= a * gradient.at(2 * k + 1);
        }
        nextEnergy = mesh.getIterativeEnergyOfRegion(localFaceRegion);
        for (int k = 0; k < mesh.V.size(); k++) {
            Vertex& v = mesh.V.at(k);
            v.x += a * gradient.at(2 * k);
            v.y += a * gradient.at(2 * k + 1);
        }
        if (a == maxStepSize && currentEnergy <= nextEnergy)
            return a;
        if (currentEnergy - nextEnergy <= a * t)
            break;
        a *= 0.5;
    }
    if (std::abs(currentEnergy - nextEnergy) <= 0.001) {
        return 0;
    }
    return a;
}

std::vector<double> Solver::estimateLocalPartialDerivativesForV(size_t vId, std::vector<int> localFaceRegion) {
    std::vector<double> changeInEnergy = {0.0, 0.0};
    Vertex& v = mesh.V.at(vId);
    
    double currentMinimumSJ = - mesh.getIterativeEnergyOfRegion(localFaceRegion);
    v.x += DELTA;
    double energyPlusX = - mesh.getIterativeEnergyOfRegion(localFaceRegion);
    v.x -= 2 * DELTA;
    double energyMinusX = - mesh.getIterativeEnergyOfRegion(localFaceRegion);
    v.x += DELTA;
    v.y += DELTA;
    double energyPlusY = - mesh.getIterativeEnergyOfRegion(localFaceRegion);
    v.y -= 2 * DELTA;
    double energyMinusY = - mesh.getIterativeEnergyOfRegion(localFaceRegion);
    v.y += DELTA;
    changeInEnergy.at(0) = energyPlusX - energyMinusX;
    changeInEnergy.at(1) = energyPlusY - energyMinusY;
    return changeInEnergy;
}
