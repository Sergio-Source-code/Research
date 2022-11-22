#include "Solver.h"
#include "MeshFileWriter.h"
#include <iostream>

using namespace std;

int getIndexOf(int num, vector<unsigned long int> vec);
double extendedACos(double value);
double barrier(double x, double s);

void Solver::getHessian() {
    cout << "In getHessian()\n";
    standardHessianCoefficients.clear();
    standardVector.setZero();
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
                // float a = 1;
                // if (f.id == 51)
                //     a = 4;
                if (difference == 0)
                    value = 3;// * a * a - 2 * a + 1;
                else if (difference == 1)
                    value = -2;// * a * a;
                else if (difference == 2)
                    value = 1;// * a * a + 2 * a - 1;
                else
                    value = -2;// * a * a;
                // float area = mesh.getAreaOfFace(f);
                // value *= std::sqrt(area);
                // if (f.id == 916) {
                //     if (difference % 2 == 0)
                //         value += 10;
                //     else
                //         value -= 10;
                // }

                addCoefficientsForIndex(v.id, u.id, value);
            }
        }
    }
    //109, 59 to 60
    // 109 + 59
    //109 + 60 - 2*59
    // int c = 3;
    // addCoefficientsForIndex(107, 107, 1 * c);
    // addCoefficientsForIndex(107, 60, -1 * c);
    // addCoefficientsForIndex(107, 109, 1 * c);
    // addCoefficientsForIndex(107, 109, -1 * c);
    // addCoefficientsForIndex(109, 109, 1 * c);
    // addCoefficientsForIndex(109, 107, -1 * c);

    standardHessian.setFromTriplets(standardHessianCoefficients.begin(), standardHessianCoefficients.end());
}

void Solver::getHessianForEdgeTerms() {
    cout << "In getHessianForEdgeTerms()\n";
    standardHessianCoefficients.clear();
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex v = mesh.V.at(i);
        for (int j = 0; j < v.N_Vids.size(); j++) {
            Vertex u = mesh.V.at(v.N_Vids.at(j));
            addCoefficientsForIndex(v.id, u.id, -1);
        }
    }
    standardHessian.setFromTriplets(standardHessianCoefficients.begin(), standardHessianCoefficients.end());
}

void Solver::getInitialTerm(std::vector<int> vIds) {
    cout << "In getInitialTerm()\n";
    getInitialTermWithLength(vIds);
}

void Solver::getInitialTermWithLength(std::vector<int> vIds) {
    Mesh previousMesh = pastMeshes.at(0);
    float averageArea = mesh.getAverageArea();
    for (int i = 0; i < vIds.size(); i++) {
        Vertex& v = mesh.V.at(vIds.at(i));
        Vertex& prevV = previousMesh.V.at(vIds.at(i));
        float length = glm::length(v - prevV);
        float w = 200 * std::pow(length, 2) / averageArea;
        addCoefficientsForIndex(v.id, v.id, w * 1);
        standardVector(2 * v.id) += - w * v.x;
        standardVector(2 * v.id + 1) += - w * v.y;
    }
    standardHessian.setFromTriplets(standardHessianCoefficients.begin(), standardHessianCoefficients.end());
}

void Solver::getInitialTermWithDepthWeight(std::vector<int> vIds) {
    // let a be the average magnitude of the value of standard terms 
    // defined by d, depth from boundary. Most are 0-2, 3 is max when selecting bad regoions using two ring from corner vertex
    // want 1 to be a lot, maybe around 5a
    // want 2 to approach inner vertex. Curve should be of shape e^-x, maybe around a
    // select 3 as 0.5 * a
    //

    // term for each v,f is (sqrt(3)v - 2v1/sqrt(3) + v2/sqrt(3) - 2v3/sqrt(3))^2
    float sumOfValues = 0;
    float sqrt3 = std::sqrt(3);
    for (int i = 0; i < mesh.F.size(); i++) {
        Face& f = mesh.F.at(i);
        for (int j = 0; j < f.Vids.size(); j++) {
            Vertex& v0 = mesh.V.at(f.Vids.at(j));
            Vertex& v1 = mesh.V.at(f.Vids.at((j + 1) % 4));
            Vertex& v2 = mesh.V.at(f.Vids.at((j + 2) % 4));
            Vertex& v3 = mesh.V.at(f.Vids.at((j + 3) % 4));
            float value = std::pow(sqrt3 * v0.x - 2 * v1.x / sqrt3 + v2.x / sqrt3 - 2 * v3.x / sqrt3, 2);
            value += std::pow(sqrt3 * v0.y - 2 * v1.y / sqrt3 + v2.y / sqrt3 - 2 * v3.y / sqrt3, 2);
            sumOfValues += value;
        }
    }
    float averageOfValues = sumOfValues / (mesh.F.size() * 4);
    double w = 1000;
    for (int i = 0; i < vIds.size(); i++) {
        Vertex& v = mesh.V.at(vIds.at(i));
        std::vector<int> currentVertices;
        currentVertices.push_back(v.id);
        float distanceFromBoundary = mesh.getDistanceFromCorner(currentVertices, 0);
        if (distanceFromBoundary == 1)
            addInitialWeightForV(v.id, w * 5 * averageOfValues);
        else if (distanceFromBoundary >= 2)
            addInitialWeightForV(v.id, w * averageOfValues / std::pow(2, distanceFromBoundary - 2));
    }
    standardHessian.setFromTriplets(standardHessianCoefficients.begin(), standardHessianCoefficients.end());
}

void Solver::getInitialTermWithConstantWeight(std::vector<int> vIds, float weight) {
    for (int i = 0; i < vIds.size(); i++) {
        Vertex& v = mesh.V.at(vIds.at(i));
        // E = (v - v_0)^2
        double w = weight; 
        addCoefficientsForIndex(v.id, v.id, w * 1);
        standardVector(2 * v.id) += - w * v.x;
        standardVector(2 * v.id + 1) += - w * v.y;
    }
    standardHessian.setFromTriplets(standardHessianCoefficients.begin(), standardHessianCoefficients.end());
}

void Solver::addInitialWeightForV(int vId, float weight) {
    Vertex& v = mesh.V.at(vId);
    addCoefficientsForIndex(v.id, v.id, weight * 1);
    standardVector(2 * v.id) += - weight * v.x;
    standardVector(2 * v.id + 1) += - weight * v.y;
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

void Solver::addAverageTermForVertex(int vId) {
    Vertex& v = mesh.V.at(vId);
    int n = v.N_Vids.size();
    float w = 5;
    //E = (v - 1\n (sum of neighbors))^2
    //partial v: (v - 1/n (sum of neighbors))
    addCoefficientsForIndex(v.id, v.id, 1.0 * w);
    for (int i = 0; i < n; i++) {
        Vertex& nV = mesh.V.at(v.N_Vids.at(i));
        addCoefficientsForIndex(v.id, nV.id, - 1.0 * w / n);
    }
    //partial nV: (-1/n v + 1/n^2 (sum of neighbots))
    for (int i = 0; i < n; i++) {
        Vertex& v1 = mesh.V.at(v.N_Vids.at(i));
        addCoefficientsForIndex(v1.id, v.id, - 1.0 * w / n);
        for (int j = 0; j < n; j++) {
            Vertex& v2 = mesh.V.at(v.N_Vids.at(j));
            addCoefficientsForIndex(v1.id, v2.id, 1.0 * w / (n * n));
        }
    }
    standardHessian.setFromTriplets(standardHessianCoefficients.begin(), standardHessianCoefficients.end());
}

void Solver::addAverageTermForFace(int fId) {
    Face& f = mesh.F.at(fId);
    float w = 1;
    for (int i = 0; i < f.Vids.size(); i++) {
        Vertex& v = mesh.V.at(f.Vids.at(i));
        //E = (v - 1/4(sum of v in f))^2
        //E = (3/4v - 1/4v1 - 1/4v2 - 1/4v3)^2
        //partial v: 3/4 * (3/4v0 - 1/4v1 - 1/4v2 - 1/4v3)
        addCoefficientsForIndex(v.id, v.id, w * 9/16);
        for (int j = 0; j < f.Vids.size(); j++) {
            if (i == j)
                continue;
            Vertex& v1 = mesh.V.at(f.Vids.at(j));
            addCoefficientsForIndex(v.id, v1.id, - w * 3 / 16);
        }
        for (int j = 0; j < f.Vids.size(); j++) {
            if (i == j)
                continue;
            //partial v1: -1/4 *  (3/4v0 - 1/4v1 - 1/4v2 - 1/4v3)^2
            Vertex& v1 = mesh.V.at(f.Vids.at(j));
            addCoefficientsForIndex(v1.id, v.id, - w * 3 / 16);
            for (int k = 0; k < f.Vids.size(); k++) {
                if (k == i)
                    continue;
                Vertex& v2 = mesh.V.at(f.Vids.at(k));
                addCoefficientsForIndex(v1.id, v2.id, w * 1 / 16);
            }
        }
    }
    standardHessian.setFromTriplets(standardHessianCoefficients.begin(), standardHessianCoefficients.end());
}

void Solver::updateVertices(std::vector<Vertex> newV) {
    if (mesh.V.size() != newV.size()) {
        std::cout << "ERROR\n";
        std::exit(0);
    }
    for (int i = 0; i < mesh.V.size(); i++) {
        mesh.V.at(i) = newV.at(i);
    }
}

std::vector<Vertex> Solver::iterativeSmoothInnerVertices() {
    std::vector<Vertex> newV;
    for (int i = 0; i < mesh.V.size(); i++) {
        Vertex v = mesh.V.at(i);
        if (v.isBoundary) {
            newV.push_back(v);
            continue;
        }

        std::vector<Vertex> A;
        std::vector<double> B;
        int size = v.N_Vids.size();
        for (int j = 0; j < size; j++) {
            Vertex& v1 = mesh.V.at(v.N_Vids.at(j));
            Vertex& v2 = mesh.V.at(v.N_Vids.at( (j + 1) % size ));
            Vertex& v4 = mesh.V.at(v.N_Vids.at( (j - 1 + size) % size ));

            glm::vec3 a = v.xyz() - v1.xyz();
            glm::vec3 b = v4.xyz() - v1.xyz();
            glm::vec3 c = v2.xyz() - v1.xyz();

            double thetaNew = (1.0 / 2) * (extendedACos((glm::dot(a, b)) / (glm::length(a) * glm::length(b)) ) - extendedACos((glm::dot(a, c)) / (glm::length(a) * glm::length(c))));
            double newX = std::cos(thetaNew) * a.x - std::sin(thetaNew) * a.y;
            double newY = std::sin(thetaNew) * a.x + std::cos(thetaNew) * a.y;
            Vertex aNew;
            aNew.x = v1.x + newX;
            aNew.y = v1.y + newY;
            A.push_back(aNew);
            B.push_back(thetaNew);
        }

        double sumTheta = 0;
        for (int j = 0; j < B.size(); j++) {
            sumTheta += std::abs(B.at(j));
        }

        Vertex newVertex;
        double wSum = 0;
        for (int j = 0; j < A.size(); j++) {
            double w = 1 - std::abs(B.at(j)) / sumTheta;
            wSum += w;
            newVertex.x += w * A.at(j).x;
            newVertex.y += w * A.at(j).y;
        }
        newVertex.x = 1 * ((newVertex.x / wSum) - v.x) + v.x;//0.8 * v.x + 0.2 * newVertex.x;
        newVertex.y = 1 * ((newVertex.y / wSum) - v.y) + v.y;//0.8 * v.y + 0.2 * newVertex.y;
        newV.push_back(newVertex);
    }
    return newV;
}

std::vector<Vertex> Solver::iterativeSmoothBoundaryVertices() {
    std::vector<Vertex> newV;

    return mesh.V;
}

void Solver::iterativeOptimization() {
    std::cout << "In iterativeOptimization()\n";
    double currentEnergy = 0;
    for (int i = 0; i < 30; i++) {
        std::vector<Vertex> newV = iterativeSmoothInnerVertices();
        updateVertices(newV);
        std::string fileName = std::string("/examples/iterative_") + std::to_string(i) + ".vtk";
        MeshFileWriter fw(mesh, fileName.c_str());
        fw.WriteFile();
        // for (int j = 0; j < 10; j++) {
        //     std::vector<Vertex> newBoundaryV = iterativeSmoothBoundaryVertices();
        //     updateVertices(newBoundaryV);
        // }
    }
    std::cout << "Out iterativeOptimization()\n";
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

double extendedACos(double value) {
    if (value >= 0.99999)
        return std::acos(0.99999);
    if (value <= -0.99999)
        return std::acos(-0.99999);
    return std::acos(value);
}

int FACE_ID = 6;
int VERTEX_ID = 486;//24;//59;//24;
double MAX_SIZE = 0.2;//1;//0.1;
double DELTA = 0.001;//0.001
double STEP_SIZE = 0.003;//2000;//20

void Solver::iterativeOptimizeMinSJ() {
    std::cout << "In iterativeOptimizeMinSJ()\n";
    double minimumEdgeLength = 1000000000;
    std::vector<std::vector<int>> localFaceRegions = mesh.getLocalFaceRegions();
    for (int i = 0; i < localFaceRegions.size(); i++) {
        // std::cout << "i: " << i << std::endl;
        std::vector<int> localFaceRegion = localFaceRegions.at(i);
        for (int j = 0; j < localFaceRegion.size(); j++) {
            std::cout << localFaceRegion.at(j) << ", ";
            Face& f = mesh.F.at(localFaceRegion.at(j));
            for (int k = 0; k < f.Eids.size(); k++) {
                Edge& e = mesh.E.at(f.Eids.at(k));
                double currentEdgeLength = glm::length(mesh.V.at(e.Vids.at(1)) - mesh.V.at(e.Vids.at(0)));
                if (currentEdgeLength < minimumEdgeLength) {
                    minimumEdgeLength = currentEdgeLength;
                }
            }
            // if (i == 4) {
            //     std::cout << f.id << ": ";
            //     for (int m = 0; m < 4; m++) {
            //         std::cout << f.Vids.at(m) << " ";
            //     }
            //     std::cout << std::endl;
            // }
        }
        std::cout << std::endl;
    }
    DELTA = minimumEdgeLength / 10;
    // std::cout << "CURRENT MIN, DELTA: " << minimumEdgeLength << " " << DELTA << std::endl;
    // return;
    int counter = 0;
    for (int l = 0; l < localFaceRegions.size(); l++) {
        std::vector<int> localFaceRegion = localFaceRegions.at(2);
        localFaceRegion.push_back(499);
        for (int i = 0; i < 100; i++) {
            std::string fileName = std::string("/examples/_iterativeMinSJ_") + std::to_string(counter++) + ".vtk";
            MeshFileWriter fw(mesh, fileName.c_str());
            fw.WriteFile();
            std::vector<double> gradient(2 * mesh.V.size());
            for (int j = 0; j < mesh.V.size(); j++) {
                Vertex& v = mesh.V.at(j);
                if (v.isBoundary)
                    continue;
                std::vector<double> estimatedPartialDerivatives = estimateMinSJGradient(v.id, localFaceRegion);
                gradient.at(2 * v.id) += estimatedPartialDerivatives.at(0);
                gradient.at(2 * v.id + 1) += estimatedPartialDerivatives.at(1);
            }

            //step size
            double t = 0;
            for (int j = 0; j < gradient.size(); j++) {
                t += std::pow(gradient.at(j), 2);
            }
            if (t <= 0) {
                std::cout << "ERROR: t (should be +) is " << t << std::endl;
                std::exit(0);
            }
            t *= -1;
            double a = 0.1;
            // a *= 10;
            double currentEnergy = mesh.getIterativeEnergyOfRegion(localFaceRegion);
            double nextEnergy;
            for (int j = 0; j < 10; j++) {
                //f(x)-f(x+ap)
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
                // std::cout << currentEnergy << " " << nextEnergy << " " << (currentEnergy - nextEnergy) << " " << (a * t) << std::endl;
                if (currentEnergy - nextEnergy <= a * t)
                    break;
                a *= 0.2;
            }
            // return;
            STEP_SIZE = a;
            //old step size
            // std::cout << "a, diff = " << a << " " << (currentEnergy - nextEnergy) << std::endl;
            if (std::abs(currentEnergy - nextEnergy) <= 0.001) {
                std::cout << "Iterative optimization done in " << i << " iterations\n";
                // break;
            }
            for (int j = 0; j < mesh.V.size(); j++) {
                Vertex& v = mesh.V.at(j);
                v.x -= STEP_SIZE * gradient.at(2 * v.id);
                v.y -= STEP_SIZE * gradient.at(2 * v.id + 1);
            }
        }
    }
    std::cout << "-------------------------\n";
}

void Solver::iterativeOptimizeCornerVertices() {
    std::cout << "In iterativeOptimizeCornerVertices()\n";
    std::vector<double> initialSJ; //stores initial sj value of all v/f
    for (int i = 0; i < mesh.F.size(); i++) {
        Face f = mesh.F.at(i);
        for (int j = 0; j < 4; j++) {
            Vertex& v = mesh.V.at(f.Vids.at(j));
        }
        initialSJ.push_back(mesh.getMinimumScaledJacobian(f));
    }

    for (int i = 0; i < 100; i++) {
        std::string fileName = std::string("/examples/corner_iterative_") + std::to_string(i) + ".vtk";
        MeshFileWriter fw(mesh, fileName.c_str());
        fw.WriteFile();
        std::vector<double> gradient(2 * mesh.V.size());
        for (int j = 0; j < mesh.F.size(); j++) {
            Face& f = mesh.F.at(j);
            for (int k = 0; k < 4; k++) {
                Vertex& v = mesh.V.at(f.Vids.at(k));
                if (v.isBoundary)
                    continue;
                // if (glm::length(v - mesh.V.at(158)) >= 30)
                //     continue;
                double localInitialSJ = initialSJ.at(f.id);
                std::vector<double> estimatedPartialDerivatives = estimateGradientForFV(f.id, v.id, localInitialSJ);
                gradient.at(2 * v.id) += estimatedPartialDerivatives.at(0);
                gradient.at(2 * v.id + 1) += estimatedPartialDerivatives.at(1);

                // if (i <= 8 && v.id == 25) {
                //     std::cout << "---gradient term---" << f.id << " " << estimatedPartialDerivatives.at(0) << " " << estimatedPartialDerivatives.at(1) << std::endl;
                // }
            }
        }

        for (int j = 0; j < mesh.V.size(); j++) {
            Vertex& v = mesh.V.at(j);
            if (glm::length(v - mesh.V.at(VERTEX_ID)) <= MAX_SIZE) {
                if (i <= 8) {
                    std::cout << i << " final term: " << gradient.at(2 * v.id) << " " << gradient.at(2 * v.id + 1) << std::endl;
                }
                v.x -= STEP_SIZE * gradient.at(2 * v.id);
                v.y -= STEP_SIZE * gradient.at(2 * v.id + 1);
            }
        }
    }
    std::cout << "-------------------------\n";
}

#include <iostream>
#include <iomanip>

std::vector<double> Solver::estimateMinSJGradient(size_t vId, std::vector<int> localFaceRegion) {
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
    // if (v.id == 431)
    //     std::cout << currentMinimumSJ << " " << energyPlusX << std::endl;
    changeInEnergy.at(0) = energyPlusX - energyMinusX;
    changeInEnergy.at(1) = energyPlusY - energyMinusY;
    return changeInEnergy;
}

std::vector<double> Solver::estimateGradientForFV(size_t fId, size_t vId, float initialSJ) {
    // only have to compute energy terms which change in energy.
    // There exists such a term for each neighboring facee, and all vertices in face other than oposite corners
    std::vector<double> changeInEnergy = {0.0, 0.0};
    Vertex& v = mesh.V.at(vId);
    Face& nF = mesh.F.at(fId);
    // // for (int i = 0; i < v.N_Fids.size(); i++) {
    // //     Face& nF = mesh.F.at(v.N_Fids.at(i));
    //     int vIndexInNF = getIndexOf(v.id, nF.Vids);
    //     for (int j = 0; j < nF.Vids.size(); j++) {
    //         Vertex& nV = mesh.V.at(nF.Vids.at(j));
    //         int nVIndexINNF = getIndexOf(nV.id, nF.Vids);
    //         if (std::abs(vIndexInNF - nVIndexINNF) == 2)
    //             continue;

    //         // here we are looking at the term for nF/nV, and addings its change in energy when changing v.x
    //         // original energy is (v1 - v) x (v2 - v) / (|v1-v| + |v2-v|)
    //         Vertex& v1 = mesh.V.at(nF.Vids.at( (nVIndexINNF + 1) % 4 ));
    //         Vertex& v2 = mesh.V.at(nF.Vids.at( (nVIndexINNF + 3) % 4 ));
    //         nV.x += DELTA;
    //         double energyPlusX = glm::cross(v1 - nV, v2 - nV).z / (glm::length(v1 - nV) * glm::length(v2 - nV));
    //         nV.x -= 2 * DELTA;
    //         double energyMinusX = glm::cross(v1 - nV, v2 - nV).z / (glm::length(v1 - nV) * glm::length(v2 - nV));
    //         nV.x += DELTA;
    //         nV.y += DELTA;
    //         double energyPlusY = glm::cross(v1 - nV, v2 - nV).z / (glm::length(v1 - nV) * glm::length(v2 - nV));
    //         nV.y -= 2 * DELTA;
    //         double energyMinusY = glm::cross(v1 - nV, v2 - nV).z / (glm::length(v1 - nV) * glm::length(v2 - nV));
    //         nV.y += DELTA;
    //         if (nF.id != FACE_ID && nF.id != 195) {
    //             energyPlusX = barrier(energyPlusX, initialSJ);
    //             energyMinusX = barrier(energyMinusX, initialSJ);
    //             energyPlusY = barrier(energyPlusY, initialSJ);
    //             energyMinusY = barrier(energyMinusY, initialSJ);
    //         } else {
    //             energyPlusX *= -10;
    //             energyMinusX *= -10;
    //             energyPlusY *= -10;
    //             energyMinusY *= -10;
    //         }
    //         if (nF.id == 421 && v.id == 435) {// && nV.id != 620) {
    //             // continue;
    //             // std::cout << "-----\n";
    //             // std::cout << energyPlusX << " " << energyMinusX << std::endl;
    //             // std::cout << energyPlusY << " " << energyMinusY << std::endl;
    //             // std::cout << nF.id << " " << nV.id << " " << v1.id << " " << v2.id << std::endl;
    //         }
    //         changeInEnergy.at(0) += energyPlusX - energyMinusX;
    //         changeInEnergy.at(1) += energyPlusY - energyMinusY;
    //     }
    // // }
    // return changeInEnergy;
    Face& f = mesh.F.at(fId);
    double minimumSJ = 1.0;
    // for (int i = 0; i < f.Vids.size(); i++) {
        Vertex& nV = v;//mesh.V.at(f.Vids.at(i));
        double currentMinimumSJ = mesh.getMinimumScaledJacobianOfRegion(mesh.V.at(VERTEX_ID), MAX_SIZE);
        nV.x += DELTA;
        double energyPlusX = mesh.getMinimumScaledJacobianOfRegion(mesh.V.at(VERTEX_ID), MAX_SIZE);
        nV.x -= 2 * DELTA;
        double energyMinusX = mesh.getMinimumScaledJacobianOfRegion(mesh.V.at(VERTEX_ID), MAX_SIZE);
        nV.x += DELTA;
        nV.y += DELTA;
        double energyPlusY = mesh.getMinimumScaledJacobianOfRegion(mesh.V.at(VERTEX_ID), MAX_SIZE);
        nV.y -= 2 * DELTA;
        double energyMinusY = mesh.getMinimumScaledJacobianOfRegion(mesh.V.at(VERTEX_ID), MAX_SIZE);
        nV.y += DELTA;
        if (f.id == 164 && v.id == 390) {
            std::cout << "164 before: " << energyPlusX << std::endl;
        }
        if (true || f.id == FACE_ID || f.id == 164) {
            energyPlusX *= -0.5;
            energyMinusX *= -0.5;
            energyPlusY *= -0.5;
            energyMinusY *= -0.5;
        } else {
            energyPlusX = -barrier(energyPlusX, initialSJ);
            energyMinusX = -barrier(energyMinusX, initialSJ);
            energyPlusY = -barrier(energyPlusY, initialSJ);
            energyMinusY = -barrier(energyMinusY, initialSJ);
            double inf = 1000000.0;
            if (energyPlusX == inf || energyMinusX == inf || energyPlusY == inf || energyMinusY == inf) {
                // DELTA /= 2;
                // std::cout << "new delta " << DELTA << std::endl;
                // return estimateGradientForFV(fId, vId, initialSJ);
                // energyPlusX = 0;
                // energyMinusX = 0;
                // energyPlusY = 0;
                // energyMinusY = 0;
                // currentMinimumSJ = -1.0;
            }
        }
        if (true || currentMinimumSJ < minimumSJ) {
            minimumSJ = currentMinimumSJ;
            changeInEnergy.at(0) = energyPlusX - energyMinusX;
            changeInEnergy.at(1) = energyPlusY - energyMinusY;
        }
        if (f.id == 164 && v.id == 390) {
            std::cout << "164 after: " << energyPlusX << std::endl;
        }
    // }
    return changeInEnergy;
}

double barrier(double x, double s) {
    if (x >= s)
        return 0.0;
    else if (x <= 0)
        return 1000000.0;
    double arg = x / s;
    double g = std::pow(arg, 3) - 3 * std::pow(arg, 2) + 3 * arg;
    return std::pow((1 / g), 10) - 1;
}