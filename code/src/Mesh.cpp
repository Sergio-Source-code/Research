/*
 * Mesh.cpp
 *
 *  Created on: Nov 6, 2016
 *      Author: cotrik
 *  Edited by: sergio
 */

#include "Mesh.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "glm/gtx/intersect.hpp"
#include <algorithm>
#include <map>
#include <iostream>
#include <queue>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkDataSet.h>
#include <vtkMeshQuality.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkProperty.h>
#include <vtkScalarBarActor.h>
#include <vtkLookupTable.h>
#include <vtkColorTransferFunction.h>
#include <vtkThreshold.h>
#include <vtkExtractEdges.h>

const size_t MAXID = 0xffffffffffffffff;

Mesh::Mesh()
: m_cellType(QUAD) {}

Mesh::Mesh(const Mesh& r)
: V(r.V)
, E(r.E)
, F(r.F)
, C(r.C)
, m_cellType(r.m_cellType)
, pointScalarFieldNames(r.pointScalarFieldNames)
, pointScalarFields(r.pointScalarFields)
, cellScalarNameFields(r.cellScalarNameFields)
, cellScalarFields(r.cellScalarFields) {}

Mesh::Mesh(const std::vector<Vertex>& V, const std::vector<Cell>& C, ElementType m_cellType)
: V(V)
, C(C)
, m_cellType(m_cellType) {}

Mesh::Mesh(const std::vector<Vertex>& V, const std::vector<Face>& F, ElementType m_cellType)
: V(V)
, F(F)
, m_cellType(m_cellType) {}

Mesh::Mesh(const Mesh& r, const std::vector<size_t>& cellIds)
: m_cellType(r.m_cellType) {
    V.resize(r.V.size());
    // Read V
    for (vtkIdType i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        const Vertex& rv = r.V.at(i);
        v.x = rv.x;
        v.y = rv.y;
        v.z = rv.z;
        v.id = i;
    }
    // Read C
    C.resize(cellIds.size());
    for (vtkIdType i = 0; i < cellIds.size(); i++)
        C.at(i).Vids = r.C.at(cellIds.at(i)).Vids;

    if (m_cellType == TRIANGLE || r.m_cellType == QUAD) {
        F.resize(C.size());
        for (vtkIdType i = 0; i < cellIds.size(); i++)
            F[i].Vids = C[i].Vids;
    }
}

Mesh::~Mesh() {
    V.clear();
    E.clear();
    F.clear();
    C.clear();

    pointScalarFieldNames.clear();
    pointScalarFields.clear();
    cellScalarNameFields.clear();
    cellScalarFields.clear();
}

const unsigned int QuadEdge[4][2] = {
    { 0, 1 },
    { 1, 2 },
    { 2, 3 },
    { 3, 0 }
};

void Mesh::process() {
    BuildAllConnectivities();
    ExtractBoundary();
    orderNeighboringVertices();
    orderVerticesInFaces();
    classifyCornerVertices();
}

void Mesh::BuildAllConnectivities() {
    BuildV_C();
    BuildE();
    BuildF();
    BuildF_F();
}

void Mesh::BuildV_C() {
    for (size_t i = 0; i < C.size(); i++)
        for (size_t j = 0; j < C[i].Vids.size(); j++)
            V[C[i].Vids[j]].N_Fids.push_back(i);
}

void Mesh::BuildE() {
    std::vector<Edge> Es;
    for (size_t i = 0; i < C.size(); i++) {
        Edge e(2);
        for (size_t j = 0; j < 4; j++) {
            e.Vids[0] = C[i].Vids[QuadEdge[j][0]];
            e.Vids[1] = C[i].Vids[QuadEdge[j][1]];
            Es.push_back(e);
        }
    }

    size_t E_N = 0;
    for (size_t i = 0; i < Es.size(); i++) {
        bool havesame = false;
        const size_t id1 = Es[i].Vids[0];
        const size_t id2 = Es[i].Vids[1];

        for (size_t j = 0; j < V[id1].N_Vids.size(); j++) {
            if (V[id1].N_Vids[j] == id2) {
                havesame = true;
                break;
            }
        }
        if (!havesame) {
            Es[i].id = E_N++;
            Es[i].isBoundary = false;

            for (size_t j = 0; j < V[Es[i].Vids[0]].N_Cids.size(); j++) {
                size_t h1 = V[Es[i].Vids[0]].N_Cids[j];
                for (size_t k = 0; k < V[Es[i].Vids[1]].N_Cids.size(); k++) {
                    size_t h2 = V[Es[i].Vids[1]].N_Cids[k];
                    if (h1 == h2)
                        Es[i].N_Cids.push_back(h1);
                }
            }

            E.push_back(Es[i]);
            V[id1].N_Eids.push_back(Es[i].id);
            V[id1].N_Vids.push_back(id2);

            V[id2].N_Eids.push_back(Es[i].id);
            V[id2].N_Vids.push_back(id1);
        }
    }

    for (size_t i = 0; i < E.size(); i++) {
        for (size_t j = 0; j < V[E[i].Vids[0]].N_Cids.size(); j++) {
            size_t nh = V[E[i].Vids[0]].N_Cids[j];
            for (size_t k = 0; k < V[E[i].Vids[1]].N_Cids.size(); k++)
                if (nh == V[E[i].Vids[1]].N_Cids[k])
                    E[i].N_Cids.push_back(nh);
        }
    }

    for (size_t i = 0; i < E.size(); i++) {
        std::vector<size_t> N_Cids;
        for (size_t j = 0; j < E[i].N_Cids.size(); j++) {
            bool have = false;
            for (size_t k = 0; k < N_Cids.size(); k++)
                if (E[i].N_Cids[j] == N_Cids[k])
                    have = true;
            if (!have)
                N_Cids.push_back(E[i].N_Cids[j]);
        }
        E[i].N_Cids.clear();
        E[i].N_Cids = N_Cids;
    }

    for (size_t i = 0; i < V.size(); i++) {
        std::vector<size_t> N_Vids;
        for (int j = 0; j < V[i].N_Vids.size(); j++) {
            bool have = false;
            for (int k = 0; k < N_Vids.size(); k++)
                if (V[i].N_Vids[j] == N_Vids[k])
                    have = true;
            if (!have)
                N_Vids.push_back(V[i].N_Vids[j]);
        }
        V[i].N_Vids.clear();
        V[i].N_Vids = N_Vids;
    }
}

void Mesh::BuildF() {
    F.resize(C.size());
    for (size_t i = 0; i < C.size(); i++) {
        F[i].Vids = C[i].Vids;
        F[i].id = C[i].id;
        F[i].isBoundary = true;
    }
    for (size_t i = 0; i < F.size(); i++) {
        F[i].Eids.resize(4);
        for (size_t j = 0; j < 4; j++) {
            bool found = false;
            for (size_t k = 0; k < V[F[i].Vids[j]].N_Eids.size(); k++) {
                size_t ne = V[F[i].Vids[j]].N_Eids[k];
                for (size_t m = 0; m < V[F[i].Vids[(j + 1) % 4]].N_Eids.size(); m++) {
                    if (ne == V[F[i].Vids[(j + 1) % 4]].N_Eids[m]) {
                        E[ne].N_Fids.push_back(i);
                        F[i].Eids[j] = ne;
                        found = true;
                        break;
                    }
                }
                if (found) break;
            }
        }
    }
    ////////////////////////////////////////////////
    std::vector<bool> Vs_flags(V.size(), false);
    for (int i = 0; i < F.size(); i++) {
        for (int j = 0; j < 3; j++)
            Vs_flags[F[i].Vids[j]] = true;
        for (int j = 0; j < F[i].N_Cids.size(); j++) {
            int nhid = F[i].N_Cids[j];
            for (int k = 0; k < 4; k++) {
                bool have_true = false;
                for (int m = 0; m < 3; m++)
                    if (Vs_flags[F[C[nhid].Fids[k]].Vids[m]]) {
                        have_true = true;
                        break;
                    }
                if (!have_true) {
                    F[i].N_Fids.push_back(C[nhid].Fids[k]);
                    break;
                }
            }
        }
    }
}

void Mesh::BuildF_F() {
    for (int i = 0; i < F.size(); i++) {
        Face& f = F[i];
        for (int j = 0; j < f.Vids.size(); j++) {
            Vertex& v = V[f.Vids[j]];
            for (int k = 0; k < v.N_Fids.size(); k++) {
                int f1Id = v.N_Fids.at(k);
                if (std::find(f.N_Fids.begin(), f.N_Fids.end(), f1Id) == f.N_Fids.end())
                    f.N_Fids.push_back(f1Id);
            }
        }
    }
}

void Mesh::ExtractBoundary() {
    for (size_t i = 0; i < V.size(); i++)
        V[i].isBoundary = false;
    for (size_t i = 0; i < E.size(); i++)
        E[i].isBoundary = false;
    for (size_t i = 0; i < F.size(); i++)
        F[i].isBoundary = false;

    for (size_t i = 0; i < E.size(); i++) {
        if (E[i].N_Fids.size() == 1) {
            V[E[i].Vids[0]].isBoundary = true;
            V[E[i].Vids[1]].isBoundary = true;
            E[i].isBoundary = true;
        }
    }
    for (size_t i = 0; i < F.size(); i++) {
        for (size_t j = 0; j < F[i].Vids.size(); j++) {
            if (V[F[i].Vids[j]].isBoundary) {
                F[i].isBoundary = true;
                break;
            }
        }
    }
}

void Mesh::orderNeighboringVertices() {
    std::vector<int> verticesToRecheck;
    // for (size_t i = 0; i < V.size(); i++) {
    //     Vertex& v = V.at(i);
    //     if (v.isBoundary)
    //         continue;

    //     bool done = false;
    //     size_t currentFaceId = v.N_Fids.at(0);
    //     int secondFaceSelected = -1;
    //     std::vector<size_t> polygon;
    //     // create ordered polygon of neighboring vertices
    //     while (!done) {
    //         Face& f = F.at(currentFaceId);
    //         int neighboringVerticesInPolygon = 0;
    //         for (int j = 0; j < f.Vids.size(); j++) {
    //             Vertex& v1 = V.at(f.Vids.at(j));
    //             if (std::find(v.N_Vids.begin(), v.N_Vids.end(), v1.id) != v.N_Vids.end()) {//  if v1 neighbors v
    //                 if (std::find(polygon.begin(), polygon.end(), v1.id) == polygon.end()) {// and is not in polygon
    //                     polygon.push_back(v1.id);
    //                     for (int k = 0; k < v1.N_Fids.size(); k++) {
    //                         Face& f1 = F.at(v1.N_Fids.at(k));
    //                         if (f1.id != currentFaceId && std::find(f1.Vids.begin(), f1.Vids.end(), v.id) != f1.Vids.end()) {
    //                             currentFaceId = f1.id;
    //                             if (secondFaceSelected == -1)
    //                                 secondFaceSelected = f1.id;
    //                             break;
    //                         }
    //                     }
    //                     break;
    //                 }
    //                 else {// v1 neighbors v and already is in polygon
    //                     neighboringVerticesInPolygon++;
    //                     if (neighboringVerticesInPolygon == 2) {
    //                         done = true;
    //                         break;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    //     if (v.N_Vids.size() != polygon.size()) {
    //         std::cout << v.N_Vids.size() << " != " << polygon.size() << " for " << v.id << std::endl;
    //         std::cout << "ERROR\n";
    //     }
       
    //     int numCounterClockwise = 0;
    //     for (int j = 0; j < polygon.size() - 1; j++) {
    //         Vertex& v1 = V.at(polygon.at(j));
    //         Vertex& v2 = V.at(polygon.at(j + 1));
    //         const glm::vec3 crossProduct = glm::cross(v1 - v, v2 - v);
    //         if (crossProduct.z >= 0)
    //             numCounterClockwise++;
    //         else
    //             numCounterClockwise--;
    //     }
    //     if (numCounterClockwise == 0) {
    //         // check later
    //         checkOrderOfVertices.push_back(v.id);
    //     }
    //     else if (numCounterClockwise < 0) {
    //         std::reverse(polygon.begin(), polygon.end());
    //     }  

    //     v.N_Vids.clear();
    //     for (int j = 0; j < polygon.size(); j++) {
    //         v.N_Vids.push_back(polygon.at(j));
    //     }
    // }
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isBoundary)
            continue;

        // create ordered polygon by traversing edges
        std::vector<int> polygon;
        for (int j = 0; j < v.N_Eids.size(); j++) {
            Edge& e = E.at(v.N_Eids.at(j));
            if (e.isBoundary) {
                int v1Id = e.Vids.at(0) != v.id ? e.Vids.at(0) : e.Vids.at(1);
                polygon.push_back(v1Id);
                break;
            }
        }
        while (polygon.size() != v.N_Vids.size()) {
            int previousSize = polygon.size();
            int previousVId = polygon.at(polygon.size() - 1);
            for (int j = 0; j < v.N_Vids.size(); j++) {
                Vertex& v2 = V.at(v.N_Vids.at(j));
                // if v2Id not in polygon and previousVId is in neighboring faces, push to polygon and break
                if (std::find(polygon.begin(), polygon.end(), v2.id) == polygon.end()) {
                    for (int k = 0; k < v2.N_Fids.size(); k++) {
                        Face& f = F.at(v2.N_Fids.at(k));
                        if (std::find(f.Vids.begin(), f.Vids.end(), previousVId) != f.Vids.end()) {
                            polygon.push_back(v2.id);
                            break;
                        }
                    }
                }
                int newSize = polygon.size();
                if (newSize != previousSize)
                    break;
            }
            int newSize = polygon.size();
            if (newSize == previousSize) {
                std::cerr << "ERROR: could not order neighboring vertices of boundary vertex " << v.id << std::endl;
                std::exit(1);
            }
        }

        int numCounterClockwise = 0;
        for (int j = 0; j < polygon.size() - 1; j++) {
            Vertex& v1 = V.at(polygon.at(j));
            Vertex& v2 = V.at(polygon.at(j + 1));
            const glm::vec3 crossProduct = glm::cross(v1 - v, v2 - v);
            if (crossProduct.z >= 0)
                numCounterClockwise++;
            else
                numCounterClockwise--;
        }
        if (numCounterClockwise == 0) {
            verticesToRecheck.push_back(v.id);
        }
        if (numCounterClockwise < 0) {
            std::reverse(polygon.begin(), polygon.end());
        }
        v.N_Vids.clear();
        for (int j = 0; j < polygon.size(); j++) {
            v.N_Vids.push_back(polygon.at(j));
        }

        // Vertex& b1 = V.at(polygon.at(0));
        // Vertex& b2 = V.at(polygon.at(polygon.size() - 1));
        // if (glm::cross(b1 - v, b2 - v).z <= -0.3 * glm::length(b1 - v) * glm::length(b2 - v)) //tolerance for straight boundary edges. float * sin of angle is tolerance
        //     v.isCorner = true;
        // else
        //     v.isCorner = false;
    }

    while (verticesToRecheck.size() != 0) {
        Vertex& v = V.at(verticesToRecheck.back());
        verticesToRecheck.pop_back();
        bool foundOrientation = false;
        for (int j = 0; j < v.N_Fids.size(); j++) {
            Face& f = F.at(v.N_Fids.at(j));
            int localId = getLocalIndexOf(v.id, f.Vids);
            Vertex& v1 = V.at(f.Vids.at( (localId + 2) % 4 ));
            if (std::find(verticesToRecheck.begin(), verticesToRecheck.end(), v1.id) == verticesToRecheck.end()) {// if v1 is not in verticesToRecheck
                // should be reversed
                int ordering = getRelativeOrientation(v.N_Vids, v1.N_Vids);
                if (ordering == 0) {
                    std::cerr << "ERROR determining order of " << v.id << std::endl;
                    std::exit(1);
                } else if (ordering == 1) {
                    std::reverse(v.N_Vids.begin(), v.N_Vids.end());
                }
                foundOrientation = true;
                break;
            }
        }
        if (!foundOrientation) {
            std::cerr << "ERROR: unable to find orientation of neighboring vertices of " << v.id << std::endl;
            std::exit(1);
        }
    }
}

void Mesh::orderVerticesInFaces() {
    std::vector<int> facesToRecheck;
    for (int i = 0; i < F.size(); i++) {
        Face& f = F.at(i);
        int numCounterClockwise = 0;
        for (int j = 0; j < f.Vids.size(); j++) {
            Vertex& v = V.at(f.Vids.at(j));
            Vertex& v1 = V.at(f.Vids.at( (j + 1) % f.Vids.size() ));
            Vertex& v2 = V.at(f.Vids.at( (f.Vids.size() + j - 1) % f.Vids.size() ));
            const glm::vec3 crossProduct = glm::cross(v1 - v, v2 - v);
            if (crossProduct.z >= 0)
                numCounterClockwise++;
            else
                numCounterClockwise--;
        }
        if (numCounterClockwise == 0) {
            facesToRecheck.push_back(f.id);
        } else if (numCounterClockwise < 0) {
            std::reverse(f.Vids.begin(), f.Vids.end());
        }
    }

    for (int i = 0; i < facesToRecheck.size(); i++) {
        Face& f = F.at(facesToRecheck.at(i));
        for (int j = 0; j < f.N_Fids.size(); j++) {
            Face& f1 = F.at(f.N_Fids.at(j));
            if (std::find(facesToRecheck.begin(), facesToRecheck.end(), f1.id) != facesToRecheck.end())
                continue;
            // f1 not in in facesToRecheck

            int ordering = getRelativeOrientation(f.Vids, f1.Vids);
            if (ordering != 0) {
                if (ordering == 1) {
                    std::reverse(f.Vids.begin(), f.Vids.end());
                }
                break;
            }
        }
    }
}

void Mesh::classifyCornerVertices() {
    for (int i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isBoundary)
            continue;
        Vertex& b1 = V.at(v.N_Vids.at(0));
        Vertex& b2 = V.at(v.N_Vids.at(v.N_Vids.size() - 1));
        if (glm::cross(b1 - v, b2 - v).z <= -0.3 * glm::length(b1 - v) * glm::length(b2 - v)) // tolerance of sin(theta) = 0.3
            v.isCorner = true;
        else
            v.isCorner = false;
    }
}

void Mesh::printQuality() {
    double minSJ = 1;
    double maxSJ = -1;
    double avgSJ = 0;
    int numInverted = 0;
    int minIndex = -1;
    std::vector<int> invertedFIds;
    for (int i = 0; i < F.size(); i++) {
        Face& f = F.at(i);
        double jacobian = 1.0;
        for (int j = 0; j < 4; j++) {
            glm::vec3 v0 = V.at(f.Vids.at(j));
            glm::vec3 v1 = V.at(f.Vids.at( (j + 1) % 4 ));
            glm::vec3 v2 = V.at(f.Vids.at( (j + 2) % 4 ));
            double newJacobian = glm::cross(glm::normalize(v2 - v1), glm::normalize(v0 - v1)).z;
            if (newJacobian < jacobian)
                jacobian = newJacobian;
        }
        if (jacobian <= 0) {
            numInverted++;
            invertedFIds.push_back(f.id);
        }
        if (jacobian < minSJ) {
            minSJ = jacobian;
            minIndex = i;
        }
        if (jacobian > maxSJ)
            maxSJ = jacobian;
        avgSJ += jacobian;
    }
    avgSJ /= F.size();
    std::cout << "-------------------------" << std::endl;
    std::cout << "minSJ:          " << minSJ << std::endl;
    std::cout << "avjSJ:          " << avgSJ << std::endl;
    std::cout << "maxSJ:          " << maxSJ << std::endl;
    std::cout << "num inverted:   " << numInverted << std::endl;
    std::cout << "min face id:    " << minIndex << std::endl;
    if (invertedFIds.size() != 0) {
        std::cout << "inverted faces:" << std::endl;
        for (int i = 0; i < invertedFIds.size(); i++) {
            std::cout << invertedFIds.at(i);
            if (i != invertedFIds.size() - 1) std::cout << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << "-------------------------" << std::endl;
}

int Mesh::getNumOfInvertedElements(std::vector<int> localFaceRegion) {
    int result = 0;
    for (int i = 0; i < localFaceRegion.size(); i++) {
        Face& f = F.at(localFaceRegion.at(i));
        if (getMinimumScaledJacobian(f) <= 0) {
            result++;
        }
    }
    return result;
}

double Mesh::getMinimumScaledJacobian(std::vector<int> faceIds) {
    double result = 1.0;
    for (int i = 0; i < faceIds.size(); i++) {
        Face& f = F.at(faceIds.at(i));
        double localMinSJ = getMinimumScaledJacobian(f);
        if (localMinSJ < result)
            result = localMinSJ;
    }
    return result;
}

double Mesh::getMinimumScaledJacobian(Face& f) {
    double result = 1.0;
    for (int j = 0; j < 4; j++) {
        glm::vec3 v0 = V.at(f.Vids.at(j));
        glm::vec3 v1 = V.at(f.Vids.at( (j + 1) % 4 ));
        glm::vec3 v2 = V.at(f.Vids.at( (j + 2) % 4 ));
        double currentScaledJacobian = glm::cross(glm::normalize(v2 - v1), glm::normalize(v0 - v1)).z;
        if (currentScaledJacobian < result)
            result = currentScaledJacobian;
    }
    return result;
}

float Mesh::getAreaOfFace(Face& f) {
    float area = 0;
    for (int i = 0; i < 4; i++) {
        Vertex& v = V.at(f.Vids.at(i));
        Vertex& v1 = V.at(f.Vids.at( (i + 1) % 4 ));
        Vertex& v2 = V.at(f.Vids.at( (i + 3) % 4 ));
        float cross = glm::length(glm::cross(v1 - v, v2 - v));
        area += cross / 4;
    }

    return area;
}

float Mesh::getMaximumEdgeLength(Face& f) {
    float maxEdgeLength = 0.0;
    for (int i = 0; i < f.Eids.size(); i++) {
        Edge& e = E.at(f.Eids.at(i));
        double currentLength = glm::length(V.at(e.Vids.at(1)) - V.at(e.Vids.at(0)));
        if (currentLength > maxEdgeLength)
            maxEdgeLength = currentLength;
    }
    return maxEdgeLength;
}

float Mesh::getMinimumEdgeLength(Vertex& v) {
    float minEdgeLength = std::pow(10, 30);
    for (int i = 0; i < v.N_Eids.size(); i++) {
        Edge& e = E.at(v.N_Eids.at(i));
        double currentLength = glm::length(V.at(e.Vids.at(1)) - V.at(e.Vids.at(0)));
        if (currentLength < minEdgeLength)
            minEdgeLength = currentLength;
    }
    return minEdgeLength;
}

std::vector<int> Mesh::getInnerVerticesOfFaces(std::vector<int> fIds) {
    std::vector<int> vertices = collapseFacesToVerticies(fIds);
    std::vector<int> boundaryVertices;

    for (int i = 0; i < fIds.size(); i++) {
        Face& f = F.at(fIds.at(i));
        for (int j = 0; j < f.Eids.size(); j++) {
            Edge& e = E.at(f.Eids.at(j));
            int nFID = -1;
            if (e.N_Fids.size() == 2)
                nFID = e.N_Fids.at(0) == f.id ? e.N_Fids.at(1) : e.N_Fids.at(0);
            if (nFID == -1 || std::find(fIds.begin(), fIds.end(), nFID) == fIds.end()) {
                boundaryVertices.push_back(e.Vids.at(0));
                boundaryVertices.push_back(e.Vids.at(1));
            }
        }
    }
    std::sort(boundaryVertices.begin(), boundaryVertices.end());
    auto last = std::unique(boundaryVertices.begin(), boundaryVertices.end());
    boundaryVertices.erase(last, boundaryVertices.end());
    std::vector<int> innerVertices;
    std::set_difference(vertices.begin(), vertices.end(), boundaryVertices.begin(), boundaryVertices.end(), std::back_inserter(innerVertices));
    return innerVertices;
}

int Mesh::getLocalIndexOf(int vId, std::vector<unsigned long int> fVIds) {
    for (int i = 0; i < fVIds.size(); i++) {
        if (fVIds.at(i) == vId)
            return i;
    }
    return -1;
}

int Mesh::getRelativeOrientation(std::vector<size_t> vec1, std::vector<size_t> vec2) {
    int index1 = -1;
    int index2 = -1;
    for (int i = 0; i < vec1.size(); i++) {
        int num1 = vec1.at(i);
        for (int j = 0; j < vec2.size(); j++) {
            int num2 = vec2.at(j);
            if (num1 == num2) {
                index1 = i;
                index2 = j;
                break;
            }
        }
        if (index1 != -1) {
            break;
        }
    }
    int vec1Prev = vec1.at( (vec1.size() + index1 - 1) % vec1.size() );
    int vec1After = vec1.at( (index1 + 1) % vec1.size() );
    int vec2Prev = vec2.at( (vec2.size() + index2 - 1) % vec2.size() );
    int vec2After = vec2.at( (index2 + 1) % vec2.size() );
    if (vec1Prev == vec2Prev || vec1After == vec2After) {
        return 1;
    } else if (vec1Prev == vec2After || vec1After == vec2Prev) {
        return -1;
    }
    return 0;
}

std::vector<std::vector<int>> Mesh::getLocalFaceRegions() {
    std::vector<std::vector<int>> result;
    for (int i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (v.isCorner && v.N_Fids.size() > 1) {
            std::vector<int> localFaceRegion;
            for (int j = 0; j < v.N_Fids.size(); j++) {
                localFaceRegion.push_back(v.N_Fids.at(j));
            }
            expandFaceRegionByOneLayer(localFaceRegion);

            if (getMinimumScaledJacobian(localFaceRegion) <= 0 && !hasUnsolvableInversion(localFaceRegion))
                result.push_back(localFaceRegion);
        }
    }
    for (int i = 0; i < F.size(); i++) {
        Face& f = F.at(i);
        if (getMinimumScaledJacobian(f) <= 0 && !hasUnsolvableInversion(f)) {
            std::vector<int> localFaceRegion;
            localFaceRegion.push_back(f.id);
            expandFaceRegionByOneLayer(localFaceRegion);
            result.push_back(localFaceRegion);
        }
    }

    result = mergeRegions(result);
    for (int i = 0; i < result.size(); i++) {
        expandFaceRegionByOneLayer(result.at(i));
    }

    for (int i = result.size() - 1; i >= 0; i--) {
        std::vector<int>& localFaceRegion = result.at(i);
        removeUnsolvableFaces(localFaceRegion);
        if (localFaceRegion.size() == 0)
            result.erase(result.begin() + i);
    }

    // std::cout << "-----\n";
    // std::cout << "Local Regions:\n";
    // for (int i = 0; i < result.size(); i++) {
    //     std::vector<int> localRegion = result.at(i);
    //     for (int j = 0; j < localRegion.size(); j++) {
    //         std::cout << localRegion.at(j) << ", ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "-----\n";

    return result;
}

double Mesh::getIterativeEnergyOfRegion(std::vector<int> faceIds) {
    double minimumSJ = 1.0;
    double sumOfInverted = 0.0;
    for (int i = 0; i < faceIds.size(); i++) {
        Face& f = F.at(faceIds.at(i));
        double localMinSJ = getMinimumScaledJacobian(f);
        if (localMinSJ < minimumSJ) {
            minimumSJ = localMinSJ;
        }
        if (localMinSJ < 0)
            sumOfInverted += localMinSJ;
    }

    return minimumSJ <= 0 ? sumOfInverted : minimumSJ;
}

std::vector<std::vector<int>> Mesh::mergeRegions(std::vector<std::vector<int>> localFaceRegions) {
    for (int i = 0; i < localFaceRegions.size(); i++) {
        std::vector<int>& localFaceRegion = localFaceRegions.at(i);
        std::sort(localFaceRegion.begin(), localFaceRegion.end());
    }
    std::vector<std::vector<int>> result;
    for (int i = 0; i < localFaceRegions.size(); i++) {
        std::vector<int>& currentFaceRegion = localFaceRegions.at(i);
        while (true) {
            bool merged = false;
            for (int j = localFaceRegions.size() - 1; j > i; j--) {
                std::vector<int> otherFaceRegion = localFaceRegions.at(j);
                std::vector<int> faceDifference;
                std::set_intersection(otherFaceRegion.begin(), otherFaceRegion.end(), currentFaceRegion.begin(), currentFaceRegion.end(), std::back_inserter(faceDifference));
                if (faceDifference.size() > 0) {
                    merged = true;
                    for (int k = 0; k < faceDifference.size(); k++) {
                        currentFaceRegion.push_back(faceDifference.at(k));
                    }
                    localFaceRegions.erase(localFaceRegions.begin() + j);
                }
            }
            if (!merged)
                break;
        }
        result.push_back(currentFaceRegion);
    }
    return result;
}

void Mesh::expandFaceRegionByOneLayer(std::vector<int>& localFaceRegion) {
    std::vector<int> newFacesToAdd;
    for (int k = 0; k < F.size(); k++) {
        Face& f = F.at(k);
        if (std::find(localFaceRegion.begin(), localFaceRegion.end(), f.id) != localFaceRegion.end())
            continue;

        bool neighborsFaceRegion = false;
        for (int l = 0; l < f.N_Fids.size(); l++) {
            Face& f1 = F.at(f.N_Fids.at(l));

            bool shareEdge = false;
            for (int m = 0; m < f.Eids.size(); m++) {
                for (int n = 0; n < f1.Eids.size(); n++) {
                    if (f.Eids.at(m) == f1.Eids.at(n)) {
                        shareEdge = true;
                        break;
                    }
                }
            }
            if (!shareEdge)
                continue;

            if (std::find(localFaceRegion.begin(), localFaceRegion.end(), f1.id) != localFaceRegion.end()) {
                neighborsFaceRegion = true;
                break;
            }
        }
        if (neighborsFaceRegion)
            newFacesToAdd.push_back(f.id);
    }
    for (int j = 0; j < newFacesToAdd.size(); j++) {
        localFaceRegion.push_back(newFacesToAdd.at(j));
    }
}

std::vector<int> Mesh::collapseFacesToVerticies(std::vector<int> fIds) {
    std::vector<int> vIds;
    for (int i = 0; i < fIds.size(); i++) {
        Face& f = F.at(fIds.at(i));
        for (int j = 0; j < f.Vids.size(); j++) {
            vIds.push_back(f.Vids.at(j));
        }
    }
    
    std::sort(vIds.begin(), vIds.end());
    auto last = std::unique(vIds.begin(), vIds.end());
    vIds.erase(last, vIds.end());
    return vIds;
}

void Mesh::removeUnsolvableFaces(std::vector<int>& localFaceRegion) {
    for (int j = localFaceRegion.size() - 1; j >= 0; j--) {
        Face& f = F.at(localFaceRegion.at(j));
        if (getMaximumEdgeLength(f) <= 0.001 || hasUnsolvableInversion(f)) {
            localFaceRegion.erase(localFaceRegion.begin() + j);
        }
    }
}

bool Mesh::hasUnsolvableInversion(std::vector<int> fIds) {
    for (int i = 0; i < fIds.size(); i++) {
        Face& f = F.at(fIds.at(i));
        if (hasUnsolvableInversion(f))
            return true;
    }
    return false;
}

bool Mesh::hasUnsolvableInversion(Face& f) {
    for (int i = 0; i < 4; i++) {
        Vertex& v0 = V.at(f.Vids.at(i));
        Vertex& v1 = V.at(f.Vids.at( (i + 1) % 4 ));
        Vertex& v2 = V.at(f.Vids.at( (i + 2) % 4 ));
        if (v0.isBoundary && v1.isBoundary && v2.isBoundary) {
            double currentScaledJacobian = glm::cross(glm::normalize(v2 - v1), glm::normalize(v0 - v1)).z;
            if (currentScaledJacobian <= 0)
                return true;
        }
    }
    return false;
}
