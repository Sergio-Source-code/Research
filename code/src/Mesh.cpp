/*
 * Mesh.cpp
 *
 *  Created on: Nov 6, 2016
 *      Author: cotrik
 */

#include "Mesh.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
//#include "MeshQuality.h"
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

//#include <eigen3/Eigen/Core>
//#include <eigen3/Eigen/Eigen>
//#include <eigen3/Eigen/Dense>
//#include <eigen3/Eigen/Sparse>
//#include <eigen3/Eigen/Cholesky>
//#include <eigen3/Eigen/LU>
//#include <eigen3/Eigen/SVD>
//#include <eigen3/Eigen/SparseCore>
//#include <eigen3/Eigen/SparseLU>
//#include <eigen3/Eigen/SparseQR>
//#include <eigen3/Eigen/SparseCholesky>
//using namespace Eigen;

const size_t MAXID = 0xffffffffffffffff;
int getIndexOf1(int num, std::vector<unsigned long int> vec);
std::vector<int> intersectVectors(std::vector<unsigned long int> vec1, std::vector<unsigned long int> vec2);
int findRelativeOrientation(std::vector<size_t> vec1, std::vector<size_t> vec2);

Mesh::Mesh()
: m_cellType(HEXAHEDRA)
, avgEdgeLength(0.0)
, numOfSharpEdges(0)
{
    // TODO Auto-generated constructor stub

}

Mesh::~Mesh()
{
    // TODO Auto-generated destructor stub
    V.clear();
    E.clear();
    F.clear();
    C.clear();

    pointScalarFieldNames.clear();
    pointScalarFields.clear();
    cellScalarNameFields.clear();
    cellScalarFields.clear();
}

Mesh::Mesh(const Mesh& r)
: V(r.V)
, E(r.E)
, F(r.F)
, C(r.C)
, m_cellType(r.m_cellType)
, pointScalarFieldNames(r.pointScalarFieldNames)
, pointScalarFields(r.pointScalarFields)
, cellScalarNameFields(r.cellScalarNameFields)
, cellScalarFields(r.cellScalarFields)
, avgEdgeLength(r.avgEdgeLength)
, numOfSharpEdges(r.numOfSharpEdges)
{
    m_refIds.resize(V.size());
    for (size_t i = 0; i < V.size(); ++i) m_refIds[i] = i;
}

Mesh::Mesh(const std::vector<Vertex>& V, const std::vector<Cell>& C, ElementType m_cellType)
: V(V)
, C(C)
, m_cellType(m_cellType)
, avgEdgeLength(0.0)
, numOfSharpEdges(0)
{
    m_refIds.resize(V.size());
    for (size_t i = 0; i < V.size(); ++i) m_refIds[i] = i;
}

Mesh::Mesh(const std::vector<Vertex>& V, const std::vector<Face>& F, ElementType m_cellType)
: V(V)
, F(F)
, m_cellType(m_cellType)
, avgEdgeLength(0.0)
, numOfSharpEdges(0)
{
    m_refIds.resize(V.size());
    for (size_t i = 0; i < V.size(); ++i) m_refIds[i] = i;
}

Mesh::Mesh(const Mesh& r, const std::vector<size_t>& cellIds)
: m_cellType(r.m_cellType)
, avgEdgeLength(0.0)
, numOfSharpEdges(0)
{
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

    m_refIds.resize(V.size());
    for (size_t i = 0; i < V.size(); ++i) m_refIds[i] = i;
}

void Mesh::BuildE()
{
    std::vector<Edge> Es;
    for (size_t i = 0; i < C.size(); i++)
    {
        Edge e(2);
        if (m_cellType == HEXAHEDRA)
            for (size_t j = 0; j < 12; j++)  // there are 12 edges in hex
            {
                e.Vids[0] = C[i].Vids[HexEdge[j][0]];
                e.Vids[1] = C[i].Vids[HexEdge[j][1]];
                Es.push_back(e);
            }
        else if (m_cellType == QUAD)
            for (size_t j = 0; j < 4; j++)  // there are 3 edges in tri
            {
                e.Vids[0] = C[i].Vids[QuadEdge[j][0]];
                e.Vids[1] = C[i].Vids[QuadEdge[j][1]];
                Es.push_back(e);
            }
        else if (m_cellType == TRIANGLE)
            for (size_t j = 0; j < 3; j++)  // there are 3 edges in tri
            {
                e.Vids[0] = C[i].Vids[TriEdge[j][0]];
                e.Vids[1] = C[i].Vids[TriEdge[j][1]];
                Es.push_back(e);
            }
        else if (m_cellType == TETRAHEDRA)
            for (size_t j = 0; j < 6; j++)  // there are 6 edges in tet
            {
                e.Vids[0] = C[i].Vids[TetEdge[j][0]];
                e.Vids[1] = C[i].Vids[TetEdge[j][1]];
                Es.push_back(e);
            }
    }

    size_t E_N = 0;
    for (size_t i = 0; i < Es.size(); i++)
    {
        bool havesame = false;
        const size_t id1 = Es[i].Vids[0];
        const size_t id2 = Es[i].Vids[1];

        for (size_t j = 0; j < V[id1].N_Vids.size(); j++)
        {
            if (V[id1].N_Vids[j] == id2)
            {
                havesame = true;
                break;
            }
        }
        if (!havesame)
        {
            Es[i].id = E_N++;
            Es[i].isBoundary = false;

            for (size_t j = 0; j < V[Es[i].Vids[0]].N_Cids.size(); j++)
            {
                size_t h1 = V[Es[i].Vids[0]].N_Cids[j];
                for (size_t k = 0; k < V[Es[i].Vids[1]].N_Cids.size(); k++)
                {
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

    for (size_t i = 0; i < E.size(); i++)
    {
        for (size_t j = 0; j < V[E[i].Vids[0]].N_Cids.size(); j++)
        {
            size_t nh = V[E[i].Vids[0]].N_Cids[j];
            for (size_t k = 0; k < V[E[i].Vids[1]].N_Cids.size(); k++)
                if (nh == V[E[i].Vids[1]].N_Cids[k])
                    E[i].N_Cids.push_back(nh);
        }
    }

    for (size_t i = 0; i < E.size(); i++)
    {
        std::vector<size_t> N_Cids;
        for (size_t j = 0; j < E[i].N_Cids.size(); j++)
        {
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

    for (size_t i = 0; i < V.size(); i++)
    {
        std::vector<size_t> N_Vids;
        for (int j = 0; j < V[i].N_Vids.size(); j++)
        {
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

void set_redundent_clearn(std::vector<size_t>& set)
{
    std::vector<size_t> set_copy;
    for (int i = 0; i < set.size(); i++)
    {
        bool have = false;
        for (size_t j = i + 1; j < set.size(); j++)
            if (set[i] == set[j])
                have = true;
        if (!have)
            set_copy.push_back(set[i]);
    }
    set = set_copy;
}

bool set_contain(std::vector<size_t>& large_set, size_t element)
{
    for (size_t j = 0; j < large_set.size(); j++)
        if (element == large_set[j])
            return true;

    return false;
}

void set_exclusion(std::vector<size_t>& large_set, std::vector<size_t>& small_set, std::vector<size_t> &result_set)
{
    result_set.clear();
    for (size_t i = 0; i < large_set.size(); i++)
    {
        bool inside = false;
        for (size_t j = 0; j < small_set.size(); j++)
        {
            if (small_set[j] == large_set[i])
            {
                inside = true;
                break;
            }
        }
        if (!inside)
            result_set.push_back(large_set[i]);
    }
}

void Mesh::BuildParallelE()
{
    for (size_t i = 0; i < F.size(); i++)
        for (size_t j = 0; j < F[i].N_Ortho_4Eids.size(); j++)
            for (size_t k = 0; k < 4; k++)
                for (size_t m = 0; m < 4; m++)
                    if (k != m)
                        E[F[i].N_Ortho_4Eids[j][k]].parallelEids.push_back(F[i].N_Ortho_4Eids[j][m]);

    for (size_t i = 0; i < E.size(); i++)
        set_redundent_clearn(E[i].parallelEids);
}

void Mesh::BuildConsecutiveE()
{
    for (size_t i = 0; i < E.size(); i++)
    {
        std::vector<size_t> fes;
        for (size_t j = 0; j < E[i].N_Fids.size(); j++)
        {
            size_t fid = E[i].N_Fids[j];
            for (size_t k = 0; k < 4; k++)
                fes.push_back(F[fid].Eids[k]);
        }
        set_redundent_clearn(fes);
        for (size_t j = 0; j < 2; j++)
        {
            int vid = E[i].Vids[j];
            std::vector<size_t> leftes;
            set_exclusion(V[vid].N_Eids, fes, leftes);
            std::copy(leftes.begin(), leftes.end(), back_inserter(E[i].consecutiveEids));
        }
        set_redundent_clearn(E[i].consecutiveEids);
    }
}

void Mesh::BuildOrthogonalE()
{
    for (size_t i = 0; i < E.size(); i++)
    {
        Edge& edge = E.at(i);
        std::vector<size_t> faceEdges;
        for (size_t j = 0; j < E[i].N_Fids.size(); j++)
        {
            size_t fid = E[i].N_Fids[j];
            for (size_t k = 0; k < 4; k++)
                faceEdges.push_back(F[fid].Eids[k]);
        }
        set_redundent_clearn(faceEdges);
        size_t vid1 = edge.Vids[0];
        size_t vid2 = edge.Vids[1];
        for (size_t j = 0; j < faceEdges.size(); j++)
        {
            bool hasVid1 = false;
            bool hasVid2 = false;
            const Edge& faceEdge = E.at(faceEdges.at(j));
            for (size_t k = 0; k < faceEdge.Vids.size(); k++){
                if (faceEdge.Vids.at(k) == vid1)
                    hasVid1 = true;
                if (faceEdge.Vids.at(k) == vid2)
                    hasVid2 = true;
            }
            if (hasVid1 ^ hasVid2)
                edge.orthogonalEids.push_back(faceEdge.id);
        }
    }
}

void Mesh::GetNormalOfSurfaceFaces()
{
    for (size_t i = 0; i < F.size(); i++) {
        Face& face = F.at(i);
        const Vertex& v0 = V[face.Vids[0]];
        const Vertex& v1 = V[face.Vids[1]];
        const Vertex& v2 = V[face.Vids[2]];

        const glm::vec3 v10 = v0.xyz() - v1.xyz();
        const glm::vec3 v12 = v2.xyz() - v1.xyz();
        const glm::vec3 n = glm::normalize(glm::cross(v12, v10));
        face.normal = n;
    }
}

void Mesh::GetNormalOfSurfaceVertices()
{
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        v.normal = glm::vec3(0.0, 0.0, 0.0);
        if (!v.isBoundary)
            continue;
        glm::vec3 sumNormal(0.0, 0.0, 0.0);
        size_t faceCount = 0;
        for (size_t j = 0; j < v.N_Fids.size(); j++) {
            const Face& face = F.at(v.N_Fids.at(j));
            if (face.isBoundary) {
                faceCount++;
                v.normal.x += face.normal.x;
                v.normal.y += face.normal.y;
                v.normal.z += face.normal.z;
            }
        }
        v.normal = glm::normalize(v.normal);
    }
}

#define CROSSVECTOR3(a,b,c)       {(a)[0]=(b)[1]*(c)[2]-(b)[2]*(c)[1]; \
    (a)[1]=(b)[2]*(c)[0]-(b)[0]*(c)[2]; \
    (a)[2]=(b)[0]*(c)[1]-(b)[1]*(c)[0];}
const double PAI = 3.1415926535898;
void Mesh::ClassifyVertexTypes()
{
    std::vector<float> face_angles;
    for (size_t i = 0; i < E.size(); i++)
    {
        const Edge& edge = E.at(i);
        size_t v1 = edge.Vids[0];
        size_t v2 = edge.Vids[1];

        std::vector<size_t> the_other_vs;
        for (size_t j = 0; j < edge.N_Fids.size(); j++) {
            const Face& face = F.at(edge.N_Fids.at(j));
            if (face.isBoundary) for (size_t k = 0; k < 3; k++) {
                const size_t vid = face.Vids[k];
                if (vid != v1 && vid != v2) the_other_vs.push_back(vid);
            }
        }
        std::vector<size_t> V1s, V2s;
        V1s.push_back(v1);
        V1s.push_back(v2);
        V1s.push_back(the_other_vs[0]);

        V2s.push_back(v1);
        V2s.push_back(v2);
        V2s.push_back(the_other_vs[1]);

        double edge1vector1[3], edge2vector1[3], edge1vector2[3], edge2vector2[3];

        for (size_t j = 0; j < 3; j++)
        {
            edge1vector1[j] = V[V1s[1]][j] - V[V1s[0]][j];
            edge2vector1[j] = V[V1s[2]][j] - V[V1s[0]][j];

            edge1vector2[j] = V[V2s[1]][j] - V[V2s[0]][j];
            edge2vector2[j] = V[V2s[2]][j] - V[V2s[0]][j];
        }

        glm::vec3 normal1, normal2;
        CROSSVECTOR3(normal1, edge1vector1, edge2vector1);
        CROSSVECTOR3(normal2, edge2vector2, edge1vector2);
        normal1 = glm::normalize(normal1);
        normal2 = glm::normalize(normal2);

        double cos_angle = normal1[0] * normal2[0] + normal1[1] * normal2[1] + normal1[2] * normal2[2];
        if (cos_angle < -1.0)
            cos_angle = -1.0;
        else if (cos_angle > 1.0)
            cos_angle = 1.0;

        double face_angle = acos(cos_angle);

        if (face_angle > PAI)
            face_angle = PAI;
        else if (face_angle < 0)
            face_angle = 0;

        face_angles.push_back(PAI - face_angle);

        E[i].face_angle = PAI - face_angle;
    }
    std::vector<std::vector<size_t> > flags(V.size());
    const double angle_threshold = 170.0 / 180 * PAI;
    for (size_t i = 0; i < E.size(); i++) {
        if (E[i].face_angle < angle_threshold) {
            flags[E[i].Vids[0]].push_back(i);
            flags[E[i].Vids[1]].push_back(i);
        }
    }
    for (size_t i = 0; i < flags.size(); i++) {
        if (flags[i].size() == 0)
            V[i].type = 0;
        else if (flags[i].size() <= 2) {
            V[i].type = 1;
            std::vector<size_t> vs_order;
            for (size_t j = 0; j < flags[i].size(); j++)
            {
                if (j == 0)
                {
                    size_t v1 = E[flags[i][j]].Vids[0];
                    size_t v2 = E[flags[i][j]].Vids[1];
                    if (v1 == i) {
                        vs_order.push_back(v2);
                        vs_order.push_back(v1);
                    }
                    else {
                        vs_order.push_back(v1);
                        vs_order.push_back(v2);
                    }
                }
                else
                {
                    int v1 = E[flags[i][j]].Vids[0];
                    int v2 = E[flags[i][j]].Vids[1];
                    if (v1 == i) {
                        vs_order.push_back(v1);
                        vs_order.push_back(v2);
                    }
                    else {
                        vs_order.push_back(v2);
                        vs_order.push_back(v1);
                    }
                }
            }
            float length_all = 0;
            for (size_t j = 1; j < vs_order.size(); j++)
            {
                const glm::vec3 dir = V[vs_order[j - 1]] - V[vs_order[j]];
                const float dis = glm::length(dir);
                V[i].tangent += glm::vec3(dis * dir.x, dis * dir.y, dis * dir.z);
            }
            //vector_normalization(tmi_sur.Vs[i].tangent);
        }
        else if (flags[i].size() > 2)
            V[i].type = 2; //tmi_sur.Vs[i].type=1;//
    }
}
void Mesh::ClearLabelOfSurface()
{
    for (size_t i = 0; i < F.size(); i++) {
        Face& face = F.at(i);
        face.label = MAXID;
    }
    for (size_t i = 0; i < E.size(); i++) {
        Edge& edge = E.at(i);
        edge.isSharpFeature = false;
    }
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        v.type = MAXID;
        v.isCorner = false;
        v.tangent = glm::vec3(0.0,0.0,0.0);
    }
}
void Mesh::LabelSurface()
{
    size_t label = 0;
    for (size_t i = 0; i < F.size(); i++) {
        Face& face = F.at(i);
        if (!face.isBoundary || face.label != MAXID)
            continue;
        LabelFace(face, label);
        label++;
    }

    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (v.isBoundary)
            v.type = REGULAR;
    }
    for (size_t i = 0; i < E.size(); i++) {
        Edge& edge = E.at(i);
        if (!edge.isBoundary)
            continue;
        Face* pFace1 = NULL;
        for (size_t j = 0; j < edge.N_Fids.size(); j++) {
            Face& face = F.at(edge.N_Fids.at(j));
            if (face.isBoundary)
                pFace1 = &face;
        }
        Face* pFace2 = NULL;
        for (size_t j = 0; j < edge.N_Fids.size(); j++) {
            Face& face = F.at(edge.N_Fids.at(j));
            if (face.isBoundary && face.id != pFace1->id)
                pFace2 = &face;
        }

        const double cos_angle = GetCosAngle(edge, *pFace1, *pFace2);
        // cos(15) = 0.9659 cos(30) = 0.866 cos(45) = 0.707

        if (pFace1->label != pFace2->label
                || fabs(cos_angle) < 0.5
            )
        {
            edge.isSharpFeature = true;
            V.at(edge.Vids[0]).type = FEATURE;
            V.at(edge.Vids[1]).type = FEATURE;
        }
    }

    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isBoundary)
            continue;
        size_t numOfSharpEdges = 0;
        for (size_t j = 0; j < v.N_Eids.size(); j++) {
            const Edge& edge = E.at(v.N_Eids.at(j));
            if (edge.isBoundary && edge.isSharpFeature)
                numOfSharpEdges++;
        }
        if (numOfSharpEdges >= 3) {
            v.isCorner = true;
            v.type = CORNER;
        }
        else if (numOfSharpEdges < 3) {
            std::vector<size_t> vs_order;
            numOfSharpEdges = 0;
            for (size_t j = 0; j < v.N_Eids.size(); j++) {
                const Edge& edge = E.at(v.N_Eids.at(j));
                if (edge.isBoundary && edge.isSharpFeature) {
                    size_t v1 = edge.Vids[0];
                    size_t v2 = edge.Vids[1];
                    if (numOfSharpEdges == 0) {
                        if (v1 == v.id) {
                            vs_order.push_back(v2);
                            vs_order.push_back(v1);
                        } else {
                            vs_order.push_back(v1);
                            vs_order.push_back(v2);
                        }
                    } else {
                        if (v1 == i) {
                            vs_order.push_back(v1);
                            vs_order.push_back(v2);
                        } else {
                            vs_order.push_back(v2);
                            vs_order.push_back(v1);
                        }
                    }
                    numOfSharpEdges++;
                }
            }
            for (size_t j = 1; j < vs_order.size(); j++) {
                const glm::vec3 dir = V[vs_order[j - 1]] - V[vs_order[j]];
                const float dis = glm::length(dir);
                v.tangent += glm::vec3(dis * dir.x, dis * dir.y, dis * dir.z);
            }
        }
    }
}

void Mesh::LabelSharpEdges(const bool breakAtConrer/* = false*/)
{
    size_t label = 0;
    for (size_t i = 0; i < E.size(); i++) {
        Edge& edge = E.at(i);
        if (!edge.isBoundary || !edge.isSharpFeature || edge.label != MAXID)
            continue;
        LabelEdge(edge, label, breakAtConrer);
        label++;
    }
    this->numOfSharpEdges = label;
}

void Mesh::LabelEdge(Edge& edge, size_t& label, const bool breakAtConrer/* = false*/)
{
    edge.label = label;
    for (size_t i = 0; i < edge.Vids.size(); i++) {
        const Vertex& v = V.at(edge.Vids.at(i));
        if (v.isCorner && breakAtConrer)
            continue;
        std::vector<Edge*> edges;
        for (size_t j = 0; j < v.N_Eids.size(); j++) {
            Edge& edge2 = E.at(v.N_Eids.at(j));
            if (edge2.id != edge.id && edge2.isBoundary && edge.isSharpFeature && edge2.label == MAXID) {
                edges.push_back(&edge2);
            }
        }
        for (size_t j = 0; j < edges.size(); j++)
            LabelEdge(*edges.at(j), label, breakAtConrer);
    }
}


//double Mesh::GetCosAngle(const Edge& edge, const Face& face1, const Face& face2)
//{
//    /*              v0
//     *             /|\
//     *            / | \
//     *      face1/  |  \face2
//     *          /   |   \
//     *         /    |    \
//     *        -------------
//     *       v2     v1     v3
//     * */
//    const Vertex& v1 = V.at(edge.Vids[0]);
//    const Vertex& v0 = V.at(edge.Vids[1]);
//    size_t v2id = MAXID;
//    for (size_t i = 0; i < face1.Eids.size(); i++) {
//        const Edge& edge1 = E.at(face1.Eids.at(i));
//        if (edge == edge1)
//            continue;
//        if (edge1.Vids[0] == v1.id)
//            v2id = edge1.Vids[1];
//        else if (edge1.Vids[1] == v1.id)
//            v2id = edge1.Vids[0];
//        if (v2id != MAXID)
//            break;
//    }
//    const Vertex& v2 = V.at(v2id);
//
//    size_t v3id = MAXID;
//    for (size_t i = 0; i < face2.Eids.size(); i++) {
//        const Edge& edge2 = E.at(face2.Eids.at(i));
//        if (edge == edge2)
//            continue;
//        if (edge2.Vids[0] == v1.id)
//            v3id = edge2.Vids[1];
//        else if (edge2.Vids[1] == v1.id)
//            v3id = edge2.Vids[0];
//        if (v3id != MAXID)
//            break;
//    }
//    const Vertex& v3 = V.at(v3id);
//    Plane plane1(v2, v1, v0);
//    const Plane plane2(v3, v0, v1);
//
//    return plane1.IntersectionAngle(plane2);
//}

double Mesh::GetCosAngle(const Edge& edge, const Face& face1, const Face& face2)
{
    /*              v0
     *             /|\
     *            / | \
     *      face1/  |  \face2
     *          /   |   \
     *         /    |    \
     *        -------------
     *       v2     v1     v3
     * */
    Vertex* pv0 = &V.at(edge.Vids[0]);
    Vertex* pv1 = &V.at(edge.Vids[1]);
    bool correct_orientation = false;
    const size_t n = face1.Vids.size();
    for (size_t i = 0; i < n; i++) {
        const Vertex& v0 = V.at(face1.Vids[i % n]);
        const Vertex& v1 = V.at(face1.Vids[(i+1) % n]);
        if (v0.id == pv0->id && v1.id == pv1->id) {
            correct_orientation = true;
            break;
        }
    }
    if (!correct_orientation) {
        pv1 = &V.at(edge.Vids[0]);
        pv0 = &V.at(edge.Vids[1]);
    }
    const Vertex& v0 = *pv0;
    const Vertex& v1 = *pv1;

    size_t v2id = MAXID;
    for (size_t i = 0; i < face1.Eids.size(); i++) {
        const Edge& edge1 = E.at(face1.Eids.at(i));
        if (edge == edge1)
            continue;
        if (edge1.Vids[0] == v1.id)
            v2id = edge1.Vids[1];
        else if (edge1.Vids[1] == v1.id)
            v2id = edge1.Vids[0];
        if (v2id != MAXID)
            break;
    }
    const Vertex& v2 = V.at(v2id);

    size_t v3id = MAXID;
    for (size_t i = 0; i < face2.Eids.size(); i++) {
        const Edge& edge2 = E.at(face2.Eids.at(i));
        if (edge == edge2)
            continue;
        if (edge2.Vids[0] == v1.id)
            v3id = edge2.Vids[1];
        else if (edge2.Vids[1] == v1.id)
            v3id = edge2.Vids[0];
        if (v3id != MAXID)
            break;
    }
    const Vertex& v3 = V.at(v3id);
    Plane plane1(v2, v1, v0);
    const Plane plane2(v3, v0, v1);

    return plane1.IntersectionAngle(plane2);
}

void Mesh::SetCosAngleThreshold(const double cos_angle/* = 0.91*/)
{
    cos_angle_threshold = cos_angle;
}

void Mesh::LabelFace(Face& face, size_t& label)
{
    face.label = label;
    for (size_t i = 0; i < face.Eids.size(); i++)
    {
        const Edge& edge = E.at(face.Eids.at(i));
        std::vector<Face*> faces;
        for (size_t j = 0; j < edge.N_Fids.size(); j++) {
            Face& face2 = F.at(edge.N_Fids.at(j));
            if (face2.isBoundary && face2.id != face.id && face2.label == MAXID) {
                const double cos_angle = GetCosAngle(edge, face, face2);
                //std::cout << "cos_angle = " << cos_angle << std::endl;
                if (cos_angle > cos_angle_threshold) // cos(15) = 0.9659 cos(30) = 0.866
                    faces.push_back(&face2);
            }
        }
        for (size_t i = 0; i < faces.size(); i++)
            LabelFace(*faces.at(i), label);
    }
}

void Mesh::RemoveUselessVertices()
{
    std::vector<size_t> v_real_index;
    size_t c_size = 0;
    for (size_t i = 0; i < C.size(); i++)
    {
        const Cell& cell = C.at(i);
        for (size_t j = 0; j < cell.Vids.size(); j++)
            v_real_index.push_back(cell.Vids.at(j));
        c_size++;
    }

    std::sort(v_real_index.begin(), v_real_index.end());
    std::vector<size_t>::iterator iter = std::unique(v_real_index.begin(), v_real_index.end());
    v_real_index.resize(std::distance(v_real_index.begin(), iter));

    std::vector<Vertex> newV(v_real_index.size());
    for (size_t i = 0; i < v_real_index.size(); i++)
    {
        const Vertex& v = V.at(v_real_index.at(i));
        Vertex& newv = newV.at(i);
        newv.id = i;
        newv.x = v.x;
        newv.y = v.y;
        newv.z = v.z;
    }
    V = newV;
    m_refIds = v_real_index;
    //////////////////////////////////////////////////////
    std::map<size_t, size_t> v_v;
    size_t index = 0;
    for (size_t i = 0; i < v_real_index.size(); i++)
    {
        v_v[v_real_index.at(i)] = i;
    }

    std::vector<Cell> newC;
    for (size_t i = 0; i < C.size(); i++)
    {
        Cell c;
        if (m_cellType == HEXAHEDRA) c.Vids.resize(8);
        else if (m_cellType == TETRAHEDRA) c.Vids.resize(4);
        else if (m_cellType == QUAD) c.Vids.resize(4);
        else if (m_cellType == TRIANGLE) c.Vids.resize(3);
        const Cell& cell = C.at(i);
        for (size_t j = 0; j < C.at(i).Vids.size(); j++) {
            c.Vids.at(j) = v_v[cell.Vids.at(j)];
            c.id = i;
        }
        newC.push_back(c);
    }
    C.resize(newC.size());
    C = newC;

    if (m_cellType == QUAD) {
        std::vector<Face> newF;
        for (size_t i = 0; i < F.size(); i++)
        {
            Face f;
            if (m_cellType == HEXAHEDRA) f.Vids.resize(8);
            else if (m_cellType == TETRAHEDRA) f.Vids.resize(4);
            else if (m_cellType == QUAD) f.Vids.resize(4);
            else if (m_cellType == TRIANGLE) f.Vids.resize(3);
            const Face& cell = F.at(i);
            for (size_t j = 0; j < F.at(i).Vids.size(); j++) {
                f.Vids.at(j) = v_v[cell.Vids.at(j)];
                f.id = i;
            }
            newF.push_back(f);
        }
        F.resize(newC.size());
        F = newF;
    }
}

bool Mesh::IsPointInTriangle(const glm::vec3& P, const glm::vec3& A, const glm::vec3& B, const glm::vec3& C) const
{
    // Compute vectors
    const glm::vec3 v0 = C - A;
    const glm::vec3 v1 = B - A;
    const glm::vec3 v2 = P - A;

    // Compute dot products
    const double dot00 = glm::dot(v0, v0);
    const double dot01 = glm::dot(v0, v1);
    const double dot02 = glm::dot(v0, v2);
    const double dot11 = glm::dot(v1, v1);
    const double dot12 = glm::dot(v1, v2);

    // Compute barycentric coordinates
    const double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    const double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    const double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Check if point is in triangle
    //return (u >= 0) && (v >= 0) && (u + v <= 1.0);
    return (u >= -2e-1) && (v >= -2e-1) && (u + v <= 1.3);
}

bool Mesh::IsPointInFace(const glm::vec3& P, const Face& face) const
{
    const glm::vec3& A = V.at(face.Vids[0]).xyz();
    const glm::vec3& B = V.at(face.Vids[1]).xyz();
    const glm::vec3& C = V.at(face.Vids[2]).xyz();

    // Compute vectors
    const glm::vec3 v0 = C - A;
    const glm::vec3 v1 = B - A;
    const glm::vec3 v2 = P - A;

    // Compute dot products
    const double dot00 = glm::dot(v0, v0);
    const double dot01 = glm::dot(v0, v1);
    const double dot02 = glm::dot(v0, v2);
    const double dot11 = glm::dot(v1, v1);
    const double dot12 = glm::dot(v1, v2);

    // Compute barycentric coordinates
    const double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    const double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    const double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    // Check if point is in triangle
    return (u >= -2e-1) && (v >= -2e-1) && (u + v <= 1.3);
    //return (u >= -1e-3) && (v >= -1e-3) && (u + v <= 1.001);
}

//bool Mesh::IsPointInTriangle(const Vertex& p, const Vertex& p0, const Vertex& p1, const Vertex& p2) const
//{
//    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3, 3);
//    Eigen::VectorXd b = Eigen::VectorXd::Zero(3);
//    A(0,0) = p0.x; A(0,1) = p1.x; A(0,2) = p2.x;
//    A(1,0) = p0.y; A(1,1) = p1.y; A(1,2) = p2.y;
//    A(2,0) = p0.z; A(2,1) = p1.z; A(2,2) = p2.z;
//    b(0) = p.x; b(1) = p.y; b(2) = p.z;
//    Eigen::VectorXd x = A.ldlt().solve(b);
//
//    if (x(0) > -1e-3 && x(1) > -1e-3 && x(2) > -1e-3 && (x(0) + x(1) + x(2)) < 1.001)
//    {
//        return true;
//    }
//
//    return false;
//}
//
//double Area(const glm::vec3& p0, const glm::vec3& p1, const glm::vec3& p2)
//{
//    const glm::vec3 v10 = p0 - p1;
//    const glm::vec3 v12 = p2 - p1;
//    const glm::vec3 normal = glm::cross(v12, v10);
//    return 0.5f * glm::length(normal);
//}
//
//bool Mesh::IsPointInTriangle(const glm::vec3& p, const glm::vec3& p0, const glm::vec3& p1, const glm::vec3& p2) const
//{
////    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(4, 3);
////    Eigen::VectorXd b = Eigen::VectorXd::Zero(4);
////    A(0,0) = p0.x; A(0,1) = p1.x; A(0,2) = p2.x;
////    A(1,0) = p0.y; A(1,1) = p1.y; A(1,2) = p2.y;
////    A(2,0) = p0.z; A(2,1) = p1.z; A(2,2) = p2.z;
////    A(3,0) = 1.0 ; A(3,1) = 1.0 ; A(3,2) = 1.0 ;
////    b(0) = p.x; b(1) = p.y; b(2) = p.z; b(3) = 1.0;
////    Eigen::MatrixXd AT = A.transpose();
////    Eigen::MatrixXd ATA = AT*A;
////    Eigen::VectorXd x = ATA.ldlt().solve(AT*b);
////
////    if (x(0) > -1e-1 && x(1) > -1e-1 && x(2) > -1e-1 && (x(0) + x(1) + x(2)) < 1.1)
////    {
////        return true;
////    }
//    const double a0 = Area(p, p0, p1);
//    const double a1 = Area(p, p1, p2);
//    const double a2 = Area(p, p2, p0);
//    const double a  = Area(p0, p1, p2);
//
//    if (fabs(a - (a0 + a1 + a2))/a < 2e-1)
//        return true;
//    return false;
//}

struct d_i {
    double d;
    glm::vec3 i;
    d_i(): d(0.0){}
    d_i(const double d, const glm::vec3 i)
    : d(d), i(i)
    {

    }
    d_i(const d_i& r)
    : d(r.d), i(r.i)
    {

    }
    bool operator < (const d_i& r) const
    {
        return d < r.d;
    }
};

glm::vec3 Mesh::GetProjectLocation(const glm::vec3& p) const
{
    for (size_t i = 0; i < F.size(); i++) {
        const Face& face = F.at(i);
        if (!face.isBoundary)
            continue;
        Plane plane(V[face.Vids[0]].xyz(), V[face.Vids[1]].xyz(), V[face.Vids[2]].xyz());
        const double d1 = glm::length(V[face.Vids[0]].xyz() - V[face.Vids[1]].xyz());
        const double d2 = glm::length(V[face.Vids[1]].xyz() - V[face.Vids[2]].xyz());
        const double d3 = glm::length(V[face.Vids[2]].xyz() - V[face.Vids[0]].xyz());
        const double d = 0.33333333*(d1 + d2 + d3);
        glm::vec3 intersection;
        const double distance = plane.DistanseFromPoint(p, intersection);
        if (distance < 1.0 * avgEdgeLength && distance < 1.0 * d)
            if (IsPointInFace(intersection, face))
            //if (IsPointInTriangle(intersection, V[face.Vids[0]].xyz(), V[face.Vids[1]].xyz(), V[face.Vids[2]].xyz()))
                return intersection;
    }
    //std::cerr << "\033[1;31mERROR in GetProjectLocation(" << p.x << "," << p.y << "," << p.z << ")\033[0m\n";
    return glm::vec3(0.0, 0.0, 0.0);
}
glm::vec3 Mesh::GetProjectLocationOnRefSurface(const glm::vec3& p, const Vertex& refV) const
{
    std::vector<d_i> d_is;
    for (size_t i = 0; i < refV.twoRingNeighborSurfaceFaceIds.size(); i++) {
        const Face& face = F.at(refV.twoRingNeighborSurfaceFaceIds[i]);
        if (!face.isBoundary)
            continue;
        for (size_t j = 0; j < face.Vids.size(); j++) {
            Plane plane(V[face.Vids[(0 + j) % 4]].xyz(), V[face.Vids[(1 + j) % 4]].xyz(), V[face.Vids[(2 + j) % 4]].xyz());
            const double d1 = glm::length(V[face.Vids[0]].xyz() - V[face.Vids[1]].xyz());
            const double d2 = glm::length(V[face.Vids[1]].xyz() - V[face.Vids[2]].xyz());
            const double d3 = glm::length(V[face.Vids[2]].xyz() - V[face.Vids[0]].xyz());
            const double d = 0.33333333 * (d1 + d2 + d3);
            glm::vec3 intersection;
            const double distance = plane.DistanseFromPoint(p, intersection);
            if (distance < 1.0 * avgEdgeLength && distance < 1.0 * d)
                if (IsPointInFace(intersection, face))
                    d_is.push_back(d_i(distance, intersection));
//            if (face.Vids.size() == 3)
//                break;
        }
    }
    if (d_is.empty()) {
        //std::cerr << "\033[1;31mERROR in GetProjectLocationOfVertex " << p.id << "(" << p.x << "," << p.y << "," << p.z << ")\033[0m\n";
        return glm::vec3(0.0, 0.0, 0.0);
    }

    std::sort(d_is.begin(), d_is.end());
    //std::vector<d_i>::iterator iter = std::unique(d_is.begin(), d_is.end());
    return d_is[0].i;
}

glm::vec3 Mesh::GetProjectLocation(const Vertex& p) const
{
    std::vector<d_i> d_is;
    for (size_t i = 0; i < F.size(); i++) {
        const Face& face = F.at(i);
        if (!face.isBoundary)
            continue;
        for (size_t j = 0; j < face.Vids.size(); j++)
        {
            Plane plane(V[face.Vids[(0 + j) % 4]].xyz(), V[face.Vids[(1 + j) % 4]].xyz(), V[face.Vids[(2 + j) % 4]].xyz());
            const double d1 = glm::length(V[face.Vids[0]].xyz() - V[face.Vids[1]].xyz());
            const double d2 = glm::length(V[face.Vids[1]].xyz() - V[face.Vids[2]].xyz());
            const double d3 = glm::length(V[face.Vids[2]].xyz() - V[face.Vids[0]].xyz());
            const double d = 0.33333333*(d1 + d2 + d3);
            glm::vec3 intersection;
            const double distance = plane.DistanseFromPoint(p, intersection);
            if (distance < 1.0 * avgEdgeLength && distance < 1.0 * d)
                if (IsPointInFace(intersection, face))
                //if (IsPointInTriangle(intersection, V[face.Vids[(0 + j) % 4]].xyz(), V[face.Vids[(1 + j) % 4]].xyz(), V[face.Vids[(2 + j) % 4]].xyz()))
                    //return intersection;
                    //if (glm::dot(p.normal, face.normal) > 0)
                    d_is.push_back(d_i(distance, intersection));

            if (face.Vids.size() == 3)
                break;
        }
    }
    if (d_is.empty()) {
        //std::cerr << "\033[1;31mERROR in GetProjectLocationOfVertex " << p.id << "(" << p.x << "," << p.y << "," << p.z << ")\033[0m\n";
        return glm::vec3(0.0, 0.0, 0.0);
    }

    std::sort(d_is.begin(), d_is.end());
    //std::vector<d_i>::iterator iter = std::unique(d_is.begin(), d_is.end());
    return d_is[0].i;
}
glm::vec3 Mesh::GetProjectLocationFast(const Vertex& p) const
{
    std::vector<d_i> d_is;
    for (size_t i = 0; i < p.twoRingNeighborSurfaceFaceIds.size(); i++) {
        const Face& face = F.at(p.twoRingNeighborSurfaceFaceIds[i]);
        if (!face.isBoundary)
            continue;
        for (size_t j = 0; j < face.Vids.size(); j++)
        {
            Plane plane(V[face.Vids[(0 + j) % 4]].xyz(), V[face.Vids[(1 + j) % 4]].xyz(), V[face.Vids[(2 + j) % 4]].xyz());
            const double d1 = glm::length(V[face.Vids[0]].xyz() - V[face.Vids[1]].xyz());
            const double d2 = glm::length(V[face.Vids[1]].xyz() - V[face.Vids[2]].xyz());
            const double d3 = glm::length(V[face.Vids[2]].xyz() - V[face.Vids[0]].xyz());
            const double d = 0.33333333*(d1 + d2 + d3);
            glm::vec3 intersection;
            const double distance = plane.DistanseFromPoint(p, intersection);
            if (distance < 1.0 * avgEdgeLength && distance < 1.0 * d)
                if (IsPointInFace(intersection, face))
                //if (IsPointInTriangle(intersection, V[face.Vids[(0 + j) % 4]].xyz(), V[face.Vids[(1 + j) % 4]].xyz(), V[face.Vids[(2 + j) % 4]].xyz()))
                    //return intersection;
                    //if (glm::dot(p.normal, face.normal) > 0)
                    d_is.push_back(d_i(distance, intersection));

            if (face.Vids.size() == 3)
                break;
        }
    }
    if (d_is.empty()) {
        //std::cerr << "\033[1;31mERROR in GetProjectLocationOfVertex " << p.id << "(" << p.x << "," << p.y << "," << p.z << ")\033[0m\n";
        return glm::vec3(0.0, 0.0, 0.0);
    }

    std::sort(d_is.begin(), d_is.end());
    //std::vector<d_i>::iterator iter = std::unique(d_is.begin(), d_is.end());
    return d_is[0].i;
}
glm::vec3 Mesh::GetProjectLocationOnTargetSurface(const glm::vec3& p, const Vertex& refV, const Mesh& targetSurfaceMesh) const
{
    std::vector<d_i> d_is;
    for (size_t i = 0; i < refV.twoRingNeighborSurfaceFaceIds.size(); i++) {
        const Face& face = targetSurfaceMesh.F.at(refV.twoRingNeighborSurfaceFaceIds[i]);
        if (!face.isBoundary)
            continue;
        for (size_t j = 0; j < face.Vids.size(); j++) {
            Plane plane(targetSurfaceMesh.V[face.Vids[(0 + j) % 4]].xyz(),
                    targetSurfaceMesh.V[face.Vids[(1 + j) % 4]].xyz(), targetSurfaceMesh.V[face.Vids[(2 + j) % 4]].xyz());
            const double d1 = glm::length(targetSurfaceMesh.V[face.Vids[0]].xyz() - targetSurfaceMesh.V[face.Vids[1]].xyz());
            const double d2 = glm::length(targetSurfaceMesh.V[face.Vids[1]].xyz() - targetSurfaceMesh.V[face.Vids[2]].xyz());
            const double d3 = glm::length(targetSurfaceMesh.V[face.Vids[2]].xyz() - targetSurfaceMesh.V[face.Vids[0]].xyz());
            const double d = 0.33333333 * (d1 + d2 + d3);
            glm::vec3 intersection;
            const double distance = plane.DistanseFromPoint(p, intersection);
            if (distance < 1.0 * avgEdgeLength && distance < 1.0 * d)
                if (targetSurfaceMesh.IsPointInFace(intersection, face))
                    d_is.push_back(d_i(distance, intersection));
//            if (face.Vids.size() == 3)
//                break;
        }
    }
    if (d_is.empty()) {
        //std::cerr << "\033[1;31mERROR in GetProjectLocationOfVertex " << p.id << "(" << p.x << "," << p.y << "," << p.z << ")\033[0m\n";
        return glm::vec3(0.0, 0.0, 0.0);
    }

    std::sort(d_is.begin(), d_is.end());
    //std::vector<d_i>::iterator iter = std::unique(d_is.begin(), d_is.end());
    return d_is[0].i;
}

void Mesh::ProjectTo(const Mesh& mesh)
{
    glm::vec3 vm(0.0, 0.0, 0.0);
    GetAvgEdgeLength();
#pragma omp parallel for
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isBoundary) continue;
        const glm::vec3 newv = mesh.GetProjectLocation(v);
        if (vm != newv) v = newv;
    }
}
void Mesh::FastProjectTo(const Mesh& mesh)
{
    glm::vec3 vm(0.0, 0.0, 0.0);
    GetAvgEdgeLength();
#pragma omp parallel for
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isBoundary) continue;
        const glm::vec3 newv = mesh.GetProjectLocationFast(v);
        if (vm != newv) v = newv;
    }
}
void Mesh::ProjectToRefMesh(const Mesh& refMesh)
{
    glm::vec3 vm(0.0, 0.0, 0.0);
    GetAvgEdgeLength();
#pragma omp parallel for
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isBoundary) continue;
        const Vertex& refV = refMesh.V.at(m_refIds[i]);
        if (!refV.isBoundary) continue;
        const glm::vec3 newv = refMesh.GetProjectLocation(refV);
        if (vm != newv) v = newv;
    }
}

void Mesh::ProjectToTargetSurface(const Mesh& refMesh, const Mesh& targetSurfaceMesh)
{
    glm::vec3 vm(0.0, 0.0, 0.0);
    GetAvgEdgeLength();
#pragma omp parallel for
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isBoundary) continue;
        const Vertex& refV = refMesh.V.at(m_refIds[i]);
        if (!refV.isBoundary) continue;
        const glm::vec3 newv = refMesh.GetProjectLocationOnTargetSurface(v.xyz(), refV, targetSurfaceMesh);
        if (vm != newv) v = newv;
    }
}

double Mesh::GetAvgEdgeLength()
{
    if (avgEdgeLength != 0)
        return avgEdgeLength;
    double sum = 0;
    size_t count = 0;
    for (size_t i = 0; i < E.size(); i++) {
        const Edge& e = E.at(i);
        if (!e.isBoundary)
            continue;
        sum += glm::length(V.at(e.Vids[0]).xyz() - V.at(e.Vids[1]).xyz());
        count++;
    }
    avgEdgeLength = sum/count;
    return avgEdgeLength;
}

glm::vec3 Mesh::GetFaceCenter(const Face& f)
{
    glm::vec3 center(0.0, 0.0, 0.0);
    for (size_t k = 0; k < f.Vids.size(); k++) {
        const Vertex& f_v = V.at(f.Vids[k]);
        center += f_v.xyz();
    }
    center.x /= f.Vids.size();
    center.y /= f.Vids.size();
    center.z /= f.Vids.size();
    return center;
}

double Mesh::SmoothVolume(const SmoothMethod smoothMethod/* = LAPLACE_EDGE*/)
{
    std::vector<Vertex> oldV = V;
    std::vector<Vertex> newV = V;
    //while (iters-- != 0)
    {
#pragma omp parallel for
        for (size_t i = 0; i < V.size(); i++) {
            const Vertex& v = V.at(i);
            if (v.isBoundary)
                continue;
            //if (v.type == REGULAR || (v.type == FEATURE && treatSharpFeatureAsRegular) || (v.type == CORNER && treatCornerAsRegular))
            {
                if (smoothMethod == LAPLACE_EDGE) {
                    newV.at(i) = LapLace(v);
                }
                else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::vec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary)
                            continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x/count;
                    newV.at(i).y = sum.y/count;
                    newV.at(i).z = sum.z/count;
                }
            }
        }

        for (size_t i = 0; i < V.size(); i++) {
            Vertex& v = V.at(i);
            if (v.isBoundary)
                continue;
            const Vertex& newv = newV.at(i);
            v = newv.xyz();
        }
    }

    double energy = 0;
    for (size_t i = 0; i < V.size(); i++) {
        const Vertex& v = V.at(i);
        if (v.isBoundary)
            continue;
        const Vertex& oldv = oldV.at(i);
        const double distance = glm::length(v.xyz() - oldv.xyz());
        energy += distance * distance;
    }

    std::cout << "Volume Energy = " << energy << std::endl;
    return energy;
}

double Mesh::SmoothSurface(size_t iters/* = 1*/, const SmoothMethod smoothMethod/* = LAPLACE_EDGE*/,
        const bool preserveSharpFeature/* = false*/, const bool treatSharpFeatureAsRegular/* = false*/, const bool treatCornerAsRegular/* = false*/)
{
    std::vector<Vertex> oldV = V;
    std::vector<Vertex> newV = V;
    while (iters-- != 0)
    {
#pragma omp parallel for
        for (size_t i = 0; i < V.size(); i++) {
            const Vertex& v = V.at(i);
            if (!v.isBoundary)
                continue;
            if (v.type == REGULAR || (v.type == FEATURE && treatSharpFeatureAsRegular) || (v.type == CORNER && treatCornerAsRegular))
            {
                if (smoothMethod == LAPLACE_EDGE) {
                    newV.at(i) = LapLace(v, treatSharpFeatureAsRegular, treatCornerAsRegular);
                }
                else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::vec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary)
                            continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x/count;
                    newV.at(i).y = sum.y/count;
                    newV.at(i).z = sum.z/count;
                }
            }
            else if (v.type == FEATURE) {
                if (preserveSharpFeature)
                    continue;
                if (smoothMethod == LAPLACE_EDGE) {
                    newV.at(i) = LapLace(v);
                }
                else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::vec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary)
                            continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x/count;
                    newV.at(i).y = sum.y/count;
                    newV.at(i).z = sum.z/count;
                }
            }
        }

        for (size_t i = 0; i < V.size(); i++) {
            Vertex& v = V.at(i);
            if (!v.isBoundary)
                continue;
            const Vertex& newv = newV.at(i);
            v = newv.xyz();
        }
    }

    double energy = 0;
    for (size_t i = 0; i < V.size(); i++) {
        const Vertex& v = V.at(i);
        if (!v.isBoundary)
            continue;
        const Vertex& oldv = oldV.at(i);
        const double distance = glm::length(v.xyz() - oldv.xyz());
        energy += distance * distance;
    }

    std::cout << "Energy = " << energy << std::endl;
    return energy;
}
double Mesh::GetMinScaledJacobian(double& avgSJ) const
{
    double minScaledJacobian = 1;
    double sum = 0;
    for (int i = 0; i < C.size(); i++) {
        const Cell& cell = C.at(i);
        float scaledJacobian = GetScaledJacobian(cell);
        minScaledJacobian = minScaledJacobian < scaledJacobian ? minScaledJacobian : scaledJacobian;
        sum += scaledJacobian;
    }
    avgSJ = sum/C.size();
    return minScaledJacobian;
}

size_t Mesh::GetQuality(const char* filename, double& minValue, double& avgValue, const double minSJ/* = 0.0*/)
{
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename);
    reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGridMesh = reader->GetOutput();
    vtkSmartPointer<vtkMeshQuality> qualityFilter = vtkSmartPointer<vtkMeshQuality>::New();
#if VTK_MAJOR_VERSION <= 5
    qualityFilter->SetInputConnection(mesh->GetProducerPort());
#else
    ///qualityFilter->SetInputData(mesh);
    qualityFilter->SetInputConnection(reader->GetOutputPort());
#endif
    minValue = 0;
    double maxValue = 0;
    avgValue = 0;
    std::vector<double> metrics;
    //qualityFilter->SetTriangleQualityMeasureToArea();
    //qualityFilter->SetHexQualityMeasureToEdgeRatio()
    //qualityFilter->SetHexQualityMeasureToMedAspectFrobenius()
    //qualityFilter->SetHexQualityMeasureToMaxAspectFrobenius()
    //qualityFilter->SetHexQualityMeasureToMaxEdgeRatios()
    //qualityFilter->SetHexQualityMeasureToSkew()
    //qualityFilter->SetHexQualityMeasureToTaper()
    //qualityFilter->SetHexQualityMeasureToVolume()
    //qualityFilter->SetHexQualityMeasureToStretch()
    //qualityFilter->SetHexQualityMeasureToDiagonal()
    //qualityFilter->SetHexQualityMeasureToDimension()
    //qualityFilter->SetHexQualityMeasureToOddy()
    //qualityFilter->SetHexQualityMeasureToCondition()
    //qualityFilter->SetHexQualityMeasureToJacobian()
    qualityFilter->SetHexQualityMeasureToScaledJacobian();
    //qualityFilter->SetHexQualityMeasureToShear()
    //qualityFilter->SetHexQualityMeasureToShape()
    //qualityFilter->SetHexQualityMeasureToRelativeSizeSquared()
    //qualityFilter->SetHexQualityMeasureToShapeAndSize()
    //qualityFilter->SetHexQualityMeasureToShearAndSize()
    //qualityFilter->SetHexQualityMeasureToDistortion()
    qualityFilter->Update();

    vtkSmartPointer<vtkDoubleArray> qualityArray = vtkDoubleArray::SafeDownCast(qualityFilter->GetOutput()->GetCellData()->GetArray("Quality"));
    minValue = qualityArray->GetValue(0);
    maxValue = qualityArray->GetValue(0);
    avgValue = 0;
    size_t numOfInvertedElements = 0;
    for (vtkIdType i = 0; i < qualityArray->GetNumberOfTuples(); i++)
    {
        double val = qualityArray->GetValue(i);
        if (minValue > val) minValue = val;
        if (maxValue < val) maxValue = val;
        avgValue += val;

        if (val < minSJ) {
            numOfInvertedElements++;
            //badCellIds.push_back(i);
        }
        else if (val > 1) numOfInvertedElements++;
        //std::cout << "value " << i << " : " << val << std::endl;
    }
    avgValue /= qualityArray->GetNumberOfTuples();
    return numOfInvertedElements;
}

#include "verdict.h"
size_t Mesh::GetQualityVerdict(double& minValue, double& avgValue, const double minSJ/* = 0.0*/)
{
    size_t numOfInvertedElements = 0;
    avgValue = 0.0;
    //badCellIds.clear();
    double minScaledJacobian = 1;
    //const std::vector<Vertex>& V = mesh.V;
    //const std::vector<Cell>& C = mesh.C;
    double coordinates[8][3];
    for (size_t i = 0; i < C.size(); i++) {
        const Cell& c = C[i];
        for (size_t j = 0; j < 8; j++) {
            const Vertex& v = V[c.Vids[j]];
            coordinates[j][0] = v.x;
            coordinates[j][1] = v.y;
            coordinates[j][2] = v.z;
        }
        double scaledJacobian = 0;// TODO: replace this v_hex_scaled_jacobian(8, coordinates);
        if (scaledJacobian < minScaledJacobian) minScaledJacobian = scaledJacobian;
        //minScaledJacobian = minScaledJacobian < scaledJacobian ? minScaledJacobian : scaledJacobian;
        if (scaledJacobian < minSJ) {
            numOfInvertedElements++;
            //badCellIds.push_back(i);
        }
        avgValue += scaledJacobian;
    }
    minValue = minScaledJacobian;
    avgValue /= C.size();

    return numOfInvertedElements;
}
// preseverQuality, so we cannot be parallel
double Mesh::SmoothAndProjectSurface(const Mesh& mesh, size_t iters/* = 1*/, const SmoothMethod smoothMethod/* = LAPLACE_EDGE*/,
        const bool preserveSharpFeature/* = false*/, const bool treatSharpFeatureAsRegular/* = false*/, const bool treatCornerAsRegular/* = false*/,
        const bool preserveQuality/* = false*/)
{
    double targetMinSJ = 0;
    double targetAvgSJ = 0;
    //double targetMinSJ = GetMinScaledJacobian(targetAvgSJ);
    //MeshFileWriter writer(*this, "temp.vtk");
    //writer.WriteFile();
    //GetQuality("temp.vtk", targetMinSJ, targetAvgSJ);
    GetQualityVerdict(targetMinSJ, targetAvgSJ);
    //GetQuality(targetMinSJ, targetAvgSJ);
    std::cout << "targetMinSJ = " << targetMinSJ << " targetAvgSJ = " << targetAvgSJ << std::endl;
    double energy = 0;
    std::vector<Vertex> oldV = V;
    std::vector<Vertex> newV = V;
    while (iters-- != 0)
    {
//#pragma omp parallel for
        for (size_t i = 0; i < V.size(); i++) {
            const Vertex& v = V.at(i);
            if (!v.isBoundary)
                continue;
            if (v.type == REGULAR || (v.type == FEATURE && treatSharpFeatureAsRegular) || (v.type == CORNER && treatCornerAsRegular))
            {
                if (smoothMethod == LAPLACE_EDGE) {
                    const glm::vec3 smoothedv = LapLace(v, treatSharpFeatureAsRegular, treatCornerAsRegular);
                    const glm::vec3 newv = mesh.GetProjectLocation(smoothedv);
                    if (!preserveQuality) {
                        V.at(i) = newv;
                    }
                    else {
                        const glm::vec3 oldv = V.at(i).xyz();
                        double oldMinSJ = 0;
                        double oldAvgSJ = 0;
                        GetQuality(V.at(i), oldMinSJ, oldAvgSJ);
                        double newMinSJ = 0;
                        double newAvgSJ = 0;
                        V.at(i) = newv;
                        GetQuality(V.at(i), newMinSJ, newAvgSJ);
                        if (newMinSJ < targetMinSJ || newAvgSJ < oldAvgSJ)
                            V.at(i) = oldv;
                        else{
                            const double distance = glm::length(newv- oldv);
                            energy += distance * distance;
                        }
                    }
                }
                else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::vec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary)
                            continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x/count;
                    newV.at(i).y = sum.y/count;
                    newV.at(i).z = sum.z/count;
                }
            }
            else if (v.type == FEATURE) {
                if (preserveSharpFeature)
                    continue;
                if (smoothMethod == LAPLACE_EDGE) {
                    newV.at(i) = LapLace(v);
                }
                else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::vec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary)
                            continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x/count;
                    newV.at(i).y = sum.y/count;
                    newV.at(i).z = sum.z/count;
                }
            }
        }
    }

    std::cout << "Energy = " << energy << std::endl;
    return energy;
}

double Mesh::ProjectSurface(const Mesh& mesh, size_t iters/* = 1*/, const SmoothMethod smoothMethod/* = LAPLACE_EDGE*/,
        const bool preserveSharpFeature/* = false*/, const bool treatSharpFeatureAsRegular/* = false*/, const bool treatCornerAsRegular/* = false*/,
        const bool preserveQuality/* = false*/)
{
    double targetAvgSJ = 0;
    double targetMinSJ = GetMinScaledJacobian(targetAvgSJ);
    std::cout << "targetMinSJ = " << targetMinSJ << " targetAvgSJ = " << targetAvgSJ << std::endl;
    double energy = 0;
    std::vector<Vertex> oldV = V;
    std::vector<Vertex> newV = V;
    while (iters-- != 0)
    {
//#pragma omp parallel for
        for (size_t i = 0; i < V.size(); i++) {
            const Vertex& v = V.at(i);
            if (!v.isBoundary)
                continue;
            if (v.type == REGULAR || (v.type == FEATURE && treatSharpFeatureAsRegular) || (v.type == CORNER && treatCornerAsRegular))
            {
                if (smoothMethod == LAPLACE_EDGE) {
                    //const glm::vec3 smoothedv = LapLace(v, treatSharpFeatureAsRegular, treatCornerAsRegular);
                    const glm::vec3 smoothedv = v.xyz();
                    const glm::vec3 newv = mesh.GetProjectLocation(smoothedv);
                    if (!preserveQuality) {
                        V.at(i) = newv;
                    }
                    else {
                        const glm::vec3 oldv = V.at(i).xyz();
                        double oldMinSJ = 0;
                        double oldAvgSJ = 0;
                        GetQuality(V.at(i), oldMinSJ, oldAvgSJ);
                        double newMinSJ = 0;
                        double newAvgSJ = 0;
                        V.at(i) = newv;
                        GetQuality(V.at(i), newMinSJ, newAvgSJ);
                        if (newMinSJ < targetMinSJ || newAvgSJ < oldAvgSJ)
                            V.at(i) = oldv;
                        else{
                            const double distance = glm::length(newv- oldv);
                            energy += distance * distance;
                        }
                    }
                }
                else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::vec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary)
                            continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x/count;
                    newV.at(i).y = sum.y/count;
                    newV.at(i).z = sum.z/count;
                }
            }
            else if (v.type == FEATURE) {
                if (preserveSharpFeature)
                    continue;
                if (smoothMethod == LAPLACE_EDGE) {
                    newV.at(i) = LapLace(v);
                }
                else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::vec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary)
                            continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x/count;
                    newV.at(i).y = sum.y/count;
                    newV.at(i).z = sum.z/count;
                }
            }
        }
    }

    std::cout << "Energy = " << energy << std::endl;
    return energy;
}

double Mesh::SmoothVolume(const Mesh& mesh, size_t iters/* = 1*/, const SmoothMethod smoothMethod/* = LAPLACE_EDGE*/,
        const bool preserveQuality/* = false*/)
{
    double targetAvgSJ = 0;
    double targetMinSJ = GetMinScaledJacobian(targetAvgSJ);
    std::cout << "targetMinSJ = " << targetMinSJ << " targetAvgSJ = " << targetAvgSJ << std::endl;
    double energy = 0;
    std::vector<Vertex> oldV = V;
    std::vector<Vertex> newV = V;
    while (iters-- != 0)
    {
//#pragma omp parallel for
        for (size_t i = 0; i < V.size(); i++) {
            const Vertex& v = V.at(i);
            if (v.isBoundary)
                continue;
            //if (v.type == REGULAR || (v.type == FEATURE && treatSharpFeatureAsRegular) || (v.type == CORNER && treatCornerAsRegular))
            {
                if (smoothMethod == LAPLACE_EDGE) {
                    const glm::vec3 smoothedv = LapLace(v);
                    //const glm::vec3 newv = mesh.GetProjectLocation(smoothedv);
                    const glm::vec3 newv = smoothedv;
                    if (!preserveQuality) {
                        V.at(i) = newv;
                    }
                    else {
                        const glm::vec3 oldv = V.at(i).xyz();
                        double oldMinSJ = 0;
                        double oldAvgSJ = 0;
                        GetQuality(V.at(i), oldMinSJ, oldAvgSJ);
                        double newMinSJ = 0;
                        double newAvgSJ = 0;
                        V.at(i) = newv;
                        GetQuality(V.at(i), newMinSJ, newAvgSJ);
                        if (newMinSJ < targetMinSJ || newAvgSJ < oldAvgSJ)
                            V.at(i) = oldv;
                        else{
                            const double distance = glm::length(newv- oldv);
                            energy += distance * distance;
                        }
                    }
                }
                else if (smoothMethod == LAPLACE_FACE_CENTER) {
                    glm::vec3 sum(0.0, 0.0, 0.0);
                    int count = 0;
                    for (size_t j = 0; j < v.N_Fids.size(); j++) {
                        const Face& f = F.at(v.N_Fids[j]);
                        if (!f.isBoundary)
                            continue;
                        sum += GetFaceCenter(f);
                        count++;
                    }
                    newV.at(i).x = sum.x/count;
                    newV.at(i).y = sum.y/count;
                    newV.at(i).z = sum.z/count;
                }
            }
        }
    }

    std::cout << "Volum Energy = " << energy << std::endl;
    return energy;
}


const int V_T[8][4] =
{
    {0, 3, 4, 1},
    {1, 0, 5, 2},
    {2, 1, 6, 3},
    {3, 2, 7, 0},
    {4, 7, 5, 0},
    {5, 4, 6, 1},
    {6, 5, 7, 2},
    {7, 6, 4, 3}
};

static float cal_volume_Tet_real(float v0[3],float v1[3],float v2[3],float v3[3])
{
    float v1v0[3], v2v0[3], v3v0[3];
    for (int i = 0; i < 3; i++) {
        v1v0[i] = v1[i] - v0[i];
        v2v0[i] = v2[i] - v0[i];
        v3v0[i] = v3[i] - v0[i];
    }

    float norm1 = sqrt(v1v0[0] * v1v0[0] + v1v0[1] * v1v0[1] + v1v0[2] * v1v0[2]);
    float norm2 = sqrt(v2v0[0] * v2v0[0] + v2v0[1] * v2v0[1] + v2v0[2] * v2v0[2]);
    float norm3 = sqrt(v3v0[0] * v3v0[0] + v3v0[1] * v3v0[1] + v3v0[2] * v3v0[2]);

    float volume = v1v0[0] * (v2v0[1] * v3v0[2] - v2v0[2] * v3v0[1]) - v1v0[1] * (v2v0[0] * v3v0[2] - v2v0[2] * v3v0[0]) + v1v0[2] * (v2v0[0] * v3v0[1] - v2v0[1] * v3v0[0]);
    return volume;
}

bool JudgeDirection(const Mesh& mesh, const Cell& c)
{
    float v[8][3];
    for (int i = 0; i < c.Vids.size(); i++) {
        v[i][0] = mesh.V.at(c.Vids.at(i)).x;
        v[i][1] = mesh.V.at(c.Vids.at(i)).y;
        v[i][2] = mesh.V.at(c.Vids.at(i)).z;
    }

    float VL[8];
    for (int i = 0; i < 8; i++) {
        const int* p = V_T[i];
        VL[i] = cal_volume_Tet_real(v[p[0]], v[p[1]], v[p[2]], v[p[3]]);
    }

    if (VL[0] + VL[1] + VL[2] + VL[3] + VL[4] + VL[5] + VL[6] + VL[7] < 0) return false;
    else
        return true;
}

static float _GetScaledJacobian(const glm::vec3& i, const glm::vec3& j, const glm::vec3& k)
{
    const glm::mat3x3 m(i, j, k);
    return glm::determinant(m);
}

const float Mesh::GetScaledJacobian(const Cell& c) const
{
    Cell c1(c);
    if (!JudgeDirection(*this, c)) {
        std::swap(c1.Vids[0], c1.Vids[3]);
        std::swap(c1.Vids[1], c1.Vids[2]);
        std::swap(c1.Vids[4], c1.Vids[7]);
        std::swap(c1.Vids[5], c1.Vids[6]);
    }
    float minScaledJacobian = 1;
    for (int n = 0; n < 7; n++) {
        const Vertex& o = V.at(c1.Vids.at(n));

        const Vertex& i = V.at(c1.Vids.at(HexPoint_Points[n][0]));
        const Vertex& j = V.at(c1.Vids.at(HexPoint_Points[n][1]));
        const Vertex& k = V.at(c1.Vids.at(HexPoint_Points[n][2]));

        const glm::vec3 ei(i.x - o.x, i.y - o.y, i.z - o.z);
        const glm::vec3 ej(j.x - o.x, j.y - o.y, j.z - o.z);
        const glm::vec3 ek(k.x - o.x, k.y - o.y, k.z - o.z);

        const float length_i = glm::length(ei);
        const float length_j = glm::length(ej);
        const float length_k = glm::length(ek);

        const glm::vec3 ni(ei.x / length_i, ei.y / length_i, ei.z / length_i);
        const glm::vec3 nj(ej.x / length_j, ej.y / length_j, ej.z / length_j);
        const glm::vec3 nk(ek.x / length_k, ek.y / length_k, ek.z / length_k);

        float scaledJacobian = _GetScaledJacobian(ni, nj, nk);
        minScaledJacobian = minScaledJacobian < scaledJacobian ? minScaledJacobian : scaledJacobian;
    }

    return minScaledJacobian;
}
//void Mesh::GetQuality(const Vertex& v, double& minSJ, double& avgSJ)
//{
//    double min = 1.0;
//    double sum = 0;
//    for (size_t i = 0; i < v.N_Cids.size(); i++) {
//        const Cell& cell = C.at(v.N_Cids.at(i));
//        const float scaledJacobian = GetScaledJacobian(cell);
//        if (scaledJacobian < min)
//            min = scaledJacobian;
//        sum += scaledJacobian;
//    }
//    minSJ = min;
//    avgSJ = sum / v.N_Cids.size();
//}

void Mesh::GetQuality(const Vertex& v, double& minSJ, double& avgSJ)
{
    double minScaledJacobian = 1;
    double coordinates[8][3];
    for (size_t i = 0; i < v.N_Cids.size(); i++) {
        const Cell& c = C.at(v.N_Cids.at(i));
        for (size_t j = 0; j < 8; j++) {
            const Vertex& v = V[c.Vids[j]];
            coordinates[j][0] = v.x;
            coordinates[j][1] = v.y;
            coordinates[j][2] = v.z;
        }
        double scaledJacobian = 0;// TODO replace this v_hex_scaled_jacobian(8, coordinates);
        if (scaledJacobian < minScaledJacobian) minScaledJacobian = scaledJacobian;
        avgSJ += scaledJacobian;
    }
    minSJ = minScaledJacobian;
    avgSJ /= v.N_Cids.size();;

//    std::vector<size_t> cellIds(v.N_Cids.size());
//    for (size_t i = 0; i < v.N_Cids.size(); i++) {
//        const Cell& cell = C.at(v.N_Cids.at(i));
//        cellIds[i] = cell.id;
//    }
    //OutputBadCells(cellIds, "temp.vtk");
    //GetQuality("temp.vtk", minSJ, avgSJ);
}

void Mesh::OutputBadCells(const std::vector<size_t>& badCellIds, const char* filename)
{
    std::vector<Cell> cells(badCellIds.size());
    for (size_t i = 0; i < badCellIds.size(); i++)
        cells.at(i) = C.at(badCellIds.at(i));
    MeshFileWriter writer(V, cells, filename);
    writer.FixMesh();
    writer.WriteFile();
}

glm::vec3 Mesh::LapLace(const Vertex& v, const bool treatSharpFeatureAsRegular/* = false*/, const bool treatCornerAsRegular/* = false*/)
{
    if (v.type == REGULAR || (v.type == FEATURE && treatSharpFeatureAsRegular) || (v.type == CORNER && treatCornerAsRegular)) {
        glm::vec3 sum(0.0, 0.0, 0.0);
        int count = 0;
        for (size_t j = 0; j < v.N_Vids.size(); j++) {
            const Vertex& n_v = V.at(v.N_Vids[j]);
            if (!n_v.isBoundary)
                continue;
            sum += n_v.xyz();
            count++;
        }
        return glm::vec3(sum.x/count, sum.y/count, sum.z/count);
    }
    else if (v.type == FEATURE) {
        glm::vec3 sum(0.0, 0.0, 0.0);
        int count = 0;
        for (size_t j = 0; j < v.N_Vids.size(); j++) {
            const Vertex& n_v = V.at(v.N_Vids[j]);
            if (!n_v.isBoundary)
                continue;
            if (n_v.type != FEATURE && n_v.type != CORNER)
                continue;
            sum += n_v.xyz();
            count++;
        }
        if (count == 1) {
            int c = 0;
            sum = glm::vec3 (0.0, 0.0, 0.0);
            for (size_t j = 0; j < v.N_Vids.size(); j++) {
                const Vertex& n_v = V.at(v.N_Vids[j]);
                if (!n_v.isBoundary)
                    continue;
                sum += n_v.xyz();
                c++;
            }
            return glm::vec3(sum.x/c, sum.y/c, sum.z/c);
        }
        else if (count == 2) {
            return glm::vec3(sum.x/count, sum.y/count, sum.z/count);
        }
        else {
            //std::cerr << "\033[1;31mERROR\033[0m in smoothing Vertex of Sharp Feature!!" << std::endl;
            return v.xyz();
        }
    }
    else if (v.type == CORNER)
        return v.xyz();
    else {
        glm::vec3 sum(0.0, 0.0, 0.0);
        int count = 0;
        for (size_t j = 0; j < v.N_Vids.size(); j++) {
            const Vertex& n_v = V.at(v.N_Vids[j]);
            //if (!n_v.isBoundary)
            //    continue;
            sum += n_v.xyz();
            count++;
        }
        return glm::vec3(sum.x/count, sum.y/count, sum.z/count);
    }
}

//glm::vec3 Mesh::LapLace(const Vertex& v)
//{
//    if (v.type == REGULAR || v.type == FEATURE) {
//        glm::vec3 sum(0.0, 0.0, 0.0);
//        int count = 0;
//        for (size_t j = 0; j < v.N_Vids.size(); j++) {
//            const Vertex& n_v = V.at(v.N_Vids[j]);
//            if (!n_v.isBoundary)
//                continue;
//            sum += n_v.xyz();
//            count++;
//        }
//        return glm::vec3(sum.x/count, sum.y/count, sum.z/count);
//    }
//    else if (v.type == FEATURE) {
//        glm::vec3 sum(0.0, 0.0, 0.0);
//        int count = 0;
//        for (size_t j = 0; j < v.N_Vids.size(); j++) {
//            const Vertex& n_v = V.at(v.N_Vids[j]);
//            if (!n_v.isBoundary)
//                continue;
//            if (n_v.type != FEATURE && n_v.type != CORNER)
//                continue;
//            sum += n_v.xyz();
//            count++;
//        }
//        if (count == 1) {
//            int c = 0;
//            sum = glm::vec3 (0.0, 0.0, 0.0);
//            for (size_t j = 0; j < v.N_Vids.size(); j++) {
//                const Vertex& n_v = V.at(v.N_Vids[j]);
//                if (!n_v.isBoundary)
//                    continue;
//                sum += n_v.xyz();
//                c++;
//            }
//            return glm::vec3(sum.x/c, sum.y/c, sum.z/c);
//        }
//        else if (count == 2) {
//            return glm::vec3(sum.x/count, sum.y/count, sum.z/count);
//        }
//        else {
//            std::cout << "\033[1;31mERROR\033[0m in smoothing Vertex of Sharp Feature!!" << std::endl;
//        }
//    }
//    else
//        return v.xyz();
//}

static void set_cross(const std::vector<size_t>& set1, const std::vector<size_t>& set2, std::vector<size_t> &result_set)
{
    result_set.clear();
    for (size_t i = 0; i < set1.size(); i++)
    {
        bool inside = false;
        for (size_t j = 0; j < set2.size(); j++)
        {
            if (set2[j] == set1[i])
            {
                inside = true;
                break;
            }
        }
        if (inside)
            result_set.push_back(set1[i]);
    }
}

void Mesh::BuildF()
{
    if (m_cellType == TRIANGLE || m_cellType == QUAD) {
        F.resize(C.size());
        for (size_t i = 0; i < C.size(); i++) {
            F[i].Vids = C[i].Vids;
            F[i].id = C[i].id;
            F[i].isBoundary = true;
        }
        if (m_cellType == TRIANGLE) {
            for (size_t i = 0; i < F.size(); i++) {
                F[i].Eids.resize(3);
                for (size_t j = 0; j < 3; j++) {
                    bool found = false;
                    for (size_t k = 0; k < V[F[i].Vids[j]].N_Eids.size(); k++) {
                        size_t ne = V[F[i].Vids[j]].N_Eids[k];
                        for (size_t m = 0; m < V[F[i].Vids[(j + 1) % 3]].N_Eids.size(); m++) {
                            if (ne == V[F[i].Vids[(j + 1) % 3]].N_Eids[m]) {
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
        else if (m_cellType == QUAD) {
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
    } else if (m_cellType == HEXAHEDRA) {
        int F_N = 0;
        for (size_t i = 0; i < C.size(); i++) {
            std::vector<Face> hf(6, Face(4, 4));
            for (size_t j = 0; j < 6; j++)
                for (size_t k = 0; k < 4; k++)
                    hf[j].Vids[k] = C[i].Vids[HexFaces[j][k]];

            for (size_t j = 0; j < 6; j++) {
                bool have = false;
                for (size_t m = 0; m < 4; m++) {
                    for (size_t n = 0; n < V[hf[j].Vids[m]].N_Fids.size(); n++) {
                        const size_t F_id = V[hf[j].Vids[m]].N_Fids[n];
                        bool all_have = true;
                        for (size_t p = 0; p < 4; p++) {
                            bool exist_v = false;
                            for (size_t q = 0; q < 4; q++)
                                if (hf[j].Vids[p] == F[F_id].Vids[q])
                                    exist_v = true;
                            if (!exist_v)
                                all_have = false;
                        }
                        if (all_have) {
                            have = true;
                            F[F_id].N_Cids.push_back(i);
                        }
                    }
                }
                if (!have) {
                    hf[j].id = F_N++;
                    for (size_t k = 0; k < 4; k++) {
                        size_t id1 = hf[j].Vids[k];
                        size_t id2 = hf[j].Vids[(k + 1) % 4];
                        bool found = false;
                        for (size_t m = 0; m < V[id1].N_Eids.size(); m++) {
                            int edge1 = V[id1].N_Eids[m];
                            for (size_t n = 0; n < V[id2].N_Eids.size(); n++) {
                                size_t edge2 = V[id2].N_Eids[n];
                                if (edge1 == edge2) {
                                    hf[j].Eids[k] = edge1;
                                    found = true;
                                }
                                if (found)
                                    break;
                            }
                            if (found)
                                break;
                        }
                    }
                    F.push_back(hf[j]);
                    V[hf[j].Vids[0]].N_Fids.push_back(hf[j].id);
                    V[hf[j].Vids[1]].N_Fids.push_back(hf[j].id);
                    V[hf[j].Vids[2]].N_Fids.push_back(hf[j].id);
                    V[hf[j].Vids[3]].N_Fids.push_back(hf[j].id);

                    F[F.size() - 1].N_Cids.push_back(i);
                }
            }
        }

        for (size_t i = 0; i < F.size(); i++) {
            std::vector<size_t> N_Cids = F[i].N_Cids;
            F[i].N_Cids.clear();
            for (size_t j = 0; j < N_Cids.size(); j++) {
                bool already = false;
                for (size_t k = 0; k < F[i].N_Cids.size(); k++) {
                    if (N_Cids[j] == F[i].N_Cids[k]) {
                        already = true;
                        break;
                    }
                }
                if (!already)
                    F[i].N_Cids.push_back(N_Cids[j]);
            }
        }

        for (size_t i = 0; i < F.size(); i++) {
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
                    if (found)
                        break;
                }
            }
        }

        for (size_t i = 0; i < F.size(); i++)
            for (size_t j = 0; j < F[i].N_Cids.size(); j++)
                C[F[i].N_Cids[j]].N_Fids.push_back(i);

        for (size_t i = 0; i < C.size(); i++) {
            C[i].Fids.resize(6);
            std::vector<size_t> f_ids = C[i].N_Fids;
            C[i].N_Fids.clear();
            for (size_t j = 0; j < f_ids.size(); j++) {
                bool havesame = false;
                for (size_t k = j + 1; k < f_ids.size(); k++)
                    if (f_ids[j] == f_ids[k])
                        havesame = true;
                if (!havesame) {
                    C[i].N_Fids.push_back(f_ids[j]);
                    C[i].Fids[C[i].N_Fids.size() - 1] = F[f_ids[j]].id;
                }
            }
        }
        ////////////////////////////////////////////////
        std::vector<bool> Vs_flags(V.size(), false);
        for (int i = 0; i < F.size(); i++) {
            for (int j = 0; j < 4; j++)
                Vs_flags[F[i].Vids[j]] = true;
            for (int j = 0; j < F[i].N_Cids.size(); j++) {
                int nhid = F[i].N_Cids[j];
                for (int k = 0; k < 6; k++) {
                    bool have_true = false;
                    for (int m = 0; m < 4; m++)
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

            std::vector<size_t> N_Ortho_4Vs1(4), N_Ortho_4Vs2(4);
            for (size_t k = 0; k < 4; k++) {
                size_t fvid = F[F[i].N_Fids[0]].Vids[k];
                N_Ortho_4Vs1[k] = fvid;
            }
            F[i].N_Ortho_4Vids.push_back(N_Ortho_4Vs1);
            for (size_t j = 0; j < F[i].N_Fids.size(); j++) {
                if (j == 1)
                    for (size_t k = 0; k < 4; k++)
                        Vs_flags[F[F[i].N_Fids[1]].Vids[k]] = true;
                for (size_t k = 0; k < 4; k++) {
                    int fvid = N_Ortho_4Vs1[k];
                    for (size_t m = 0; m < V[fvid].N_Vids.size(); m++) {
                        if (Vs_flags[V[fvid].N_Vids[m]]) {
                            N_Ortho_4Vs2[k] = V[fvid].N_Vids[m];
                            break;
                        }
                    }
                }
                F[i].N_Ortho_4Vids.push_back(N_Ortho_4Vs2);
                N_Ortho_4Vs1 = N_Ortho_4Vs2;
                for (size_t k = 0; k < 4; k++)
                    Vs_flags[N_Ortho_4Vs1[k]] = false;
            }

            std::vector<size_t> N_4Eids(4);
            for (size_t j = 1; j < F[i].N_Ortho_4Vids.size(); j++) {
                for (size_t k = 0; k < 4; k++) {
                    std::vector<size_t> sharedEids;
                    set_cross(V[F[i].N_Ortho_4Vids[j - 1][k]].N_Eids, V[F[i].N_Ortho_4Vids[j][k]].N_Eids, sharedEids);
                    N_4Eids[k] = sharedEids[0];
                }
                F[i].N_Ortho_4Eids.push_back(N_4Eids);
            }
        }
    } else if (m_cellType == TETRAHEDRA) {
        int F_N = 0;
        for (size_t i = 0; i < C.size(); i++) {
            std::vector<Face> hf(4, Face(3, 3));
            for (size_t j = 0; j < 4; j++)
                for (size_t k = 0; k < 3; k++)
                    hf[j].Vids[k] = C[i].Vids[TetFaces[j][k]];

            for (size_t j = 0; j < 4; j++) {
                bool have = false;
                for (size_t m = 0; m < 3; m++) {
                    for (size_t n = 0; n < V[hf[j].Vids[m]].N_Fids.size(); n++) {
                        const size_t F_id = V[hf[j].Vids[m]].N_Fids[n];
                        bool all_have = true;
                        for (size_t p = 0; p < 3; p++) {
                            bool exist_v = false;
                            for (size_t q = 0; q < 3; q++)
                                if (hf[j].Vids[p] == F[F_id].Vids[q])
                                    exist_v = true;
                            if (!exist_v)
                                all_have = false;
                        }
                        if (all_have) {
                            have = true;
                            F[F_id].N_Cids.push_back(i);
                        }
                    }
                }
                if (!have) {
                    hf[j].id = F_N++;
                    for (size_t k = 0; k < 3; k++) {
                        size_t id1 = hf[j].Vids[k];
                        size_t id2 = hf[j].Vids[(k + 1) % 3];
                        bool found = false;
                        for (size_t m = 0; m < V[id1].N_Eids.size(); m++) {
                            int edge1 = V[id1].N_Eids[m];
                            for (size_t n = 0; n < V[id2].N_Eids.size(); n++) {
                                size_t edge2 = V[id2].N_Eids[n];
                                if (edge1 == edge2) {
                                    hf[j].Eids[k] = edge1;
                                    found = true;
                                }
                                if (found)
                                    break;
                            }
                            if (found)
                                break;
                        }
                    }
                    F.push_back(hf[j]);
                    V[hf[j].Vids[0]].N_Fids.push_back(hf[j].id);
                    V[hf[j].Vids[1]].N_Fids.push_back(hf[j].id);
                    V[hf[j].Vids[2]].N_Fids.push_back(hf[j].id);

                    F[F.size() - 1].N_Cids.push_back(i);
                }
            }
        }

        for (size_t i = 0; i < F.size(); i++) {
            std::vector<size_t> N_Cids = F[i].N_Cids;
            F[i].N_Cids.clear();
            for (size_t j = 0; j < N_Cids.size(); j++) {
                bool already = false;
                for (size_t k = 0; k < F[i].N_Cids.size(); k++) {
                    if (N_Cids[j] == F[i].N_Cids[k]) {
                        already = true;
                        break;
                    }
                }
                if (!already)
                    F[i].N_Cids.push_back(N_Cids[j]);
            }
        }

        for (size_t i = 0; i < F.size(); i++) {
            for (size_t j = 0; j < 3; j++) {
                bool found = false;
                for (size_t k = 0; k < V[F[i].Vids[j]].N_Eids.size(); k++) {
                    size_t ne = V[F[i].Vids[j]].N_Eids[k];
                    for (size_t m = 0; m < V[F[i].Vids[(j + 1) % 3]].N_Eids.size(); m++) {
                        if (ne == V[F[i].Vids[(j + 1) % 3]].N_Eids[m]) {
                            E[ne].N_Fids.push_back(i);
                            F[i].Eids[j] = ne;
                            found = true;
                            break;
                        }
                    }
                    if (found)
                        break;
                }
            }
        }

        for (size_t i = 0; i < F.size(); i++)
            for (size_t j = 0; j < F[i].N_Cids.size(); j++)
                C[F[i].N_Cids[j]].N_Fids.push_back(i);

        for (size_t i = 0; i < C.size(); i++) {
            C[i].Fids.resize(4);
            std::vector<size_t> f_ids = C[i].N_Fids;
            C[i].N_Fids.clear();
            for (size_t j = 0; j < f_ids.size(); j++) {
                bool havesame = false;
                for (size_t k = j + 1; k < f_ids.size(); k++)
                    if (f_ids[j] == f_ids[k])
                        havesame = true;
                if (!havesame) {
                    C[i].N_Fids.push_back(f_ids[j]);
                    C[i].Fids[C[i].N_Fids.size() - 1] = F[f_ids[j]].id;
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
}

void Mesh::BuildV_V()
{

}
void Mesh::BuildV_E()
{

}
void Mesh::BuildV_F()
{

}
void Mesh::BuildV_C()
{
    if (m_cellType == HEXAHEDRA || m_cellType == TETRAHEDRA)
        for (size_t i = 0; i < C.size(); i++)
            for (size_t j = 0; j < C[i].Vids.size(); j++)
                V[C[i].Vids[j]].N_Cids.push_back(i);
    else if (m_cellType == TRIANGLE || m_cellType == QUAD)
        for (size_t i = 0; i < C.size(); i++)
            for (size_t j = 0; j < C[i].Vids.size(); j++)
                V[C[i].Vids[j]].N_Fids.push_back(i);
}
void Mesh::BuildE_V()
{
    for (size_t i = 0; i < E.size(); i++){
        for (size_t j = 0; j < E[i].Vids.size(); j++){
            for (size_t k = 0; k < V[E[i].Vids[j]].N_Vids.size(); k++) {
                E[i].N_Vids.push_back(V[E[i].Vids[j]].N_Vids[k]);
            }
        }
        for (size_t j = 0; j < E[i].Vids.size(); j++) {
            for (std::vector<size_t>::iterator iter = E[i].N_Vids.begin(); iter != E[i].N_Vids.end();) {
                if (*iter == E[i].Vids[j])
                    iter = E[i].N_Vids.erase(iter);
                else
                    ++iter;
            }
        }
    }
}
void Mesh::BuildE_E()
{
    for (size_t i = 0; i < E.size(); i++){
        for (size_t j = 0; j < E[i].Vids.size(); j++){
            for (size_t k = 0; k < V[E[i].Vids[j]].N_Eids.size(); k++) {
                E[i].N_Eids.push_back(V[E[i].Vids[j]].N_Eids[k]);
            }
        }
        for (size_t j = 0; j < E[i].Vids.size(); j++) {
            for (std::vector<size_t>::iterator iter = E[i].N_Eids.begin(); iter != E[i].N_Eids.end();) {
                if (*iter == E[i].id)
                    iter = E[i].N_Eids.erase(iter);
                else
                    ++iter;
            }
        }
    }
}
void Mesh::BuildE_F()
{
    for (size_t i = 0; i < E.size(); i++){
        for (size_t j = 0; j < E[i].Vids.size(); j++){
            for (size_t k = 0; k < V[E[i].Vids[j]].N_Fids.size(); k++) {
                const size_t N_Fid = V[E[i].Vids[j]].N_Fids[k];
                if (std::find(E[i].N_Fids.begin(), E[i].N_Fids.end(), N_Fid) == E[i].N_Fids.end())
                    E[i].N_Fids.push_back(N_Fid);
            }
        }
    }
}
void Mesh::BuildE_C()
{

}
void Mesh::BuildF_V()
{

}
void Mesh::BuildF_E()
{

}
void Mesh::BuildF_F()
{
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
void Mesh::BuildF_C()
{

}
void Mesh::BuildC_V()
{

}
void Mesh::BuildC_E()
{
    if (m_cellType == HEXAHEDRA)
        for (size_t i = 0; i < C.size(); i++) {
            Cell& cell = C.at(i);
            cell.Eids.resize(12);
            for (size_t k = 0; k < 12; k++) {
                Edge e1(2);
                e1.Vids[0] = cell.Vids[HexEdge[k][0]];
                e1.Vids[1] = cell.Vids[HexEdge[k][1]];
                for (size_t j = 0; j < cell.Fids.size(); j++) {
                    bool found = false;
                    const Face& face = F.at(cell.Fids.at(j));
                    for (size_t n = 0; n < face.Eids.size(); n++) {
                        const Edge& e2 = E[face.Eids[n]];
                        if (e1 == e2) {
                            cell.Eids[k] = e2.id;
                            found = true;
                            break;
                        }
                    }
                    if (found)
                        break;
                }
            }
        }
    else if (m_cellType == TETRAHEDRA)
        for (size_t i = 0; i < C.size(); i++) {
            Cell& cell = C.at(i);
            cell.Eids.resize(6);
            for (size_t k = 0; k < 6; k++) {
                Edge e1(2);
                e1.Vids[0] = cell.Vids[TetEdge[k][0]];
                e1.Vids[1] = cell.Vids[TetEdge[k][1]];
                for (size_t j = 0; j < cell.Fids.size(); j++) {
                    bool found = false;
                    const Face& face = F.at(cell.Fids.at(j));
                    for (size_t n = 0; n < face.Eids.size(); n++) {
                        const Edge& e2 = E[face.Eids[n]];
                        if (e1 == e2) {
                            cell.Eids[k] = e2.id;
                            found = true;
                            break;
                        }
                    }
                    if (found)
                        break;
                }
            }
        }
}
void Mesh::BuildC_F()
{

}
void Mesh::BuildC_C()
{

}

void Mesh::BuildAllConnectivities()
{
    BuildV_C();
    BuildE(); // BuildV_V(); BuildV_E(); BuildV_F();
    BuildF(); // BuildF_C();
    BuildC_E();
//    if (m_cellType != HEXAHEDRA){
//        BuildE_V();
//        BuildE_E();
//        BuildE_F();
//    }
    BuildF_F();
}

//static void set_cross(const std::vector<size_t> set1, const std::vector<size_t> set2, std::vector<size_t> &result_set)
//{
//    result_set.clear();
//    for (int i = 0; i < set1.size(); i++)
//    {
//        bool inside = false;
//        for (int j = 0; j < set2.size(); j++)
//        {
//            if (set2[j] == set1[i])
//            {
//                inside = true;
//                break;
//            }
//        }
//        if (inside)
//            result_set.push_back(set1[i]);
//    }
//}

void Mesh::ExtractBoundary()
{
    if (m_cellType == TRIANGLE || m_cellType == QUAD) {
        for (size_t i = 0; i < V.size(); i++)
            V[i].isBoundary = false;
        for (size_t i = 0; i < E.size(); i++)
            E[i].isBoundary = false;
        for (size_t i = 0; i < F.size(); i++)
            F[i].isBoundary = false;

        if (m_cellType == TRIANGLE)
            return;

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
        return;
    }

    for (size_t i = 0; i < E.size(); i++) {
        std::vector<size_t> neighborhs = E[i].N_Cids;
        E[i].N_Cids.clear();
        size_t pre_ = neighborhs[0];
        size_t aft_;
        E[i].N_Cids.push_back(pre_);
        std::vector<bool> tags(neighborhs.size(), false);
        tags[0] = true;
        for (size_t m = 1; m < neighborhs.size(); m++) {
            for (size_t j = 0; j < neighborhs.size(); j++) {
                if (tags[j])
                    continue;

                std::vector<size_t> fs_pre, fs_cur;
                if (m_cellType == HEXAHEDRA)
                    for (int k = 0; k < 6; k++) {
                        fs_pre.push_back(C[pre_].Fids[k]);
                        fs_cur.push_back(C[neighborhs[j]].Fids[k]);
                    }
                else if (m_cellType == TETRAHEDRA)
                    for (int k = 0; k < 4; k++) {
                        fs_pre.push_back(C[pre_].Fids[k]);
                        fs_cur.push_back(C[neighborhs[j]].Fids[k]);
                    }
                std::vector<size_t> sharedf;
                set_cross(fs_pre, fs_cur, sharedf);
                if (sharedf.size()) {
                    E[i].N_Cids.push_back(neighborhs[j]);
                    tags[j] = true;
                    pre_ = neighborhs[j];
                }
            }
        }
        if (E[i].N_Cids.size() != neighborhs.size()) {
            //cout<<"ERROR here"<<endl;//boundary edge have not been handled yet
            E[i].N_Cids = neighborhs;
        }
    }
    for (size_t i = 0; i < F.size(); i++) {
        F[i].isBoundary = false;
        if (F[i].N_Cids.size() == 1) {
            F[i].isBoundary = true;
            if (m_cellType == HEXAHEDRA)
                for (size_t j = 0; j < 4; j++) {
                    E[F[i].Eids[j]].isBoundary = true;
                    V[F[i].Vids[j]].isBoundary = true;
                }
            else if (m_cellType == TETRAHEDRA)
                for (size_t j = 0; j < 3; j++) {
                    E[F[i].Eids[j]].isBoundary = true;
                    V[F[i].Vids[j]].isBoundary = true;
                }
        }
    }
    for (size_t i = 0; i < C.size(); i++) {
        C[i].isBoundary = false;
        for (size_t j = 0; j < C[i].N_Fids.size(); j++) {
            size_t fid = C[i].N_Fids[j];
            if (F[fid].isBoundary)
                C[i].isBoundary = true;
        }
    }
}

size_t Mesh::ExtractLayers()
{
    std::vector<bool> savedCellBoundary(C.size(), false);
    for (size_t i = 0; i < C.size(); i++)
        savedCellBoundary.at(i) = C[i].isBoundary;

    std::vector<Cell> surfaceCells;
    std::vector<Cell> innerCells;

    int layerCount = 1;
    for (size_t i = 0; i < C.size(); i++)
    {
        if (C[i].isBoundary)
            surfaceCells.push_back(C[i]);
        else
            innerCells.push_back(C[i]);
    }
    while (!innerCells.empty()) {
        // -------------------------------------------------------
        // build outLayer
        {
            Layer layer;
            for (size_t i = 0; i < surfaceCells.size(); i++)
                layer.Cids.push_back(surfaceCells[i].id);

            for (size_t i = 0; i < layer.Cids.size(); i++){
                const Cell& cell = C[layer.Cids[i]];
                for (size_t j = 0; j < cell.Fids.size(); j++)
                    layer.Fids.push_back(cell.Fids[j]);
                for (size_t j = 0; j < cell.Eids.size(); j++)
                    layer.Eids.push_back(cell.Eids[j]);
                for (size_t j = 0; j < cell.Vids.size(); j++)
                    layer.Vids.push_back(cell.Vids[j]);
            }
            std::sort(layer.Fids.begin(), layer.Fids.end());
            std::vector<size_t>::iterator iterF = std::unique(layer.Fids.begin(), layer.Fids.end());
            layer.Fids.resize(std::distance(layer.Fids.begin(), iterF));

            std::sort(layer.Eids.begin(), layer.Eids.end());
            std::vector<size_t>::iterator iterE = std::unique(layer.Eids.begin(), layer.Eids.end());
            layer.Eids.resize(std::distance(layer.Eids.begin(), iterE));

            std::sort(layer.Vids.begin(), layer.Vids.end());
            std::vector<size_t>::iterator iterV = std::unique(layer.Vids.begin(), layer.Vids.end());
            layer.Vids.resize(std::distance(layer.Vids.begin(), iterV));

            layer.fixed.resize(layer.Vids.size(), false);
            for (size_t i = 0; i < layer.Vids.size(); i++)
                if (V[layer.Vids[i]].isBoundary)
                    layer.fixed.at(i) = true;
            layers.push_back(layer);
        }
        // build innerLayer
        {
            Layer layer;
            for (size_t i = 0; i < innerCells.size(); i++)
                layer.Cids.push_back(innerCells[i].id);

            for (size_t i = 0; i < layer.Cids.size(); i++){
                const Cell& cell = C[layer.Cids[i]];
                for (size_t j = 0; j < cell.Fids.size(); j++)
                    layer.Fids.push_back(cell.Fids[j]);
                for (size_t j = 0; j < cell.Eids.size(); j++)
                    layer.Eids.push_back(cell.Eids[j]);
                for (size_t j = 0; j < cell.Vids.size(); j++)
                    layer.Vids.push_back(cell.Vids[j]);
            }
            std::sort(layer.Fids.begin(), layer.Fids.end());
            std::vector<size_t>::iterator iterF = std::unique(layer.Fids.begin(), layer.Fids.end());
            layer.Fids.resize(std::distance(layer.Fids.begin(), iterF));

            std::sort(layer.Eids.begin(), layer.Eids.end());
            std::vector<size_t>::iterator iterE = std::unique(layer.Eids.begin(), layer.Eids.end());
            layer.Eids.resize(std::distance(layer.Eids.begin(), iterE));

            std::sort(layer.Vids.begin(), layer.Vids.end());
            std::vector<size_t>::iterator iterV = std::unique(layer.Vids.begin(), layer.Vids.end());
            layer.Vids.resize(std::distance(layer.Vids.begin(), iterV));

            layer.fixed.resize(layer.Vids.size(), false);
            for (size_t i = 0; i < layer.Vids.size(); i++)
                if (V[layer.Vids[i]].isBoundary)
                    layer.fixed.at(i) = true;
            innerLayers.push_back(layer);
        }
        // -------------------------------------------------------
        std::string surfaceCellsFileName = std::string("OutLayer") + std::to_string(layerCount) + ".vtk";
        MeshFileWriter surfaceCellsWriter(V, surfaceCells, surfaceCellsFileName.c_str(), HEXAHEDRA);
        surfaceCellsWriter.SetFixFlag(true);
        surfaceCellsWriter.WriteFile();

        std::string innerCellsFileName = std::string("InnerLayer") + std::to_string(layerCount) + ".vtk";
        MeshFileWriter innerCellsWriter(V, innerCells, innerCellsFileName.c_str(), HEXAHEDRA);
        innerCellsWriter.SetFixFlag(true);
        innerCellsWriter.WriteFile();

        MeshFileReader innerCellsReader(innerCellsFileName.c_str());
        Mesh& mesh = (Mesh&)innerCellsReader.GetMesh();
        mesh.BuildAllConnectivities();
        mesh.ExtractBoundary();

        std::vector<Cell> newC = innerCells;
        surfaceCells.clear();
        innerCells.clear();
        for (size_t i = 0; i < newC.size(); i++)
        {
            if (mesh.C[i].isBoundary)
                surfaceCells.push_back(newC[i]);
            else
                innerCells.push_back(newC[i]);
        }

        layerCount++;
    }

    for (size_t i = 0; i < C.size(); i++)
        C[i].isBoundary = savedCellBoundary.at(i);

    return layerCount;
}

//size_t Mesh::ExtractLayers()
//{
//    std::vector<Cell> surfaceCells;
//    std::vector<Cell> innerCells;
//
//    int layerCount = 1;
//    for (size_t i = 0; i < C.size(); i++)
//    {
//        if (C[i].isBoundary)
//            surfaceCells.push_back(Cell(C[i].Vids));
//        else
//            innerCells.push_back(Cell(C[i].Vids));
//    }
//    while (!innerCells.empty()) {
//        std::string surfaceCellsFileName = std::string("OutLayer") + std::to_string(layerCount) + ".vtk";
//        MeshFileWriter surfaceCellsWriter(V, surfaceCells, surfaceCellsFileName.c_str(), HEXAHEDRA);
//        surfaceCellsWriter.SetFixFlag(true);
//        surfaceCellsWriter.WriteFile();
//
//        std::string innerCellsFileName = std::string("InnerLayer") + std::to_string(layerCount) + ".vtk";
//        MeshFileWriter innerCellsWriter(V, innerCells, innerCellsFileName.c_str(), HEXAHEDRA);
//        innerCellsWriter.WriteFile();
//
//        MeshFileReader innerCellsReader(innerCellsFileName.c_str());
//        Mesh& mesh = (Mesh&)innerCellsReader.GetMesh();
//        mesh.BuildAllConnectivities();
//        mesh.ExtractBoundary();
//
//        surfaceCells.clear();
//        innerCells.clear();
//        for (size_t i = 0; i < mesh.C.size(); i++)
//        {
//            if (mesh.C[i].isBoundary)
//                surfaceCells.push_back(Cell(mesh.C[i].Vids));
//            else
//                innerCells.push_back(Cell(mesh.C[i].Vids));
//        }
//
//        layerCount++;
//    }
//
//    {
//        Layer layer;
//        for (size_t i = 0; i < C.size(); i++)
//            if (C[i].isBoundary)
//                layer.Cids.push_back(C[i].id);
//        for (size_t i = 0; i < layer.Cids.size(); i++){
//            const Cell& cell = C[layer.Cids[i]];
//            for (size_t j = 0; j < cell.Fids.size(); j++)
//                layer.Fids.push_back(cell.Fids[j]);
//            for (size_t j = 0; j < cell.Eids.size(); j++)
//                layer.Eids.push_back(cell.Eids[j]);
//            for (size_t j = 0; j < cell.Vids.size(); j++)
//                layer.Vids.push_back(cell.Vids[j]);
//        }
//        std::sort(layer.Fids.begin(), layer.Fids.end());
//        std::vector<size_t>::iterator iterF = std::unique(layer.Fids.begin(), layer.Fids.end());
//        layer.Fids.resize(std::distance(layer.Fids.begin(), iterF));
//
//        std::sort(layer.Eids.begin(), layer.Eids.end());
//        std::vector<size_t>::iterator iterE = std::unique(layer.Eids.begin(), layer.Eids.end());
//        layer.Eids.resize(std::distance(layer.Eids.begin(), iterE));
//
//        std::sort(layer.Vids.begin(), layer.Vids.end());
//        std::vector<size_t>::iterator iterV = std::unique(layer.Vids.begin(), layer.Vids.end());
//        layer.Vids.resize(std::distance(layer.Vids.begin(), iterV));
//
//        layer.fixed.resize(layer.Vids.size(), false);
//        for (size_t i = 0; i < layer.Vids.size(); i++)
//            if (V[layer.Vids[i]].isBoundary)
//                layer.fixed.at(i) = true;
//        layers.push_back(layer);
//    }
//
//    return layerCount;
//}

void Mesh::ExtractSingularities()
{
    if (m_cellType == HEXAHEDRA) {
        for (size_t i = 0; i < E.size(); i++) {
            Edge& edge = E.at(i);
            if (edge.isBoundary) {
                if (edge.N_Cids.size() == 1 || edge.N_Cids.size() == 3) edge.isSingularity = true;
            }
            else {
                if (edge.N_Cids.size() != 4) edge.isSingularity = true;
            }

            if (edge.isSingularity) {
                Vertex& v1 = V.at(edge.Vids.at(0));
                Vertex& v2 = V.at(edge.Vids.at(1));
                v1.isSingularity = true;
                v2.isSingularity = true;
            }
        }
    }
    else if (m_cellType == QUAD) {
        for (size_t i = 0; i < V.size(); i++) {
            Vertex& v = V.at(i);
            if (v.N_Fids.size() == 3 || v.N_Cids.size() > 4) v.isSingularity = true;
        }
    }
}

void Mesh::ExtractTwoRingNeighborSurfaceFaceIdsForEachVertex(int N/* = 2*/)
{
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V[i];
        if (!v.isBoundary) continue;
        v.twoRingNeighborSurfaceFaceIds.clear();
        int n = N;
        std::vector<size_t> surfaceFaceVids(1, v.id);
        while (n-- != 0) {
            for (size_t j = 0; j < surfaceFaceVids.size(); j++) {
                const Vertex& vv = V[surfaceFaceVids[j]];
                for (size_t k = 0; k < vv.N_Fids.size(); k++) {
                    const Face& neighboringF = F[vv.N_Fids[k]];
                    if (neighboringF.isBoundary) v.twoRingNeighborSurfaceFaceIds.push_back(neighboringF.id);
                }
            }
            for (size_t j = 0; j < v.twoRingNeighborSurfaceFaceIds.size(); j++) {
                const Face& f = F[v.twoRingNeighborSurfaceFaceIds[j]];
                for (size_t k = 0; k < f.Vids.size(); k++)
                    surfaceFaceVids.push_back(f.Vids[k]);
            }
            std::sort(surfaceFaceVids.begin(), surfaceFaceVids.end());
            std::vector<size_t>::iterator iter = std::unique(surfaceFaceVids.begin(), surfaceFaceVids.end());
            surfaceFaceVids.resize(std::distance(surfaceFaceVids.begin(), iter));
        }
        std::sort(v.twoRingNeighborSurfaceFaceIds.begin(), v.twoRingNeighborSurfaceFaceIds.end());
        std::vector<size_t>::iterator iter = std::unique(v.twoRingNeighborSurfaceFaceIds.begin(), v.twoRingNeighborSurfaceFaceIds.end());
        v.twoRingNeighborSurfaceFaceIds.resize(std::distance(v.twoRingNeighborSurfaceFaceIds.begin(), iter));
    }
//    for (size_t i = 0; i < V.size(); i++) {
//        Vertex& v = V[i];
//        if (!v.isBoundary) continue;
//        v.twoRingNeighborSurfaceFaceIds.clear();
//        for (size_t j = 0; j < v.N_Vids.size(); j++) {
//            const Vertex& neighboringV = V[v.N_Vids[j]];
//            for (size_t k = 0; k < neighboringV.N_Fids.size(); k++) {
//                const Face& neighboringF = F[neighboringV.N_Fids[k]];
//                if (neighboringF.isBoundary) v.twoRingNeighborSurfaceFaceIds.push_back(neighboringF.id);
//            }
//        }
//        std::sort(v.twoRingNeighborSurfaceFaceIds.begin(), v.twoRingNeighborSurfaceFaceIds.end());
//        std::vector<size_t>::iterator iter = std::unique(v.twoRingNeighborSurfaceFaceIds.begin(), v.twoRingNeighborSurfaceFaceIds.end());
//        v.twoRingNeighborSurfaceFaceIds.resize(std::distance(v.twoRingNeighborSurfaceFaceIds.begin(), iter));
//    }
}
glm::vec3 GetCenter(const std::vector<Vertex>& V)
{
    glm::vec3 sum(0.0, 0.0, 0.0);
    for (size_t i = 0; i < V.size(); i++)
        sum += V.at(i).xyz();
    return glm::vec3(sum.x/V.size(), sum.y/V.size(), sum.z/V.size());
}

void Mesh::Zoom(const glm::vec3& ref, const double scale/* = 1*/)
{
    if (scale == 1.0)
        return;
//    if (glm::length(m_center) == 0)
//        m_center = GetCenter(V);
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        const glm::vec3 dir = v - ref;
        v.x = ref.x + scale * dir.x;
        v.y = ref.y + scale * dir.y;
        v.z = ref.z + scale * dir.z;
    }
}

const unsigned int HexEdge[12][2] =
{
    { 0, 1 },
    { 1, 2 },
    { 2, 3 },
    { 3, 0 },
    { 4, 5 },
    { 5, 6 },
    { 6, 7 },
    { 7, 4 },
    { 0, 4 },
    { 1, 5 },
    { 2, 6 },
    { 3, 7 },
};

const unsigned int TetEdge[6][2] =
{
    { 0, 1 },
    { 0, 2 },
    { 0, 3 },
    { 1, 2 },
    { 1, 3 },
    { 2, 3 }
};

const unsigned int QuadEdge[4][2] =
{
    { 0, 1 },
    { 1, 2 },
    { 2, 3 },
    { 3, 0 }
};

const unsigned int TriEdge[3][2] =
{
    { 0, 1 },
    { 1, 2 },
    { 2, 0 },
};

const unsigned int HexPoint_Points[8][3] =
{
    { 1, 3, 4 },
    { 2, 0, 5 },
    { 3, 1, 6 },
    { 0, 2, 7 },
    { 7, 5, 0 },
    { 4, 6, 1 },                                  //const unsigned int HexEdge[12][2] =
    { 5, 7, 2 },                                  //{
    { 6, 4, 3 }                                   //    { 0, 1 }, // 0
};                                                //    { 1, 2 }, // 1
                                                  //    { 2, 3 }, // 2
const unsigned int HexPoint_Edges[8][3] =         //    { 3, 0 }, // 3
{                                                 //    { 4, 5 }, // 4
    { 0, 3, 8 },  //{ (0,1), (0,3), (0,4) },      //    { 5, 6 }, // 5
    { 1, 0, 9 },  //{ (1,2), (1,0), (1,5) },      //    { 6, 7 }, // 6
    { 2, 1, 10},  //{ (2,3), (2,1), (2,6) },      //    { 7, 4 }, // 7
    { 3, 2, 11},  //{ (3,0), (3,2), (3,7) },      //    { 0, 4 }, // 8
    { 7, 4, 8 },  //{ (4,7), (4,5), (4,0) },      //    { 1, 5 }, // 9
    { 4, 5, 9 },  //{ (5,4), (5,6), (5,1) },      //    { 2, 6 }, // 10
    { 5, 6, 10},  //{ (6,5), (6,7), (6,2) },      //    { 3, 7 }, // 11
    { 6, 7, 11}   //{ (7,6), (7,4), (7,3) },      //};
};

const unsigned int TetPoint_Points[4][3] =
{
    { 1, 2, 3 },
    { 0, 2, 3 },
    { 0, 1, 3 },
    { 0, 1, 2 }
};

const unsigned int QuadPoint_Points[4][2] =
{
    { 1, 3 },
    { 2, 0 },
    { 3, 1 },
    { 0, 2 }
};

const unsigned int TriPoint_Points[3][2] =
{
    { 1, 2 },
    { 2, 0 },
    { 0, 1 }
};

const unsigned int HexFaces1[6][4] =
{
    {0, 3, 2, 1},
    {0, 1, 5, 4},
    {0, 4, 7, 3},
    {1, 2, 6, 5},
    {3, 7, 6, 2},
    {4, 5, 6, 7}
};

const unsigned int HexFaces[6][4] =
{
    {0, 3, 2, 1},
    {4, 5, 6, 7},
    {0, 4, 7, 3},
    {1, 2, 6, 5},
    {0, 1, 5, 4},
    {3, 7, 6, 2}

//    {0, 3, 2, 1},
//    {0, 1, 5, 4},
//    {0, 4, 7, 3},
//    {1, 2, 6, 5},
//    {2, 3, 7, 6},
//    {4, 5, 6, 7}
};

const unsigned int TetFaces[4][3] =
{
    {0, 2, 1},
    {0, 3, 2},
    {0, 1, 3},
    {1, 2, 3}
};

const unsigned int HexPoint_Faces[8][3][4] =
{
    {
    {0, 3, 2, 1},
    {0, 4, 7, 3},
    {0, 1, 5, 4}},

    {
    {0, 3, 2, 1},
    {0, 1, 5, 4},
    {1, 2, 6, 5}},

    {
    {0, 3, 2, 1},
    {1, 2, 6, 5},
    {3, 7, 6, 2}},

    {
    {0, 3, 2, 1},
    {3, 7, 6, 2},
    {0, 4, 7, 3}},
    //
    {
    {4, 5, 6, 7},
    {0, 4, 7, 3},
    {0, 1, 5, 4}},

    {
    {4, 5, 6, 7},
    {0, 1, 5, 4},
    {1, 2, 6, 5}},

    {
    {4, 5, 6, 7},
    {1, 2, 6, 5},
    {3, 7, 6, 2}},

    {
    {4, 5, 6, 7},
    {3, 7, 6, 2},
    {0, 4, 7, 3}}
};

const unsigned int TetPoint_Faces[4][3][3] =
{
    {{0, 1, 2},
    {0, 2, 3},
    {0, 3, 1}},

    {{0, 1, 2},
    {1, 3, 2},
    {0, 3, 1}},

    {{0, 1, 2},
    {0, 2, 3},
    {1, 3, 2}},

    {{1, 3, 2},
    {0, 2, 3},
    {0, 3, 1}}
};

const unsigned long hexTet[5][4] =
{
    {0, 4, 5, 7},
    {2, 5, 6, 7},
    {0, 2, 3, 7},
    {0, 1, 2, 5},
    {0, 2, 5, 7}
};

///////////////////////////////////////////////
const unsigned int HexPoint_Points_CW[8][3] =
{
    { 1, 3, 4 },
    { 0, 2, 5 },
    { 1, 3, 6 },
    { 0, 2, 7 },
    { 0, 5, 7 },
    { 1, 4, 6 },
    { 2, 5, 7 },
    { 3, 4, 6 }
};

const unsigned int TetPoint_Points_CW[4][3] =
{
    { 1, 2, 3 },
    { 0, 2, 3 },
    { 0, 1, 3 },
    { 0, 1, 2 }
};

const unsigned int HexFaces_CW[6][4] =
{
    {0, 1, 2, 3},
    {0, 4, 5, 1},
    {0, 3, 7, 4},
    {1, 5, 6, 2},
    {2, 6, 7, 3},
    {4, 5, 6, 7}
};

const unsigned int TetFaces_CW[4][3] =
{
    {0, 1, 2},
    {0, 2, 3},
    {0, 3, 1},
    {1, 3, 2}
};

const unsigned int HexPoint_Faces_CW[8][3][4] =
{
    {{0, 1, 2, 3},
    {0, 3, 7, 4},
    {0, 4, 5, 1}},

    {{0, 1, 2, 3},
    {0, 4, 5, 1},
    {1, 5, 6, 2}},

    {{0, 1, 2, 3},
    {1, 5, 6, 2},
    {2, 6, 7, 3}},

    {{0, 1, 2, 3},
    {2, 6, 7, 3},
    {0, 3, 7, 4}},
    //
    {{4, 5, 6, 7},
    {0, 3, 7, 4},
    {0, 4, 5, 1}},

    {{4, 5, 6, 7},
    {0, 4, 5, 1},
    {1, 5, 6, 2}},

    {{4, 5, 6, 7},
    {1, 5, 6, 2},
    {2, 6, 7, 3}},

    {{4, 5, 6, 7},
    {2, 6, 7, 3},
    {0, 3, 7, 4}}
};

const unsigned int TetPoint_Faces_CW[4][3][3] =
{
    {{0, 1, 2},
    {0, 2, 3},
    {0, 3, 1}},

    {{0, 1, 2},
    {1, 3, 2},
    {0, 3, 1}},

    {{0, 1, 2},
    {0, 2, 3},
    {1, 3, 2}},

    {{1, 3, 2},
    {0, 2, 3},
    {0, 3, 1}}
};

const unsigned long INVALID_NUM = 3277;

const unsigned int HexTrippleEdge[12][6] =
{
    { 0, 1, 4, 3, 5, 2},
    { 1, 2, 5, 0, 6, 3},
    { 2, 3, 6, 1, 7, 0},
    { 3, 0, 7, 2, 4, 1},
    { 0, 4, 3, 1, 7, 5},
    { 1, 5, 0, 2, 4, 6},
    { 2, 6, 1, 3, 5, 7},
    { 3, 7, 2, 0, 6, 4},
    { 4, 5, 7, 0, 6, 1},
    { 5, 6, 4, 1, 7, 2},
    { 6, 7, 5, 2, 4, 3},
    { 7, 4, 6, 3, 5, 0},
};


bool IsOverlap(const Face& f1, const Face& f2)
{
    bool bRet = false;
    for (size_t i = 0; i < f1.Vids.size(); i++){
        for (size_t j = 0; j < f2.Vids.size(); j++){
            if (f1.Vids.at(i) == f2.Vids.at(j)){
                bRet = true;
                break;
            }
        }
    }
    return bRet;
}


bool Find(const std::vector<size_t>& Ids, const size_t targetId)
{
    bool isFound = false;
    for (size_t i = 0; i < Ids.size(); i++){
        if (targetId == Ids.at(i)){
            isFound = true;
            break;
        }
    }
    return isFound;
}

size_t GetoppositeFaceId(const Mesh& mesh, const size_t cellId, const size_t faceId)
{
    size_t oppositeFaceId = MAXID;
    const Face& face = mesh.F.at(faceId);
    const Cell& cell = mesh.C.at(cellId);
    for (size_t i = 0; i < cell.Fids.size(); i++) {
        const Face& cellFace = mesh.F.at(cell.Fids.at(i));
        if (!IsOverlap(face, cellFace)) {
            oppositeFaceId = cell.Fids.at(i);
            break;
        }
    }
    return oppositeFaceId;
}

bool IsEdgeInCell(const Mesh& mesh, const size_t cellId, const size_t edgeId)
{
    const Cell& cell = mesh.C.at(cellId);
    for (size_t i = 0; i < cell.Eids.size(); i++)
        if (edgeId == cell.Eids.at(i))
            return true;
    return false;
}

bool Mesh::IsPointInside(const glm::vec3& orig, const glm::vec3 dir) const
{
    bool bInside = false;
    for (size_t i = 0; i < F.size(); i++)
    {
        const Face& tri = F.at(i);
        if (!tri.isBoundary)
            continue;
        const glm::vec3& v0 = V.at(tri.Vids.at(0));
        const glm::vec3& v1 = V.at(tri.Vids.at(1));
        const glm::vec3& v2 = V.at(tri.Vids.at(2));
        glm::vec3 position;
        if (glm::intersectRayTriangle(orig, dir, v0, v1, v2, position)
                || glm::intersectRayTriangle(orig, dir, v0, v2, v1, position))
        {
            bInside = !bInside;
        }
    }
    return bInside;
}

double Mesh::getSumOfMinimumScaledJacobian(std::vector<std::vector<int>> regions) {
    // std::vector<int> d = {245, 420, 456, 274, 388, 445, 455, 469};
    // return getMinimumScaledJacobian(d);
    double result = 0;
    for (int i = 0; i < regions.size(); i++) {
        std::vector<int> region = regions.at(i);
        result += getMinimumScaledJacobian(region);
    }
    return result;
}

double Mesh::getIterativeEnergyOfRegion(std::vector<int> faceIds) {
    double minimumSJ = 1.0;
    double sumOfInverted = 0.0;
    double s = 1.0;
    for (int i = 0; i < faceIds.size(); i++) {
        Face& f = F.at(faceIds.at(i));
        double localMinSJ = getMinimumScaledJacobian(f);
        if (localMinSJ < minimumSJ) {
            s = minimumSJ;
            minimumSJ = localMinSJ;
        }
        if (localMinSJ < 0)
            sumOfInverted += localMinSJ;
        
        if (localMinSJ > minimumSJ && localMinSJ < s)
            s = localMinSJ;
    }
    // return minimumSJ <= 0 ? sumOfInverted : (minimumSJ + s) / 2;
    return minimumSJ <= 0 ? sumOfInverted : minimumSJ;
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

double Mesh::getMinimumScaledJacobianOfRegion(Vertex& v, double length) {
    double result = 1.0;
    for (int i = 0; i < F.size(); i++) {
        Face& f = F.at(i);
        double lengthToV = glm::length(v - V.at(f.Vids.at(0)));
        for (int j = 1; j < 4; j++) {
            double newL = glm::length(v - V.at(f.Vids.at(j)));
            if (newL < lengthToV)
                lengthToV = newL;
        }
        if (lengthToV > length)
            continue;

        double minSJOfF = getMinimumScaledJacobian(f);
        if (minSJOfF < result)
            result = minSJOfF;
    }
    return result;
}

double Mesh::getMinimumScaledJacobianAtV(Face& f, Vertex& v) {
    int localIndex = getIndexOf1(v.id, f.Vids);
    glm::vec3 v0 = V.at(f.Vids.at(localIndex));
    glm::vec3 v1 = V.at(f.Vids.at( (localIndex + 1) % 4 ));
    glm::vec3 v2 = V.at(f.Vids.at( (localIndex + 3) % 4 ));
    return glm::cross(glm::normalize(v1 - v0), glm::normalize(v2 - v0)).z;
}

void Mesh::printQuality() {
    double min = 1;
    double max = -1;
    double avg = 0;
    int num = 0;
    int minIndex = -1;
    std::vector<int> invertedFIds;
    for (int i = 0; i < F.size(); i++) {
        Face& f = F.at(i);
        double jacobian = 1.0;
        for (int j = 0; j < 4; j++) {
            glm::vec3 v0 = V.at(f.Vids.at(j));
            glm::vec3 v1 = V.at(f.Vids.at( (j + 1) % 4 ));
            glm::vec3 v2 = V.at(f.Vids.at( (j + 2) % 4 ));
            double newJ = glm::cross(glm::normalize(v2 - v1), glm::normalize(v0 - v1)).z;
            if (newJ < jacobian)
                jacobian = newJ;
        }
        if (jacobian <= 0) {
            num++;
            invertedFIds.push_back(f.id);
        }
        if (jacobian < min) {
            min = jacobian;
            minIndex = i;
        }
        if (jacobian > max) max = jacobian;
        avg += jacobian;
    }
    avg /= F.size();
    std::cout << "-------------------------" << std::endl;
    std::cout << "minSJ: " << min << std::endl;
    std::cout << "avjSJ: " << avg << std::endl;
    std::cout << "maxSJ: " << max << std::endl;
    std::cout << "num inverted: " << num << std::endl;
    std::cout << "minID: " << minIndex << std::endl;
    for (int i = 0; i < invertedFIds.size(); i++) {
        std::cout << invertedFIds.at(i) << ", ";
    }
    std::cout << std::endl;
    std::cout << "-------------------------" << std::endl;
}

void Mesh::orderVerticesInFaces() {
    std::vector<int> checkLater;
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
            checkLater.push_back(f.id);
        } else if (numCounterClockwise < 0) {
            std::reverse(f.Vids.begin(), f.Vids.end());
        }
    }
    for (int i = 0; i < checkLater.size(); i++) {
        Face& f = F.at(checkLater.at(i));
        for (int j = 0; j < f.N_Fids.size(); j++) {
            Face& f1 = F.at(f.N_Fids.at(j));
            // f1 cannot be in checkLater, later add more guarantees
            if (std::find(checkLater.begin(), checkLater.end(), f1.id) != checkLater.end())
                continue;

            int ordering = findRelativeOrientation(f.Vids, f1.Vids);
            if (ordering != 0) {
                if (ordering == 1) {
                    std::reverse(f.Vids.begin(), f.Vids.end());
                }
                break;
            }
        }
    }
}

void Mesh::orderNeighboringVertices() {
    std::vector<int> checkOrderOfVertices;
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (v.isBoundary || v.N_Vids.size() == 0) {
            continue;
        }
        bool done = false;
        size_t currentFaceId = v.N_Fids.at(0);
        int secondFaceSelected = -1;
        std::vector<size_t> polygon;
        while (!done) {
            Face& f = F.at(currentFaceId);
            int neighboringVerticesInPolygon = 0;
            for (int j = 0; j < f.Vids.size(); j++) {
                Vertex& v1 = V.at(f.Vids.at(j));
                if (std::find(v.N_Vids.begin(), v.N_Vids.end(), v1.id) != v.N_Vids.end()) {//if in v.N_Vids
                    if (std::find(polygon.begin(), polygon.end(), v1.id) == polygon.end()) {//if not in polygon
                        polygon.push_back(v1.id);
                        for (int k = 0; k < v1.N_Fids.size(); k++) {//set next currentFaceId
                            Face& f1 = F.at(v1.N_Fids.at(k));
                            if (f1.id != currentFaceId && std::find(f1.Vids.begin(), f1.Vids.end(), v.id) != f1.Vids.end()) {
                                currentFaceId = f1.id;
                                if (secondFaceSelected == -1)
                                    secondFaceSelected = f1.id;
                                break;
                            }
                        }
                        break;
                    }
                    else {//if in polygon
                        neighboringVerticesInPolygon++;
                        if (neighboringVerticesInPolygon == 2) {
                            done = true;
                            break;
                        }
                    }
                }
            }
        }
        if (v.N_Vids.size() != polygon.size()) {
            std::cout << v.N_Vids.size() << " != " << polygon.size() << " for " << v.id << std::endl;
            std::cout << "ERROR\n";
        }
       
        //checking ordering
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
            // check later
            checkOrderOfVertices.push_back(v.id);
        }
        else if (numCounterClockwise < 0) {
            std::reverse(polygon.begin(), polygon.end());
        }  

        v.N_Vids.clear();
        for (int j = 0; j < polygon.size(); j++) {
            v.N_Vids.push_back(polygon.at(j));
        }
    }
    for (size_t i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isBoundary || v.N_Vids.size() == 0) {
            continue;
        }

        //vertices could be in random order
        //most general solution to order: find boundary edge, traverse through neighboring vertices, selecting one that is boundary with latest selected vertex
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
                //if v2Id not in polygon and previousVId is in vertices of neighboring faces, push to polygon and break
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
                std::cout << "ERROR: could not order neighboring vertices of boundary vertex " << v.id << std::endl;
                std::exit(0);
            }
        }

        // now that vertices are ordered, how to determine if vertex is corner and if polygon should be reversed
        // count number of counter-clockwise/clockwise vertices to select ordering. Then use ends of polygon to determine if vertex is corner
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
            // check later
            checkOrderOfVertices.push_back(v.id);
        }
        if (numCounterClockwise < 0) {
            std::reverse(polygon.begin(), polygon.end());
        }
        v.N_Vids.clear();
        for (int j = 0; j < polygon.size(); j++) {
            v.N_Vids.push_back(polygon.at(j));
        }

        Vertex& b1 = V.at(polygon.at(0));
        Vertex& b2 = V.at(polygon.at(polygon.size() - 1));
        if (glm::cross(b1 - v, b2 - v).z <= -0.3 * glm::length(b1 - v) * glm::length(b2 - v)) //tolerance for straight boundary edges. float * sin of angle is tolerance
            v.isCorner = true;
        else
            v.isCorner = false;
    }

    while (checkOrderOfVertices.size() != 0) {
        Vertex& v = V.at(checkOrderOfVertices.back());
        checkOrderOfVertices.pop_back();
        bool found = false;
        for (int j = 0; j < v.N_Fids.size(); j++) {
            Face& f = F.at(v.N_Fids.at(j));
            int localId = getIndexOf1(v.id, f.Vids);
            Vertex& v1 = V.at(f.Vids.at( (localId + 2) % 4 ));
            // if v1 is not in checkOrderOfVertices
            if (std::find(checkOrderOfVertices.begin(), checkOrderOfVertices.end(), v1.id) == checkOrderOfVertices.end()) {
                // returns 1 if ordered consistently, -1 if reversed, 0 if not.
                // should be reversed
                int ordering = findRelativeOrientation(v.N_Vids, v1.N_Vids);
                if (ordering == 0) {
                    std::cout << "ERROR: ordering == 0\n";
                    std::exit(0);
                }
                if (ordering == 1) {
                    std::cout << "Reversing N_V for vertex " << v.id << std::endl;
                    std::reverse(v.N_Vids.begin(), v.N_Vids.end());
                }
                found = true;
                break;
            }
        }
        if (!found) {
            std::cout << "ERROR: checkOrderOfVertices\n";
            std::exit(0);
            checkOrderOfVertices.push_back(v.id);
        }
    }
    classifyCornerVertices();

    // std::cout << "CORNER VERTICES\n";
    // for (int i = 0; i < V.size(); i++) {
    //     Vertex& v = V.at(i);
    //     if (v.isCorner)
    //         std::cout << v.id << ", ";
    // }
    // std::cout << std::endl;
}

void Mesh::classifyCornerVertices() {
    for (int i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isBoundary)
            continue;
        Vertex& b1 = V.at(v.N_Vids.at(0));
        Vertex& b2 = V.at(v.N_Vids.at(v.N_Vids.size() - 1));
        if (glm::cross(b1 - v, b2 - v).z <= -0.3 * glm::length(b1 - v) * glm::length(b2 - v)) //tolerance for straight boundary edges. float * sin of angle is tolerance
            v.isCorner = true;
        else
            v.isCorner = false;
    }
}

float Mesh::getAverageArea() {
    float sumOfAreas = 0;
    for (int i = 0; i < F.size(); i++) {
        Face& f = F.at(i);
        sumOfAreas += getAreaOfFace(f);
    }
    return sumOfAreas / F.size();
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

int getIndexOf1(int num, std::vector<unsigned long int> vec) {
    for (int i = 0; i < vec.size(); i++) {
        if (vec.at(i) == num)
            return i;
    }
    return -1;
}

Vertex Mesh::getAverageOfVIds(std::vector<size_t> points) {
    Vertex result;
    for (int i = 0; i < points.size(); i++) {
        result.x += V.at(points.at(i)).x / points.size();
        result.y += V.at(points.at(i)).y / points.size();
    }
    return result;
}

void Mesh::splitFace(int fId) {
    Face& f = F.at(fId);
    int boundaryVId = -1;
    int oppositeVId = -1;
    int localIndex = -1;
    for (int i = 0; i < f.Vids.size(); i++) {
        Vertex& v = V.at(f.Vids.at(i));
        if (v.type == 2 || (f.id == 51 && v.id == 59)) {
            if (boundaryVId != -1) {
                std::cout << "ERROR: splitFace called with face containing more than 1 type 2 vertex\n";
                std::exit(0);
            }
            boundaryVId = v.id;
            localIndex = i;
            oppositeVId = f.Vids.at((i + 2) % 4);
        }
    }
    if (boundaryVId == -1) {
        std::cout << "ERROR: splitFace called with face containing no type 2 vertices\n";
        std::exit(0);
    }

    Vertex& cornerV = V.at(boundaryVId);
    Vertex& oppositeV = V.at(oppositeVId);

    Vertex midpointV;
    midpointV.x = (cornerV.x + oppositeV.x) / 2;
    midpointV.y = (cornerV.y + oppositeV.y) / 2;
    midpointV.id = V.size();
    V.push_back(midpointV);

    //face will have points localIndex, localIndex + 1, opposite, midpoint
    //new face will have points localIndex, midpoint, opposite, localIndex + 3
    Cell& c = C.at(fId);
    c.Vids = { (size_t)boundaryVId, f.Vids.at( (localIndex + 1) % 4 ), oppositeV.id, midpointV.id };
    Cell newC;
    newC.id = C.size();
    newC.Vids = { (size_t)boundaryVId, midpointV.id, oppositeV.id, f.Vids.at( (localIndex + 3) % 4 ) };
    C.push_back(newC);
}

void Mesh::expandFace(int fId) {
    //find average, find midpoints to average, add correct points, add correct cells
    Face& f = F.at(fId);
    Vertex average = getAverageOfVIds(f.Vids);
    size_t baseIndex = V.size();
    for (int i = 0; i < f.Vids.size(); i++) {
        Vertex& v = V.at(f.Vids.at(i));
        Vertex midPoint;
        midPoint.id = baseIndex + i;
        midPoint.x = (v.x + average.x) / 2;
        midPoint.y = (v.y + average.y) / 2;
        V.push_back(midPoint);
    }
    for (int i = 0; i < f.Vids.size(); i++) {
        Cell newC;
        newC.id = C.size();
        newC.Vids = {f.Vids.at(i), f.Vids.at( (i + 1) % 4), baseIndex + ((i + 1) % 4), baseIndex + i};
        C.push_back(newC);
    }

    Cell& c = C.at(fId);
    c.Vids = {baseIndex, baseIndex + 1, baseIndex + 2, baseIndex + 3};
}

std::vector<int> intersectVectors(std::vector<unsigned long int> vec1, std::vector<unsigned long int> vec2) {
    std::vector<int> result;
    for (int i = 0; i < vec1.size(); i++) {
        for (int j = 0; j < vec2.size(); j++) {
            if (vec1.at(i) == vec2.at(j)) {
                result.push_back(vec1.at(i));
                break;
            }
        }
    }
    return result;
}

void Mesh::collapseEdge(int eId) {
    std::vector<size_t> visitedEdges;
    std::vector<size_t> visitedFaces;
    std::map<int, int> vertexMap;
    std::queue<int> edgeQueue;
    edgeQueue.push(eId);
    while (!edgeQueue.empty()) {
        Edge& e = E.at(edgeQueue.front());
        visitedEdges.push_back(e.id);
        edgeQueue.pop();
        for (int i = 0; i < e.N_Fids.size(); i++) {
            Face& f = F.at(e.N_Fids.at(i));
            if (getIndexOf1(f.id, visitedFaces) == -1)
                visitedFaces.push_back(f.id);
            for (int j = 0; j < f.Eids.size(); j++) {
                Edge& newEdge = E.at(f.Eids.at(j));
                if (intersectVectors(e.Vids, newEdge.Vids).size() == 0) {
                    if (getIndexOf1(newEdge.id, visitedEdges) == -1) {
                        edgeQueue.push(newEdge.id);
                        break;
                    }
                }
            }
        }
    }

    std::sort(visitedFaces.begin(), visitedFaces.end());
    for (int i = visitedFaces.size() - 1; i >= 0; i--) {
        Face& f = F.at(visitedFaces.at(i));
        std::cout << "Visited face " << f.id << std::endl;
        C.erase(C.begin() + f.id);
    }
    for (int i = 0; i < visitedEdges.size(); i++) {
        Edge& e = E.at(visitedEdges.at(i));
        vertexMap[e.Vids.at(1)] = e.Vids.at(0);
    }
    for (int i = 0; i < C.size(); i++) {
        Cell& c = C.at(i);
        for (int j = 0; j < c.Vids.size(); j++) {
            Vertex& v = V.at(c.Vids.at(j));
            if (vertexMap.find(v.id) != vertexMap.end()) {
                c.Vids.at(j) = vertexMap[v.id];
            }
        }
    }
}

std::vector<int> Mesh::selectBoundaryForIterativeOptimization() {
    return selectBoundaryAsNFacesFromCornerVertices(2);
}

std::vector<int> Mesh::selectBoundrayAsNFacesFromInversion(int n) {
    std::vector<int> facesToFree;
    for (int i = 0; i < F.size(); i++) {
        Face& f = F.at(i);
        if (isFaceInverted(f.id)) {
            facesToFree.push_back(f.id);
        }
    }

    //expand faces to free by one ring
    for (int i = 0; i < n - 1; i++) {
        std::vector<int> newFacesToFree;
        for (int j = 0; j < F.size(); j++) {
            Face& f = F.at(j);
            if (std::find(facesToFree.begin(), facesToFree.end(), f.id) != facesToFree.end())
                continue;
            // f is not in facesToFree
            bool neighborsFreeFace = false;
            for (int k = 0; k < f.N_Fids.size(); k++) {
                Face& f1 = F.at(f.N_Fids.at(k));

                //make sure f and f1 share an edge
                bool shareEdge = false;
                for (int m = 0; m < f.Eids.size(); m++) {
                    for (int n = 0; n < f1.Eids.size(); n++) {
                        if (f.Eids.at(m) == f1.Eids.at(n)) {
                            shareEdge = true;
                            break;
                        }
                    }
                    if (shareEdge)
                        break;
                }
                if (!shareEdge)
                    continue;

                if (std::find(facesToFree.begin(), facesToFree.end(), f1.id) != facesToFree.end()) {
                    // f1 is in facesToFree
                    neighborsFreeFace = true;
                    break;
                }
            }
            if (neighborsFreeFace)
                newFacesToFree.push_back(f.id);
        }
        for (int j = 0; j < newFacesToFree.size(); j++) {
            facesToFree.push_back(newFacesToFree.at(j));
        }
    }

    std::vector<int> verticesToFree;
    for (int i = 0; i < facesToFree.size(); i++) {
        Face& f = F.at(facesToFree.at(i));
        for (int j = 0; j < f.Vids.size(); j++) {
            Vertex& v = V.at(f.Vids.at(j));
            if (v.isBoundary)
                continue;
            if (std::find(verticesToFree.begin(), verticesToFree.end(), v.id) == verticesToFree.end()) {
                verticesToFree.push_back(v.id);
            }
        }
    }

    // reset boundary/free vertices
    for (int i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        v.isBoundary = true;
    }
    // std::cout << "Freeing vertices: " << std::endl;
    for (int i = 0; i < verticesToFree.size(); i++) {
        Vertex& v = V.at(verticesToFree.at(i));
        // std::cout << v.id << ", ";
        v.isBoundary = false;
    }
    // std::cout << std::endl;
    return verticesToFree;
}

std::vector<int> Mesh::selectBoundaryAsNeighboringFacesFromInversion() {
    std::vector<int> verticesToFree;
    for (int i = 0; i < F.size(); i++) {
        Face& f = F.at(i);
        if (isFaceInverted(f.id)) {
            for (int j = 0; j < f.Vids.size(); j++) {
                Vertex& v = V.at(f.Vids.at(j));
                if (v.isBoundary)
                    continue;
                if (std::find(verticesToFree.begin(), verticesToFree.end(), v.id) == verticesToFree.end()) {
                    verticesToFree.push_back(v.id);
                }
            }
        }
    }

    // reset boundary/free vertices
    for (int i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        v.isBoundary = true;
    }
    // std::cout << "Freeing vertices: " << std::endl;
    for (int i = 0; i < verticesToFree.size(); i++) {
        Vertex& v = V.at(verticesToFree.at(i));
        // std::cout << v.id << ", ";
        v.isBoundary = false;
    }
    // std::cout << std::endl;
    return verticesToFree;
}

bool Mesh::isFaceInverted(int fId) {
    Face& f = F.at(fId);
    for (int i = 0; i < f.Vids.size(); i++) {
        Vertex& v = V.at(f.Vids.at(i));
        Vertex& v1 = V.at(f.Vids.at( (i + 1) % f.Vids.size() ));
        Vertex& v2 = V.at(f.Vids.at( (f.Vids.size() + i - 1) % f.Vids.size() ));
        if (glm::cross(v1 - v, v2 - v).z <= 0) {
            return true;
        }
    }
    return false;
}

std::vector<int> Mesh::selectBoundaryAsTwoFacesFromCornerVertices() {
    return selectBoundaryAsNFacesFromCornerVertices(2);
}

std::vector<std::vector<int>> Mesh::getLocalFaceRegions() {
    std::vector<std::vector<int>> result;
    for (int i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isCorner)
            continue;
        std::vector<int> localFaceRegion;
        for (int j = 0; j < v.N_Fids.size(); j++) {
            localFaceRegion.push_back(v.N_Fids.at(j));
        }

        for (int j = 0; j < 1; j++) {
            expandFaceRegionByOneLayer(localFaceRegion);
        }

        double minSJOfRegion = getMinimumScaledJacobian(localFaceRegion);
        if (minSJOfRegion <= 0)
            result.push_back(localFaceRegion);
    }
    for (int i = 0; i < F.size(); i++) {
        Face& f = F.at(i);
        if (getMinimumScaledJacobian(f) > 0)
            continue;
        std::vector<int> localFaceRegion;
        localFaceRegion.push_back(f.id);
        expandFaceRegionByOneLayer(localFaceRegion);
        // expandFaceRegionByOneLayer(localFaceRegion);
        result.push_back(localFaceRegion);
    }

    result = mergeRegions(result);

    //remove unsolvable faces
    for (int i = result.size() - 1; i >= 0; i--) {
        std::vector<int>& localRegion = result.at(i);
        for (int j = localRegion.size() - 1; j >= 0; j--) {
            Face& f = F.at(localRegion.at(j));
            double faceArea = getAreaOfFace(f);
            if (std::abs(faceArea) <= 0.00001) {
                localRegion.erase(localRegion.begin() + j);
            }
        }
        if (localRegion.size() == 0)
            result.erase(result.begin() + i);
    }
    return result;
}

// 5, 16, 140, 170, 227, 232, 278, 440, 14, 194, 242, 261, 276, 290, 294, 419, 474, 178, 194, 276, 290, 294, 88, 158, 211, 214, 473, 20, 81, 88, 152, 158, 211, 248, 473,
// 19, 161, 184, 305, 353, 13, 143, 154, 393, 395, 17, 154, 391,
// 12, 36, 82, 235, 329, 342, 351, 358, 386, 84, 174, 231, 296, 432, 475, 14, 112, 190, 233, 242, 282, 434, 112, 190, 221, 434,
// 80, 98, 148, 182, 200, 293, 347, 455, 245, 274, 388, 420, 445, 456, 469,
// 159, 172, 186, 195, 266, 421

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
                std::vector<int> verticiesInCurrent = collapseFacesToVerticies(currentFaceRegion);
                std::vector<int> otherFaceRegion = localFaceRegions.at(j);
                std::vector<int> verticiesInOther = collapseFacesToVerticies(otherFaceRegion);
                std::vector<int> vertexIntersection;
                std::set_intersection(verticiesInCurrent.begin(), verticiesInCurrent.end(), verticiesInOther.begin(), verticiesInOther.end(), std::back_inserter(vertexIntersection));
                if (vertexIntersection.size() != 0) {
                    merged = true;
                    // verticies overlap, add elements of otherFaceRegion not in currentFaceRegion to currentFaceRegion
                    std::vector<int> faceDifference;
                    std::set_difference(otherFaceRegion.begin(), otherFaceRegion.end(), currentFaceRegion.begin(), currentFaceRegion.end(), std::back_inserter(faceDifference));
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

void Mesh::expandFaceRegionByOneLayer(std::vector<int>& localFaceRegion) {
    std::vector<int> newFacesToAdd;
    for (int k = 0; k < F.size(); k++) {
        Face& f = F.at(k);
        if (std::find(localFaceRegion.begin(), localFaceRegion.end(), f.id) != localFaceRegion.end())
            continue;
        // f is not in localFaceRegion
        bool neighborsFaceRegion = false;
        for (int l = 0; l < f.N_Fids.size(); l++) {
            Face& f1 = F.at(f.N_Fids.at(l));

            //make sure f and f1 share an edge
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
                // f1 is in localFaceRegion
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

std::vector<int> Mesh::selectBoundaryAsNFacesFromCornerVertices(int n) {
    std::vector<int> facesToFree;
    for (int i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        if (!v.isBoundary)
            continue;
        int v1Id = -1;
        int v2Id = -1;
        for (int j = 0; j < v.N_Vids.size(); j++) {
            Vertex& v1 = V.at(v.N_Vids.at(j));
            if (!v1.isBoundary)
                continue;
            if (v1Id == -1)
                v1Id = v1.id;
            else if (v2Id == -1)
                v2Id = v1.id;
            else {
                v1Id = -1;
                v2Id = -1;
                break;
            }
        }
        if (v1Id == -1 || v2Id == -1)
            continue;
        // for each boundary vertex with two neighboring boundary vertices
        Vertex& boundaryV1 = V.at(v1Id);
        Vertex& boundaryV2 = V.at(v2Id);
        if (v.isCorner) {
            for (int j = 0; j < v.N_Fids.size(); j++) {
                facesToFree.push_back(v.N_Fids.at(j));
            }
        }
    }

    //expand faces to free by one ring
    for (int i = 0; i < n - 1; i++) {
        std::vector<int> newFacesToFree;
        for (int j = 0; j < F.size(); j++) {
            Face& f = F.at(j);
            if (std::find(facesToFree.begin(), facesToFree.end(), f.id) != facesToFree.end())
                continue;
            // f is not in facesToFree
            bool neighborsFreeFace = false;
            for (int k = 0; k < f.N_Fids.size(); k++) {
                Face& f1 = F.at(f.N_Fids.at(k));

                //make sure f and f1 share an edge
                bool shareEdge = false;
                for (int m = 0; m < f.Eids.size(); m++) {
                    for (int n = 0; n < f1.Eids.size(); n++) {
                        if (f.Eids.at(m) == f1.Eids.at(n)) {
                            shareEdge = true;
                            break;
                        }
                    }
                    if (shareEdge)
                        break;
                }
                if (!shareEdge)
                    continue;

                if (std::find(facesToFree.begin(), facesToFree.end(), f1.id) != facesToFree.end()) {
                    // f1 is in facesToFree
                    neighborsFreeFace = true;
                    break;
                }
            }
            if (neighborsFreeFace)
                newFacesToFree.push_back(f.id);
        }
        for (int j = 0; j < newFacesToFree.size(); j++) {
            facesToFree.push_back(newFacesToFree.at(j));
        }
    }

    // select all inner vertices of facesToFree as free vertices
    std::vector<int> verticesToFree;
    std::vector<int> verticesToIgnoreDELETE = {496, 488, 487, 505, 507, 480, 481, 467, 468};
    for (int i = 0; i < facesToFree.size(); i++) {
        Face& f = F.at(facesToFree.at(i));
        for (int j = 0; j < f.Vids.size(); j++) {
            Vertex& v = V.at(f.Vids.at(j));
            if (v.isBoundary)
                continue;
            if (std::find(verticesToFree.begin(), verticesToFree.end(), v.id) == verticesToFree.end()) {
                if (std::find(verticesToIgnoreDELETE.begin(), verticesToIgnoreDELETE.end(), v.id) == verticesToIgnoreDELETE.end())
                    verticesToFree.push_back(v.id);
            }
        }
    }

    // reset boundary/free vertices
    for (int i = 0; i < V.size(); i++) {
        Vertex& v = V.at(i);
        v.isBoundary = true;
    }
    std::cout << "Freeing vertices: " << std::endl;
    for (int i = 0; i < verticesToFree.size(); i++) {
        Vertex& v = V.at(verticesToFree.at(i));
        std::cout << v.id << ", ";
        v.isBoundary = false;
    }
    std::cout << std::endl;
    return verticesToFree;
}

int Mesh::getDistanceFromCorner(std::vector<int> vIds, int currentDistance) {
    std::vector<int> newVerticesToAdd;
    for (int i = 0; i < vIds.size(); i++) {
        Vertex& v = V.at(vIds.at(i));
        for (int j = 0; j < v.N_Vids.size(); j++) {
            Vertex& nV = V.at(v.N_Vids.at(j));
            newVerticesToAdd.push_back(nV.id); // inefficiency is okay for now
            if (nV.isCorner) {
                return currentDistance + 1;
            }
        }
    }

    for (int i = 0; i < newVerticesToAdd.size(); i++) {
        vIds.push_back(newVerticesToAdd.at(i));
    }
    return getDistanceFromCorner(vIds, currentDistance + 1);
}

int findRelativeOrientation(std::vector<size_t> vec1, std::vector<size_t> vec2) {
    // first find matching element
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

void Mesh::fixInvertedFaces() {
    std::cout << "-----------------\n";
    std::cout << "in fixInvertedFaces()\n";
    // find list of critical vertices
    std::vector<int> criticalVertexIds;
    for (int i = 0; i < F.size(); i++) {
        Face& f = F.at(i);
        int cornerVId = -1;
        int secondBoundaryVId = -1;
        int numOfBoundaryVertices = 0;
        for (int j = 0; j < f.Vids.size(); j++) {
            Vertex& v = V.at(f.Vids.at(j));
            if (v.isBoundary) {
                numOfBoundaryVertices++;
                if (v.isCorner) {
                    cornerVId = v.id;
                }
                else {
                    secondBoundaryVId = v.id;
                }
            }
        }
        if (cornerVId == -1)
            continue;
        
        if (numOfBoundaryVertices == 2 && cornerVId != -1 && secondBoundaryVId != -1) {
            // for half pull, one vertex past straight edge
            Vertex& v = V.at(cornerVId);
            Vertex& v1 = V.at(secondBoundaryVId);
            int v2Id = f.Vids.at((getIndexOf1(v1.id, f.Vids) + 2) % 4);
            Vertex& v2 = V.at(v2Id);

            std::vector<size_t> vertices = {v2.id, v.id, v1.id};
            bool isOrientedCounterClockwise = (findRelativeOrientation(f.Vids, vertices) == 1);
            if ((glm::cross(v1 - v, v2 - v).z <= 0) == isOrientedCounterClockwise) {
                std::cout << v2.id << ", ";
                // find new optimal position:
                float length = glm::length(v2 - v);
                glm::vec3 delta = (v - v1) * length / glm::length(v - v1);
                glm::vec3 optimalV;
                float angle = 10.0 * (3.14159265 / 180);
                if (isOrientedCounterClockwise)
                    angle *= -1;
                optimalV.x = v.x + std::cos(angle) * delta.x - std::sin(angle) * delta.y;
                optimalV.y = v.y + std::sin(angle) * delta.x + std::cos(angle) * delta.y;
                v2.x = optimalV.x;
                v2.y = optimalV.y;
                v2.isBoundary = true;
            }
        }
        if (numOfBoundaryVertices == 1) {
            // for full pull, angle at corner is inverted
        }
    }
    std::cout << "\nout of fixInvertedFaces()\n";
    std::cout << "-----------------\n";
}

void Mesh::getiterativeCornerEnergy() {
}