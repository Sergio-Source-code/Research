/*
 * Mesh.h
 *
 *  Created on: Nov 6, 2016
 *      Author: cotrik
 *  Edited by: sergio
 */

#ifndef MESH_H_
#define MESH_H_

#include <vector>
#include <string>
#include <glm/glm.hpp>

extern const size_t MAXID;
extern const unsigned int QuadEdge[4][2];

enum ElementType {
    POLYGON,
    TRIANGLE,
    QUAD,
    TETRAHEDRA,
    HEXAHEDRA
};

class NeighborInfo {
public:
    NeighborInfo(){}
    NeighborInfo(const NeighborInfo& rhs)
    : N_Vids(rhs.N_Vids)
    , N_Eids(rhs.N_Eids)
    , N_Fids(rhs.N_Fids)
    , N_Cids(rhs.N_Cids)
    {}
    NeighborInfo& operator = (const NeighborInfo& rhs) {
        N_Vids = rhs.N_Vids;
        N_Eids = rhs.N_Eids;
        N_Fids = rhs.N_Fids;
        N_Cids = rhs.N_Cids;
        return *this;
    }
    ~NeighborInfo(){}
    std::vector<size_t> N_Vids;
    std::vector<size_t> N_Eids;
    std::vector<size_t> N_Fids;
    std::vector<size_t> N_Cids;
};

struct GeoInfo {
public:
    GeoInfo()
    : id(MAXID)
    , isBoundary(false)
    , isSingularity(false)
    {}
    GeoInfo(const GeoInfo& r)
    : id(r.id)
    , isBoundary(r.isBoundary)
    , isSingularity(r.isSingularity)
    {}
    GeoInfo& operator = (const GeoInfo& rhs) {
        id = rhs.id;
        isBoundary = rhs.isBoundary;
        isSingularity = rhs.isSingularity;
        return *this;
    }
    virtual ~GeoInfo()
    {}
public:
    size_t id;
    bool isBoundary;
    bool isSingularity;
};

enum SurfaceVertexType {
    REGULAR = 0,
    FEATURE,
    CORNER
};

class Vertex : public glm::vec3, public GeoInfo, public NeighborInfo {
public:
    Vertex()
    : hvid(MAXID)
    , triVid(MAXID)
    , type(MAXID)
    , isCorner(false)
    {}
    Vertex(const Vertex& r)
    : glm::vec3(r)
    , GeoInfo(r)
    , NeighborInfo(r)
    , normal(r.normal)
    , tangent(r.tangent)
    , hvid(r.hvid)
    , triVid(r.triVid)
    , type(r.type)
    , isCorner(r.isCorner)
    {}
    Vertex(const glm::vec3& v)
    : glm::vec3(v)
    , hvid(MAXID)
    , triVid(MAXID)
    , type(MAXID)
    , isCorner(false)
    {}
    virtual ~Vertex()
    {}
    Vertex& operator = (const Vertex& r) {
        if (r == *this)
            return *this;
        x = r.x;
        y = r.y;
        z = r.z;

        return *this;
    }
    Vertex& operator = (const glm::vec3& r) {
        if (r == *this)
            return *this;
        x = r.x;
        y = r.y;
        z = r.z;

        return *this;
    }
    glm::vec3 xyz() const {
        return glm::vec3(x,y,z);
    }
    glm::vec3 normal;

    glm::vec3 tangent;
    size_t hvid;
    size_t triVid;
    size_t type;
    bool isCorner;
    std::vector<size_t> twoRingNeighborSurfaceFaceIds;
};

class Edge : public GeoInfo, public NeighborInfo {
public:
    Edge()
    : length(0.0)
    , energySingularity(0.0)
    , energyOrthogonality(0.0)
    , energyStraightness(0.0)
    , face_angle(0.0)
    , isSharpFeature(false)
    , label(MAXID)
    {}
    Edge(const Edge& r)
    : GeoInfo(r)
    , NeighborInfo(r)
    , Vids(r.Vids)
    , parallelEids(r.parallelEids)
    , consecutiveEids(r.consecutiveEids)
    , orthogonalEids(r.orthogonalEids)
    , length(r.length)
    , energySingularity(r.energySingularity)
    , energyOrthogonality(r.energyOrthogonality)
    , energyStraightness(r.energyStraightness)
    , face_angle(r.face_angle)
    , isSharpFeature(r.isSharpFeature)
    , label(r.label)
    {}
    Edge(size_t vnum)
    : length(0.0)
    , energySingularity(0.0)
    , energyOrthogonality(0.0)
    , energyStraightness(0.0)
    , face_angle(0.0)
    , isSharpFeature(false)
    , label(MAXID) {
        Vids.resize(vnum);
    }
    Edge(const std::vector<size_t>& Vids)
    : Vids(Vids)
    , length(0.0)
    , energySingularity(0.0)
    , energyOrthogonality(0.0)
    , energyStraightness(0.0)
    , face_angle(0.0)
    , isSharpFeature(false)
    , label(MAXID)
    {}
    virtual ~Edge()
    {}
public:
    bool operator == (const Edge& e) const {
        return ((Vids[0] == e.Vids[0] && Vids[1] == e.Vids[1]) ||
                (Vids[0] == e.Vids[1] && Vids[1] == e.Vids[0]) );
    }
public:
    std::vector<size_t> Vids;
    std::vector<size_t> parallelEids;
    std::vector<size_t> consecutiveEids;
    std::vector<size_t> orthogonalEids;
    double length;
    double energySingularity;
    double energyOrthogonality;
    double energyStraightness;

    double face_angle;

    bool isSharpFeature;
    size_t label;
};

class Face : public GeoInfo, public NeighborInfo {
public:
    Face()
    : label(MAXID)
    {}
    Face(const Face& r)
    : GeoInfo(GeoInfo(r))
    , NeighborInfo(NeighborInfo(r))
    , Vids(r.Vids)
    , Eids(r.Eids)
    , normal(r.normal)
    , label(r.label)
    {}
    Face(size_t vnum)
    : label(MAXID) {
        Vids.resize(vnum);
    }
    Face(size_t vnum, size_t eNum)
    : label(MAXID) {
        Vids.resize(vnum);
        Eids.resize(eNum);
    }
    Face(const std::vector<size_t>& Vids)
    : Vids(Vids)
    , label(MAXID)
    {}

    Face(const std::vector<size_t>& Vids, const std::vector<size_t> Eids)
    : Vids(Vids)
    , Eids(Eids)
    , label(MAXID)
    {}
    virtual ~Face()
    {}

public:
    std::vector<size_t> Vids;
    std::vector<size_t> Eids;
    std::vector<std::vector<size_t> > N_Ortho_4Vids;
    std::vector<std::vector<size_t> > N_Ortho_4Eids;

    glm::vec3 normal;
    size_t label;
};

class Cell : public GeoInfo, public NeighborInfo {
public:
    Cell()
    {}
    Cell(const Cell& r)
    : GeoInfo(GeoInfo(r))
    , NeighborInfo(NeighborInfo(r))
    , Vids(r.Vids)
    , Eids(r.Eids)
    , Fids(r.Fids)
    {}
    Cell(size_t vnum) {
        Vids.resize(vnum);
    }
    Cell(size_t vnum, size_t eNum, size_t fNum) {
        Vids.resize(vnum);
        Eids.resize(eNum);
        Fids.resize(fNum);
    }
    Cell(const std::vector<size_t>& Vids)
    : Vids(Vids)
    {}

    Cell(const std::vector<size_t>& Vids, const std::vector<size_t> Eids, const std::vector<size_t> Fids)
    : Vids(Vids)
    , Eids(Eids)
    , Fids(Fids)
    {}

    virtual ~Cell()
    {}

public:
    std::vector<size_t> Vids;
    std::vector<size_t> Eids;
    std::vector<size_t> Fids;
};

class Mesh {
public:
    Mesh();
    Mesh(const Mesh& r);
    Mesh(const std::vector<Vertex>& V, const std::vector<Cell>& C, ElementType m_cellType);
    Mesh(const std::vector<Vertex>& V, const std::vector<Face>& F, ElementType m_cellType = QUAD);
    Mesh(const Mesh& r, const std::vector<size_t>& cellIds);
    virtual ~Mesh();

    void process();
    void printQuality();
    int getNumOfInvertedElements(std::vector<int> localFaceRegion);
    double getMinimumScaledJacobian(std::vector<int> faceIds);
    float getMaximumEdgeLength(Face& f);
    float getMinimumEdgeLength(Vertex& v);
    std::vector<int> getInnerVerticesOfFaces(std::vector<int> fIds);
    int getLocalIndexOf(int vId, std::vector<size_t> fVIds);
    int getRelativeOrientation(std::vector<size_t> vec1, std::vector<size_t> vec2);

    std::vector<std::vector<int>> getLocalFaceRegions();
    double getIterativeEnergyOfRegion(std::vector<int> faceIds);

private:
    void BuildAllConnectivities();
    void BuildV_C();
    void BuildE();
    void BuildF();
    void BuildF_F();
    void ExtractBoundary();
    void orderNeighboringVertices();
    void orderVerticesInFaces();
    void classifyCornerVertices();
    float getAreaOfFace(Face& f);

    std::vector<std::vector<int>> mergeRegions(std::vector<std::vector<int>> localFaceRegions);
    void expandFaceRegionByOneLayer(std::vector<int>& localFaceRegion);
    std::vector<int> collapseFacesToVerticies(std::vector<int> fIds);
    void removeUnsolvableFaces(std::vector<int>& localFaceRegion);
    bool hasUnsolvableInversion(std::vector<int> fIds);
    bool hasUnsolvableInversion(Face& f);
    double getMinimumScaledJacobian(Face& f);

public:
    std::vector<Vertex> V;
    std::vector<Edge> E;
    std::vector<Face> F;
    std::vector<Cell> C;
    ElementType m_cellType;

    std::vector<std::string> pointScalarFieldNames;
    std::vector<std::string> cellScalarNameFields;
    std::vector<std::vector<double> > pointScalarFields;
    std::vector<std::vector<double> > cellScalarFields;
};

#endif /* MESH_H_ */
