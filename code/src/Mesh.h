/*
 * Mesh.h
 *
 *  Created on: Nov 6, 2016
 *      Author: cotrik
 */

#ifndef MESH_H_
#define MESH_H_

#include <vector>
#include <string>
#include <glm/glm.hpp>
class FeatureLine;
enum ElementType
{
    POLYGON,
    TRIANGLE,
    QUAD,
    TETRAHEDRA,
    HEXAHEDRA
};

extern const size_t MAXID;
extern const unsigned int HexEdge[12][2];
extern const unsigned int HexFaces[6][4];

extern const unsigned int TetEdge[6][2];
extern const unsigned int TetFaces[4][3];

extern const unsigned int TriEdge[3][2];
extern const unsigned int QuadEdge[4][2];

extern const unsigned int HexPoint_Points[8][3];
extern const unsigned int HexPoint_Edges[8][3];

class NeighborInfo
{
public:
    NeighborInfo(){}
    NeighborInfo(const NeighborInfo& rhs)
    : N_Vids(rhs.N_Vids)
    , N_Eids(rhs.N_Eids)
    , N_Fids(rhs.N_Fids)
    , N_Cids(rhs.N_Cids)
    {}
    NeighborInfo& operator = (const NeighborInfo& rhs)
    {
        N_Vids = rhs.N_Vids;
        N_Eids = rhs.N_Eids;
        N_Fids = rhs.N_Fids;
        N_Cids = rhs.N_Cids;
        return *this;
    }
    ~NeighborInfo(){}
    std::vector<size_t> N_Vids;  // neighboring vertices ids
    std::vector<size_t> N_Eids;  // neighboring edges ids
    std::vector<size_t> N_Fids;  // neighboring faces ids
    std::vector<size_t> N_Cids;  // neighboring cells ids
};

struct GeoInfo
{
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
    GeoInfo& operator = (const GeoInfo& rhs)
    {
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

enum SurfaceVertexType
{
    REGULAR = 0,
    FEATURE,
    CORNER
};

enum SmoothMethod
{
    LAPLACE_EDGE = 0,
    LAPLACE_FACE_CENTER,
    LAPLACE_FACE_CENTER_ANGLE,
    LAPLACE_BELTRAMI
};

class Vertex : public glm::vec3, public GeoInfo, public NeighborInfo
{
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
    Vertex& operator = (const Vertex& r)
    {
        if (r == *this)
            return *this;
        x = r.x;
        y = r.y;
        z = r.z;

        return *this;
    }
    Vertex& operator = (const glm::vec3& r)
    {
        if (r == *this)
            return *this;
        x = r.x;
        y = r.y;
        z = r.z;

        return *this;
    }
    glm::vec3 xyz() const
    {
        return glm::vec3(x,y,z);
    }
    glm::vec3 normal;

    // for edge cone triangle surface
    glm::vec3 tangent;
    size_t hvid;
    size_t triVid;
    size_t type; //0 regular, 1 feature, 2 corner.
    bool isCorner;
    std::vector<size_t> twoRingNeighborSurfaceFaceIds; // for surface projecting;
};

class Edge : public GeoInfo, public NeighborInfo
{
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
    , label(MAXID)
    {
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
    {
    }
    virtual ~Edge()
    {}
public:
    bool operator == (const Edge& e) const
    {
        return ((Vids[0] == e.Vids[0] && Vids[1] == e.Vids[1]) ||
                (Vids[0] == e.Vids[1] && Vids[1] == e.Vids[0]) );
    }
public:
    std::vector<size_t> Vids;
    std::vector<size_t> parallelEids;
    std::vector<size_t> consecutiveEids;
    std::vector<size_t> orthogonalEids;
    double length;                             // target length that is used in frameOpt
    double energySingularity;
    double energyOrthogonality;
    double energyStraightness;

    // for edge cone triangle surface
    double face_angle;

    bool isSharpFeature;
    size_t label;     // patch number starting from 0, MAXID is invalid
};

class Face : public GeoInfo, public NeighborInfo
{
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
    : label(MAXID)
    {
        Vids.resize(vnum);
    }
    Face(size_t vnum, size_t eNum)
    : label(MAXID)
    {
        Vids.resize(vnum);
        Eids.resize(eNum);
    }
    Face(const std::vector<size_t>& Vids)
    : Vids(Vids)
    , label(MAXID)
    {
    }

    Face(const std::vector<size_t>& Vids, const std::vector<size_t> Eids)
    : Vids(Vids)
    , Eids(Eids)
    , label(MAXID)
    {
    }
    virtual ~Face()
    {}

public:
    std::vector<size_t> Vids;
    std::vector<size_t> Eids;
    std::vector<std::vector<size_t> > N_Ortho_4Vids;   // Neighboring Orthogonal 4 Vids of Edges
    std::vector<std::vector<size_t> > N_Ortho_4Eids;   // Neighboring Orthogonal 4 Eids of Edges

    glm::vec3 normal;
    size_t label;     // patch number starting from 0, MAXID is invalid
};

class Cell : public GeoInfo, public NeighborInfo
{
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
    Cell(size_t vnum)
    {
        Vids.resize(vnum);
    }
    Cell(size_t vnum, size_t eNum, size_t fNum)
    {
        Vids.resize(vnum);
        Eids.resize(eNum);
        Fids.resize(fNum);
    }
    Cell(const std::vector<size_t>& Vids)
    : Vids(Vids)
    {
    }

    Cell(const std::vector<size_t>& Vids, const std::vector<size_t> Eids, const std::vector<size_t> Fids)
    : Vids(Vids)
    , Eids(Eids)
    , Fids(Fids)
    {
    }

    virtual ~Cell()
    {}

public:
    std::vector<size_t> Vids;
    std::vector<size_t> Eids;
    std::vector<size_t> Fids;
    //glm::vec3 cv;                 // center vertex's x,y,z coordinate of the Cell;  For Node of Frame
};

class Layer
{
public:
    Layer(){}
    Layer(const Layer& layer)
    : Vids(layer.Vids)
    , Eids(layer.Eids)
    , Fids(layer.Fids)
    , Cids(layer.Cids)
    , fixed(layer.fixed)
    {}
    ~Layer(){}
    std::vector<size_t> Vids;
    std::vector<size_t> Eids;
    std::vector<size_t> Fids;
    std::vector<size_t> Cids;
    std::vector<bool>  fixed;  // refer to this->Vids;
};

class Plane
{
public:
    Plane(const double a, const double b, const double c, const double d)
        : a(a), b(b), c(c), d(d)
    {

    }

    Plane(const Vertex& p1, const Vertex& p2, const Vertex& p3)
    {
        a = ((p2.y - p1.y) * (p3.z - p1.z) - (p2.z - p1.z) * (p3.y - p1.y));
        b = ((p2.z - p1.z) * (p3.x - p1.x) - (p2.x - p1.x) * (p3.z - p1.z));
        c = ((p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x));
        d = (0 - (a * p1.x + b * p1.y + c * p1.z));
    }

    Plane(const Plane& plane)
    : a(plane.a), b(plane.b), c(plane.c), d(plane.d)
    {
    }

    double IntersectionAngle(const Plane& p) const
    {
        double cosangle = (a*p.a + b*p.b + c*p.c) / (sqrt(a*a + b*b + c*c)* sqrt(p.a*p.a + p.b*p.b + p.c*p.c));
        return cosangle;
    }

    double DistanseFromPoint(const glm::vec3& p, glm::vec3& intersection) const
    {
        const double distance = fabs(a*p.x + b*p.y + c*p.z + d) / sqrt(a*a + b*b + c*c);
        const double t = (a*p.x + b*p.y + c*p.z + d) / (a*a + b*b + c*c);
        intersection = glm::vec3(p.x - a*t, p.y - b*t, p.z - c*t);
        return distance;
    }

    ~Plane(){}
    bool operator < (const Plane& right) const
    {
        return d < right.d;
    }

public:
    double a;
    double b;
    double c;
    double d;
};

class Line
{
public:
    Line(const glm::vec3& o, const glm::vec3& d)
        : o(o), dir(d)
    {
    }

    Line(const Vertex& p1, const Vertex& p2)
        : o(p1.xyz()), dir(p2.xyz() - p1.xyz())
    {

    }

    Line(const Line& r)
    : o(r.o), dir(r.dir)
    {
    }

    double Perpendicular(const glm::vec3& a, glm::vec3& intersection) const
    {
        const double t = ((a.x - o.x) * dir.x + (a.y - o.y) * dir.y + (a.z - o.z) * dir.z) / (dir.x * dir.x + dir.y * dir.y + dir.z * dir.z);
        intersection.x = o.x + dir.x * t;
        intersection.y = o.y + dir.y * t;
        intersection.z = o.z + dir.z * t;
        return glm::length(intersection - a);
    }

    ~Line(){}
//    bool operator < (const Plane& right) const
//    {
//        return d < right.d;
//    }

public:
    glm::vec3 o;
    glm::vec3 dir;
};


class Mesh {
public:
    Mesh();// : m_cellType(QUAD) {};
    Mesh(const Mesh& r);
    Mesh(const std::vector<Vertex>& V, const std::vector<Cell>& C, ElementType m_cellType);
    Mesh(const std::vector<Vertex>& V, const std::vector<Face>& F, ElementType m_cellType = QUAD);
    Mesh(const Mesh& r, const std::vector<size_t>& cellIds);
    virtual ~Mesh();
    void process();


    std::vector<std::vector<int>> getLocalFaceRegions();
    double getIterativeEnergyOfRegion(std::vector<int> faceIds);
    void printQuality();
    int getNumOfInvertedElements(std::vector<int> localFaceRegion);
private:
    void BuildAllConnectivities();
    void BuildE();
    void BuildF();
    void BuildV_C();
    void BuildE_V();
    void BuildE_E();
    void BuildE_F();
    void BuildF_F();
    void ExtractBoundary();
    void orderNeighboringVertices();
    void classifyCornerVertices();
    void orderVerticesInFaces();

    double getMinimumScaledJacobian(Face& f);
    std::vector<std::vector<int>> mergeRegions(std::vector<std::vector<int>> localFaceRegions);
    std::vector<int> collapseFacesToVerticies(std::vector<int> fIds);
    void expandFaceRegionByOneLayer(std::vector<int>& localFaceRegion);
    void removeUnsolvableFaces(std::vector<int>& localFaceRegion);
    bool hasUnsolvableInversion(Face& f);
    bool hasUnsolvableInversion(std::vector<int> fIds);

    // maybe
    float getAverageArea();
    float getAreaOfFace(Face& f);

public:
    double getMinimumScaledJacobian(std::vector<int> faceIds);
    float getMaximumEdgeLength(Face& f);
    float getMinimumEdgeLength(Vertex& v);
    size_t GetQuality(const char* filename, double& minValue, double& avgValue, const double minSJ = 0.0);
    size_t GetQualityVerdict(double& minValue, double& avgValue, const double minSJ = 0.0);
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

    std::vector<Layer> layers;
    std::vector<Layer> innerLayers;

    double avgEdgeLength;
    size_t numOfSharpEdges;

    glm::vec3 m_center;
    std::vector<size_t> m_refIds;
    double cos_angle_threshold = 0.984807753;
};

#endif /* MESH_H_ */
