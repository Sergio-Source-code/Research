#ifndef Mesh_File_Writer_H_
#define Mesh_File_Writer_H_

#include <vector>
#include <string>
#include <fstream>

#include "Mesh.h"

class MeshFileWriter
{
public:
	MeshFileWriter(const Mesh& mesh, const char* pFileName = "tempfile");
	MeshFileWriter(const std::vector<Vertex>& v, const std::vector<Cell>& c, const char* pFileName, const ElementType cellType = HEXAHEDRA);
	MeshFileWriter(const std::vector<Vertex>& v, const std::vector<Face>& f, const char* pFileName, const ElementType cellType = QUAD);
	MeshFileWriter();
	~MeshFileWriter();

	void WriteMeshBoundaryInfo(const char* pFileName = NULL);
	void WriteCellData(const std::vector<int>& cellData, const char* dataName = "cell_scalar");
	void WritePointData(const std::vector<int>& pointData, const char* dataName = "point_scalar");
	void WriteEdgesVtk();
	void WriteSharpEdgesVtk();
	void WriteCornersVtk();
	void WriteEdgesVtk(const std::vector<size_t>& edgeIds);
	void WriteFacesVtk();
	void WriteVerticesVtk();
	void WriteSurfaceOff();

public:
	void WriteFile();
	void WriteMeshFile();
	void WriteVtkFile();
	void WriteVtkPolyDataFile();
	void WriteOffFile();
	void WriteObjFile();
	void WriteStlFile();
	void WritePlyFile();
	void WriteSingularitiesVtkFile();
    void SetFixFlag(bool bFixed = true);
	void FixMesh();

private:
	std::string m_strFileName;
	Mesh m_mesh;
	bool m_bFixed;
};

#endif // Mesh_File_Writer_H_
