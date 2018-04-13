#pragma once

#include "Object3D.h"
#include <stack>

class Terrain
{
public:
	Terrain(const string& pTerrLoc);
	~Terrain();

	// Overridden intersect function
	void draw( );

	void get_Point_Pos( float fPosX, float fPosZ );
	void toggleHeightMap();
	void lockPoint();
	void saveSelection();
	void clearSelection();
	void applyTerrain(const Terrain* vMesh);
	void growTerrain();
	void reduceTerrain();

private:
	Terrain( const Terrain* pNewPlane );  // Protected Copy Constructor

	// Mesh Structure for Storing unique attributes to a Terrain Mesh
	struct tMesh
	{
		// Vertex Data
		vector< vec3 > m_vVertices, m_vNormals;
		vector< unsigned int > m_vIndices, m_vBarries;
		vec3 m_vStartPos, m_vEndPos;
		vector< vec2 > m_vUVs;

		// Dimensions of Mesh
		unsigned int m_iUSize, m_iVSize;
		float m_fWidth, m_fDepth, m_fTileWidth, m_fTileDepth;

		// Stack for Reverse Subdivision
		stack<pair<bool, unsigned int>> m_MrMap;

		tMesh()
		{
			clear();
		}

		void clear()
		{
			m_vVertices.clear();
			m_vIndices.clear();
			m_vBarries.clear();
			m_vUVs.clear();
			m_vNormals.clear();
			m_vEndPos = m_vStartPos = vec3(0.0f);
			m_fWidth = m_fDepth = -1.0f;
			m_iUSize = m_iVSize = 0;
			m_fTileDepth = m_fTileWidth = -1.0f;
		}

		void outputMesh()
		{
			cout << "Mesh Details:\n\t";
			cout << "Num Verts: " << m_vVertices.size() << "; Indices: " << m_vIndices.size()
				<< "; Barries: " << m_vBarries.size() << "; UVs: " << m_vUVs.size() << "; Normals: " << m_vNormals.size() << endl;
			cout << "\tStartPos: {" << m_vStartPos.x << ", " << m_vStartPos.y << ", " << m_vStartPos.z << "}; EndPos: {"
				<< m_vEndPos.x << ", " << m_vEndPos.y << ", " << m_vEndPos.z << "}\n\t";
			cout << "Width: " << m_fWidth << "; Depth: " << m_fDepth << endl;
			cout << "\tUSize: " << m_iUSize << "; VSize: " << m_iVSize << endl;
			cout << "\tTile Width: " << m_fTileWidth << "; Tile Depth: " << m_fTileDepth << endl;

		}
	};

	// Main TerrainMesh
	tMesh m_defaultTerrain;
	tMesh m_vSavedSubset, m_vCurrentSubset, m_vApplicationMesh;
	
	//vector< vec3 > m_vSavedSubset, m_vCurrentSubset;

	vec3 m_vSelector, m_vLockedStart;
	unsigned int m_iSelector;
	int m_iLockedStart;

	// Height Map for Main Terrain, helps for selections
	vector< float > m_vHeightMap;
	bool m_bHeightMapStored;
	bool m_bApplicationMode;

	// Calculation functions
	void initMesh(tMesh& pTerrain);
	void calculateDimensions( tMesh& pTerrain );
	void generateIndices( tMesh& pTerrain );
	void generateNormals(tMesh& pTerrain);
	void calculateBarries(tMesh& pTerrain);
	void getSelectedArea(unsigned int &iStartU, unsigned int &iStartV, unsigned int& iEndU, unsigned int& iEndV);
	void loadMeshData(const tMesh& pTerrain);
	void blendMesh();
	int blendCoarsePoints(unsigned int iStartU, unsigned int iStartV, unsigned int iEndU, unsigned int iEndV);
	int replaceDetails(unsigned int iStartU, unsigned int iStartV, unsigned int iEndU, unsigned int iEndV, unsigned int iTerrRowSize, unsigned int iAppRowSize );
	int replaceSquaredDetails(unsigned int iStartU, unsigned int iStartV, unsigned int iEndU, unsigned int iEndV, unsigned int& iTerrIndex,
								unsigned int& iTerrUSize, unsigned int& iTerrVSize, unsigned int& iAppUSize, unsigned int& iAppVSize);

	GLuint m_iVertexArray, m_iBaryCentric, m_iVertexBuffer, m_iNormalBuffer, m_iTextureBuffer, m_iIndicesBuffer;

	void get_Quad_Points( float fPosX, float fPosZ, int &iIndex1, int &iIndex2, int &iIndex3, int &iIndex4 );
	void get_Triangle_Points(float fPosX, float fPosZ, int &index1, int &index2, int &index3);
	void flip(tMesh&);
	void reduce(tMesh&);
	void reduceU(tMesh&);
	void grow(tMesh&);	
	void growU(tMesh&);


};
