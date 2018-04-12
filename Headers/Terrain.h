#pragma once

#include "Object3D.h"

class Terrain
{
public:
	Terrain(const string& pTerrLoc);
	~Terrain();

	// Overridden intersect function
	void draw( );

	void get_Triangle_Points(float fPosX, float fPosZ, int &index1, int &index2, int &index3);
	void get_Point_Pos( float fPosX, float fPosZ );
	void toggleHeightMap();
	void lockPoint();
	void saveSelection();
	void grow();	
	void reduce(vector<vec3>& out_Mesh, unsigned int& uSize, unsigned int& vSize);
	void reduceTerrain();
private:
	Terrain( const Terrain* pNewPlane );  // Protected Copy Constructor

	// Normal of the Terrain.
	vector< vec3 > m_vVertices, m_vNormals;
	vector< vec3 > m_vSavedSubset, m_vCurrentSubset;
	vector< float > m_vHeightMap;
	vector< vec2 > m_vUVs;
	vector< unsigned int > m_vIndices;
	vec3 m_vTempSelectedQuad[4];
	vec3 m_vStartPos, m_vEndPos;
	vec3 m_vSelector, m_vLockedStart;
	unsigned int m_iSelector;
	int m_iLockedStart;
	float m_fWidth, m_fDepth, m_fTileWidth, m_fTileDepth;
	unsigned int m_iUSize, m_iVSize;
	
	// Calculation functions
	void calculateDimensions();
	void generateIndices();
	void generateNormals();
	void getSelectedArea(unsigned int &iStartU, unsigned int &iStartV, unsigned int& iEndU, unsigned int& iEndV);

	GLuint m_iVertexArray, m_iBaryCentric, m_iVertexBuffer, m_iNormalBuffer, m_iTextureBuffer, m_iIndicesBuffer;

	void get_Quad_Points( float fPosX, float fPosZ, int &iIndex1, int &iIndex2, int &iIndex3, int &iIndex4 );
	void flip(vector<vec3>&, unsigned int&, unsigned int&);
	void reduceU(vector<vec3>& meshV, unsigned int& uSize, unsigned int& vSize);
};
