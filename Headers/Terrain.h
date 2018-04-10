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

private:
	Terrain( const Terrain* pNewPlane );  // Protected Copy Constructor

	// Normal of the Terrain.
	vector< vec3 > m_vVertices, m_vNormals;
	vector< vec2 > m_vUVs;
	vector< unsigned int > m_vIndices;
	vec3 m_vTempSelectedQuad[4];
	vec3 m_vStartPos, m_vEndPos;
	float m_fWidth, m_fDepth, m_fTileWidth, m_fTileDepth;
	unsigned int m_iUSize, m_iVSize;

	GLuint m_iVertexArray, m_iBaryCentric, m_iVertexBuffer, m_iNormalBuffer, m_iTextureBuffer, m_iIndicesBuffer;

	void get_Quad_Points( float fPosX, float fPosZ, int &iIndex1, int &iIndex2, int &iIndex3, int &iIndex4 );
};
