#pragma once
#include "Object3D.h"
class Terrain 
{
public:
	Terrain();
	~Terrain();

	// Overridden intersect function
	void draw( );

	void get_Triangle_Points(float fPosX, float fPosZ, int &index1, int &index2, int &index3);

private:
	Terrain( const Terrain* pNewPlane );  // Protected Copy Constructor

	// Normal of the Terrain.
	vector< vec3 > m_vVertices, m_vNormals;
	vector< vec2 > m_vUVs;
	vector< int > m_vIndices;
	vec3 m_vTempSelectedQuad[4];
	vec3 m_vStartPos, m_vEndPos;


	GLuint m_iVertexArray, m_iVertexBuffer, m_iNormalBuffer, m_iTextureBuffer, m_iIndexBuffer;


};

