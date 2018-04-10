#include "Terrain.h"
#include "ObjLoader.h"

#define NUM_CORNERS 4
#define DEFAULT_U 10
#define DEFAULT_V 10
#define DEFAULT_DEPTH 50.0f
#define DEFAULT_WIDTH 50.0f
#define TILE_WIDTH (DEFAULT_WIDTH / (float)DEFAULT_V)
#define TILE_DEPTH (DEFAULT_DEPTH / (float)DEFAULT_U)
#define PI					3.14159265f

// Constructor
Terrain::Terrain()
{
	// Set up generation for Terrain
	// Length = U
	// Width = V
	vector<int> BCverts;
	vector<vec2> uvs; 

	objLoader::loadOBJ("models/terrain2.obj",m_vVertices, m_vIndices, uvs, m_vNormals);

	m_vStartPos = m_vVertices.front();
	m_vEndPos = m_vVertices.back();

	m_fWidth = abs( m_vEndPos.x - m_vStartPos.x);
	m_fDepth = abs( m_vEndPos.z - m_vStartPos.z);

	m_fTileWidth = abs( m_vVertices[1].x - m_vStartPos.x);
	m_iUSize = (unsigned int)ceil(m_fWidth / m_fTileWidth) + 1;
	m_fTileDepth = abs(m_vVertices[m_iUSize].z - m_vStartPos.z);
	m_iVSize = (unsigned int)ceil(m_fDepth / m_fTileDepth);


	unsigned int iArry[3] = {0, 1, 2};
	unsigned int iX = 0;

	BCverts.reserve( m_vIndices.size());

	for( unsigned int v = 0; v < m_iVSize; ++v )
	{
		for( unsigned int u = 0; u < m_iUSize; u += 3)
		{
			BCverts.push_back(iArry[iX]);
			if( u + 1 < m_iUSize )
				BCverts.push_back(iArry[(iX + 1) % 3]);
			if( u + 2 < m_iUSize )
				BCverts.push_back(iArry[(iX + 2) % 3]);
		}
		iX = (iX + 2) % 3;
	}

	glGenVertexArrays( 1, &m_iVertexArray );

	m_iVertexBuffer = ShaderManager::getInstance()->genVertexBuffer( m_iVertexArray,
									    							 0, 3, m_vVertices.data(),
																	 m_vVertices.size() * sizeof( glm::vec3 ), GL_STATIC_DRAW );

	m_iBaryCentric = ShaderManager::getInstance()->genVertexBuffer( m_iVertexArray,
									    							 1, 1, BCverts.data(),
																	 BCverts.size() * sizeof( int ), GL_STATIC_DRAW );

	m_iNormalBuffer = ShaderManager::getInstance()->genVertexBuffer( m_iVertexArray,
																	 2, 3, m_vNormals.data(),
																	 m_vNormals.size() * sizeof( glm::vec3 ), GL_STATIC_DRAW );

	m_iTextureBuffer = ShaderManager::getInstance()->genVertexBuffer(m_iVertexArray,
																	 3, 2, m_vUVs.data(),
																	 m_vUVs.size() * sizeof(vec2), GL_STATIC_DRAW);

  	m_iIndicesBuffer = ShaderManager::getInstance()->genIndicesBuffer( m_iVertexArray,
																	 m_vIndices.data(),
																	 m_vIndices.size() * sizeof( unsigned int ),
																	 GL_STATIC_DRAW );

}

// Destructor
Terrain::~Terrain()
{
	glDeleteBuffers( 1, &m_iNormalBuffer );
	glDeleteBuffers( 1, &m_iVertexBuffer );
	glDeleteBuffers( 1, &m_iBaryCentric );
	glDeleteBuffers( 1, &m_iTextureBuffer );
	glDeleteBuffers( 1, &m_iIndicesBuffer );
	glDeleteVertexArrays( 1, &m_iVertexArray );
}

// Setup OpenGl to draw the Terrain using the Plane Shader.
void Terrain::draw(  )
{


	// Colors for drawing
	vec3 RED = vec3(1.0, 0.0, 0.0);
	vec3 BLACK = vec3(0.0);

	glBindVertexArray( m_iVertexArray );
	/* Not yet Implemented
	if ( nullptr != m_pTexture )
	{
		m_pTexture->bindTexture( ShaderManager::eShaderType::PLANE_SHDR, "gSampler" );
	} //*/
	glBindBuffer(GL_ARRAY_BUFFER, m_iVertexBuffer);
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, m_iIndicesBuffer );
	glBufferData(GL_ARRAY_BUFFER, m_vVertices.size() * sizeof(vec3), m_vVertices.data(), GL_STATIC_DRAW);

	// Draw Points
	ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::WORLD_SHDR, "vColor", &BLACK);
	glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::TERRAIN_GRID_SHDR));
	glPointSize( 5.0f );
	glDrawArrays( GL_POINTS, 0, m_vVertices.size() );
	glPointSize(1.0f);
	//glDrawElements( GL_LINE_STRIP, m_vIndices.size(), GL_UNSIGNED_INT, 0 );

	// Draw Mesh
	glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::TERRAIN_SHDR ) );
	glDrawElements( GL_TRIANGLES, m_vIndices.size(), GL_UNSIGNED_INT, 0 );

	glBufferData(GL_ARRAY_BUFFER, 4 * sizeof(vec3), &m_vTempSelectedQuad, GL_STATIC_DRAW);

	//glDrawArrays(GL_TRIANGLES, 0, 4);


	/* Not Yet Implemented
	if ( nullptr != m_pTexture )
		m_pTexture->unbindTexture();
	//*/

	glUseProgram(0);
	glBindVertexArray( 0 );
}

void Terrain::get_Triangle_Points(float fPosX, float fPosZ, int &index1, int &index2, int &index3)
{
	int iIndex4;

	// Get Indices that position lies within.
	get_Quad_Points( fPosX, fPosZ, index1, index2, index3, iIndex4 );

	// Position is within Terrain?
	if( index1 > -1 )
	{
			vec3 v2, v3, v4;
			vec3 vPosition = vec3( fPosX, 0.0f, fPosZ );

			// Compute Barycentric Coordinates
			v4 = normalize(vPosition - m_vVertices[iIndex4]);
			v2 = normalize(vPosition - m_vVertices[index2]);
			v3 = normalize(vPosition - m_vVertices[index3]);

			float fTheta1;

			// Check Tri v2,v4,v3
			fTheta1 = acos(dot(v2, v4))
					+ acos(dot(v4, v3))
					+ acos(dot(v3, v2));

			// Intersected through 2nd? Otherwise, we know it intersected with first.
			if (fabs(fTheta1 - (2 * PI) < 0.1))
				index1 = iIndex4;
	}

}

void Terrain::get_Quad_Points( float fPosX, float fPosZ, int &iIndex1, int &iIndex2, int &iIndex3, int &iIndex4 )
{
	// Locals/Initialization
	vec3 vOffset;
	iIndex1 = iIndex2 = iIndex3 = iIndex4 = -1;

	// Check that it's within the terrain
	if (fPosX >= m_vStartPos.x && fPosX < (m_vEndPos.x - TILE_WIDTH) &&
		fPosZ >= m_vStartPos.z && fPosZ < (m_vEndPos.z - TILE_DEPTH))
	{
		// Get step position
		vOffset = vec3(fPosX, 0.0, fPosZ) - vec3(m_vStartPos.x, 0.0, m_vStartPos.z);
		int u = (int)(vOffset.x / (float)TILE_WIDTH);
		int v = (int)(vOffset.z / (float)TILE_DEPTH);

		/*
			iV1 --------- iV2
			|							|
			|							|
			|							|
			iV3 --------- iV4
		*/
		// Return Indices
		iIndex1 = (u * DEFAULT_U) + v;
		iIndex2 = iIndex1 + 1;
		iIndex3 = iIndex1 + DEFAULT_U;
		iIndex4 = iIndex3 + 1;
	}
}
