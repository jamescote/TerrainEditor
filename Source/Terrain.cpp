#include "Terrain.h"

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
	float fUV_Step_U = (float)(1 / DEFAULT_U);
	float fUV_Step_V = (float)(1 / DEFAULT_V);
	vec3 vStartingPos = vec3(-DEFAULT_WIDTH / 2.0f, 0.0, -DEFAULT_DEPTH / 2.0f);
	m_vStartPos = vStartingPos;
	m_vEndPos = vStartingPos + vec3(DEFAULT_WIDTH, 0.0, DEFAULT_DEPTH);
	vec2 vStartingUV = vec2(0.0);
	vector< vec3 > vWorkingVerts, vWorkingNormals;
	vector< vec2 > vWorkingUVs;
	for (unsigned int i = 0; i < 4; ++i)
		m_vTempSelectedQuad[i] = vec3(0.0);

	// Generate Vertices for the Terrain
	for (unsigned int u = 0; u < DEFAULT_U; ++u)
	{
		for (unsigned int v = 0; v < DEFAULT_V; ++v)
		{
			// Add Normal and Vertex
			m_vVertices.push_back( vStartingPos );
			m_vNormals.push_back(vec3(0.0, 1.0, 0.0));
			m_vUVs.push_back(vStartingUV);
			vStartingUV.y += fUV_Step_V;
			vStartingPos.z += TILE_DEPTH;	// Increment to next vPosition
		}
		// Store the vector and reset for next u-strip
		vStartingPos.z -= DEFAULT_DEPTH;
		vStartingPos.x += TILE_WIDTH;
		vStartingUV.y = 0.0f;
		vStartingUV.x += fUV_Step_U;
	}

	// Generate Indices for drawing
	vStartingPos.x -= DEFAULT_WIDTH;
	vStartingPos += vec3( TILE_WIDTH / 2.0f, 0.0, TILE_DEPTH / 2.0f );
	int iIndex0, iIndex1, iIndex2, iIndex3;
	for (unsigned int u = 0; u < DEFAULT_U; ++u)
	{
		for (unsigned int v = 0; v < DEFAULT_V; ++v)
		{
			get_Quad_Points( vStartingPos.x, vStartingPos.z, iIndex0, iIndex1, iIndex2, iIndex3 );
			if( iIndex0 > -1 )
			{
				// Push Indices for first Triangle
				m_vIndices.push_back(iIndex0);
				m_vIndices.push_back(iIndex1);
				m_vIndices.push_back(iIndex2);
				// Push Indices for Second Triangle
				m_vIndices.push_back(iIndex1);
				m_vIndices.push_back(iIndex2);
				m_vIndices.push_back(iIndex3);
			}
			vStartingPos.z += TILE_DEPTH;
		}
		vStartingPos.z -= DEFAULT_DEPTH;
		vStartingPos.x += TILE_WIDTH;
	}

	glGenVertexArrays( 1, &m_iVertexArray );

	m_iVertexBuffer = ShaderManager::getInstance()->genVertexBuffer( m_iVertexArray,
									    							 0, 3, m_vVertices.data(),
																	 m_vVertices.size() * sizeof( glm::vec3 ), GL_STATIC_DRAW );
	m_iNormalBuffer = ShaderManager::getInstance()->genVertexBuffer( m_iVertexArray,
																	 1, 3, m_vNormals.data(),
																	 m_vNormals.size() * sizeof( glm::vec3 ), GL_STATIC_DRAW );
	m_iTextureBuffer = ShaderManager::getInstance()->genVertexBuffer(m_iVertexArray,
																	 2, 2, m_vUVs.data(),
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
	glDrawElements( GL_LINE_STRIP, m_vIndices.size(), GL_UNSIGNED_INT, 0 );

	// Draw Mesh
	glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::TERRAIN_SHDR ) );
	glDrawElements( GL_TRIANGLES, m_vIndices.size(), GL_UNSIGNED_INT, 0 );


	glBufferData(GL_ARRAY_BUFFER, 4 * sizeof(vec3), &m_vTempSelectedQuad, GL_STATIC_DRAW);

	glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);


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
