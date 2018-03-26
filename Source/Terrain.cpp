#include "Terrain.h"

#define NUM_CORNERS 4
#define DEFAULT_U 10
#define DEFAULT_V 10
#define DEFAULT_DEPTH 50.0f
#define DEFAULT_WIDTH 50.0f
#define TILE_WIDTH (DEFAULT_WIDTH / (float)DEFAULT_V)
#define TILE_DEPTH (DEFAULT_DEPTH / (float)DEFAULT_U)

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
}

// Destructor
Terrain::~Terrain()
{
	glDeleteBuffers( 1, &m_iNormalBuffer );
	glDeleteBuffers( 1, &m_iVertexBuffer );
	glDeleteBuffers( 1, &m_iTextureBuffer );
	glDeleteVertexArrays( 1, &m_iVertexArray );
}

// Setup OpenGl to draw the Terrain using the Plane Shader.
void Terrain::draw(  )
{
	glBindVertexArray( m_iVertexArray );
	/* Not yet Implemented
	if ( nullptr != m_pTexture )
	{
		m_pTexture->bindTexture( ShaderManager::eShaderType::PLANE_SHDR, "gSampler" );
	} //*/
	glBindBuffer(GL_ARRAY_BUFFER, m_iVertexBuffer);
	glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::PLANE_SHDR ) );
	glBufferData(GL_ARRAY_BUFFER, m_vVertices.size() * sizeof(vec3), m_vVertices.data(), GL_STATIC_DRAW);
	glDrawArrays(GL_LINE_STRIP, 0, m_vVertices.size());
	glPointSize(4.0f);
	glDrawArrays(GL_POINTS, 0, m_vVertices.size());
	glPointSize(1.0f);

	vec3 RED = vec3(1.0, 0.0, 0.0);

	glUseProgram(ShaderManager::getInstance()->getProgram(ShaderManager::eShaderType::WORLD_SHDR));
	ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::WORLD_SHDR, "vColor", &RED);
	glBufferData(GL_ARRAY_BUFFER, 4 * sizeof(vec3), &m_vTempSelectedQuad, GL_STATIC_DRAW);

	glDrawArrays(GL_LINE_STRIP, 0, 4);
	glPointSize(6.0f);
	glDrawArrays(GL_POINTS, 0, 4);
	glPointSize(1.0f);


	/* Not Yet Implemented
	if ( nullptr != m_pTexture )
		m_pTexture->unbindTexture();
	//*/

	glUseProgram(0);
	glBindVertexArray( 0 );
}

void Terrain::get_Triangle_Points(float fPosX, float fPosZ, int &index1, int &index2, int &index3)
{
	vec3 vOffset;
	fPosZ = -fPosZ;
	index1 = index2 = index3 = -1;
	int iV1, iV2, iV3, iV4;

	// Check that it's within the terrain
	if (fPosX >= m_vStartPos.x && fPosX < m_vEndPos.x &&
		fPosZ >= m_vStartPos.z && fPosZ < m_vEndPos.z)
	{
		// Get step position
		vOffset = vec3(fPosX, 0.0, fPosZ) - vec3(m_vStartPos.x, 0.0, m_vStartPos.z);
		cout << "TILE_WIDTH: " << TILE_WIDTH << endl;
		cout << "TILE_DEPTH: " << TILE_DEPTH << endl;
		cout << "float division: " << vOffset.x / (float)TILE_WIDTH << "|" << vOffset.z / (float)TILE_DEPTH << "\n";

		int u = vOffset.x / (float)TILE_WIDTH;
		int v = vOffset.z / (float)TILE_DEPTH;
		float fXWithinQuad = fmodf(vOffset.x, TILE_WIDTH);
		float fZWithinQuad = fmodf(vOffset.z, TILE_DEPTH);

		/*
			iV1 --------- iV2
			|				|
			|				|
			|				|
			iV3 --------- iV4
		*/
		iV1 = (v * DEFAULT_V) + u;
		iV2 = iV1 + 1;
		iV3 = iV1 + DEFAULT_V;
		iV4 = iV3 + 1;

		m_vTempSelectedQuad[0] = m_vVertices[iV1]; 
		m_vTempSelectedQuad[1] = m_vVertices[iV2];
		m_vTempSelectedQuad[2] = m_vVertices[iV3];
		m_vTempSelectedQuad[3] = m_vVertices[iV4];
	}
}
