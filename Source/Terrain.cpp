#include "Terrain.h"
#include "ObjLoader.h"

#define PI					3.14159265f

// Constructor
Terrain::Terrain(const string& pTerrLoc)
{
	// Set up generation for Terrain
	// Length = U
	// Width = V
	vector<int> BCverts;
	vector<vec2> uvs; 

	m_vEndPos = m_vStartPos = vec3(0.0f);
	m_fWidth = m_fDepth = -1.0f;
	m_iUSize = m_iVSize = 0;
	m_fTileDepth = m_fTileWidth = -1.0f;

	objLoader::loadOBJ(pTerrLoc.data(),m_vVertices, m_vIndices, uvs, m_vNormals);


	calculateDimensions();
	generateIndices();
	generateNormals();


	


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

	/* Draw Points
	ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::WORLD_SHDR, "vColor", &BLACK);
	glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::TERRAIN_GRID_SHDR));
	glPointSize( 5.0f );  // divide by 4th positon in the projected vector
	glDrawArrays( GL_POINTS, 0, m_vVertices.size() );
//*/

	// Draw Mesh
	glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::TERRAIN_SHDR ) );
	glDrawElements( GL_TRIANGLES, m_vIndices.size(), GL_UNSIGNED_INT, 0 );

	glBufferData(GL_ARRAY_BUFFER, 4 * sizeof(vec3), &m_vTempSelectedQuad, GL_STATIC_DRAW);
	
	if( vec3(-1.0) != m_vSelectedPoint )
	{
		ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::WORLD_SHDR, "vColor", &RED);
		glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::WORLD_SHDR));
		glBufferData(GL_ARRAY_BUFFER, sizeof(vec3), &m_vSelectedPoint, GL_STATIC_DRAW);
		glPointSize( 5.0f );  // divide by 4th positon in the projected vector
		glDrawArrays( GL_POINTS, 0, 1 );
		glPointSize( 1.0f );
	}

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
	if (fPosX >= m_vStartPos.x && fPosX < (m_vEndPos.x - m_fTileWidth) &&
		fPosZ >= m_vStartPos.z && fPosZ < (m_vEndPos.z - m_fTileDepth))
	{
		// Get step position
		vOffset = vec3(fPosX, 0.0, fPosZ) - vec3(m_vStartPos.x, 0.0, m_vStartPos.z);
		int u = (int)(vOffset.x / (float)m_fTileWidth);
		int v = (int)(vOffset.z / (float)m_fTileDepth);

		/*
			iV1 --------- iV2
			|							|
			|							|
			|							|
			iV3 --------- iV4
		*/
		// Return Indices
		iIndex1 = (u * m_iUSize) + v;
		iIndex2 = iIndex1 + 1;
		iIndex3 = iIndex1 + m_iUSize;
		iIndex4 = iIndex3 + 1;
	}
}

void Terrain::get_Point_Pos(float fPosX, float fPosZ)
{
	int iIndex1, iIndex2, iIndex3;
	get_Triangle_Points(fPosX, fPosZ, iIndex1, iIndex2, iIndex3);
	m_vSelectedPoint = vec3(-1.0f);

	if (iIndex1 > -1) // Triangle Found
	{
		vec2 vCompVert = vec2(fPosX, fPosZ);
		vec2 vVert1 = vec2(m_vVertices[iIndex1].x, m_vVertices[iIndex1].z);
		vec2 vVert2 = vec2(m_vVertices[iIndex2].x, m_vVertices[iIndex2].z);
		vec2 vVert3 = vec2(m_vVertices[iIndex3].x, m_vVertices[iIndex3].z);
		float fDist = length(vCompVert - vVert1);
		float fCompDist;
		vec3 m_vSelectedPoint = m_vVertices[iIndex1];

		// Vert2 is closer
		if ((fCompDist = length(vCompVert - vVert2)) < fDist)
		{
			fDist = fCompDist;
			m_vSelectedPoint = m_vVertices[iIndex2];
		}

		// Vert3 is closer
		if ((fCompDist = length(vCompVert - vVert3)) < fDist)
			m_vSelectedPoint = m_vVertices[iIndex3];
	}
}

// General function to compute Terrain Dimensions
void Terrain::calculateDimensions()
{
	if (!m_vVertices.empty())
	{
		m_vStartPos = m_vVertices.front();
		m_vEndPos = m_vVertices.back();

		m_fWidth = abs(m_vEndPos.x - m_vStartPos.x);
		m_fDepth = abs(m_vEndPos.z - m_vStartPos.z);

		m_fTileWidth = abs(m_vVertices[1].x - m_vStartPos.x);
		m_iUSize = (unsigned int)round(m_fWidth / m_fTileWidth) + 1;
		m_fTileDepth = abs(m_vVertices[m_iUSize].z - m_vStartPos.z);
		m_iVSize = (unsigned int)round(m_fDepth / m_fTileDepth) + 1;
	}
	else
		cout << "Unable to compute dimensions; no vertices loaded.\n";
}

void Terrain::generateIndices()
{
	if (!m_vVertices.empty())
	{
		if (!m_vIndices.empty())
			m_vIndices.clear();

		for (unsigned int v = 0; v < m_iVSize - 1; ++v)
		{
			for (unsigned int u = 0; u < m_iUSize - 1; ++u)
			{
				unsigned int iA = u + (v*m_iUSize);
				unsigned int iB = iA + 1;
				unsigned int iC = iA + m_iUSize;
				unsigned int iD = iC + 1;

				m_vIndices.push_back(iA);
				m_vIndices.push_back(iB);
				m_vIndices.push_back(iC);

				m_vIndices.push_back(iC);
				m_vIndices.push_back(iD);
				m_vIndices.push_back(iB);
			}
		}
	}
	else
		cout << "Unable to generate Indices; no vertices loaded.\n";
}

// Computes and Generates Normals per vertex.
void Terrain::generateNormals()
{
	// Only attempt to compute if Vertices are loaded.
	if (!m_vVertices.empty())
	{
		// Compute Indices if not calculated yet.
		if (m_vIndices.empty())
			generateIndices();

		// Vertices and Indices Loaded, need to compute Normals
		m_vNormals.resize(m_vVertices.size(), vec3(0.0));
		for (unsigned int i = 0; i < m_vIndices.size(); i += 3)
		{
			assert((i + 2) < m_vIndices.size());
			vec3 vTri[3] = { m_vVertices[m_vIndices[i]], m_vVertices[m_vIndices[i + 1]], m_vVertices[m_vIndices[i + 2]] };

			// Calculate Triangle Normal: (v1 - v0) x (v2 - v0)
			vec3 vTriNormal = cross((vTri[1] - vTri[0]), (vTri[2] - vTri[0]));

			// Accumulate Normals Per Vertex;
			m_vNormals[m_vIndices[i]] += vTriNormal;
			m_vNormals[m_vIndices[i + 1]] += vTriNormal;
			m_vNormals[m_vIndices[i + 2]] += vTriNormal;
		}

		// Normalize all Accumulated Normals
		for (vector< vec3 >::iterator vNormIter = m_vNormals.begin();
			vNormIter != m_vNormals.end();
			++vNormIter)
			(*vNormIter) = normalize((*vNormIter));
	}
	else
		cout << "Unable to compute Normals; Vertices aren't loaded.\n";
}
