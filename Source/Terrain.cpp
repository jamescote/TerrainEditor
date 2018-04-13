#include "Terrain.h"
#include "ObjLoader.h"

#define PI					3.14159265f
#define TWOPI				6.283185307f
#define MRLEVEL_MIN			5

#define HALF 0.5f
#define QUARTER 0.25f
#define THREEQUARTER .75f

// Constructor
Terrain::Terrain(const string& pTerrLoc)
{
	// Set up generation for Terrain
	// Length = U
	// Width = V
	vector<vec2> uvs; 

	m_vEndPos = m_vStartPos = vec3(0.0f);
	m_fWidth = m_fDepth = -1.0f;
	m_iUSize = m_iVSize = 0;
	m_fTileDepth = m_fTileWidth = -1.0f;
	m_iLockedStart = -1;
	m_iSelector = 0;

	objLoader::loadOBJ(pTerrLoc.data(),m_vVertices, uvs, m_vNormals);



	calculateDimensions();
	generateIndices();
	generateNormals();

	glGenVertexArrays( 1, &m_iVertexArray );

	m_iVertexBuffer = ShaderManager::getInstance()->genVertexBuffer( m_iVertexArray,
									    							 0, 3, m_vVertices.data(),
																	 m_vVertices.size() * sizeof( glm::vec3 ), GL_STATIC_DRAW );

	m_iBaryCentric = ShaderManager::getInstance()->genVertexBuffer( m_iVertexArray,
									    							 1, 1, nullptr,
																	 0, GL_STATIC_DRAW );

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

	calculateBarries();
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
	vec3 BLUE = vec3(0.0, 0.0, 1.0);
	vec3 WHITE = vec3(1.0);

	glBindVertexArray( m_iVertexArray );
	/* Not yet Implemented
	if ( nullptr != m_pTexture )
	{
		m_pTexture->bindTexture( ShaderManager::eShaderType::PLANE_SHDR, "gSampler" );
	} //*/
	glBindBuffer(GL_ARRAY_BUFFER, m_iVertexBuffer);
	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, m_iIndicesBuffer );
	glBufferData(GL_ARRAY_BUFFER, m_iUSize*m_iVSize * sizeof(vec3), m_vVertices.data(), GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_vIndices.size() * sizeof( unsigned int ), m_vIndices.data(), GL_DYNAMIC_DRAW);

	//* Draw Points
	ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::WORLD_SHDR, "vColor", &BLACK);
	glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::WORLD_SHDR));
	glPointSize( 5.0f );  // divide by 4th positon in the projected vector
	glDrawArrays( GL_POINTS, 0, m_iUSize*m_iVSize );
//*/


	// Draw Mesh
	glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::TERRAIN_SHDR ) );
	glDrawElements( GL_TRIANGLES, m_vIndices.size(), GL_UNSIGNED_INT, 0 );
	
	ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::WORLD_SHDR, "vColor", &RED);
	glUseProgram(ShaderManager::getInstance()->getProgram(ShaderManager::eShaderType::WORLD_SHDR));
	
	vector< vec3 > vTempVerts;
	if (-1 != m_iLockedStart)
	{
		unsigned int iStartU, iStartV, iEndU, iEndV;
		getSelectedArea(iStartU, iStartV, iEndU, iEndV);;
		for (unsigned int u = iStartU; u <= iEndU; ++u)
			for (unsigned int v = iStartV; v <= iEndV; ++v)
				vTempVerts.push_back(m_vVertices[u + (v*m_iUSize)]);
	}
	else if( !m_vCurrentSubset.empty() )
	{
		glBufferData(GL_ARRAY_BUFFER, m_vCurrentSubset.size() * sizeof(vec3), m_vCurrentSubset.data(), GL_STATIC_DRAW);
		glPointSize(7.0f);  // divide by 4th positon in the projected vector
		glDrawArrays(GL_POINTS, 0, m_vCurrentSubset.size());
		glPointSize(1.0f);
		vTempVerts.push_back(m_vVertices[m_iSelector]);
		ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::WORLD_SHDR, "vColor", &WHITE);
	}
	else
	{
		vTempVerts.push_back(m_vVertices[m_iSelector]);
	}
	
	glBufferData(GL_ARRAY_BUFFER, vTempVerts.size() * sizeof(vec3), vTempVerts.data(), GL_STATIC_DRAW);
	glPointSize(7.0f);  // divide by 4th positon in the projected vector
	glDrawArrays(GL_POINTS, 0, vTempVerts.size());
	glPointSize(1.0f);

	// Draw Saved Selection
	if (!m_vSavedSubset.empty())
	{
		ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::WORLD_SHDR, "vColor", &BLUE);
		glBufferData(GL_ARRAY_BUFFER, m_vSavedSubset.size() * sizeof(vec3), m_vSavedSubset.data(), GL_STATIC_DRAW);
		glPointSize(7.0f);  // divide by 4th positon in the projected vector
		glDrawArrays(GL_POINTS, 0, m_vSavedSubset.size());
		glPointSize(1.0f);
	}


	// glDrawArrays(GL_POINTS, 0, m_vVertices.size());


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
			/*
				Quad:
				1 - 2
				| / |
				3 - 4
			*/
			v4 = normalize(vPosition - m_vVertices[iIndex4]);
			v2 = normalize(vPosition - m_vVertices[index2]);
			v3 = normalize(vPosition - m_vVertices[index3]);

			float fTheta1;
			float fDiagonal = acos(dot(v3, v2));

			// Check Tri v2,v4,v3
			fTheta1 = acos(dot(v2, v4))
					+ acos(dot(v4, v3))
					+ fDiagonal;

			// Intersected through 2nd? Otherwise, we know it intersected with first.
			if (fabs(fTheta1 - (TWOPI)) < 0.1)
				index1 = iIndex4;
			else
			{
				vec3 v1 = normalize(vPosition - m_vVertices[index1]);

				// Check Tri v1,v2,v3
				fTheta1 = acos(dot(v1, v2)) +
						  fDiagonal +
						  acos(dot(v1, v3));

				if (fabs(fTheta1 - TWOPI) >= 0.1f)
					index1 = -1;
			}
	}

}

void Terrain::toggleHeightMap()
{
	if (!m_vVertices.empty())
	{
		if (m_vHeightMap.empty())
			m_vHeightMap.resize(m_iUSize * m_iVSize, 0.0f);

		float fTemp;
		for (unsigned int u = 0; u < m_iUSize; ++u)
		{
			for (unsigned int v = 0; v < m_iVSize; ++v)
			{
				unsigned int iIndex = u + (v*m_iUSize);
				fTemp = m_vVertices[iIndex].y;
				m_vVertices[iIndex].y = m_vHeightMap[iIndex];
				m_vHeightMap[iIndex] = fTemp;

			}
		}
	}

}

void Terrain::grow()
{
	// if (m_vVertices.size == m_iUSize*m_iVSize)
	// {
	// 	//add empty details
	// }
	// vector< vec3 > newMeshVerts.reserve((2*m_iUSize-1) * (2*m_iVSize-1));
	// float C1, C2, C3, C4;

	// C1 = m_vVertices.at(0);
	// C2 = m_vVertices.at(1);
	// C3 = m_vVertices.at(2);
	// C4 = m_vVertices.at(3);

	// newMeshVerts.push_back((-1/2) * C1) + (1 * C2) - ((3/4) * C3) + ((1/4) * C4);
	// newMeshVerts.push_back(((-1/4) * C1) + ((3/4) * C2) - ((3/4) * C3) + ((1/4) * C4));

	// int i = 3;

	// for (; i < m_iUSize- 5; i+=2)
	// {
	// 	C1 = m_vVertices.at(i);
	// 	C2 = m_vVertices.at(i+1);
	// 	C3 = m_vVertices.at(i+2);
	// 	C4 = m_vVertices.at(i+3);

	// 	newMeshVerts.push_back(((-3/4)*C1) + ((3/4)*C2) + ((3/4)C3) - ((1/4)C4));
	// }

	// C1 = m_vVertices.at(i);
	// C2 = m_vVertices.at(i+1);
	// C3 = m_vVertices.at(i+2);
	// C4 = m_vVertices.at(i+3);


	// newMeshVerts.push_back(())

	// newMeshVerts.push_back()
}

void Terrain::reduceTerrain()
{
	reduce(m_vVertices, m_iUSize, m_iVSize);
	calculateDimensions();
	generateIndices();
	calculateBarries();
}

void Terrain::reduce(vector<vec3>& out_Mesh, unsigned int& uSize, unsigned int& vSize)
{
	if (uSize >= MRLEVEL_MIN && vSize >= MRLEVEL_MIN)
	{
		flip(out_Mesh, uSize, vSize);
		reduceU(out_Mesh, uSize, vSize);


		flip(out_Mesh, uSize, vSize);
		reduceU(out_Mesh, uSize, vSize);
	}


}

void Terrain::flip(vector<vec3>& out_meshV, unsigned int& uSize, unsigned int& vSize)
{
	vector<vec3> flippedMesh;
	for (unsigned int u = 0; u < uSize; u++)
	{
		for (unsigned int v = 0; v < vSize; v++)
		{
			flippedMesh.push_back(out_meshV.at(v*uSize + u));
		}
	}

	flippedMesh.insert(flippedMesh.end(), out_meshV.begin()+flippedMesh.size(), out_meshV.end());
	out_meshV = flippedMesh;

	uSize = uSize ^ vSize;
	vSize = vSize ^ uSize;
	uSize = uSize ^ vSize;
}

void Terrain::reduceU(vector<vec3>& out_meshV, unsigned int& uSize, unsigned int& vSize)
{
	vec3 C1, C2, C3, C4;
	vector<vec3> meshV;
	vector<vec3> meshD;
	vector< vec3 > vApplicationCurve;
	bool isOdd = uSize%2 != 0;
	unsigned int tmpUSize;

	int i;

	for (unsigned int v = 0; v < vSize; v++)
	{
		unsigned int vIndex = v*uSize;
		for(int j = 0; j < uSize; ++j)
			vApplicationCurve.push_back(out_meshV[j+vIndex]);

		if (isOdd)
		{	
			vec3 vTranslateVector = vApplicationCurve.back() - *(vApplicationCurve.end() - 2);
			vApplicationCurve.push_back(vApplicationCurve.back() + vec3(vTranslateVector.x, 0.0f, vTranslateVector.z));
		}

		i = 0;

		C1 = vApplicationCurve.at( i ); 
		C2 = vApplicationCurve.at(i+1);
		C3 = vApplicationCurve.at(i+2);
		C4 = vApplicationCurve.at(i+3);

		meshV.push_back(C1);
		meshV.push_back(((-HALF)*C1) + (1*C2) + ((THREEQUARTER)*C3) + ((-QUARTER)*C4));

		meshD.push_back(((-HALF)*C1) + (1*C2) + ((-THREEQUARTER)*C3) + ((-QUARTER)*C4));

		C1 = vApplicationCurve.at( i+2 ); 
		C2 = vApplicationCurve.at(i+3);
		C3 = vApplicationCurve.at(i+4);
		C4 = vApplicationCurve.at(i+5);
		meshD.push_back(((-QUARTER)*C1) + ((THREEQUARTER)*C2) + ((-THREEQUARTER)*C3) + ((QUARTER)*C4));


		for (i=2; i < vApplicationCurve.size() - 5; i+=2)
		{
			C1 = vApplicationCurve.at( i );
			C2 = vApplicationCurve.at(i+1);
			C3 = vApplicationCurve.at(i+2);
			C4 = vApplicationCurve.at(i+3);

			meshV.push_back(((-QUARTER)*C1) + ((THREEQUARTER)*C2) + ((THREEQUARTER)*C3) + ((-QUARTER)*C4));

			if (i >= 4)
			{
				meshD.push_back(((QUARTER)*C1) - ((THREEQUARTER)*C2) + ((THREEQUARTER)*C3) - ((QUARTER)*C4));
			}
		}
	
		C1 = vApplicationCurve.at( i );
		C2 = vApplicationCurve.at(i+1);
		C3 = vApplicationCurve.at(i+2);
		C4 = vApplicationCurve.at(i+3);

		meshV.push_back(((-QUARTER)*C1) + ((THREEQUARTER)*C2) + ((1)*C3) + ((-HALF)*C4));
		meshD.push_back(((QUARTER)*C1) - ((THREEQUARTER)*C2) + ((1)*C3) - ((HALF)*C4));
			
		meshV.push_back(C4);

		vApplicationCurve.clear();

		if( v == 0 )
			tmpUSize = meshV.size();

	}

	meshV.insert( meshV.end(), meshD.begin(), meshD.end());
	meshV.insert( meshV.end(), out_meshV.begin()+(uSize*vSize), out_meshV.end());
	uSize = tmpUSize;
	out_meshV = meshV;
}

void Terrain::get_Quad_Points( float fPosX, float fPosZ, int &iIndex1, int &iIndex2, int &iIndex3, int &iIndex4 )
{
	// Locals/Initialization
	vec3 vOffset;
	iIndex1 = iIndex2 = iIndex3 = iIndex4 = -1;

	// Check that it's within the terrain
	// Get step position
	float fTrueX = (fPosX > m_vStartPos.x + m_fWidth) ? m_vStartPos.x + m_fWidth - (0.5f * m_fTileWidth) : fPosX;
	float fTrueZ = (fPosZ > m_vStartPos.z + m_fDepth) ? m_vStartPos.z + m_fDepth - (0.5f * m_fTileDepth) : fPosZ;
	vOffset = vec3(fTrueX, 0.0, fTrueZ) - vec3(m_vStartPos.x, 0.0, m_vStartPos.z);

	int u = (int)(vOffset.x / (float)m_fTileWidth);
	int v = (int)(vOffset.z / (float)m_fTileDepth);

	if (u < 0)
		u = 0;
	if (v < 0)
		v = 0;
		

	/*
		iV1 --------- iV2
		|				|
		|				|
		|				|
		iV3 --------- iV4
	*/
	// Return Indices
	iIndex1 = (v * m_iUSize) + u;
	iIndex2 = iIndex1 + 1;
	iIndex3 = iIndex1 + m_iUSize;
	iIndex4 = iIndex3 + 1;

}

void Terrain::get_Point_Pos(float fPosX, float fPosZ)
{
	int iIndex1, iIndex2, iIndex3, iIndex4;
	get_Quad_Points(fPosX, fPosZ, iIndex1, iIndex2, iIndex3, iIndex4);

	if (iIndex1 > -1) // Triangle Found
	{
		vec2 vCompVert = vec2(fPosX, fPosZ);
		vec2 vVert1 = vec2(m_vVertices[iIndex1].x, m_vVertices[iIndex1].z);
		vec2 vVert2 = vec2(m_vVertices[iIndex2].x, m_vVertices[iIndex2].z);
		vec2 vVert3 = vec2(m_vVertices[iIndex3].x, m_vVertices[iIndex3].z);
		vec2 vVert4 = vec2(m_vVertices[iIndex4].x, m_vVertices[iIndex4].z);
		float fDist = length(vCompVert - vVert1);
		float fCompDist;
		//m_vSelector = m_vVertices[iIndex1];
		m_iSelector = iIndex1;
		// Vert2 is closer
		if ((fCompDist = length(vCompVert - vVert2)) < fDist)
		{
			fDist = fCompDist;
			//m_vSelector = m_vVertices[iIndex2];
			m_iSelector = iIndex2;
		}

		// Vert3 is closer
		if ((fCompDist = length(vCompVert - vVert3)) < fDist)
		{
			fDist = fCompDist;
			m_iSelector = iIndex3;
		}
		
		// Vert4 is closer
		if ((fCompDist = length(vCompVert - vVert4)) < fDist)
			m_iSelector = iIndex4;
	}
}

void Terrain::lockPoint()
{
	if (-1 == m_iLockedStart)
	{
		m_iLockedStart = m_iSelector;
		m_vCurrentSubset.clear();
	}
	else
	{
		unsigned int iStartU, iStartV, iEndU, iEndV;
		getSelectedArea(iStartU, iStartV, iEndU, iEndV);;
		for (unsigned int u = iStartU; u <= iEndU; ++u)
			for (unsigned int v = iStartV; v <= iEndV; ++v)
				m_vCurrentSubset.push_back(m_vVertices[u + (v*m_iUSize)]);

		m_iLockedStart = -1;
	}
}

void Terrain::saveSelection()
{
	if (!m_vCurrentSubset.empty())
	{
		m_vSavedSubset = m_vCurrentSubset;
		m_vCurrentSubset.clear();
	}
}

void Terrain::getSelectedArea(unsigned int &iStartU, unsigned int &iStartV, unsigned int& iEndU, unsigned int& iEndV)
{
	unsigned int iTemp;

	iStartV = m_iLockedStart / m_iUSize;
	iStartU = m_iLockedStart % m_iUSize;
	iEndV = m_iSelector / m_iUSize;
	iEndU = m_iSelector % m_iUSize;

	if (iEndU < iStartU)
	{
		iTemp = iEndU;
		iEndU = iStartU;
		iStartU = iTemp;
	}
	if (iEndV < iStartV)
	{
		iTemp = iEndV;
		iEndV = iStartV;
		iStartV = iTemp;
	}
}

// General function to compute Terrain Dimensions
void Terrain::calculateDimensions()
{
	if (!m_vVertices.empty())
	{
		m_vStartPos = m_vVertices.front();
		m_vEndPos = (m_iUSize == 0 || m_iVSize == 0) ? m_vVertices.back() : m_vVertices[m_iUSize*m_iVSize - 1];

		m_fWidth = abs(m_vEndPos.x - m_vStartPos.x);
		m_fDepth = abs(m_vEndPos.z - m_vStartPos.z);

		m_fTileWidth = abs(m_vVertices[1].x - m_vStartPos.x);
		
		m_fTileDepth = abs(m_vVertices[m_iUSize].z - m_vStartPos.z);

		if (m_iUSize == 0 || m_iVSize == 0)
		{
			for (vector<vec3>::const_iterator iter = m_vVertices.begin();
				iter->z == m_vStartPos.z;
				++iter)
				m_iUSize++;
			for (vector<vec3>::const_iterator iter = m_vVertices.begin();
				iter != m_vVertices.end() && iter->x == m_vStartPos.x;
				iter += m_iUSize)
				m_iVSize++;
		}
	}
	else
		cout << "Unable to compute dimensions; no vertices loaded.\n";
}

void Terrain::generateIndices( )
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

void Terrain::calculateBarries()
{
	if (!m_vIndices.empty())
	{
		unsigned int iArry[3] = { 0, 1, 2 };
		unsigned int iX = 0;

		vector<unsigned int> BCverts;
		BCverts.reserve(m_vIndices.size());

		for (unsigned int v = 0; v < m_iVSize; ++v)
		{
			for (unsigned int u = 0; u < m_iUSize; u += 3)
			{
				BCverts.push_back(iArry[iX]);
				if (u + 1 < m_iUSize)
					BCverts.push_back(iArry[(iX + 1) % 3]);
				if (u + 2 < m_iUSize)
					BCverts.push_back(iArry[(iX + 2) % 3]);
			}
			iX = (iX + 2) % 3;
		}

		glBindVertexArray(m_iVertexArray);
		glBindBuffer(GL_ARRAY_BUFFER, m_iBaryCentric);
		glBufferData(GL_ARRAY_BUFFER, BCverts.size() * sizeof(unsigned int), BCverts.data(), GL_STATIC_DRAW);
		glBindVertexArray(0);
	}
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
