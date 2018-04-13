#include "Terrain.h"
#include "ObjLoader.h"

#define PI					3.14159265f
#define TWOPI				6.283185307f
#define MRLIMIT_MIN			5

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

	vector<vec3> mVerts, mNorms;

	m_defaultTerrain.m_vEndPos = m_defaultTerrain.m_vStartPos = vec3(0.0f);
	m_defaultTerrain.m_fWidth = m_defaultTerrain.m_fDepth = -1.0f;
	m_defaultTerrain.m_iUSize = m_defaultTerrain.m_iVSize = 0;
	m_defaultTerrain.m_fTileDepth = m_defaultTerrain.m_fTileWidth = -1.0f;

	m_iLockedStart = -1;
	m_iSelector = 0;

	objLoader::loadOBJ(pTerrLoc.data(),mVerts, uvs, mNorms);

	m_defaultTerrain.m_vVertices = mVerts;
	m_defaultTerrain.m_vNormals = mNorms;
	m_defaultTerrain.m_vUVs = uvs;
	//m_defaultTerrain.m_vNormals;

	calculateDimensions(); //TODO: Update for tMesh
	generateIndices(); //TODO: Update for tMesh
	generateNormals(); //TODO: Update for tMesh

	glGenVertexArrays( 1, &m_iVertexArray );

	m_iVertexBuffer = ShaderManager::getInstance()->genVertexBuffer( m_iVertexArray,
									    							 0, 3, m_defaultTerrain.m_vVertices.data(),
																	 m_defaultTerrain.m_vVertices.size() * sizeof( glm::vec3 ), GL_STATIC_DRAW );

	m_iBaryCentric = ShaderManager::getInstance()->genVertexBuffer( m_iVertexArray,
									    							 1, 1, nullptr,
																	 0, GL_STATIC_DRAW );

	m_iNormalBuffer = ShaderManager::getInstance()->genVertexBuffer( m_iVertexArray,
																	 2, 3, m_defaultTerrain.m_vNormals.data(),
																	 m_defaultTerrain.m_vNormals.size() * sizeof( glm::vec3 ), GL_STATIC_DRAW );

	m_iTextureBuffer = ShaderManager::getInstance()->genVertexBuffer(m_iVertexArray,
																	 3, 2, m_defaultTerrain.m_vUVs.data(),
																	 m_defaultTerrain.m_vUVs.size() * sizeof(vec2), GL_STATIC_DRAW);

  	m_iIndicesBuffer = ShaderManager::getInstance()->genIndicesBuffer( m_iVertexArray,
																	 m_defaultTerrain.m_vIndices.data(),
																	 m_defaultTerrain.m_vIndices.size() * sizeof( unsigned int ),
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
	glBufferData(GL_ARRAY_BUFFER, m_defaultTerrain.m_iUSize*m_defaultTerrain.m_iVSize * sizeof(vec3), m_defaultTerrain.m_vVertices.data(), GL_STATIC_DRAW);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, m_defaultTerrain.m_vIndices.size() * sizeof(unsigned int), m_defaultTerrain.m_vIndices.data(), GL_DYNAMIC_DRAW );

	//* Draw Points
	ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::WORLD_SHDR, "vColor", &BLACK);
	glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::WORLD_SHDR));
	glPointSize( 5.0f );  // divide by 4th positon in the projected vector
	glDrawArrays( GL_POINTS, 0, m_defaultTerrain.m_iUSize*m_defaultTerrain.m_iVSize );
//*/


	// Draw Mesh
	glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::TERRAIN_SHDR ) );
	glDrawElements( GL_TRIANGLES, m_defaultTerrain.m_vIndices.size(), GL_UNSIGNED_INT, 0 );

	ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::WORLD_SHDR, "vColor", &RED);
	glUseProgram(ShaderManager::getInstance()->getProgram(ShaderManager::eShaderType::WORLD_SHDR));
	
	vector< vec3 > vTempVerts;
	if (-1 != m_iLockedStart)
	{
		unsigned int iStartU, iStartV, iEndU, iEndV;
		getSelectedArea(iStartU, iStartV, iEndU, iEndV);;
		for (unsigned int u = iStartU; u <= iEndU; ++u)
			for (unsigned int v = iStartV; v <= iEndV; ++v)
				vTempVerts.push_back(m_defaultTerrain.m_vVertices[u + (v*m_defaultTerrain.m_iUSize)]);
	}
	else if( !m_vCurrentSubset.empty() )
	{
		glBufferData(GL_ARRAY_BUFFER, m_vCurrentSubset.size() * sizeof(vec3), m_vCurrentSubset.data(), GL_STATIC_DRAW);
		glPointSize(7.0f);  // divide by 4th positon in the projected vector
		glDrawArrays(GL_POINTS, 0, m_vCurrentSubset.size());
		glPointSize(1.0f);
		vTempVerts.push_back(m_defaultTerrain.m_vVertices[m_iSelector]);
		ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::WORLD_SHDR, "vColor", &WHITE);
	}
	else
	{
		vTempVerts.push_back(m_defaultTerrain.m_vVertices[m_iSelector]);
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
			v4 = normalize(vPosition - m_defaultTerrain.m_vVertices[iIndex4]);
			v2 = normalize(vPosition - m_defaultTerrain.m_vVertices[index2]);
			v3 = normalize(vPosition - m_defaultTerrain.m_vVertices[index3]);

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
				vec3 v1 = normalize(vPosition - m_defaultTerrain.m_vVertices[index1]);

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
	if (!m_defaultTerrain.m_vVertices.empty())
	{
		if (m_defaultTerrain.m_vHeightMap.empty())
			m_defaultTerrain.m_vHeightMap.resize(m_defaultTerrain.m_iUSize * m_defaultTerrain.m_iVSize, 0.0f);

		float fTemp;
		for (unsigned int u = 0; u < m_defaultTerrain.m_iUSize; ++u)
		{
			for (unsigned int v = 0; v < m_defaultTerrain.m_iVSize; ++v)
			{
				unsigned int iIndex = u + (v*m_defaultTerrain.m_iUSize);
				fTemp = m_defaultTerrain.m_vVertices[iIndex].y;
				m_defaultTerrain.m_vVertices[iIndex].y = m_defaultTerrain.m_vHeightMap[iIndex];
				m_defaultTerrain.m_vHeightMap[iIndex] = fTemp;

			}
		}
	}

}

void Terrain::growTerrain()
{
	grow(m_defaultTerrain);
	calculateDimensions();
	generateIndices();
}

void Terrain::grow(tMesh& terrain)
{

		flip(terrain);
		growU(terrain);

		flip(terrain);
		growU(terrain);
}

void Terrain::growU(tMesh& terrain)
{	
	vector<vec3> E;
	vector<vec3> meshV;
	vector<vec3> meshD;

	vec3 D1,D2,D3;
	

	unsigned int detailOffset = terrain.m_MrMap.empty() ? 0 : terrain.m_MrMap.top().second;
	cout << "detail Offset: " << detailOffset << endl; cin.get();
	unsigned int remainSize = 0;
	bool bExtraPoint = terrain.m_MrMap.empty() ? false : terrain.m_MrMap.top().first;

	cout << "bExtraPoint: " << bExtraPoint << endl; cin.get();

	if( !terrain.m_MrMap.empty() )
		terrain.m_MrMap.pop();

	// Extrude Details From Mesh
	for (int v = 0; v < terrain.m_iVSize - 1; v++)
	{

		if( detailOffset != 0 )
		{
			cout << "LOOKING FOR DAT D" << endl; cin.get();

			D1 = terrain.m_vVertices.at(detailOffset   + terrain.m_iUSize*v);
			D2 = terrain.m_vVertices.at(detailOffset+1 + terrain.m_iUSize*v);
			D3 = terrain.m_vVertices.at(detailOffset+2 + terrain.m_iUSize*v);

			// Apply Starting Computation
			E.push_back(vec3(0));	// E1 = 0 D1
			E.push_back(HALF * D1);	// E2 = 1/2 D1
			E.push_back((-THREEQUARTER * D1) + (QUARTER * D2));	// E3 = -3/4 D1 + 1/4 D2
			E.push_back((-QUARTER * D1) + (THREEQUARTER * D2));	// E4 = -1/4 D1 + 3/4 D2
			E.push_back((-THREEQUARTER * D2) - (QUARTER * D3));	// E5 = -3/4 D2 - 1/4 D3
			E.push_back((-QUARTER * D2) + (-THREEQUARTER * D3));	// E6 = -1/4 D2 - 3/4 D3

			unsigned int i;
			for (i = 3; i < terrain.m_vVertices.size() - detailOffset - 1; i++)
			{
				vec3 DI, DII;
				DI = terrain.m_vVertices.at(detailOffset + i);
				DII = terrain.m_vVertices.at(detailOffset + i+1);
				E.push_back((THREEQUARTER * DI) + (-QUARTER * DII));
				E.push_back((QUARTER * DI) + (-THREEQUARTER * DII));
			}


			// Final E Calculations
			vec3 DS = terrain.m_vVertices.at(i);

			E.push_back(HALF * DS);
			E.push_back(vec3(0));
			remainSize = i; // details size
			cout << "e size: " << E.size() << endl; cin.get();
		}
		else
			E.resize(terrain.m_vVertices.size()+1, vec3(0)); //FIXME: +1 was a hack, remove if we can figure out why



		cout << "e size: " << E.size() << endl; cin.get();


		meshV.push_back(terrain.m_vVertices.at(v) + E.at(v));
				cout << "vert 0 " << " {" << meshV.back().x << ", " << meshV.back().y << ", " << meshV.back().z << "}" << endl; cin.get();
		meshV.push_back((HALF * terrain.m_vVertices.at(v)) + (HALF * terrain.m_vVertices.at(v+1)) + (E.at(v+1)));
				cout << "vert 1 " << " {" << meshV.back().x << ", " << meshV.back().y << ", " << meshV.back().z << "}" <<endl; cin.get();

		unsigned int i;
		unsigned int j;
		vec3 CI,CII;
		j = 3;

		cout << "this " << terrain.m_iUSize << " should be less than " << terrain.m_vVertices.size() << endl;
		for (i = 2; i < terrain.m_iUSize - 2; i+=2)
		{	

			CI = terrain.m_vVertices.at(i);
			CII = terrain.m_vVertices.at(i+1);
			meshV.push_back((THREEQUARTER * CI) + (QUARTER * CII) + E.at(j));
			cout << "vert " << i << " {" << meshV.back().x << ", " << meshV.back().y << ", " << meshV.back().z << "}" << endl; cin.get();
			meshV.push_back((QUARTER * CI) + (THREEQUARTER * CII) + E.at(j+1));
			cout << "vert " << i+1 << " {" << meshV.back().x << ", " << meshV.back().y << ", " << meshV.back().z << "}" << endl; cin.get();


			j+=2;
		}

		CI = terrain.m_vVertices.at(i);
		CII = terrain.m_vVertices.at(i+1);
		
		meshV.push_back((HALF * CI) + (HALF * CII) + E.at(j));
		cout << "vert " << i << " {" << meshV.back().x << ", " << meshV.back().y << ", " << meshV.back().z << "}" << endl; cin.get();
		meshV.push_back(CII + E.at(j+1));
		cout << "vert " << i+1 << " {" << meshV.back().x << ", " << meshV.back().y << ", " << meshV.back().z << "}" << endl; cin.get();
	}	
	terrain.m_vVertices = meshV;


}

void Terrain::reduceTerrain()
{
	reduce(m_defaultTerrain);
	calculateDimensions();
	generateIndices();
	calculateBarries();
}

void Terrain::reduce(tMesh& terrain)
{
	if( m_defaultTerrain.m_iUSize >= MRLIMIT_MIN && m_defaultTerrain.m_iVSize >= MRLIMIT_MIN )
	{
		reduceU(terrain);
		flip(terrain);

		reduceU(terrain);
		flip(terrain);
	}


}

void Terrain::flip(tMesh& terrain)
{
	vector<vec3> flippedMesh;
	for (unsigned int u = 0; u < terrain.m_iUSize; u++)
	{
		for (unsigned int v = 0; v < terrain.m_iVSize; v++)
		{
			flippedMesh.push_back(terrain.m_vVertices.at(v*terrain.m_iUSize + u));
		}
	}

	flippedMesh.insert(flippedMesh.end(), terrain.m_vVertices.begin()+flippedMesh.size(), terrain.m_vVertices.end());
	terrain.m_vVertices = flippedMesh;

	terrain.m_iUSize = terrain.m_iUSize ^ terrain.m_iVSize;
	terrain.m_iVSize = terrain.m_iVSize ^ terrain.m_iUSize;
	terrain.m_iUSize = terrain.m_iUSize ^ terrain.m_iVSize;
}

void Terrain::reduceU(tMesh& terrain)
{
	vec3 C1, C2, C3, C4;
	vector<vec3> meshV;
	vector<vec3> meshD;
	vector< vec3 > vApplicationCurve;
	bool isOdd = terrain.m_iUSize%2 != 0;
	unsigned int tmpUSize;

	int i;

	pair <bool, unsigned int> newMrMap;

	for (unsigned int v = 0; v < terrain.m_iVSize; v++)
	{
		unsigned int vIndex = v*terrain.m_iUSize;
		for(int j = 0; j < terrain.m_iUSize; ++j)
			vApplicationCurve.push_back(terrain.m_vVertices[j+vIndex]);

		newMrMap.first = false;

		if (isOdd)
		{	
			vec3 vTranslateVector = vApplicationCurve.back() - *(vApplicationCurve.end() - 2);
			vApplicationCurve.push_back(vApplicationCurve.back() + vec3(vTranslateVector.x, 0.0f, vTranslateVector.z));
			newMrMap.first = true;
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
	newMrMap.second = meshD.size();
	m_defaultTerrain.m_MrMap.push(newMrMap);
	meshV.insert( meshV.end(), meshD.begin(), meshD.end());
	meshV.insert( meshV.end(), terrain.m_vVertices.begin()+(terrain.m_iUSize*terrain.m_iVSize), terrain.m_vVertices.end());
	terrain.m_iUSize = tmpUSize;
	terrain.m_vVertices = meshV;
}

void Terrain::get_Quad_Points( float fPosX, float fPosZ, int &iIndex1, int &iIndex2, int &iIndex3, int &iIndex4 )
{
	// Locals/Initialization
	vec3 vOffset;
	iIndex1 = iIndex2 = iIndex3 = iIndex4 = -1;

	// Check that it's within the terrain
	// Get step position
	float fTrueX = (fPosX > m_defaultTerrain.m_vStartPos.x + m_defaultTerrain.m_fWidth) ? m_defaultTerrain.m_vStartPos.x + m_defaultTerrain.m_fWidth - (0.5f * m_defaultTerrain.m_fTileWidth) : fPosX;
	float fTrueZ = (fPosZ > m_defaultTerrain.m_vStartPos.z + m_defaultTerrain.m_fDepth) ? m_defaultTerrain.m_vStartPos.z + m_defaultTerrain.m_fDepth - (0.5f * m_defaultTerrain.m_fTileDepth) : fPosZ;
	vOffset = vec3(fTrueX, 0.0, fTrueZ) - vec3(m_defaultTerrain.m_vStartPos.x, 0.0, m_defaultTerrain.m_vStartPos.z);

	int u = (int)(vOffset.x / (float)m_defaultTerrain.m_fTileWidth);
	int v = (int)(vOffset.z / (float)m_defaultTerrain.m_fTileDepth);

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
	iIndex1 = (v * m_defaultTerrain.m_iUSize) + u;
	iIndex2 = iIndex1 + 1;
	iIndex3 = iIndex1 + m_defaultTerrain.m_iUSize;
	iIndex4 = iIndex3 + 1;

}

void Terrain::get_Point_Pos(float fPosX, float fPosZ)
{
	int iIndex1, iIndex2, iIndex3, iIndex4;
	get_Quad_Points(fPosX, fPosZ, iIndex1, iIndex2, iIndex3, iIndex4);

	if (iIndex1 > -1) // Triangle Found
	{
		vec2 vCompVert = vec2(fPosX, fPosZ);
		vec2 vVert1 = vec2(m_defaultTerrain.m_vVertices[iIndex1].x, m_defaultTerrain.m_vVertices[iIndex1].z);
		vec2 vVert2 = vec2(m_defaultTerrain.m_vVertices[iIndex2].x, m_defaultTerrain.m_vVertices[iIndex2].z);
		vec2 vVert3 = vec2(m_defaultTerrain.m_vVertices[iIndex3].x, m_defaultTerrain.m_vVertices[iIndex3].z);
		vec2 vVert4 = vec2(m_defaultTerrain.m_vVertices[iIndex4].x, m_defaultTerrain.m_vVertices[iIndex4].z);
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
				m_vCurrentSubset.push_back(m_defaultTerrain.m_vVertices[u + (v*m_defaultTerrain.m_iUSize)]);

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

	iStartV = m_iLockedStart / m_defaultTerrain.m_iUSize;
	iStartU = m_iLockedStart % m_defaultTerrain.m_iUSize;
	iEndV = m_iSelector / m_defaultTerrain.m_iUSize;
	iEndU = m_iSelector % m_defaultTerrain.m_iUSize;

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
	if (!m_defaultTerrain.m_vVertices.empty())
	{
		m_defaultTerrain.m_vStartPos = m_defaultTerrain.m_vVertices.front();
		m_defaultTerrain.m_vEndPos = (m_defaultTerrain.m_iUSize == 0 || m_defaultTerrain.m_iVSize == 0) ? m_defaultTerrain.m_vVertices.back() : m_defaultTerrain.m_vVertices[m_defaultTerrain.m_iUSize * m_defaultTerrain.m_iVSize - 1];

		m_defaultTerrain.m_fWidth = abs(m_defaultTerrain.m_vEndPos.x - m_defaultTerrain.m_vStartPos.x);
		m_defaultTerrain.m_fDepth = abs(m_defaultTerrain.m_vEndPos.z - m_defaultTerrain.m_vStartPos.z);

		m_defaultTerrain.m_fTileWidth = abs(m_defaultTerrain.m_vVertices[1].x - m_defaultTerrain.m_vStartPos.x);
		

		if (m_defaultTerrain.m_iUSize == 0 || m_defaultTerrain.m_iVSize == 0)
		{
			for (vector<vec3>::const_iterator iter = m_defaultTerrain.m_vVertices.begin();
				iter->z == m_defaultTerrain.m_vStartPos.z;
				++iter)
				m_defaultTerrain.m_iUSize++;
			for (vector<vec3>::const_iterator iter = m_defaultTerrain.m_vVertices.begin();
				iter != m_defaultTerrain.m_vVertices.end() && iter->x == m_defaultTerrain.m_vStartPos.x;
				iter += m_defaultTerrain.m_iUSize)
				m_defaultTerrain.m_iVSize++;
		}
		m_defaultTerrain.m_fTileDepth = abs(m_defaultTerrain.m_vVertices[m_defaultTerrain.m_iUSize].z - m_defaultTerrain.m_vStartPos.z);
	}
	else
		cout << "Unable to compute dimensions; no vertices loaded.\n";
}

void Terrain::generateIndices( )
{
	if (!m_defaultTerrain.m_vVertices.empty())
	{
		if (!m_defaultTerrain.m_vIndices.empty())
			m_defaultTerrain.m_vIndices.clear();

		for (unsigned int v = 0; v < m_defaultTerrain.m_iVSize - 1; ++v)
		{
			for (unsigned int u = 0; u < m_defaultTerrain.m_iUSize - 1; ++u)
			{
				unsigned int iA = u + (v*m_defaultTerrain.m_iUSize);
				unsigned int iB = iA + 1;
				unsigned int iC = iA + m_defaultTerrain.m_iUSize;
				unsigned int iD = iC + 1;

				m_defaultTerrain.m_vIndices.push_back(iA);
				m_defaultTerrain.m_vIndices.push_back(iB);
				m_defaultTerrain.m_vIndices.push_back(iC);

				m_defaultTerrain.m_vIndices.push_back(iC);
				m_defaultTerrain.m_vIndices.push_back(iD);
				m_defaultTerrain.m_vIndices.push_back(iB);
			}
		}
	}
	else
		cout << "Unable to generate Indices; no vertices loaded.\n";
}

void Terrain::calculateBarries()
{
	if (!m_defaultTerrain.m_vIndices.empty())
	{
		unsigned int iArry[3] = { 0, 1, 2 };
		unsigned int iX = 0;

		vector<unsigned int> BCverts;
		BCverts.reserve(m_defaultTerrain.m_vIndices.size());

		for (unsigned int v = 0; v < m_defaultTerrain.m_iVSize; ++v)
		{
			for (unsigned int u = 0; u < m_defaultTerrain.m_iUSize; u += 3)
			{
				BCverts.push_back(iArry[iX]);
				if (u + 1 < m_defaultTerrain.m_iUSize)
					BCverts.push_back(iArry[(iX + 1) % 3]);
				if (u + 2 < m_defaultTerrain.m_iUSize)
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
	if (!m_defaultTerrain.m_vVertices.empty())
	{
		// Compute Indices if not calculated yet.
		if (m_defaultTerrain.m_vIndices.empty())
			generateIndices();

		// Vertices and Indices Loaded, need to compute Normals
		m_defaultTerrain.m_vNormals.resize(m_defaultTerrain.m_vVertices.size(), vec3(0.0));
		for (unsigned int i = 0; i < m_defaultTerrain.m_vIndices.size(); i += 3)
		{
			assert((i + 2) < m_defaultTerrain.m_vIndices.size());
			vec3 vTri[3] = { m_defaultTerrain.m_vVertices[m_defaultTerrain.m_vIndices[i]], m_defaultTerrain.m_vVertices[m_defaultTerrain.m_vIndices[i + 1]], m_defaultTerrain.m_vVertices[m_defaultTerrain.m_vIndices[i + 2]] };

			// Calculate Triangle Normal: (v1 - v0) x (v2 - v0)
			vec3 vTriNormal = cross((vTri[1] - vTri[0]), (vTri[2] - vTri[0]));

			// Accumulate Normals Per Vertex;
			m_defaultTerrain.m_vNormals[m_defaultTerrain.m_vIndices[i]] += vTriNormal;
			m_defaultTerrain.m_vNormals[m_defaultTerrain.m_vIndices[i + 1]] += vTriNormal;
			m_defaultTerrain.m_vNormals[m_defaultTerrain.m_vIndices[i + 2]] += vTriNormal;
		}

		// Normalize all Accumulated Normals
		for (vector< vec3 >::iterator vNormIter = m_defaultTerrain.m_vNormals.begin();
			vNormIter != m_defaultTerrain.m_vNormals.end();
			++vNormIter)
			(*vNormIter) = normalize((*vNormIter));
	}
	else
		cout << "Unable to compute Normals; Vertices aren't loaded.\n";
}
