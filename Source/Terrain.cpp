#include "Terrain.h"
#include "ObjLoader.h"

#define PI					3.14159265f
#define TWOPI				6.283185307f
#define MRLIMIT_MIN			5
#define MIN_BLEND_RS		3

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

	// Initialize Data.
	m_iLockedStart = -1;
	m_iSelector = 0;
	m_bHeightMapStored = m_bApplicationMode = false;

	// Load the OBJ File; pulls Verts, UVS and Normals
	objLoader::loadOBJ(pTerrLoc.data(), m_defaultTerrain.m_vVertices, m_defaultTerrain.m_vUVs, m_defaultTerrain.m_vNormals);

	// Initialize Mesh Details
	initMesh(m_defaultTerrain);

	// Generate Open GL Buffers
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
	vec3 TERRAIN_COLOR = vec3(0.586f, 0.551f, 0.598f);
	vec3 BLACK = vec3(0.0);
	vec3 BLUE = vec3(0.0, 0.0, 1.0);
	vec3 WHITE = vec3(1.0);
	vec3 GREEN = vec3(0.0, 1.0f, 0.0f);

	glBindVertexArray( m_iVertexArray );
	/* Not yet Implemented
	if ( nullptr != m_pTexture )
	{
		m_pTexture->bindTexture( ShaderManager::eShaderType::PLANE_SHDR, "gSampler" );
	} //*/
	

	/* Draw Points
	ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::WORLD_SHDR, "vColor", &BLACK);
	glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::WORLD_SHDR));
	glPointSize( 5.0f );  // divide by 4th positon in the projected vector
	glDrawArrays( GL_POINTS, 0, m_defaultTerrain.m_iUSize*m_defaultTerrain.m_iVSize );
//*/

	// Draw Additional Selections
	ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::WORLD_SHDR, "vColor", &RED );
	vector< vec3 > vTempVerts;
	vTempVerts.push_back(m_defaultTerrain.m_vVertices[m_iSelector]);

	if (-1 != m_iLockedStart || m_bApplicationMode) // Generate a list of Selected Vertices
	{
		unsigned int iStartU, iStartV, iEndU, iEndV;
		getSelectedArea(iStartU, iStartV, iEndU, iEndV);
		vTempVerts.clear();
		for (unsigned int u = iStartU; (m_bApplicationMode ? u < iEndU : u <= iEndU); ++u)
			for (unsigned int v = iStartV; (m_bApplicationMode ? v < iEndV : v <= iEndV); ++v)
				vTempVerts.push_back(m_defaultTerrain.m_vVertices[u + (v*m_defaultTerrain.m_iUSize)]);
	}
	else if( !m_vCurrentSubset.m_vVertices.empty() ) // Draw the Currently Selected Subset
	{
		ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::TERRAIN_SHDR, "vColor", &RED);
		glUseProgram(ShaderManager::getInstance()->getProgram(ShaderManager::eShaderType::TERRAIN_SHDR));
		loadMeshData(m_vCurrentSubset);
		glDrawElements(GL_TRIANGLES, m_vCurrentSubset.m_vIndices.size(), GL_UNSIGNED_INT, 0);

		// Set the Cursor Color to be White
		ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::WORLD_SHDR, "vColor", &WHITE);
	}

	if (m_bApplicationMode)
	{
		// Load Application Mesh Draw Information
		loadMeshData(m_vApplicationMesh);
		glBindBuffer(GL_ARRAY_BUFFER, m_iVertexBuffer);
		glBufferData(GL_ARRAY_BUFFER, vTempVerts.size() * sizeof(vec3), vTempVerts.data(), GL_STATIC_DRAW);
		glUseProgram(ShaderManager::getInstance()->getProgram(ShaderManager::eShaderType::TERRAIN_SHDR));
		ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::TERRAIN_SHDR, "vColor", &GREEN);
		glDrawElements(GL_TRIANGLES, m_vApplicationMesh.m_vIndices.size(), GL_UNSIGNED_INT, 0);
	}
	else
	{
		glBindBuffer(GL_ARRAY_BUFFER, m_iVertexBuffer);
		glBufferData(GL_ARRAY_BUFFER, vTempVerts.size() * sizeof(vec3), vTempVerts.data(), GL_STATIC_DRAW);
		glUseProgram(ShaderManager::getInstance()->getProgram(ShaderManager::eShaderType::WORLD_SHDR));
		glPointSize(10.0f);  // divide by 4th positon in the projected vector
		glDrawArrays(GL_POINTS, 0, vTempVerts.size());
		glPointSize(1.0f);
	}

	// Draw Saved Selection
	if (!m_vSavedSubset.m_vVertices.empty() && !m_bApplicationMode)
	{
		ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::TERRAIN_SHDR, "vColor", &BLUE);
		glUseProgram(ShaderManager::getInstance()->getProgram(ShaderManager::eShaderType::TERRAIN_SHDR));
		loadMeshData(m_vSavedSubset);
		
		glDrawElements(GL_TRIANGLES, m_vSavedSubset.m_vIndices.size(), GL_UNSIGNED_INT, 0);
		
	}

	// Draw Main Terrain
	loadMeshData(m_defaultTerrain);
	ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::TERRAIN_SHDR, "vColor", &TERRAIN_COLOR);
	glUseProgram(ShaderManager::getInstance()->getProgram(ShaderManager::eShaderType::TERRAIN_SHDR));
	glDrawElements(GL_TRIANGLES, m_defaultTerrain.m_vIndices.size(), GL_UNSIGNED_INT, 0);
	/* Not Yet Implemented
	if ( nullptr != m_pTexture )
		m_pTexture->unbindTexture();
	//*/

	glUseProgram(0);
	glBindVertexArray( 0 );
}

void Terrain::loadMeshData(const tMesh& pTerrain)
{
	// Reload Vertices
	glBindBuffer(GL_ARRAY_BUFFER, m_iVertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, pTerrain.m_iUSize*pTerrain.m_iVSize * sizeof(vec3), pTerrain.m_vVertices.data(), GL_STATIC_DRAW);

	// Reload Indices
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_iIndicesBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, pTerrain.m_vIndices.size() * sizeof(unsigned int), pTerrain.m_vIndices.data(), GL_STATIC_DRAW);

	// Reload BaryCentricCoords
	glBindBuffer(GL_ARRAY_BUFFER, m_iBaryCentric);
	glBufferData(GL_ARRAY_BUFFER, pTerrain.m_vBarries.size() * sizeof(unsigned int), pTerrain.m_vBarries.data(), GL_STATIC_DRAW);
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

// Flattens Terrain to xz-plane. Will restore the y-values if they've already been stored.
void Terrain::toggleHeightMap()
{
	// Only apply if Vertices are loaded.
	if (!m_defaultTerrain.m_vVertices.empty())
	{
		if (m_vHeightMap.empty())
			m_vHeightMap.resize(m_defaultTerrain.m_iUSize * m_defaultTerrain.m_iVSize, 0.0f);

		// Loop Through and swap HeightMap Values
		float fTemp;
		for (unsigned int u = 0; u < m_defaultTerrain.m_iUSize; ++u)
		{
			for (unsigned int v = 0; v < m_defaultTerrain.m_iVSize; ++v)
			{
				unsigned int iIndex = u + (v*m_defaultTerrain.m_iUSize);
				fTemp = m_defaultTerrain.m_vVertices[iIndex].y;
				m_defaultTerrain.m_vVertices[iIndex].y = m_vHeightMap[iIndex];
				m_vHeightMap[iIndex] = fTemp;

			}
		}

		// Toggle Height Map Check.
		m_bHeightMapStored = !m_bHeightMapStored;
	}

}

void Terrain::growTerrain()
{
	grow(m_defaultTerrain);
	calculateDimensions( m_defaultTerrain);
	generateIndices(m_defaultTerrain);
}

void Terrain::grow(tMesh& terrain)
{

		growU(terrain);
		flip(terrain);

		growU(terrain);
		flip(terrain);
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
	for (int v = 0; v < terrain.m_iVSize; v++)
	{
		unsigned int vOffset = v*terrain.m_iUSize;
		cout << "vOffset " << vOffset << endl;

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


		meshV.push_back(terrain.m_vVertices.at(vOffset) + E.at(vOffset));
				cout << "vert 0 " << " {" << meshV.back().x << ", " << meshV.back().y << ", " << meshV.back().z << "}" << endl; cin.get();
		meshV.push_back((HALF * terrain.m_vVertices.at(vOffset)) + (HALF * terrain.m_vVertices.at(vOffset+1)) + (E.at(vOffset+1)));
				cout << "vert 1 " << " {" << meshV.back().x << ", " << meshV.back().y << ", " << meshV.back().z << "}" <<endl; cin.get();
				

		unsigned int i;
		unsigned int j;
		vec3 CI,CII;
		j = 2;

		// *************E is incorrect for now since they are all zeros ****************
		cout << "loop " << endl;
		for (i = 1; i < terrain.m_iUSize - 2; i++)
		{	

			CI = terrain.m_vVertices.at(i + vOffset);
			CII = terrain.m_vVertices.at(i+1 + vOffset);
			meshV.push_back((THREEQUARTER * CI) + (QUARTER * CII) + E.at(j));
			cout << "vert " << meshV.size() << " {" << meshV.back().x << ", " << meshV.back().y << ", " << meshV.back().z << "}" << endl; cin.get();
		
			meshV.push_back((QUARTER * CI) + (THREEQUARTER * CII) + E.at(j+1));
			cout << "vert " << meshV.size() << " {" << meshV.back().x << ", " << meshV.back().y << ", " << meshV.back().z << "}" << endl; cin.get();

			j+=2;
		}
		cout << "end loop " << endl;
		CI = terrain.m_vVertices.at(i-1);
		cout << "i " << i << "/" << terrain.m_iUSize << endl;
		CII = terrain.m_vVertices.at(i);
		
		//meshV.push_back((HALF * CI) + (HALF * CII) + E.at(j)); // skip if mesh is size 2
		//cout << "vert " << meshV.size() << " {" << meshV.back().x << ", " << meshV.back().y << ", " << meshV.back().z << "}" << endl; cin.get();
		meshV.push_back(CII + E.at(j+1));
		cout << "LAST vert " << meshV.size() << " {" << meshV.back().x << ", " << meshV.back().y << ", " << meshV.back().z << "}" << endl; cin.get();

	}	


	terrain.m_vVertices = meshV;
}


/********************************************************************\
 * Reverse Subdivision Section                                      *
\********************************************************************/

// Reduces Main Terrain
void Terrain::reduceTerrain()
{
	// Turn on HeightMap before Reverse Subdividing.
	bool bReturn = m_bHeightMapStored;
	if (m_bHeightMapStored)
		toggleHeightMap();

	// Reduce
	reduce(m_defaultTerrain);

	if (bReturn)
		toggleHeightMap();
}

// Given a Terrain Mesh, Reverse Subdivide it. Set our min limit to MRLIMIT_MIN
void Terrain::reduce(tMesh& terrain)
{
	if( terrain.m_iUSize >= MRLIMIT_MIN && terrain.m_iVSize >= MRLIMIT_MIN )
	{
		// Flip to Apply U-Reverse Subdivision on V-Direction first
		flip(terrain);
		reduceU(terrain);

		// Flip again to apply U-Rverse Subdivision on U-Direction second
		flip(terrain);
		reduceU(terrain);

	}

	initMesh(terrain);
}

// Transposes the Vertices of a Mesh, preserving the Details Information.
void Terrain::flip(tMesh& terrain)
{
	//Flip Vertices
	vector<vec3> flippedMesh;
	for (unsigned int u = 0; u < terrain.m_iUSize; u++)
	{
		for (unsigned int v = 0; v < terrain.m_iVSize; v++)
		{
			flippedMesh.push_back(terrain.m_vVertices.at(v*terrain.m_iUSize + u));
		}
	}

	// Append on Details Information
	flippedMesh.insert(flippedMesh.end(), terrain.m_vVertices.begin()+flippedMesh.size(), terrain.m_vVertices.end());
	terrain.m_vVertices = flippedMesh;

	// Swap U and V sizes
	terrain.m_iUSize = terrain.m_iUSize ^ terrain.m_iVSize;
	terrain.m_iVSize = terrain.m_iVSize ^ terrain.m_iUSize;
	terrain.m_iUSize = terrain.m_iUSize ^ terrain.m_iVSize;
}

// Apply Reverse Subdivision on terrain (in the U-direction)
void Terrain::reduceU(tMesh& terrain)
{
	// Local Variables
	vec3 C1, C2, C3, C4;
	vector<vec3> meshV;
	vector<vec3> meshD;
	vector< vec3 > vApplicationCurve;
	bool isOdd = terrain.m_iUSize%2 != 0;
	unsigned int tmpUSize;
	unsigned int i;

	pair <bool, unsigned int> newMrMap;

	// Loop for Each U-Curve on Mesh
	for (unsigned int v = 0; v < terrain.m_iVSize; v++)
	{
		// Index into Row of Mesh
		unsigned int vIndex = v*terrain.m_iUSize;

		// Loop and store each relevant control point into curve.
		for(unsigned int j = 0; j < terrain.m_iUSize; ++j)
			vApplicationCurve.push_back(terrain.m_vVertices[j+vIndex]);

		newMrMap.first = false;

		// Add an additional point if Dealing with Odd Vertices. Translate it out by Tile Distance for more uniform division
		if (isOdd)
		{	
			vec3 vTranslateVector = vApplicationCurve.back() - *(vApplicationCurve.end() - 2);
			vApplicationCurve.push_back(vApplicationCurve.back() + vec3(vTranslateVector.x, 0.0f, vTranslateVector.z));
			newMrMap.first = true;
		}

		// Initialize Index
		i = 0;

		// First 4 Points
		C1 = vApplicationCurve.at( i ); 
		C2 = vApplicationCurve.at(i+1);
		C3 = vApplicationCurve.at(i+2);
		C4 = vApplicationCurve.at(i+3);

		// Apply Boundary Conditions for Vertices
		meshV.push_back(C1);
		meshV.push_back(((-HALF)*C1) + (1*C2) + ((THREEQUARTER)*C3) + ((-QUARTER)*C4));

		// Apply Boundary Conditions for Details
		meshD.push_back(((-HALF)*C1) + (1*C2) + ((-THREEQUARTER)*C3) + ((-QUARTER)*C4));

		C1 = vApplicationCurve.at( i+2 ); 
		C2 = vApplicationCurve.at(i+3);
		C3 = vApplicationCurve.at(i+4);
		C4 = vApplicationCurve.at(i+5);
		meshD.push_back(((-QUARTER)*C1) + ((THREEQUARTER)*C2) + ((-THREEQUARTER)*C3) + ((QUARTER)*C4));


		// Recursively Compute Information For middle curve
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
	
		// Apply End Boundary Filters
		C1 = vApplicationCurve.at( i );
		C2 = vApplicationCurve.at(i+1);
		C3 = vApplicationCurve.at(i+2);
		C4 = vApplicationCurve.at(i+3);

		meshV.push_back(((-QUARTER)*C1) + ((THREEQUARTER)*C2) + ((1)*C3) + ((-HALF)*C4));
		meshD.push_back(((QUARTER)*C1) - ((THREEQUARTER)*C2) + ((1)*C3) - ((HALF)*C4));
			
		meshV.push_back(C4);

		// Reset Curve.
		vApplicationCurve.clear();

		// Store first size for new Usize.
		if( v == 0 )
			tmpUSize = meshV.size();
	}
	newMrMap.second = meshD.size();
	m_defaultTerrain.m_MrMap.push(newMrMap);

	// Append Details onto Coarse Points
	meshV.insert( meshV.end(), meshD.begin(), meshD.end());

	// Apply any additional Details from past Reverse Subdivision
	meshV.insert( meshV.end(), terrain.m_vVertices.begin()+(terrain.m_iUSize*terrain.m_iVSize), terrain.m_vVertices.end());

	// New USize
	terrain.m_iUSize = tmpUSize;

	// Store Vertices.
	terrain.m_vVertices = meshV;
}

void Terrain::get_Quad_Points( float fPosX, float fPosZ, int &iIndex1, int &iIndex2, int &iIndex3, int &iIndex4 )
{
	// Locals/Initialization
	vec3 vOffset;
	iIndex1 = iIndex2 = iIndex3 = iIndex4 = -1;
	int u, v;

	// Check that it's within the terrain
	// Get step position
	if (!m_bApplicationMode)
	{
		u = (fPosX > m_defaultTerrain.m_vEndPos.x) ? m_defaultTerrain.m_iUSize - 2 : -1;
		v = (fPosZ < m_defaultTerrain.m_vStartPos.z) ? 0 : -1;
		u = (fPosX < m_defaultTerrain.m_vStartPos.x) ? 0 : u;
		v = (fPosZ > m_defaultTerrain.m_vEndPos.z) ? m_defaultTerrain.m_iVSize - 2 : v;
	}
	else // Limit Boundaries based on Application Mesh Size.
	{
		unsigned int uOffset = m_vApplicationMesh.m_iUSize >> 1;
		unsigned int vOffset = m_vApplicationMesh.m_iVSize >> 1;
		float fHalfAppWidth = uOffset * m_defaultTerrain.m_fTileWidth;
		float fHalfAppDepth = vOffset * m_defaultTerrain.m_fTileDepth;
		u = (fPosX > (m_defaultTerrain.m_vEndPos.x - fHalfAppWidth)) ? (m_defaultTerrain.m_iUSize - 2) - uOffset : -1;
		v = (fPosZ < m_defaultTerrain.m_vStartPos.z + fHalfAppDepth) ? vOffset : -1;
		u = (fPosX < m_defaultTerrain.m_vStartPos.x + fHalfAppWidth) ? uOffset : u;
		v = (fPosZ > m_defaultTerrain.m_vEndPos.z - fHalfAppDepth) ? m_defaultTerrain.m_iVSize - 2 - vOffset : v;
	}

	vOffset = vec3(fPosX, 0.0, fPosZ) - vec3(m_defaultTerrain.m_vStartPos.x, 0.0, m_defaultTerrain.m_vStartPos.z);

	if( -1 == u )
		u = (int)(vOffset.x / (float)m_defaultTerrain.m_fTileWidth);

	if( -1 == v )
		v = (int)(vOffset.z / (float)m_defaultTerrain.m_fTileDepth);

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
		m_iSelector = iIndex1;
		// Vert2 is closer
		if ((fCompDist = length(vCompVert - vVert2)) < fDist)
		{
			fDist = fCompDist;
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

// Locks the Lock Point for Selecting a Subsection of the mesh.
void Terrain::lockPoint()
{
	if (!m_bApplicationMode)
	{
		// Nothing locked yet? -> Lock The current selection
		if (-1 == m_iLockedStart || m_iLockedStart == m_iSelector)
		{
			m_iLockedStart = m_iSelector;
			m_vCurrentSubset.clear();
		}
		else
		{	// Second Lock -> Store the set.
			bool bHeightActive = m_bHeightMapStored;

			if (m_bHeightMapStored)
				toggleHeightMap();

			unsigned int iStartU, iStartV, iEndU, iEndV;
			getSelectedArea(iStartU, iStartV, iEndU, iEndV);
			m_vCurrentSubset.m_iUSize = iEndU - iStartU + 1;
			m_vCurrentSubset.m_iVSize = iEndV - iStartV + 1;

			if (m_vCurrentSubset.m_iUSize > 1 && m_vCurrentSubset.m_iVSize > 1)
			{
				for (unsigned int v = iStartV; v <= iEndV; ++v)
				{
					for (unsigned int u = iStartU; u <= iEndU; ++u)
					{
						m_vCurrentSubset.m_vVertices.push_back(m_defaultTerrain.m_vVertices[u + (v*m_defaultTerrain.m_iUSize)]);

					}
				}

				// Compute Values for Current Subset.
				initMesh(m_vCurrentSubset);
			}

			if (bHeightActive)
				toggleHeightMap();

			// Reset Lock for new selection
			m_iLockedStart = -1;
		}
	}
	else
		blendMesh();
}

/***********************************************************************\
 * Terrain Blending Logic                                              *
\***********************************************************************/
void Terrain::applyTerrain(const Terrain* pTerrain)
{
	if (nullptr != pTerrain)
	{
		// Grab one of the stored Meshes in the referenced Terrain
		if (!pTerrain->m_vSavedSubset.m_vVertices.empty())
			m_vApplicationMesh = pTerrain->m_vSavedSubset;
		else if (!pTerrain->m_vCurrentSubset.m_vVertices.empty())
			m_vApplicationMesh = pTerrain->m_vCurrentSubset;

		// Evaluate Pulled Application Mesh
		if (m_vApplicationMesh.m_vVertices.empty())
			cout << "No Mesh to Apply.\n";
		else if (m_vApplicationMesh.m_fWidth > m_defaultTerrain.m_fWidth || m_vApplicationMesh.m_fDepth > m_defaultTerrain.m_fDepth)
			cout << "Mesh too Large to blend.\n";
		else	// Attempt to Subdivide Mesh into blendable area.
		{
			// Apply Reverse Subdivision a minimum number of times.
			for (unsigned int i = 0; i < MIN_BLEND_RS; ++i)
				reduce(m_vApplicationMesh);

			// Further Reverse Subdivide until able to apply Mesh
			while (m_vApplicationMesh.m_iUSize > m_defaultTerrain.m_iUSize && m_vApplicationMesh.m_iVSize > m_defaultTerrain.m_iVSize)
			{
				reduce(m_vApplicationMesh);

				if (m_vApplicationMesh.m_iUSize < MRLIMIT_MIN || m_vApplicationMesh.m_iVSize < MRLIMIT_MIN)
					break;
			}

			// Check to make sure it's possible to subdivide
			if (m_vApplicationMesh.m_iUSize > m_defaultTerrain.m_iUSize || m_vApplicationMesh.m_iVSize > m_defaultTerrain.m_iVSize)
			{
				cout << "Unable to Subdivide Enough to fit Base Mesh.\n";
				m_vApplicationMesh.clear();
			}
			else  // Set Environment for Applying Mesh
			{
				// Force Height Map active
				if (m_bHeightMapStored)
					toggleHeightMap();

				// Swap to Application Mode.
				m_bApplicationMode = true;

				m_vCurrentSubset.clear();
				m_iLockedStart = -1;
				m_vSavedSubset.clear();
			}
				   
		}
		
	}
}

void Terrain::blendMesh()
{
	unsigned int iStartU, iStartV, iEndU, iEndV;
	getSelectedArea(iStartU, iStartV, iEndU, iEndV);
	unsigned int iTerrIndex = m_defaultTerrain.m_iUSize * m_defaultTerrain.m_iVSize;
	unsigned int iAppIndex = m_vApplicationMesh.m_iUSize * m_vApplicationMesh.m_iVSize;
	unsigned int iAppUIndex = m_vApplicationMesh.m_iUSize;
	unsigned int iAppVIndex = m_vApplicationMesh.m_iVSize;
	unsigned int iTerrUIndex = m_defaultTerrain.m_iUSize;
	unsigned int iTerrVIndex = m_defaultTerrain.m_iVSize;
	stack<pair<bool, unsigned int>> pMrMapCollector;
	stack<pair<bool, unsigned int>> pAppMapCollector = m_vApplicationMesh.m_MrMap;
	pair<bool, unsigned int> pNextMap;

	blendCoarsePoints(iStartU, iStartV, iEndU, iEndV);

	int iCount = 0;

	// Replace existing Details
	while (!m_defaultTerrain.m_MrMap.empty() && !pAppMapCollector.empty())
	{
		iCount += replaceSquaredDetails(iStartU, iStartV, iEndU, iEndV, iTerrIndex, iTerrUIndex, iTerrVIndex, iAppUIndex, iAppVIndex);

		// Something went seriously wrong if the MapStacks aren't even-sized.
		assert(0 != ((pAppMapCollector.size() + m_defaultTerrain.m_MrMap.size()) % 2));
		for (int i = 0; i < 2; ++i)
		{
			// Pull Mapping for First Subdivision
			pNextMap = m_defaultTerrain.m_MrMap.top();
			m_defaultTerrain.m_MrMap.pop();

			// Modify Vertex Addition
			pNextMap.first ^= pAppMapCollector.top().first;
			pAppMapCollector.pop();

			// Collect an Ordered Division Map
			pMrMapCollector.push(pNextMap);
		}
	
	}
	
	// Completed all pending Subdivision Details for the main mesh, need to apply any remaining Subdivision Calls for the Target Mesh
	while (!pAppMapCollector.empty())
	{
		m_defaultTerrain.m_vVertices.insert(m_defaultTerrain.m_vVertices.end(), m_defaultTerrain.m_vVertices.size() * 3, vec3(0.0f)); // GOOD!
		iCount += replaceSquaredDetails(iStartU, iStartV, iEndU, iEndV, iTerrIndex, iTerrUIndex, iTerrVIndex, iAppUIndex, iAppVIndex);

		// Pop Twice for a Squared Subdivision
		pMrMapCollector.push(pAppMapCollector.top());
		pAppMapCollector.pop();
		pMrMapCollector.push(pAppMapCollector.top());
		pAppMapCollector.pop();
	}

	//. Populate the Correct MR Map for the Main Terrain.
	while (!pMrMapCollector.empty())
	{
		m_defaultTerrain.m_MrMap.push(pMrMapCollector.top());
		pMrMapCollector.pop();
	}

	cout << "Applied Blend, Touched " << iCount << " Vertices.\n";
	m_bApplicationMode = false;
}

int Terrain::blendCoarsePoints(unsigned int iStartU, unsigned int iStartV, unsigned int iEndU, unsigned int iEndV)
{
	int iCount = 0;
	// Average y values
	for (unsigned int u = iStartU; u < iEndU; ++u)
		for (unsigned int v = iStartV; v < iEndV; ++v)
		{
			m_defaultTerrain.m_vVertices[u + (v*m_defaultTerrain.m_iUSize)].y *= m_vApplicationMesh.m_vVertices[(u - iStartU) + ((v - iStartV) * m_vApplicationMesh.m_iUSize)].y;
			m_defaultTerrain.m_vVertices[u + (v*m_defaultTerrain.m_iUSize)].y *= 0.5f;
			iCount++;
		}

	return iCount;
}

int Terrain::replaceSquaredDetails(unsigned int iStartU, unsigned int iStartV, unsigned int iEndU, unsigned int iEndV, unsigned int& iTerrIndex, 
	                                unsigned int& iTerrUSize, unsigned int& iTerrVSize, unsigned int& iAppUSize, unsigned int& iAppVSize)
{
	int iCount = 0;
	// Apply 1 Squared Level of Details
	// U then V (Reverse Order from Reverse Subdivision
	iTerrUSize += m_defaultTerrain.m_iUSize;
	iAppUSize += m_vApplicationMesh.m_iUSize;
	iCount += replaceDetails(iStartU + iTerrIndex, iStartV + iTerrIndex, iEndU + iTerrIndex, iEndV + iTerrIndex, iTerrUSize, iAppUSize);
	iTerrVSize += m_defaultTerrain.m_iVSize;
	iAppVSize += m_vApplicationMesh.m_iVSize;
	iCount += replaceDetails(iStartU + iTerrIndex, iStartV + iTerrIndex, iEndU + iTerrIndex, iEndV + iTerrIndex, iTerrVSize, iAppVSize);
	iTerrIndex += iTerrIndex;

	return iCount;
}

int Terrain::replaceDetails(unsigned int iStartU, unsigned int iStartV, unsigned int iEndU, unsigned int iEndV, unsigned int iTerrRowSize, unsigned int iAppRowSize )
{
	int iCount = 0;
	for (unsigned int u = iStartU; u < iEndU; ++u)
		for (unsigned int v = iStartV; v < iEndV; ++v)
		{
			m_defaultTerrain.m_vVertices[u + (v*iTerrRowSize)] = m_vApplicationMesh.m_vVertices[(u - iStartU) + ((v - iStartV) * iAppRowSize)];
			iCount++;
		}

	return iCount;
}

// Saves the Selection into a buffer.
void Terrain::saveSelection()
{
	// Make sure we have a Subset to store.
	if (!m_vCurrentSubset.m_vVertices.empty())
	{
		m_vSavedSubset = m_vCurrentSubset;
		m_vCurrentSubset.clear();
	}
}

void Terrain::clearSelection()
{
	if (!m_vCurrentSubset.m_vVertices.empty())
		m_vCurrentSubset.clear();
	else if (!m_vSavedSubset.m_vVertices.empty())
		m_vSavedSubset.clear();
}

// Get A UV Start and End that is in the correct order.
void Terrain::getSelectedArea(unsigned int &iStartU, unsigned int &iStartV, unsigned int& iEndU, unsigned int& iEndV)
{
	unsigned int iTemp;

	iEndV = m_iSelector / m_defaultTerrain.m_iUSize;
	iEndU = m_iSelector % m_defaultTerrain.m_iUSize;

	if (!m_bApplicationMode)
	{
		// Compute Us and Vs for start and end positions
		iStartV = m_iLockedStart / m_defaultTerrain.m_iUSize;
		iStartU = m_iLockedStart % m_defaultTerrain.m_iUSize;

		// Order them in correct order.
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
	else
	{
		unsigned int iHalfU = m_vApplicationMesh.m_iUSize >> 1;
		unsigned int iHalfV = m_vApplicationMesh.m_iVSize >> 1;

		

		iStartU = iEndU - iHalfU;
		iStartV = iEndV - iHalfV;
		iEndU += iHalfU - (m_vApplicationMesh.m_iUSize % 2);
		iEndV += iHalfV - (m_vApplicationMesh.m_iVSize % 2);
	}
}

void Terrain::initMesh(tMesh& pTerrain)
{
	calculateDimensions(pTerrain);
	generateIndices(pTerrain);
	calculateBarries(pTerrain);
	generateNormals(pTerrain);
}

// General function to compute Terrain Dimensions
void Terrain::calculateDimensions(tMesh& pTerrain)
{
	// Ensure that Vertices are loaded.
	if (!pTerrain.m_vVertices.empty())
	{
		// Get Start and End Position for Terrain Intersections
		pTerrain.m_vStartPos = pTerrain.m_vVertices.front();
		pTerrain.m_vEndPos = (pTerrain.m_iUSize == 0 || pTerrain.m_iVSize == 0) ? pTerrain.m_vVertices.back() : pTerrain.m_vVertices[pTerrain.m_iUSize * pTerrain.m_iVSize - 1];

		// Compute Width and Depth for Terrain Intersections
		pTerrain.m_fWidth = abs(pTerrain.m_vEndPos.x - pTerrain.m_vStartPos.x);
		pTerrain.m_fDepth = abs(pTerrain.m_vEndPos.z - pTerrain.m_vStartPos.z);

		// Compute U and V sizes if they haven't been initialized yet.
		if (pTerrain.m_iUSize == 0 || pTerrain.m_iVSize == 0)
		{
			for (vector<vec3>::const_iterator iter = pTerrain.m_vVertices.begin();
				iter->z == pTerrain.m_vStartPos.z;
				++iter)
				pTerrain.m_iUSize++;
			for (vector<vec3>::const_iterator iter = pTerrain.m_vVertices.begin();
				iter != pTerrain.m_vVertices.end() && iter->x == pTerrain.m_vStartPos.x;
				iter += pTerrain.m_iUSize)
				pTerrain.m_iVSize++;
		}

		// Compute TileDimensions
		pTerrain.m_fTileWidth = abs(pTerrain.m_vVertices[1].x - pTerrain.m_vStartPos.x);
		pTerrain.m_fTileDepth = abs(pTerrain.m_vVertices[pTerrain.m_iUSize].z - pTerrain.m_vStartPos.z);
	}
	else
		cout << "Unable to compute dimensions; no vertices loaded.\n";
}

// Computes the Indices for a Mesh when an ordered set of vertices is loaded.
void Terrain::generateIndices( tMesh& pTerrain )
{
	if (!pTerrain.m_vVertices.empty())
	{
		if (!pTerrain.m_vIndices.empty())
			pTerrain.m_vIndices.clear();

		/* Compute as:
			1----2,6
			|   /|
			|  / |
			| /  |
			|/   | 
		  3,4----5
		*/
		for (unsigned int v = 0; v < pTerrain.m_iVSize - 1; ++v)
		{
			for (unsigned int u = 0; u < pTerrain.m_iUSize - 1; ++u)
			{
				unsigned int iA = u + (v*pTerrain.m_iUSize);
				unsigned int iB = iA + 1;
				unsigned int iC = iA + pTerrain.m_iUSize;
				unsigned int iD = iC + 1;

				pTerrain.m_vIndices.push_back(iA);
				pTerrain.m_vIndices.push_back(iB);
				pTerrain.m_vIndices.push_back(iC);

				pTerrain.m_vIndices.push_back(iC);
				pTerrain.m_vIndices.push_back(iD);
				pTerrain.m_vIndices.push_back(iB);
			}
		}
	}
	else
		cout << "Unable to generate Indices; no vertices loaded.\n";
}

// Compute Barycentric Coordinates for rendering Wireframe on Terrain.
void Terrain::calculateBarries(tMesh& pTerrain)
{
	// Ensure Indices are loaded first.
	if (!pTerrain.m_vIndices.empty())
	{
		if (!pTerrain.m_vBarries.empty())
			pTerrain.m_vBarries.clear();

		// Local Variables
		unsigned int iArry[3] = { 0, 1, 2 };
		unsigned int iX = 0;
		pTerrain.m_vBarries.reserve(m_defaultTerrain.m_vIndices.size());

		// Loop through Each Point and apply the Barycentric coordinates
		for (unsigned int v = 0; v < pTerrain.m_iVSize; ++v)
		{
			for (unsigned int u = 0; u < pTerrain.m_iUSize; u += 3)
			{
				pTerrain.m_vBarries.push_back(iArry[iX]);
				if (u + 1 < pTerrain.m_iUSize)
					pTerrain.m_vBarries.push_back(iArry[(iX + 1) % 3]);
				if (u + 2 < pTerrain.m_iUSize)
					pTerrain.m_vBarries.push_back(iArry[(iX + 2) % 3]);
			}
			iX = (iX + 2) % 3;
		}
	}
}

// Computes and Generates Normals per vertex.
void Terrain::generateNormals(tMesh& pTerrain)
{
	// Only attempt to compute if Vertices are loaded.
	if (!pTerrain.m_vVertices.empty())
	{
		// Compute Indices if not calculated yet.
		if (pTerrain.m_vIndices.empty())
			generateIndices(pTerrain);

		// Vertices and Indices Loaded, need to compute Normals
		pTerrain.m_vNormals.resize(pTerrain.m_vVertices.size(), vec3(0.0));
		for (unsigned int i = 0; i < pTerrain.m_vIndices.size(); i += 3)
		{
			assert((i + 2) < pTerrain.m_vIndices.size());
			vec3 vTri[3] = { pTerrain.m_vVertices[pTerrain.m_vIndices[i]], pTerrain.m_vVertices[pTerrain.m_vIndices[i + 1]], pTerrain.m_vVertices[pTerrain.m_vIndices[i + 2]] };

			// Calculate Triangle Normal: (v1 - v0) x (v2 - v0)
			vec3 vTriNormal = cross((vTri[1] - vTri[0]), (vTri[2] - vTri[0]));

			// Accumulate Normals Per Vertex;
			pTerrain.m_vNormals[pTerrain.m_vIndices[i]] += vTriNormal;
			pTerrain.m_vNormals[pTerrain.m_vIndices[i + 1]] += vTriNormal;
			pTerrain.m_vNormals[pTerrain.m_vIndices[i + 2]] += vTriNormal;
		}

		// Normalize all Accumulated Normals
		for (vector< vec3 >::iterator vNormIter = pTerrain.m_vNormals.begin();
			vNormIter != pTerrain.m_vNormals.end();
			++vNormIter)
			(*vNormIter) = normalize((*vNormIter));
	}
	else
		cout << "Unable to compute Normals; Vertices aren't loaded.\n";
}
