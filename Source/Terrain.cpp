#include "ObjLoader.h"
#include "GraphicsManager.h" // #include "Terrain.h"
#include "Mouse_Handler.h"

#define PI					3.14159265f
#define TWOPI				6.283185307f
#define MRLIMIT_MIN			5
#define MIN_BLEND_RS		3
#define MIN_SELECTION_SIZE	16

#define HALF 0.5f
#define QUARTER 0.25f
#define THREEQUARTER .75f

// Constructor
Terrain::Terrain(const string& pTerrLoc)
{
	// Set up generation for Terrain
	// Depth = V
	// Width = U
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
	vTempVerts.push_back(m_defaultTerrain.m_vVertices[m_defaultTerrain.coarseToIndexSpace(m_iSelector)]);

	if (-1 != m_iLockedStart || m_bApplicationMode) // Generate a list of Selected Vertices
	{
		unsigned int iStartU, iStartV, iEndU, iEndV;
		bool bValidArea = getSelectedArea(iStartU, iStartV, iEndU, iEndV, false);
		if (bValidArea)
		{
			vTempVerts.clear();
			for (unsigned int u = iStartU; u <= iEndU; ++u)
				for (unsigned int v = iStartV; v <= iEndV; ++v)
				{
					unsigned int iIndexU, iIndexV;
					iIndexU = u;
					iIndexV = v;
					m_defaultTerrain.coarseToIndexSpace(iIndexU, iIndexV);
					vTempVerts.push_back(m_defaultTerrain.m_vVertices[iIndexU + (iIndexV*m_defaultTerrain.m_iUSize)]);
				}
		}
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
	glBufferData(GL_ARRAY_BUFFER, pTerrain.m_vVertices.size() * sizeof(vec3), pTerrain.m_vVertices.data(), GL_STATIC_DRAW);

	// Reload Indices
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_iIndicesBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, pTerrain.m_vIndices.size() * sizeof(unsigned int), pTerrain.m_vIndices.data(), GL_STATIC_DRAW);

	// Reload BaryCentricCoords
	glBindBuffer(GL_ARRAY_BUFFER, m_iBaryCentric);
	glBufferData(GL_ARRAY_BUFFER, pTerrain.m_vBarries.size() * sizeof(unsigned int), pTerrain.m_vBarries.data(), GL_STATIC_DRAW);
	
	// Reload Normals
	glBindBuffer(GL_ARRAY_BUFFER, m_iNormalBuffer);
	glBufferData(GL_ARRAY_BUFFER, pTerrain.m_vNormals.size() * sizeof(vec3), pTerrain.m_vNormals.data(), GL_STATIC_DRAW);
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

		flip(terrain);
		growU(terrain);

		flip(terrain);
		growU(terrain);
}

void Terrain::growU(tMesh& terrain)
{	
	/*
	vector<vec3> E;
	vector<vec3> meshV;
	vector<vec3> meshD;

	vec3 D1,D2,D3;
	
	unsigned int detailOffset = terrain.m_MrMap.empty() ? 0 : terrain.m_MrMap.top().second;
	unsigned int remainSize = 0;
	bool bExtraPoint = terrain.m_MrMap.empty() ? false : terrain.m_MrMap.top().first;

	if( !terrain.m_MrMap.empty() )
		terrain.m_MrMap.pop();

	// Extrude Details From Mesh
	if( detailOffset != 0 )
	{
		D1 = terrain.m_vVertices.at(detailOffset  );
		D2 = terrain.m_vVertices.at(detailOffset+1);
		D3 = terrain.m_vVertices.at(detailOffset+2);

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
	}
	else
		E.resize(detailOffset, vec3(0));


	meshV.push_back(terrain.m_vVertices.at(0) + E.at(0));
	meshV.push_back((HALF * terrain.m_vVertices.at(0)) + (HALF * terrain.m_vVertices.at(1)) + (E.at(1)));

	
	unsigned int i;
	unsigned int j;
	vec3 CI,CII;
	j = 3;
	for (i = 2; i < terrain.m_vVertices.size(); i++)
	{	
		CI = terrain.m_vVertices.at(i);
		CII = terrain.m_vVertices.at(i+1);
		meshV.push_back((THREEQUARTER * CI) + (QUARTER * CII) + E.at(j));
		meshV.push_back((QUARTER * CI) + (THREEQUARTER * CII) + E.at(j+1));

		j+=2;
	}
	CI = terrain.m_vVertices.at(i);
	CII = terrain.m_vVertices.at(i+1);
	
	meshV.push_back((HALF * CI) + (HALF * CII) + E.at(j));
	meshV.push_back(CII + E.at(j+1));

	terrain.m_vVertices = meshV;
//*/

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

	m_iSelector = 0;
}

// Given a Terrain Mesh, Reverse Subdivide it. Set our min limit to MRLIMIT_MIN
void Terrain::reduce(tMesh& terrain)
{
	// Get CoarseSize to see if Reverse Subdivision is possible.
	unsigned int iCoarseUSize = terrain.getCoarseUSize();

	if( iCoarseUSize >= MRLIMIT_MIN && iCoarseUSize >= MRLIMIT_MIN )
	{
		// Apply Reverse Subdivision on Us
		reduceU(terrain);
		// Flip terrain, apply Reverse Subdivision on Vs
		terrain.m_iUSize += (iCoarseUSize & 1);
		flip(terrain);
		reduceV(terrain);
		terrain.m_iUSize += (iCoarseUSize & 1);
		flip(terrain);

		// Adjust for New Details
		terrain.m_iStep <<= 1;
		terrain.m_iExp++;
		terrain.m_bAddedPointFlags.push_back((iCoarseUSize & 1) != 0);

		initMesh(terrain);
	}
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
	terrain.m_vVertices = flippedMesh;

	// Swap U and V sizes
	terrain.m_iUSize = terrain.m_iUSize ^ terrain.m_iVSize;
	terrain.m_iVSize = terrain.m_iVSize ^ terrain.m_iUSize;
	terrain.m_iUSize = terrain.m_iUSize ^ terrain.m_iVSize;
}

void Terrain::applyUReverseSubdivision(const vector< vec3 >& vApplicationCurve, vector<vec3>::iterator& vMeshInserter, unsigned int iStepSize )
{

	// Initialize Index
	unsigned int i = 0;

	// First 4 Points
	vec3 F1 = vApplicationCurve.at(i);
	vec3 F2 = vApplicationCurve.at(i + 1);
	vec3 F3 = vApplicationCurve.at(i + 2);
	vec3 F4 = vApplicationCurve.at(i + 3);

	// Apply Boundary Conditions for Vertices
	//meshV.push_back(C1);
	//*vMeshInserter = C1;																		// C1 = F1
	vMeshInserter += iStepSize;
	// Apply Boundary Conditions for Details
	//meshD.push_back(((-HALF)*C1) + (1 * C2) + ((-THREEQUARTER)*C3) + ((-QUARTER)*C4));
	*vMeshInserter = ((-HALF)*F1) + F2 - ((THREEQUARTER)*F3) + ((QUARTER)*F4);			// D1 = 1/2F1 + F2 - 3/4F2 + 1/4F4
	vMeshInserter += iStepSize;

	//meshV.push_back(((-HALF)*C1) + (1*C2) + ((THREEQUARTER)*C3) + ((-QUARTER)*C4));
	*vMeshInserter = ((-HALF)*F1) + F2 + ((THREEQUARTER)*F3) - ((QUARTER)*F4);			// C2 = -1/2F1 + F2 + 3/4F3 - 1/4F4
	vMeshInserter += iStepSize;

	F1 = vApplicationCurve.at(i + 2);
	F2 = vApplicationCurve.at(i + 3);
	F3 = vApplicationCurve.at(i + 4);
	F4 = vApplicationCurve.at(i + 5);
	//meshD.push_back(((-QUARTER)*C1) + ((THREEQUARTER)*C2) + ((-THREEQUARTER)*C3) + ((QUARTER)*C4));
	*vMeshInserter = ((-QUARTER)*F1) + ((THREEQUARTER)*F2) - ((THREEQUARTER)*F3) + ((QUARTER)*F4); // D2 = -1/4F3 + 3/4F4 - 3/4F5 + 1/4F6
	vMeshInserter += iStepSize;

	// Recursively Compute Information For middle curve
	for (i = 2; i < vApplicationCurve.size() - 5; i += 2)
	{
		F1 = vApplicationCurve.at(i);
		F2 = vApplicationCurve.at(i + 1);
		F3 = vApplicationCurve.at(i + 2);
		F4 = vApplicationCurve.at(i + 3);

		// Details
		if (i >= 4)
		{
			//meshD.push_back(((QUARTER)*C1) - ((THREEQUARTER)*C2) + ((THREEQUARTER)*C3) - ((QUARTER)*C4));
			*vMeshInserter = ((QUARTER)*F1) - ((THREEQUARTER)*F2) + ((THREEQUARTER)*F3) - ((QUARTER)*F4); // Dj = 1/4Fi - 3/4Fi+1 + 3/4Fi+2 - 1/4Fi+3
			vMeshInserter += iStepSize;
		}

		// Coarse Point
		//meshV.push_back(((-QUARTER)*C1) + ((THREEQUARTER)*C2) + ((THREEQUARTER)*C3) + ((-QUARTER)*C4));
		*vMeshInserter = ((-QUARTER)*F1) + ((THREEQUARTER)*F2) + ((THREEQUARTER)*F3) - ((QUARTER)*F4); // Cj = -1/4Fi + 3/4Fi+1 + 3/4Fi+2 - 1/4Fi+3
		vMeshInserter += iStepSize;
	}

	// m-3 zero based
	i = vApplicationCurve.size() - 4;
	// Apply End Boundary Filters
	F1 = vApplicationCurve.at(i);
	F2 = vApplicationCurve.at(i + 1);
	F3 = vApplicationCurve.at(i + 2);
	F4 = vApplicationCurve.at(i + 3);

	//meshD.push_back(((QUARTER)*C1) - ((THREEQUARTER)*C2) + ((1)*C3) - ((HALF)*C4));
	*vMeshInserter = ((QUARTER)*F1) - ((THREEQUARTER)*F2) + F3 - ((HALF)*F4); // Dj = 1/4Fm-3 - 3/4Fm-2 + Fm-1 - 1/2Fm
	++vMeshInserter;

	//meshV.push_back(((-QUARTER)*C1) + ((THREEQUARTER)*C2) + ((1)*C3) + ((-HALF)*C4));
	*vMeshInserter = ((-QUARTER)*F1) + ((THREEQUARTER)*F2) + F3 - ((HALF)*F4); // Cj = -1/4Fm-3 + 3/4Fm-2 + Fm-1 -1/2Fm
	++vMeshInserter;
}

// Apply Reverse Subdivision on terrain (in the U-direction)
void Terrain::reduceU(tMesh& terrain)
{
	// Local Variables
	vec3 F1, F2, F3, F4;
	vector< vec3 > vApplicationCurve;
	vector< vec3 >::iterator vMeshInserter = terrain.m_vVertices.begin();
	unsigned int iCoarseUSize = terrain.getCoarseUSize();
	bool isOdd = (iCoarseUSize & 1) != 0;
	unsigned int vIndex = 0;

	// Loop for Each U-Curve on Mesh
	for (unsigned int v = 0; v < terrain.m_iVSize; v++)
	{
		// Loop and store each relevant control point into curve.
		for (unsigned int j = 0; j < (iCoarseUSize >> terrain.m_iExp); ++j )
			vApplicationCurve.push_back(terrain.m_vVertices[(j*terrain.m_iStep) + vIndex]);

		// Add the last two points as per the data structure.
		if (vApplicationCurve.size() != iCoarseUSize)
		{
			vApplicationCurve.push_back(terrain.m_vVertices[terrain.m_iUSize - 2 + vIndex]);
			vApplicationCurve.push_back(terrain.m_vVertices[terrain.m_iUSize - 1 + vIndex]);
		}

		assert(vApplicationCurve.size() == iCoarseUSize);

		// Add an additional point if Dealing with Odd Vertices. Translate it out by Tile Distance for more uniform division
		if (isOdd)
		{	
			vec3 vTranslateVector = vApplicationCurve.back() - *(vApplicationCurve.end() - 2);
			vApplicationCurve.push_back(vApplicationCurve.back() + vec3(vTranslateVector.x, 0.0f, vTranslateVector.z));
		}

		// Apply Algorithm on U-Curve
		applyUReverseSubdivision(vApplicationCurve, vMeshInserter, terrain.m_iStep);

		if (isOdd) // Add the new Vertex in Odd Case
			terrain.m_vVertices.insert(vMeshInserter, F4);
		else
			++vMeshInserter;

		// Reset Curve.
		vApplicationCurve.clear();

		// Index into Row of Mesh
		vIndex += terrain.m_iUSize;
	}
}  

void Terrain::reduceV(tMesh& terrain)
{
	// Local Variables
	vec3 C1, C2, C3, C4;
	vector< vec3 > vApplicationCurve;
	vector< vec3 >::iterator vMeshInserter = terrain.m_vVertices.begin();
	unsigned int iCoarseUSize = terrain.getCoarseUSize();
	bool isOdd = (iCoarseUSize & 1) != 0;
	unsigned int vIndex = 0;

	// Loop for Each U-Curve on Mesh
	for (unsigned int v = 0; v < terrain.m_iUSize; v += (terrain.m_iStep << 1))
	{
		// Loop and store each relevant control point into curve.
		for (unsigned int j = 0; j < (iCoarseUSize >> terrain.m_iExp); ++j)
			vApplicationCurve.push_back(terrain.m_vVertices[(j*terrain.m_iStep) + vIndex]);

		// Add the last two points as per the data structure.
		if (vApplicationCurve.size() != iCoarseUSize)
		{
			vApplicationCurve.push_back(terrain.m_vVertices[terrain.m_iUSize - 2 + vIndex]);
			vApplicationCurve.push_back(terrain.m_vVertices[terrain.m_iUSize - 1 + vIndex]);
		}

		assert(vApplicationCurve.size() == iCoarseUSize);

		// Add an additional point if Dealing with Odd Vertices. Translate it out by Tile Distance for more uniform division
		if (isOdd)
		{
			vec3 vTranslateVector = vApplicationCurve.back() - *(vApplicationCurve.end() - 2);
			vApplicationCurve.push_back(vApplicationCurve.back() + vec3(vTranslateVector.x, 0.0f, vTranslateVector.z));
		}

		applyUReverseSubdivision(vApplicationCurve, vMeshInserter, terrain.m_iStep);

		if (isOdd) // Add the new Vertex in Odd Case
			terrain.m_vVertices.insert(vMeshInserter, C4);
		else
			++vMeshInserter;

		vMeshInserter += terrain.m_iUSize;
		// Reset Curve.
		vApplicationCurve.clear();

		// Index into Row of Mesh
		vIndex += 2 * terrain.m_iUSize;
	}
}

// Given a World x and z position, returns the 4 indices of the Quad the position is situated in
// These indices are in Coarse Space.
void Terrain::get_Quad_Points( float fPosX, float fPosZ, int &iIndex1, int &iIndex2, int &iIndex3, int &iIndex4 )
{
	// Locals/Initialization
	vec3 vOffset;
	iIndex1 = iIndex2 = iIndex3 = iIndex4 = -1;
	int u, v;
	unsigned int iCoarseUSize = m_defaultTerrain.getCoarseUSize(); // Get USize wrt Coarse Space

	// Check that it's within the terrain
	// Get step position
	u = (fPosX > m_defaultTerrain.m_vEndPos.x) ? iCoarseUSize - 2 : -1;
	v = (fPosZ < m_defaultTerrain.m_vStartPos.z) ? 0 : -1;
	u = (fPosX < m_defaultTerrain.m_vStartPos.x) ? 0 : u;
	v = (fPosZ > m_defaultTerrain.m_vEndPos.z) ? iCoarseUSize/*Assume Square-sized Geometry*/ - 2 : v;
	vOffset = vec3(fPosX, 0.0, fPosZ) - vec3(m_defaultTerrain.m_vStartPos.x, 0.0, m_defaultTerrain.m_vStartPos.z);

	// Calculate Within Area of Terrain
	if( -1 == u )
		u = (int)(vOffset.x / (float)m_defaultTerrain.m_fTileWidth);
	if( -1 == v )
		v = (int)(vOffset.z / (float)m_defaultTerrain.m_fTileDepth);

	// Safety check
	if (u < 0)
		u = 0;
	if (v < 0)
		v = 0;
	if (u >= iCoarseUSize - 1)
		u = iCoarseUSize - 2;
	if (v >= iCoarseUSize - 1)
		v = iCoarseUSize - 2;

	/*
		iV1 --------- iV2
		|				|
		|				|
		|				|
		iV3 --------- iV4
	*/
	// Return Indices in Coarse Space
	iIndex1 = (v * iCoarseUSize)+ u;
	iIndex2 = iIndex1 + 1;
	iIndex3 = iIndex1 + iCoarseUSize;
	iIndex4 = iIndex3 + 1;

}

// Stores a Coarse Index for the nearest Control Point to fPosX and fPosZ
void Terrain::get_Point_Pos(float fPosX, float fPosZ)
{
	// Local Variables
	int iCoarseIndex1, iCoarseIndex2, iCoarseIndex3, iCoarseIndex4;

	// Fetches relevant indices in Coarse Space
	get_Quad_Points(fPosX, fPosZ, iCoarseIndex1, iCoarseIndex2, iCoarseIndex3, iCoarseIndex4);

	if (iCoarseIndex1 > -1) // Triangle Found
	{
		// Convert Indices to Index Space from Coarse Space
		unsigned int iIndex1 = m_defaultTerrain.coarseToIndexSpace((unsigned int)iCoarseIndex1);
		unsigned int iIndex2 = m_defaultTerrain.coarseToIndexSpace((unsigned int)iCoarseIndex2);
		unsigned int iIndex3 = m_defaultTerrain.coarseToIndexSpace((unsigned int)iCoarseIndex3);
		unsigned int iIndex4 = m_defaultTerrain.coarseToIndexSpace((unsigned int)iCoarseIndex4);

		// Generate Vertex Coordinates (in xz-plane)
		vec2 vCompVert = vec2(fPosX, fPosZ);
		vec2 vVert1 = vec2(m_defaultTerrain.m_vVertices[iIndex1].x, m_defaultTerrain.m_vVertices[iIndex1].z);
		vec2 vVert2 = vec2(m_defaultTerrain.m_vVertices[iIndex2].x, m_defaultTerrain.m_vVertices[iIndex2].z);
		vec2 vVert3 = vec2(m_defaultTerrain.m_vVertices[iIndex3].x, m_defaultTerrain.m_vVertices[iIndex3].z);
		vec2 vVert4 = vec2(m_defaultTerrain.m_vVertices[iIndex4].x, m_defaultTerrain.m_vVertices[iIndex4].z);

		// Distance from first vertex of Quad to Intersection Position
		float fDist = length(vCompVert - vVert1);
		float fCompDist; // Comparative Distance
		m_iSelector = iCoarseIndex1; // Initialize as closest to first index

		// Find closest Vertex to Position.
		// Vert2 is closer
		if ((fCompDist = length(vCompVert - vVert2)) < fDist)
		{
			fDist = fCompDist;
			m_iSelector = iCoarseIndex2;
		}

		// Vert3 is closer
		if ((fCompDist = length(vCompVert - vVert3)) < fDist)
		{
			fDist = fCompDist;
			m_iSelector = iCoarseIndex3;
		}
		
		// Vert4 is closer
		if ((fCompDist = length(vCompVert - vVert4)) < fDist)
			m_iSelector = iCoarseIndex4;
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

			// Get Full area, with Details
			if (getSelectedArea(iStartU, iStartV, iEndU, iEndV))
			{
				m_vCurrentSubset.m_iUSize = iEndU - iStartU + 1;
				m_vCurrentSubset.m_iVSize = iEndV - iStartV + 1;
				m_vCurrentSubset.m_iExp = m_defaultTerrain.m_iExp;
				m_vCurrentSubset.m_iStep = m_defaultTerrain.m_iStep;

				for (unsigned int v = iStartV; v <= iEndV; ++v)
				{
					for (unsigned int u = iStartU; u <= iEndU; ++u)
					{
						int iIndex = u + (v*m_defaultTerrain.m_iUSize);
						m_vCurrentSubset.m_vVertices.push_back(m_defaultTerrain.m_vVertices[iIndex]);
						m_vCurrentSubset.m_vNormals.push_back(m_defaultTerrain.m_vNormals[iIndex]);
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
	//unsigned int iStartU, iStartV, iEndU, iEndV;
	//getSelectedArea(iStartU, iStartV, iEndU, iEndV);
	//unsigned int iTerrIndex = m_defaultTerrain.m_iUSize * m_defaultTerrain.m_iVSize;
	//unsigned int iAppIndex = m_vApplicationMesh.m_iUSize * m_vApplicationMesh.m_iVSize;
	//unsigned int iAppUIndex = m_vApplicationMesh.m_iUSize;
	//unsigned int iAppVIndex = m_vApplicationMesh.m_iVSize;
	//unsigned int iTerrUIndex = m_defaultTerrain.m_iUSize;
	//unsigned int iTerrVIndex = m_defaultTerrain.m_iVSize;
	//stack<pair<bool, unsigned int>> pMrMapCollector;
	//stack<pair<bool, unsigned int>> pAppMapCollector = m_vApplicationMesh.m_MrMap;
	//pair<bool, unsigned int> pNextMap;
	//
	//blendCoarsePoints(iStartU, iStartV, iEndU, iEndV);
	//
	//int iCount = 0;
	//
	//// Replace existing Details
	//while (!m_defaultTerrain.m_MrMap.empty() && !pAppMapCollector.empty())
	//{
	//	iCount += replaceSquaredDetails(iStartU, iStartV, iEndU, iEndV, iTerrIndex, iTerrUIndex, iTerrVIndex, iAppUIndex, iAppVIndex);
	//
	//	// Something went seriously wrong if the MapStacks aren't even-sized.
	//	assert(0 != ((pAppMapCollector.size() + m_defaultTerrain.m_MrMap.size()) % 2));
	//	for (int i = 0; i < 2; ++i)
	//	{
	//		// Pull Mapping for First Subdivision
	//		pNextMap = m_defaultTerrain.m_MrMap.top();
	//		m_defaultTerrain.m_MrMap.pop();
	//
	//		// Modify Vertex Addition
	//		pNextMap.first ^= pAppMapCollector.top().first;
	//		pAppMapCollector.pop();
	//
	//		// Collect an Ordered Division Map
	//		pMrMapCollector.push(pNextMap);
	//	}
	//
	//}
	//
	//// Completed all pending Subdivision Details for the main mesh, need to apply any remaining Subdivision Calls for the Target Mesh
	//while (!pAppMapCollector.empty())
	//{
	//	m_defaultTerrain.m_vVertices.insert(m_defaultTerrain.m_vVertices.end(), m_defaultTerrain.m_vVertices.size() * 3, vec3(0.0f)); // GOOD!
	//	iCount += replaceSquaredDetails(iStartU, iStartV, iEndU, iEndV, iTerrIndex, iTerrUIndex, iTerrVIndex, iAppUIndex, iAppVIndex);
	//
	//	// Pop Twice for a Squared Subdivision
	//	pMrMapCollector.push(pAppMapCollector.top());
	//	pAppMapCollector.pop();
	//	pMrMapCollector.push(pAppMapCollector.top());
	//	pAppMapCollector.pop();
	//}
	//
	////. Populate the Correct MR Map for the Main Terrain.
	//while (!pMrMapCollector.empty())
	//{
	//	m_defaultTerrain.m_MrMap.push(pMrMapCollector.top());
	//	pMrMapCollector.pop();
	//}
	//
	//cout << "Applied Blend, Touched " << iCount << " Vertices.\n";
	//m_bApplicationMode = false;
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
// Returns a Start and end UV position in Index Space 
// Returns a boolean that states whether the selected space is valid or not.
bool Terrain::getSelectedArea(unsigned int &iStartU, unsigned int &iStartV, unsigned int& iEndU, unsigned int& iEndV, bool bIndexSpace)
{
	// Local Variables
	unsigned int iCoarseUSize = m_defaultTerrain.getCoarseUSize();	// Get U wrt Coarse Points
	int iWorkingSU, iWorkingSV, iWorkingEU, iWorkingEV; // Working Points
	bool bReturnValue = true;

	// Calculate End U and V in Coarse-Point Space
	iWorkingEV = m_iSelector / iCoarseUSize;
	iWorkingEU = m_iSelector % iCoarseUSize;

	if (!m_bApplicationMode)
	{
		// Calculate Start U and V in Coarse-Point Space
		iWorkingSV = m_iLockedStart / iCoarseUSize;
		iWorkingSU = m_iLockedStart % iCoarseUSize;

		// Square off the selection
		bReturnValue = squareArea(iWorkingSU, iWorkingSV, iWorkingEU, iWorkingEV);

		// Order them in correct order.
		if (iWorkingEU < iWorkingSU)
		{
			iWorkingEU = iWorkingEU ^ iWorkingSU;
			iWorkingSU = iWorkingSU ^ iWorkingEU;
			iWorkingEU = iWorkingEU ^ iWorkingSU;
		}
		if (iWorkingEV < iWorkingSV)
		{
			iWorkingEV = iWorkingEV ^ iWorkingSV;
			iWorkingSV = iWorkingSV ^ iWorkingEV;
			iWorkingEV = iWorkingEV ^ iWorkingSV;
		}
	}
	else // Computing Application Mesh Size
	{
		unsigned int iHalfU = m_vApplicationMesh.getCoarseUSize() >> 1;
		unsigned int iHalfV = iHalfU * iCoarseUSize;

		if (iCoarseUSize & 1) // Odd Number of Control Points
		{
			iWorkingSU = iEndU - iHalfU;	// iStartU
			iWorkingSV = iWorkingSU - iHalfV;	// iStartV
			iWorkingEU = iEndU + iHalfU;	// iEndU
			iWorkingEV = iWorkingEU + iHalfV;	// iEndV
		}
		else // Even Number of Control Points
		{
			vec2 vPoint = Mouse_Handler::getInstance(nullptr)->getPosition();
			vec3 vIntersection = GraphicsManager::getInstance(nullptr)->getIntersection(vPoint.x, vPoint.y);
			get_Quad_Points(vIntersection.x, vIntersection.z, iWorkingSU, iWorkingSV, iWorkingEU, iWorkingEV);

			// Break if No valid Quad Found, should be capped to Terrain Limits
			assert(iWorkingSU > -1);

			/*
	  iStartU,iStartV-------------*
			  |					  |
			  |	   1 ------- 2	  |
			  |	   |         |    |
			  |    |         |    |
			  |    3 ------- 4    |
			  |					  |
			  *-------------------iEndU,iEndV
			*/
			iWorkingSU = iWorkingSU - iHalfU;		// iStartU
			iWorkingSV = iWorkingSU - iHalfV;		// iStartV
			iWorkingEU = iWorkingEV + iHalfU;		// iEndU
			iWorkingEV = iWorkingEU + iHalfV;		// iEndV
		}
	}
	
	// Bound the Expanded area to within Default Terrain Boundaries.
	bReturnValue &= boundArea(iWorkingSU, iWorkingSV, iWorkingEU, iWorkingEV);

	// Store adjusted Coarse Indices.
	iStartU = (unsigned int)iWorkingSU;
	iStartV = (unsigned int)iWorkingSV;
	iEndU = (unsigned int)iWorkingEU;
	iEndV = (unsigned int)iWorkingEV;

	// Switch from Coarse Space to Index Space.
	if (bIndexSpace)
	{
		m_defaultTerrain.coarseToIndexSpace(iStartU, iStartV);
		m_defaultTerrain.coarseToIndexSpace(iEndU, iEndV);
	}

	return bReturnValue;
}

// Returns a New Area within iCoarseUSize and starting parameters that is squared to the nearest power of 2 >= 16
bool Terrain::squareArea(int& iStartU, int& iStartV, int& iEndU, int& iEndV)
{
	unsigned int iCoarseUSize = m_defaultTerrain.getCoarseUSize();
	int iWidth = iEndU - iStartU;
	int iDepth = iEndV - iStartV;
	int iSquareCounter = MIN_SELECTION_SIZE;
	int iWidthDirection = (0 < iWidth) - (iWidth < 0);
	int iDepthDirection = (0 < iDepth) - (iDepth < 0);

	if (MIN_SELECTION_SIZE > iCoarseUSize)
		return false;

	// Half the Difference to get an average, use max(that, Minimum size);
	iSquareCounter = abs((iWidthDirection * iWidth) + (iDepthDirection * iDepth)) >> 1;
	iSquareCounter = (MIN_SELECTION_SIZE > iSquareCounter ) ? MIN_SELECTION_SIZE : iSquareCounter;

	iEndU = iStartU + iSquareCounter * iWidthDirection;
	iEndV = iStartV + iSquareCounter * iDepthDirection;

	return true;

	/* Finds the Nearest Power of 2
	unsigned int iNearestP2 = iWidth;
	unsigned int i = 0;
	for (i; iNearestP2 > 1; i++) {
		iNearestP2 = iNearestP2 >> 1;
	}
	// Round to nearest power
	if (iNearestP2 & 1 << i - 1) { i++; }
	iSquareCounter = 1 << i;
	//*/
}

// Bounds the Specified Area into the Default Terrain's area in Coarse Space.
// Returns false if area is larger than Default Terrain's area.
bool Terrain::boundArea(int& iStartU, int& iStartV, int& iEndU, int& iEndV)
{
	/*
  iStartU,iStartV-----------*
		|	 			    |
		|	 1 ------- 2	|
		|	 |         |    |
		|    |         |    |
		|    3 ------- 4    |
		|					|
		*-------------------iEndU,iEndV
	*/
	int iAdjustor;
	unsigned int iCoarseUSize = m_defaultTerrain.getCoarseUSize();

	// Adjust Bounds
	if (iStartU < 0)	// Shift Right
	{
		iAdjustor = iStartU;
		iStartU -= iAdjustor;
		iEndU -= iAdjustor;
	}
	if (iStartV < 0)	// Shift Down
	{
		iAdjustor = iStartV;
		iStartV -= iAdjustor;
		iEndV -= iAdjustor;
	}
	if (iEndU >= (int)iCoarseUSize)	// Shift Left
	{
		iAdjustor = iEndU - iCoarseUSize;
		iEndU -= iAdjustor + 1;
		iStartU -= iAdjustor;
	}
	if (iEndV >= (int)iCoarseUSize) // Shift Up
	{
		iAdjustor = iEndV - iCoarseUSize; 
		iEndV -= iAdjustor + 1;
		iStartV -= iAdjustor;
	}

	return (iStartU >= 0 && iStartV >= 0);
}

void Terrain::initMesh(tMesh& pTerrain)
{
	calculateDimensions(pTerrain);
	generateIndices(pTerrain);
	calculateBarries(pTerrain);
}

// General function to compute Terrain Dimensions
void Terrain::calculateDimensions(tMesh& pTerrain)
{
	// Ensure that Vertices are loaded.
	if (!pTerrain.m_vVertices.empty())
	{
		// Get Start and End Position for Terrain Intersections
		pTerrain.m_vStartPos = pTerrain.m_vVertices.front();	// Will always be a true vertex, by RS algorithm.
		pTerrain.m_vEndPos = pTerrain.m_vVertices.back();		// Will always be a true vertex, by property of n/2 + 1

		// Compute Width and Depth for Terrain Intersections
		pTerrain.m_fWidth = abs(pTerrain.m_vEndPos.x - pTerrain.m_vStartPos.x);
		pTerrain.m_fDepth = abs(pTerrain.m_vEndPos.z - pTerrain.m_vStartPos.z);

		// Compute U and V sizes and Tile Dimensions
		if (pTerrain.m_iUSize == 0 || pTerrain.m_iVSize == 0)
		{
			pTerrain.m_fTileWidth = pTerrain.m_fTileDepth = 0.0f;
			assert(pTerrain.m_iExp == 0); // If subdivision Occurred, this would grab the Coarse U Size with no way to get the Total Usize.
			for (vector<vec3>::const_iterator iter = pTerrain.m_vVertices.begin();
				iter->z == pTerrain.m_vStartPos.z;
				iter++)
			{
				pTerrain.m_iUSize++;

				if (iter != pTerrain.m_vVertices.begin())
					pTerrain.m_fTileWidth += abs(iter->x - (iter - pTerrain.m_iStep)->x);
			}
			pTerrain.m_fTileWidth /= (pTerrain.m_iUSize - 1); // Get average Tile Width

			// Depth/V Direction
			for (vector<vec3>::const_iterator iter = pTerrain.m_vVertices.begin();
				iter != pTerrain.m_vVertices.end() && iter->x == pTerrain.m_vStartPos.x;
				iter += pTerrain.m_iUSize)
			{
				pTerrain.m_iVSize++;

				if (iter != pTerrain.m_vVertices.begin()) // Summate the Total Tile Depth
					pTerrain.m_fTileDepth += abs(iter->z - (iter - pTerrain.m_iUSize)->z);
			}
			pTerrain.m_fTileDepth /= (pTerrain.m_iVSize - 1); // Average Tile Depth
		}
		else
		{
			
			pTerrain.m_fTileWidth = pTerrain.m_vVertices[pTerrain.m_iStep].x - pTerrain.m_vStartPos.x;
			pTerrain.m_fTileDepth = pTerrain.m_vVertices[pTerrain.m_iStep * pTerrain.m_iUSize].z - pTerrain.m_vStartPos.z;
		}
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
		// Compute in Coarse Space, store in Index Space.
		unsigned int iCoarseUSize = pTerrain.getCoarseUSize();
		for (unsigned int v = 0; v < iCoarseUSize - 1; ++v)
		{
			for (unsigned int u = 0; u < iCoarseUSize - 1; ++u)
			{
				unsigned int iA = u + (v*iCoarseUSize);
				unsigned int iC = iA + iCoarseUSize;
				unsigned int iD = pTerrain.coarseToIndexSpace(iC + 1);
				unsigned int iB = pTerrain.coarseToIndexSpace(iA + 1);
				iA = pTerrain.coarseToIndexSpace(iA);
				iC = pTerrain.coarseToIndexSpace(iC);
				

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
		pTerrain.m_vBarries.reserve(pTerrain.m_vIndices.size());
		unsigned int iCoarseUSize = pTerrain.getCoarseUSize();

		assert(pTerrain.m_iUSize == pTerrain.m_iVSize); // Ensure that We're dealing with Square Terrains
		// Loop through Each Point and apply the Barycentric coordinates (in Coarse Space)
		for (unsigned int v = 0; v < iCoarseUSize; ++v)
		{
			for (unsigned int u = 0; u < iCoarseUSize; u += 3)
			{
				pTerrain.m_vBarries.push_back(iArry[iX]);
				if (u + 1 < iCoarseUSize)
					pTerrain.m_vBarries.push_back(iArry[(iX + 1) % 3]);
				if (u + 2 < iCoarseUSize)
					pTerrain.m_vBarries.push_back(iArry[(iX + 2) % 3]);
			}
			iX = (iX + 2) % 3;
		}
	}
}
