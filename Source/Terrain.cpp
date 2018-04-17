#include "ObjLoader.h"
#include "GraphicsManager.h" // #include "Terrain.h"
#include "Mouse_Handler.h"

#define PI					3.14159265f
#define TWOPI				6.283185307f
#define MRLIMIT_MIN			8
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
	objLoader::loadOBJ(pTerrLoc.data(), m_defaultTerrain.m_vVertices, m_defaultTerrain.m_vUVs);

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
		loadMeshData(m_vApplicationMask);
		glBindBuffer(GL_ARRAY_BUFFER, m_iVertexBuffer);
		glBufferData(GL_ARRAY_BUFFER, vTempVerts.size() * sizeof(vec3), vTempVerts.data(), GL_DYNAMIC_DRAW);
		glUseProgram(ShaderManager::getInstance()->getProgram(ShaderManager::eShaderType::TERRAIN_SHDR));
		ShaderManager::getInstance()->setUniformVec3(ShaderManager::eShaderType::TERRAIN_SHDR, "vColor", &GREEN);
		glDrawElements(GL_TRIANGLES, m_vApplicationMask.m_vIndices.size(), GL_UNSIGNED_INT, 0);
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

// Adds zeroed details for further subdivisions, only works for un-reverse subdivided meshes.
void Terrain::addDetails(tMesh& pTerrain, unsigned int iLevelOfDetail)
{
	if (!pTerrain.m_iExp)
	{
		pTerrain.m_iStep = (1 << iLevelOfDetail);
		pTerrain.m_iExp = iLevelOfDetail;

		for (unsigned int i = 0; i < 2; ++i)
		{
			flip(pTerrain);
			// Local Variables
			vector<vec3>::iterator vMeshIterator = pTerrain.m_vVertices.begin() + 1;
			unsigned int iCounter = 0;
			const unsigned int iNumNewDetailsPerRow = pTerrain.m_iUSize - 2;
			unsigned int iNumDetails = ((pTerrain.m_iUSize - 2) << 1) - 1;
			pTerrain.m_iUSize += pTerrain.m_iUSize - 2;
			pTerrain.m_bAddedPointFlagsU.push_back(false);
			unsigned int iInsertCount = pTerrain.m_iStep - 1;
			unsigned int iLastInsertCount = (iInsertCount >> 1) + 1;

			// compute the additional verts that need to be removed at each step in the subdivision.
			for (unsigned int i = 1; i < iLevelOfDetail; ++i)
			{
				//pTerrain.m_iUSize += pTerrain.m_iUSize - 2;
				pTerrain.m_bAddedPointFlagsU.push_back(false);
				if ((iNumDetails != (pTerrain.m_iUSize - 2)))
					pTerrain.m_iUSize++;
				
				pTerrain.m_iUSize += iNumDetails;
				iNumDetails <<= 1;
				//pTerrain.m_bAddedPointFlagsU.push_back(bAddedPoint);
				if (i + 1 == iLevelOfDetail)
					pTerrain.m_iUSize++; // Usize is Num of points in a row, add one since base zero.
			}

			// Add the additional Points
			for (vMeshIterator; vMeshIterator != pTerrain.m_vVertices.end() - 2; vMeshIterator += pTerrain.m_iStep)
			{
				if ((++iCounter) != iNumNewDetailsPerRow)
					vMeshIterator = pTerrain.m_vVertices.insert(vMeshIterator, iInsertCount, vec3(0.0));
				else
				{
					vMeshIterator = pTerrain.m_vVertices.insert(vMeshIterator, iLastInsertCount, vec3(0.0));
					vMeshIterator -= ((int)pTerrain.m_iStep - (int)iLastInsertCount - 3);
					iCounter = 0;
				}
			}
			pTerrain.m_vVertices.insert(pTerrain.m_vVertices.end() - 2, iLastInsertCount, vec3(0.0)); // Insert Last Element @ end - 2

			// Setup Terrain to reflect the new details
			pTerrain.m_bEvenSplitU = false;
		}

		initMesh(pTerrain);
	}
	else
		"Unsure how to add details to a mesh that already has details.";
}

void Terrain::growTerrain()
{
	// Turn on HeightMap before Reverse Subdividing.
	bool bReturn = m_bHeightMapStored;
	if (m_bHeightMapStored)
		toggleHeightMap();

	// Reduce
	grow(m_defaultTerrain);

	if (bReturn)
		toggleHeightMap();

	m_iSelector = 0;
}

void Terrain::grow(tMesh& terrain)
{

	/**********************If Empty*********************/
	if( terrain.m_iStep == 1 )
	{
		addDetails(terrain, 1);
	}
	/***************************************************/

	// V-Columns
	flip(terrain);
	growU(terrain);

	// U-Rows
	flip(terrain);
	if (terrain.m_bAddedPointFlagsV.back())
	{
		terrain.m_vVertices.erase(terrain.m_vVertices.end() - terrain.m_iUSize, terrain.m_vVertices.end());
		terrain.m_iVSize--;
	}
	growU(terrain);
	
	terrain.m_iStep >>= 1;
	terrain.m_iExp--;
	
	if (terrain.m_bAddedPointFlagsU.back())
	{
		flip(terrain);
		terrain.m_vVertices.erase(terrain.m_vVertices.end() - terrain.m_iUSize, terrain.m_vVertices.end());
		terrain.m_iVSize--; // Flip Swaps this to Usize
		flip(terrain);
	}
	terrain.m_bAddedPointFlagsU.pop_back();
	terrain.m_bAddedPointFlagsV.pop_back();

	initMesh(terrain);
}

void Terrain::growU(tMesh& terrain)
{	
	// Local Variables
	unsigned int vIndex = 0;
	vector<vec3> E;
	vector<vec3> C;
	unsigned int fullStep = terrain.m_iStep;
	unsigned int halfStep = terrain.m_iStep >> 1;
	unsigned int iNumDetails = terrain.getCoarseUSize() - 2;

	// Loop Through Every Row and apply Subdivision
	for (unsigned int v = 0; v < terrain.m_iVSize; ++v)
	{
		/*************************E*************************/
		vec3 D1, D2, D3;

		D1 = terrain.m_vVertices.at(vIndex + halfStep);
		D2 = terrain.m_vVertices.at(vIndex + halfStep + (1 * fullStep));
		D3 = terrain.m_vVertices.at(vIndex + halfStep + (2 * fullStep));

		// Apply Starting Computation
		E.clear();
		E.push_back(vec3(0));	// E1 = 0 D1
		E.push_back(HALF * D1);	// E2 = 1/2 D1
		E.push_back((-THREEQUARTER * D1) + (QUARTER * D2));		// E3 = -3/4 D1 + 1/4 D2
		E.push_back((-QUARTER * D1) + (THREEQUARTER * D2));		// E4 = -1/4 D1 + 3/4 D2
		E.push_back((-THREEQUARTER * D2) - (QUARTER * D3));		// E5 = -3/4 D2 - 1/4 D3
		E.push_back((-QUARTER * D2) + (-THREEQUARTER * D3));	// E6 = -1/4 D2 - 3/4 D3

		unsigned int i;
		for (i = 2; i < iNumDetails - 1; ++i)
		{
			vec3 DI, DII;
			DI = terrain.m_vVertices.at(halfStep + (i  *fullStep) + vIndex);
			DII = terrain.m_vVertices.at(halfStep + ((i + 1)*fullStep) + vIndex);
			E.push_back((THREEQUARTER * DI) + (-QUARTER * DII));
			E.push_back((QUARTER * DI) + (-THREEQUARTER * DII));
		}

		// Final E Calculations
		vec3 DS = terrain.m_vVertices.at(vIndex + terrain.m_iUSize - 3);
		E.push_back(HALF * DS);
		E.push_back(vec3(0));	

		/*************************E*************************/

		/***********************Extract C**********************/

		for( unsigned int i = 0; i < iNumDetails; ++i )
			C.push_back(terrain.m_vVertices[(i*fullStep) + vIndex]);

		C.push_back(terrain.m_vVertices[(terrain.m_iUSize - 2) + vIndex]);
		C.push_back(terrain.m_vVertices[(terrain.m_iUSize - 1) + vIndex]);

		/***********************C TO F**********************/
		vec3 CI, CII, EI, EII;

		CI = C[0];
		CII = C[1];

		EI = E[0];
		EII = E[1];

		terrain.m_vVertices[vIndex] = CI + EI;
		terrain.m_vVertices[vIndex + halfStep] = ((HALF*CI) + (HALF*CII) + (EII));

		// Loop through coarse points and apply Efficient Subdivision Algorithm.
		unsigned int j = 2;
		for (i = 1; i < C.size() - 2; i++)
		{
			CI = C[i];
			CII = C[i+1];

			EI = E[j];
			EII = E[j + 1];

			terrain.m_vVertices[vIndex + (j * halfStep)] = ((THREEQUARTER*CI) + (QUARTER*CII) + EI);
			terrain.m_vVertices[vIndex + ((j + 1)* halfStep)] = ((QUARTER*CI) + (THREEQUARTER*CII) + EII);

			j += 2;
		}

		// Grab last two in row and apply subdivision on them.
		CI = *(C.end() - 2);
		CII = *(C.end() - 1);

		// Set them using boundary end condition
		terrain.m_vVertices[(terrain.m_iUSize - 2) + vIndex] = ((HALF * CI) + (HALF * CII) + E.at(j));
		terrain.m_vVertices[(terrain.m_iUSize - 1) + vIndex] = (CII + E.at(j + 1));

		// Increment row index
		vIndex += terrain.m_iUSize;
		C.clear();
	}
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
	unsigned int iCoarseVSize = terrain.getCoarseVSize();

	if( iCoarseUSize >= MRLIMIT_MIN && iCoarseVSize >= MRLIMIT_MIN )
	{
		if (iCoarseUSize & 1)
		{
			flip(terrain);
			terrain.m_vVertices.insert(terrain.m_vVertices.end(), terrain.m_iUSize, vec3(0.0));
			terrain.m_iVSize++; // Flip Swaps this to Usize
			flip(terrain);
			
		}
		// Apply Reverse Subdivision on Us
		reduceU(terrain);

		if (iCoarseVSize & 1)
		{
			terrain.m_vVertices.insert(terrain.m_vVertices.end(), terrain.m_iUSize, vec3(0.0));
			terrain.m_iVSize++;
		}
		// Flip terrain, apply Reverse Subdivision on Vs
		flip(terrain);
		reduceU(terrain);
		flip(terrain);

		// Adjust for New Details
		terrain.m_iStep <<= 1;
		terrain.m_iExp++;
		terrain.m_bAddedPointFlagsU.push_back((iCoarseUSize & 1) != 0);
		terrain.m_bAddedPointFlagsV.push_back((iCoarseVSize & 1) != 0);
		terrain.m_bEvenSplitU = terrain.m_bEvenSplitV = false;

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
			unsigned int iIndex = v*terrain.m_iUSize + u;
			flippedMesh.push_back(terrain.m_vVertices.at(iIndex));
		}
	}

	// Append on Details Information
	terrain.m_vVertices = flippedMesh;

	// Swap Boolean vectors
	vector< bool > bTemp = terrain.m_bAddedPointFlagsU;
	terrain.m_bAddedPointFlagsU = terrain.m_bAddedPointFlagsV;
	terrain.m_bAddedPointFlagsV = bTemp;

	// Swap Even Identifiers
	bool bEvenTemp = terrain.m_bEvenSplitU;
	terrain.m_bEvenSplitU = terrain.m_bEvenSplitV;
	terrain.m_bEvenSplitV = bEvenTemp;

	// Swap U and V sizes
	terrain.m_iUSize = terrain.m_iUSize ^ terrain.m_iVSize;
	terrain.m_iVSize = terrain.m_iVSize ^ terrain.m_iUSize;
	terrain.m_iUSize = terrain.m_iUSize ^ terrain.m_iVSize;
}

// Apply Reverse Subdivision on terrain (in the U-direction)
void Terrain::reduceU(tMesh& terrain)
{
	// Local Variables
	vector< vec3 > vApplicationCurve;
	vector< vec3 >::iterator vMeshInserter = terrain.m_vVertices.begin();
	unsigned int iCoarseUSize = terrain.getCoarseUSize();
	bool isOdd = (iCoarseUSize & 1) != 0;
	//iCoarseUSize += isOdd;
	unsigned int vIndex = 0;
	unsigned int iStepSize = terrain.m_iStep;

	// Loop for Each U-Curve on Mesh
	for (unsigned int v = 0; v < terrain.m_iVSize; v++)
	{
		// Loop and store each relevant control point into curve.
		for (unsigned int j = 0; j < terrain.m_iUSize - 2; j+=terrain.m_iStep )
			vApplicationCurve.push_back(terrain.m_vVertices[j + vIndex]);

		// Add the last point as per the data structure.
		vApplicationCurve.push_back(terrain.m_vVertices[terrain.m_iUSize - 2 + vIndex]);
		vApplicationCurve.push_back(terrain.m_vVertices[terrain.m_iUSize - 1 + vIndex]);

		// Crash if not all the required points were collected.
		//assert(vApplicationCurve.size() == iCoarseUSize);

		// Add an additional point if Dealing with Odd Vertices. Translate it out by Tile Distance for more uniform division
		if (isOdd)
		{	
			vec3 vPushOut = *(vApplicationCurve.end() - 2);
			vec3 vTranslateVector = vPushOut - *(vApplicationCurve.end() - 3);
			vApplicationCurve.back() = vPushOut + vec3(vTranslateVector.x, 0.0f, vTranslateVector.z);
		}

		// First 4 Points
		vec3 F1 = vApplicationCurve.at(0);
		vec3 F2 = vApplicationCurve.at(1);
		vec3 F3 = vApplicationCurve.at(2);
		vec3 F4 = vApplicationCurve.at(3);

		// Apply Boundary Conditions for Vertices
		*vMeshInserter = F1;																		// C1 = F1
		vMeshInserter += iStepSize;
		// Apply Boundary Conditions for Details
		*vMeshInserter = ((-HALF)*F1) + F2 - ((THREEQUARTER)*F3) + ((QUARTER)*F4);			// D1 = -1/2F1 + F2 - 3/4F2 + 1/4F4
		vMeshInserter += iStepSize;

		*vMeshInserter = ((-HALF)*F1) + F2 + ((THREEQUARTER)*F3) - ((QUARTER)*F4);			// C2 = -1/2F1 + F2 + 3/4F3 - 1/4F4
		vMeshInserter += iStepSize;

		F1 = vApplicationCurve.at(2);
		F2 = vApplicationCurve.at(3);
		F3 = vApplicationCurve.at(4);
		F4 = vApplicationCurve.at(5);
		
		*vMeshInserter = ((-QUARTER)*F1) + ((THREEQUARTER)*F2) - ((THREEQUARTER)*F3) + ((QUARTER)*F4); // D2 = -1/4F3 + 3/4F4 - 3/4F5 + 1/4F6
		vMeshInserter += iStepSize;

		// Recursively Compute Information For middle curve
		unsigned int i = 2;
		for (i; i < vApplicationCurve.size() - 5; i += 2)
		{
			F1 = vApplicationCurve.at(i);
			F2 = vApplicationCurve.at(i + 1);
			F3 = vApplicationCurve.at(i + 2);
			F4 = vApplicationCurve.at(i + 3);

			// Details
			if (i >= 4)
			{
				*vMeshInserter = ((QUARTER)*F1) - ((THREEQUARTER)*F2) + ((THREEQUARTER)*F3) - ((QUARTER)*F4); // Dj = 1/4Fi - 3/4Fi+1 + 3/4Fi+2 - 1/4Fi+3
				vMeshInserter += iStepSize;
			}

			// Coarse Point
			*vMeshInserter = ((-QUARTER)*F1) + ((THREEQUARTER)*F2) + ((THREEQUARTER)*F3) - ((QUARTER)*F4); // Cj = -1/4Fi + 3/4Fi+1 + 3/4Fi+2 - 1/4Fi+3
			vMeshInserter += iStepSize;
		}

		// m-3 zero based
		//i = vApplicationCurve.size() - 4;
		// Apply End Boundary Filters
		F1 = vApplicationCurve.at(i);
		F2 = vApplicationCurve.at(i + 1);
		F3 = vApplicationCurve.at(i + 2);
		F4 = vApplicationCurve.at(i + 3);

		vMeshInserter = terrain.m_vVertices.begin() + terrain.m_iUSize + vIndex - 3;
		*vMeshInserter = ((QUARTER)*F1) - ((THREEQUARTER)*F2) + F3 - ((HALF)*F4); // Dj = 1/4Fm-3 - 3/4Fm-2 + Fm-1 - 1/2Fm
		++vMeshInserter;

		//meshV.push_back(((-QUARTER)*C1) + ((THREEQUARTER)*C2) + ((1)*C3) + ((-HALF)*C4));
		*vMeshInserter = ((-QUARTER)*F1) + ((THREEQUARTER)*F2) + F3 - ((HALF)*F4); // Cj = -1/4Fm-3 + 3/4Fm-2 + Fm-1 -1/2Fm
		++vMeshInserter;
		*vMeshInserter = F4;
		++vMeshInserter; // Cj+1 == Fm

		// Reset Curve.
		vApplicationCurve.clear();

		// Index into Row of Mesh
		vIndex += terrain.m_iUSize;
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
	unsigned int iCoarseVSize = m_defaultTerrain.getCoarseVSize(); // Get VSize wrt Coarse Space

	// Check that it's within the terrain
	// Get step position
	u = (fPosX > m_defaultTerrain.m_vEndPos.x) ? iCoarseUSize - 2 : -1;
	v = (fPosZ < m_defaultTerrain.m_vStartPos.z) ? 0 : -1;
	u = (fPosX < m_defaultTerrain.m_vStartPos.x) ? 0 : u;
	v = (fPosZ > m_defaultTerrain.m_vEndPos.z) ? iCoarseVSize - 2 : v;
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
	//if (u >= iCoarseUSize - 1)
	//	u = iCoarseUSize - 2;
	//if (v >= iCoarseUSize - 1)
	//	v = iCoarseUSize - 2;

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
				m_vCurrentSubset.m_bEvenSplitU = (iEndU < m_defaultTerrain.m_iUSize - 1);
				m_vCurrentSubset.m_bEvenSplitV = (iEndV < m_defaultTerrain.m_iVSize - 1);

				for (unsigned int v = iStartV; v <= iEndV; ++v)
				{
					for (unsigned int u = iStartU; u <= iEndU; ++u)
					{
						int iIndex = u + (v*m_defaultTerrain.m_iUSize);
						m_vCurrentSubset.m_vVertices.push_back(m_defaultTerrain.m_vVertices[iIndex]);
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

			unsigned int iCurrUCoarseSize = m_defaultTerrain.getCoarseUSize();
			unsigned int iCurrVCoarseSize = m_defaultTerrain.getCoarseVSize();

			unsigned int iAppUCoarseSize = m_vApplicationMesh.getCoarseUSize();
			unsigned int iAppVCoarseSize = m_vApplicationMesh.getCoarseVSize();

			// Further Reverse Subdivide until able to apply Mesh
			while (iAppUCoarseSize > iCurrUCoarseSize && iAppVCoarseSize > iCurrVCoarseSize)
			{
				reduce(m_vApplicationMesh);
				
				iAppUCoarseSize = m_vApplicationMesh.getCoarseUSize();
				iAppVCoarseSize = m_vApplicationMesh.getCoarseVSize();

				if (iAppUCoarseSize < MRLIMIT_MIN || iAppVCoarseSize < MRLIMIT_MIN)
					break;
			}

			if (iAppUCoarseSize > iCurrUCoarseSize || iAppVCoarseSize > iCurrVCoarseSize)
			{
				cout << "Unable to blend mesh, cannot subdivide far enough.";
				m_vApplicationMesh.clear();
			}
			else
			{
				// Refine Default terrain to default resolution
				//while (m_defaultTerrain.m_iExp != 0)
				//	grow(m_defaultTerrain);

				addDetails(m_defaultTerrain, m_vApplicationMesh.m_iExp);

				// Force Height Map active
				if (m_bHeightMapStored)
					toggleHeightMap();

				// Swap to Application Mode.
				m_bApplicationMode = true;

				m_vCurrentSubset.clear();
				m_iLockedStart = -1;
				m_vSavedSubset.clear();

				setupApplicationMask();
			}
		}
	}
}

void Terrain::setupApplicationMask()
{
	m_vApplicationMask.m_iUSize = m_vApplicationMesh.getCoarseUSize();
	m_vApplicationMask.m_iVSize = m_vApplicationMesh.getCoarseVSize();
	m_vApplicationMask.m_vVertices.resize(m_vApplicationMask.m_iUSize * m_vApplicationMask.m_iVSize, vec3(0.0));
	m_vApplicationMask.m_vVertices[0] = m_vApplicationMesh.m_vVertices[0];
	m_vApplicationMask.m_vVertices.back() = m_vApplicationMask.m_vVertices.back();

	initMesh(m_vApplicationMask);
}

void Terrain::blendMesh()
{
	unsigned int iStartU, iStartV, iEndU, iEndV;
	if (getSelectedArea(iStartU, iStartV, iEndU, iEndV, false))
	{
		iEndU-=2;
		iEndV-=2;

		m_defaultTerrain.coarseToIndexSpace(iStartU, iStartV);
		m_defaultTerrain.coarseToIndexSpace(iEndU, iEndV);
		unsigned int u = iStartU;
		unsigned int v = iStartV;
		unsigned int iAppU, iAppV;

		unsigned int iAppCoarseU = m_vApplicationMesh.getCoarseUSize();
		unsigned int iAppCoarseV = m_vApplicationMesh.getCoarseVSize();

		for (iAppV = 0; iAppV <= (m_vApplicationMesh.m_iVSize - 2 - ((m_vApplicationMesh.m_iStep >> 1) + 1)); ++iAppV)
		{
			for (iAppU = 0; iAppU <= (m_vApplicationMesh.m_iUSize - 2 - ((m_vApplicationMesh.m_iStep >> 1) + 1)); ++iAppU)
			{
				if (iAppU % m_vApplicationMesh.m_iStep && iAppV % m_vApplicationMesh.m_iStep )
					m_defaultTerrain.m_vVertices[u + v*m_defaultTerrain.m_iUSize] = m_vApplicationMesh.m_vVertices[iAppU + (iAppV * m_vApplicationMesh.m_iUSize)];
				u++;
				//iAppU = u - iStartU;
				//iAppV = v - iStartV;
				//if (u % m_defaultTerrain.m_iStep || v % m_defaultTerrain.m_iStep)
				//{
				//	cout << "Replaced: " << m_defaultTerrain.m_vVertices[u + v*m_defaultTerrain.m_iUSize].x << ", " << m_defaultTerrain.m_vVertices[u + v*m_defaultTerrain.m_iUSize].y << ", " << m_defaultTerrain.m_vVertices[u + v*m_defaultTerrain.m_iUSize].z <<
				//		" with: {" << m_vApplicationMesh.m_vVertices[iAppU + (iAppV * m_vApplicationMesh.m_iUSize)].x << ", " << m_vApplicationMesh.m_vVertices[iAppU + (iAppV * m_vApplicationMesh.m_iUSize)].y << ", " << m_vApplicationMesh.m_vVertices[iAppU + (iAppV * m_vApplicationMesh.m_iUSize)].z << "}\n";
				//	
				//	m_defaultTerrain.m_vVertices[u + v*m_defaultTerrain.m_iUSize] = m_vApplicationMesh.m_vVertices[iAppU + (iAppV * m_vApplicationMesh.m_iUSize)];
				//}
				//else
				//{
				//	cout << "skipped: " << m_defaultTerrain.m_vVertices[u + v*m_defaultTerrain.m_iUSize].x << ", " << m_defaultTerrain.m_vVertices[u + v*m_defaultTerrain.m_iUSize].y << ", " << m_defaultTerrain.m_vVertices[u + v*m_defaultTerrain.m_iUSize].z <<
				//		"vs. {"  << m_vApplicationMesh.m_vVertices[iAppU + (iAppV * m_vApplicationMesh.m_iUSize)].x << ", " << m_vApplicationMesh.m_vVertices[iAppU + (iAppV * m_vApplicationMesh.m_iUSize)].y << ", " << m_vApplicationMesh.m_vVertices[iAppU + (iAppV * m_vApplicationMesh.m_iUSize)].z << "}\n";
				//}
			}
			u = iStartU;
			v++;
		}

		m_bApplicationMode = false;
		m_vApplicationMask.clear();
		m_vApplicationMesh.clear();
	}
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

	if (!m_bApplicationMode)
	{
		// Calculate Start U and V in Coarse-Point Space
		iWorkingSV = m_iLockedStart / iCoarseUSize;
		iWorkingSU = m_iLockedStart % iCoarseUSize;

		// Calculate End U and V in Coarse-Point Space
		iWorkingEV = m_iSelector / iCoarseUSize;
		iWorkingEU = m_iSelector % iCoarseUSize;

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
		unsigned int iHalfV = m_vApplicationMesh.getCoarseVSize() >> 1;
		int iIndex1, iIndex2, iIndex3, iIndex4;
		int iQuadSU, iQuadSV, iQuadEU, iQuadEV;

		vec2 vPoint = Mouse_Handler::getInstance(nullptr)->getPosition();
		vec3 vIntersection = GraphicsManager::getInstance(nullptr)->getIntersection(vPoint.x, vPoint.y);
		get_Quad_Points(vIntersection.x, vIntersection.z, iIndex1, iIndex2, iIndex3, iIndex4);
		bool bOddU = (m_vApplicationMesh.getCoarseUSize() & 1);
		bool bOddV = (m_vApplicationMesh.getCoarseVSize() & 1);

		bReturnValue &= (iIndex1 >= 0);

		// Calculate Start U and V in Coarse-Point Space
		iQuadSV = iIndex1 / iCoarseUSize;
		iQuadSU = iIndex1 % iCoarseUSize;

		// Calculate End U and V in Coarse-Point Space
		iQuadEV = iIndex4 / iCoarseUSize;
		iQuadEU = iIndex4 % iCoarseUSize;

		iWorkingSU = (bOddU) ? iQuadSU - iHalfU : iQuadSU - iHalfU + 1;
		iWorkingSV = (bOddV) ? iQuadSV - iHalfV : iQuadSV - iHalfV + 1;
		iWorkingEU = (bOddU) ? iQuadSU + iHalfU : iQuadEU + iHalfU - 1;
		iWorkingEV = (bOddV) ? iQuadSV + iHalfV : iQuadEV + iHalfV - 1;
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
	int iWidth = iEndU - iStartU;
	int iDepth = iEndV - iStartV;
	int iSquareCounter = MIN_SELECTION_SIZE;
	int iWidthDirection = (0 < iWidth) - (iWidth < 0);
	int iDepthDirection = (0 < iDepth) - (iDepth < 0);

	if (MIN_SELECTION_SIZE > m_defaultTerrain.getCoarseUSize() ||
		MIN_SELECTION_SIZE > m_defaultTerrain.getCoarseVSize() )
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
	unsigned int iCoarseVSize = m_defaultTerrain.getCoarseVSize();

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
		iAdjustor = iEndU - iCoarseUSize + 1; 
		iEndU -= iAdjustor;
		iStartU -= iAdjustor;
	}
	if (iEndV >= (int)iCoarseVSize) // Shift Up
	{
		iAdjustor = iEndV - iCoarseVSize + 1; 
		iEndV -= iAdjustor;
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
			pTerrain.m_fTileWidth = pTerrain.m_fWidth / (pTerrain.getCoarseUSize() - 1);
			pTerrain.m_fTileDepth = pTerrain.m_fDepth / (pTerrain.getCoarseVSize() - 1);
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
		unsigned int iCoarseVSize = pTerrain.getCoarseVSize();
		for (unsigned int v = 0; v < iCoarseVSize - 1; ++v)
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
		pTerrain.m_vBarries.resize(pTerrain.m_vVertices.size(), 0);

		// Loop through Each Point and apply the Barycentric coordinates (in Coarse Space)
		int uSize = pTerrain.getCoarseUSize();

		for (unsigned int i = 0; i < pTerrain.m_vIndices.size(); i+=3) // Loop through entire quad, ignoring duplicates
		{

			if (i/6 % (uSize-1) == 0 && i != 0) // 6 Quad size
			{
				iX = ((iX + abs((uSize%3)-3)) % 3);
			}

			iX = (iX + 1) % 3;
			pTerrain.m_vBarries[pTerrain.m_vIndices[i+1]] = iArry[(iX)];

			iX = (iX + 1) % 3;
			pTerrain.m_vBarries[pTerrain.m_vIndices[i+2]] = iArry[(iX)];

			i+=3;

			iX = (iX + 1) % 3;
			pTerrain.m_vBarries[pTerrain.m_vIndices[i+1]] = iArry[(iX)];

			iX = (iX + 1) % 3; 											 
			pTerrain.m_vBarries[pTerrain.m_vIndices[i+2]] = iArry[(iX)];
		}
	}
}
