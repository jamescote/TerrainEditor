#pragma once

#include "Object3D.h"
#include <stack>

class Terrain
{
public:
	Terrain(const string& pTerrLoc);
	~Terrain();

	// Overridden intersect function
	void draw( );

	void get_Point_Pos( float fPosX, float fPosZ );
	void toggleHeightMap();
	void lockPoint();
	void saveSelection();
	void clearSelection();
	void applyTerrain(const Terrain* vMesh);
	void growTerrain();
	void reduceTerrain();

private:
	Terrain( const Terrain* pNewPlane );  // Protected Copy Constructor

	// Mesh Structure for Storing unique attributes to a Terrain Mesh
	struct tMesh
	{
		// Vertex Data
		vector< vec3 > m_vVertices, m_vNormals;
		vector< unsigned int > m_vIndices, m_vBarries;
		vec3 m_vStartPos, m_vEndPos;
		vector< vec2 > m_vUVs;

		// Dimensions of Mesh
		unsigned int m_iUSize, m_iVSize;
		float m_fWidth, m_fDepth, m_fTileWidth, m_fTileDepth;	
		unsigned int m_iStep;
		unsigned int m_iExp;
		bool m_bEvenSplitU, m_bEvenSplitV;
		vector< bool > m_bAddedPointFlagsU;
		vector< bool > m_bAddedPointFlagsV;

		tMesh()
		{
			clear();
		}

		void clear()
		{
			m_vVertices.clear();
			m_vIndices.clear();
			m_vBarries.clear();
			m_vUVs.clear();
			m_vNormals.clear();
			m_vEndPos = m_vStartPos = vec3(0.0f);
			m_fWidth = m_fDepth = -1.0f;
			m_iUSize = m_iVSize = 0;
			m_fTileDepth = m_fTileWidth = -1.0f;
			m_iStep = 1;
			m_iExp = 0;
			m_bEvenSplitU = m_bEvenSplitV = false;
		}

		unsigned int getCoarseUSize()
		{

			// Adjust for Reverse Subdivision. Total size is divided by the number of steps between points.
			//								   There is the possibility of an additional point added (odd-> add 1)
			//								   If there was some reverse subdivision, there's an additional point added at the end.
			return (m_iUSize / m_iStep) + didAddPointU() + (m_iExp > 0 ? 1 : 0);
		}

		unsigned int getCoarseVSize()
		{

			// Adjust for Reverse Subdivision. Total size is divided by the number of steps between points.
			//								   There is the possibility of an additional point added (odd-> add 1)
			//								   If there was some reverse subdivision, there's an additional point added at the end.
			return (m_iVSize / m_iStep) + didAddPointV() + (m_iExp > 0 ? 1 : 0);
		}


		// Converts a given pair of Coarse U and V coordinates into Index space within the Vertex array.
		void coarseToIndexSpace(unsigned int& iCoarseU, unsigned int& iCoarseV)
		{
			// If No Reverse Subdivision has taken place, this function is onto and one to one
			if (m_iExp)
			{
				// Else, get the Coarse Size and set into an array for processing
				unsigned int iCoarseArry[] = { getCoarseUSize(), getCoarseVSize() };
				unsigned int* iArry[] = { &iCoarseU, &iCoarseV };
				bool bArry[] = {m_bEvenSplitU, m_bEvenSplitV};
				unsigned int iSizeArry[] = { m_iUSize, m_iVSize };

				// for each one, apply the same algorithm
				for (int i = 0; i < 2; ++i)
				{
					if (*iArry[i] == iCoarseArry[i] - 1) // The Last Point is always the same Size - 1
						*iArry[i] = iSizeArry[i] - 1;
					else if (!bArry[i] && (*iArry[i] == iCoarseArry[i] - 2) )
						*iArry[i] = iSizeArry[i] - 2; // If the Full size is odd, there's an added point. This fullfills the same condition as the first condition for this case.
					else // Otherwise, The resulting index is a multiple of the stepsize.
						*iArry[i] *= m_iStep;

				}
			}
		}

		// Overloaded Function, converts an index from Coarse Space to Index Space
		unsigned int coarseToIndexSpace(unsigned int iCoarseIndex) 
		{
			if (!m_iExp)
				return iCoarseIndex; // Nothing to compute if every point is a control point.

			unsigned int iCoarseUSize = getCoarseUSize();
			unsigned int iU = iCoarseIndex % iCoarseUSize;
			unsigned int iV = iCoarseIndex / iCoarseUSize;
			coarseToIndexSpace(iU, iV); // Convert U and V to Index Space
			// compute the new Index
			return iU + (iV * m_iUSize);
		}

		void outputMesh()
		{
			cout << "Mesh Details:\n\t";
			cout << "Num Verts: " << m_vVertices.size() << "; Indices: " << m_vIndices.size()
				<< "; Barries: " << m_vBarries.size() << "; UVs: " << m_vUVs.size() << "; Normals: " << m_vNormals.size() << endl;
			cout << "\tStartPos: {" << m_vStartPos.x << ", " << m_vStartPos.y << ", " << m_vStartPos.z << "}; EndPos: {"
				<< m_vEndPos.x << ", " << m_vEndPos.y << ", " << m_vEndPos.z << "}\n\t";
			cout << "Width: " << m_fWidth << "; Depth: " << m_fDepth << endl;
			cout << "\tUSize: " << m_iUSize << "; VSize: " << m_iVSize << endl;
			cout << "\tTile Width: " << m_fTileWidth << "; Tile Depth: " << m_fTileDepth << endl;

		}

		bool didAddPointU()
		{
			bool bAddedPoint = false;
			for (vector< bool >::const_iterator iter = m_bAddedPointFlagsU.begin();
				iter != m_bAddedPointFlagsU.end();
				++iter)
				bAddedPoint ^= (*iter);
			return bAddedPoint;
		}
		bool didAddPointV()
		{
			bool bAddedPoint = false;
			for (vector< bool >::const_iterator iter = m_bAddedPointFlagsV.begin();
				iter != m_bAddedPointFlagsV.end();
				++iter)
				bAddedPoint ^= (*iter);
			return bAddedPoint;
		}
	};

	// Main TerrainMesh
	tMesh m_defaultTerrain;
	tMesh m_vSavedSubset, m_vCurrentSubset, m_vApplicationMesh;
	
	//vector< vec3 > m_vSavedSubset, m_vCurrentSubset;

	vec3 m_vSelector, m_vLockedStart;
	unsigned int m_iSelector;
	int m_iLockedStart;

	// Height Map for Main Terrain, helps for selections
	vector< float > m_vHeightMap;
	bool m_bHeightMapStored;
	bool m_bApplicationMode;

	// Calculation functions
	void initMesh(tMesh& pTerrain);
	void calculateDimensions( tMesh& pTerrain );
	void generateIndices( tMesh& pTerrain );
	void calculateBarries(tMesh& pTerrain);
	bool getSelectedArea(unsigned int &iStartU, unsigned int &iStartV, unsigned int& iEndU, unsigned int& iEndV, bool bIndexSpace = true);
	bool squareArea(int& iStartU, int& iStartV, int& iEndU, int& iEndV);
	bool boundArea(int& iStartU, int& iStartV, int& iEndU, int& iEndV);
	void loadMeshData(const tMesh& pTerrain);
	void blendMesh();
	int blendCoarsePoints(unsigned int iStartU, unsigned int iStartV, unsigned int iEndU, unsigned int iEndV);
	int replaceDetails(unsigned int iStartU, unsigned int iStartV, unsigned int iEndU, unsigned int iEndV, unsigned int iTerrRowSize, unsigned int iAppRowSize );
	int replaceSquaredDetails(unsigned int iStartU, unsigned int iStartV, unsigned int iEndU, unsigned int iEndV, unsigned int& iTerrIndex,
								unsigned int& iTerrUSize, unsigned int& iTerrVSize, unsigned int& iAppUSize, unsigned int& iAppVSize);

	GLuint m_iVertexArray, m_iBaryCentric, m_iVertexBuffer, m_iNormalBuffer, m_iTextureBuffer, m_iIndicesBuffer;

	void get_Quad_Points( float fPosX, float fPosZ, int &iIndex1, int &iIndex2, int &iIndex3, int &iIndex4 );
	void flip(tMesh&);
	void reduce(tMesh&);
	void reduceU(tMesh&);
	void grow(tMesh&);	
	void growU(tMesh&);
	void addDetails(tMesh& m_defaultTerrain, unsigned int iLevelOfDetail);


};
