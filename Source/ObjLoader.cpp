#include "ObjLoader.h"
#include <sstream>

namespace objLoader
{
	bool loadOBJ(const char * path, vector<vec3>& out_vertices, vector<vec2>& out_uvs, vector<vec3>& out_normals)
	{
		// Locals
		ifstream in(path);
		bool bReturnValue = in.is_open();

		// Failed to load file.
		if (!bReturnValue)
			cerr << "Could not Load obj File: " << path << ".\n";
		else
		{
			// Local for Extraction
			string sInput;
			vec3 vTempVec;
			vector<vec3> vPulledNormals;
			out_vertices.clear();
			out_uvs.clear();
			out_normals.clear();
			
			// Read through file.
			while (getline(in, sInput))
			{
				stringstream ssLine(sInput);
				string sToken;

				// Read first token from line
				ssLine >> sToken;

				if ("v" == sToken)	// Vertex Data
				{
					ssLine >> vTempVec.x >> vTempVec.y >> vTempVec.z;
					out_vertices.push_back(vTempVec);
				}
				else if ("vt" == sToken) // uv-coords
				{
					ssLine >> vTempVec.x >> vTempVec.y;
					out_uvs.push_back(vec2(vTempVec));
				}
				else if ("vn" == sToken) // Normals
				{
					ssLine >> vTempVec.x >> vTempVec.y >> vTempVec.z;
					vPulledNormals.push_back(vTempVec);
				}
				else if ("g" == sToken) // Group
				{/* Ignored */
				}
				else if ("o" == sToken) // Object Name
				{/* Ignored */
				}
				else if ("f" == sToken) // Face
				{
					// Local Variables for reading face.
					int iIndices[3] = { -1, -1, -1 };
					vector< int > vFaceVerts;
					if (out_normals.empty())
						out_normals.resize(out_vertices.size());

					// Read through rest of line
					while (ssLine >> sInput)
					{
						// Set up parser and reset previous indices
						stringstream ssItem(sInput);
						iIndices[0] = iIndices[1] = iIndices[2] = -1;

						// Read in face line
						for (int j = 0; getline(ssItem, sInput, '/') && j < 3; ++j)
						{
							stringstream ssIndex(sInput);
							ssIndex >> iIndices[j];
						}

						// Convert to 0-based indices
						iIndices[0] = (-1 == iIndices[0] ? iIndices[0] : iIndices[0] - 1); // Vertex
						iIndices[1] = (-1 == iIndices[1] ? iIndices[1] : iIndices[1] - 1); // UV-Coord; Not Used
						iIndices[2] = (-1 == iIndices[2] ? iIndices[2] : iIndices[2] - 1); // Normal

						// Only Store the Vertex Indices
						vFaceVerts.push_back(iIndices[0]);
						// Store The Normal with the Corresponding Vertex.
						out_normals[iIndices[0]] += vPulledNormals[iIndices[2]];
					}
				}// END Face
			} // END While

			  // Close file when finished.
			in.close();

			for (vector< vec3 >::iterator iter = out_normals.begin();
				iter != out_normals.end();
				++iter)
				(*iter) = normalize(*iter);

		}

		return bReturnValue;
	}
	
} 
