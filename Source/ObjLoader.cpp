#include "ObjLoader.h"
#include <sstream>

namespace objLoader
{
	bool loadOBJ(const char * path, vector<vec3>& out_vertices, vector<vec2>& out_uvs)
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
					break;
			} // END While

			  // Close file when finished.
			in.close();
		}

		return bReturnValue;
	}
	
} 
