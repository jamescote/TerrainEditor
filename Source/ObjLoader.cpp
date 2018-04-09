#include "ObjLoader.h"

namespace objLoader
{
	bool loadOBJ(const char * path, vector<vec3>& out_vertices, vector<unsigned int>& out_vertIndicies, vector<vec2>& out_uvs, vector<vec3>& out_normals)
	{
		FILE * file = nullptr;
		fopen_s(&file, path, "r");
		if( file == nullptr ){
			printf("Unable to open \"%s\"!\n", path);
			return false;
		}

		while( 1 )
		{
			char lineHeader[128];

			// read the first word of the line
			int res = fscanf_s(file, "%s", lineHeader, 128);

			if (res == EOF)
				break; // EOF = End Of File. Quit the loop
				
			if (strcmp(lineHeader, "v") == 0) {
				vec3 vertex;
				fscanf_s(file, "%f %f %f\n", &vertex.x, &vertex.y, &vertex.z);

				// printf("vertex: (%f, %f, %f)\n",vertex.x, vertex.y, vertex.z);

				out_vertices.push_back(vertex);
			}
			else if (strcmp(lineHeader, "vn") == 0) {
				vec3 normal;
				fscanf_s(file, "%f %f %f\n", &normal.x, &normal.y, &normal.z);

				// printf("normal: (%f, %f, %f)\n",normal.x, normal.y, normal.z);

				out_normals.push_back(normal);
			}
			else if ( strcmp( lineHeader, "vt" ) == 0 )
			{
				vec2 uv;
				fscanf_s(file, "%f %f\n", &uv.x, &uv.y );

				// printf("uv: (%f, %f)\n",uv.x, uv.y);

				out_uvs.push_back(uv);	
			}
			else if ( strcmp( lineHeader, "f" ) == 0 )
			{
				std::string vertex1, vertex2, vertex3;
				unsigned int vertexIndex[3], uvIndex[3], normalIndex[3];
				//int matches = fscanf(file, "%d/%d/%d %d/%d/%d %d/%d/%d\n", &vertexIndex[0], &uvIndex[0], &normalIndex[0], &vertexIndex[1], &uvIndex[1], &normalIndex[1], &vertexIndex[2], &uvIndex[2], &normalIndex[2] );
				int matches = fscanf_s(file, "%d//%d %d//%d %d//%d\n", &vertexIndex[0], &normalIndex[0], &vertexIndex[1], &normalIndex[1], &vertexIndex[2], &normalIndex[2] );
	
				// printf("F: %d//%d %d//%d %d//%d\n matches: %d\n",vertexIndex[0], normalIndex[0], vertexIndex[1], normalIndex[1], vertexIndex[2], normalIndex[2],matches);

				if (matches != 6)
				{
					printf("File can't be read, Try exporting with other options or ask Brad\n");
					return false;
				}
				out_vertIndicies.push_back(vertexIndex[0]-1);
				out_vertIndicies.push_back(vertexIndex[1]-1);
				out_vertIndicies.push_back(vertexIndex[2]-1);
				// uvIndices    .push_back(uvIndex[0]-1);
				// uvIndices    .push_back(uvIndex[1]-1);
				// uvIndices    .push_back(uvIndex[2]-1);
				// normalIndices.push_back(normalIndex[0]-1);
				// normalIndices.push_back(normalIndex[1]-1);
				// normalIndices.push_back(normalIndex[2]-1);
				//cin.get();
			}
		}

		// Close the File once finished
		fclose(file);
		return true;
	}
}
