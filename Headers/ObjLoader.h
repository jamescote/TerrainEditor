#pragma once

#include <vector>
#include <glm/glm.hpp>
#include <fstream>
#include <iostream>

using namespace std;
using namespace glm;

namespace objLoader
{
	bool loadOBJ(
		const char * path,
		vector < vec3 > & out_vertices,
		vector < vec2 > & out_uvs
	);
}
