#version 430 core

uniform mat4 modelview;
uniform mat4 projection;
uniform vec3 lightPosition;

layout (location = 0) in vec3 vertex;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec2 uv;

//attributes in camera coordinates
out vec3 N;
out vec3 L;
out vec3 V;
out vec3 P;
out vec3 W;

void main(void)
{
	vec4 lightCameraSpace = modelview * vec4(lightPosition, 1.0);
	lightCameraSpace /= lightCameraSpace.w;
	W = vertex;

    vec4 positionCameraSpace = modelview * vec4(vertex, 1.0);
	P = positionCameraSpace.xyz/positionCameraSpace.w;

	// create the Normal Matrix to correct Normal into camera space
	mat3 normalMatrix = transpose(inverse(mat3(modelview)));
	N = normalize( normalMatrix * normal );

	L = normalize(lightCameraSpace.xyz - P);

	V = normalize( - P);

    gl_Position = projection * positionCameraSpace;
}
