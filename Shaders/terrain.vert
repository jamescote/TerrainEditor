#version 430 core

uniform mat4 modelview;
uniform mat4 projection;
uniform vec3 lightPosition;


layout (location = 0) in vec3 vertex;
layout (location = 1) in int bc;
layout (location = 2) in vec3 normal;
layout (location = 3) in vec2 uv;

//attributes in camera coordinates
out vec3 N;
out vec3 L;
out vec3 P;
out vec3 vBC;

vec3[] bcVerts = {vec3(1,0,0),
				vec3(0,1,0),
				vec3(0,0,1)};


void main(void)
{
	vec4 lightCameraSpace = modelview * vec4(lightPosition, 1.0);
	lightCameraSpace /= lightCameraSpace.w;

    vec4 positionCameraSpace = modelview * vec4(vertex, 1.0);
	P = positionCameraSpace.xyz/positionCameraSpace.w;

	// create the Normal Matrix to correct Normal into camera space
	mat3 normalMatrix = transpose(inverse(mat3(modelview)));
	N = normalize( normalMatrix * normal );

	L = normalize(lightCameraSpace.xyz - P);
	
	gl_PointSize = gl_Position.w * .5;
    gl_Position = projection * positionCameraSpace;
    vBC=bcVerts[bc];
    
    
    vec3 ndc = gl_Position.xyz / gl_Position.w ; // perspective divide.

	float zDist = 1.0-ndc.z ; // 1 is close (right up in your face,)
	// 0 is far (at the far plane)
	gl_PointSize = 500.0*zDist ; // between 0 and 50 now.
}
