#version 430 core

// Gooch Variables
uniform float b = 0.1;		// Determine maxmimum "blueness" in temperature shift
uniform float y = 0.0;		// Determine maximmum "yellowness" in temperature shift
uniform float alpha = 0.5;	// Amount of which the Color is visible in the Cool Temperature
uniform float beta = 0.5;	// Amount of which the Color is visible in the Warm Temperature
uniform float rc = 0.5;		// Red Value of the Object's Color
uniform float gc = 0.5;		// Green Value of the Object's Color
uniform float bc = 0.5;		// Blue Value of the Object's Color
uniform vec3 vEdgeColor = vec3(0.0);
uniform vec3 vColor = vec3(0.5);



uniform bool bUsingLinux = true;

out vec4 color;

in vec3 N;
in vec3 L;
in vec3 P;
in vec3 vBC;

uniform sampler2D mySampler;

void main(void)
{
	//vec2 UV;
	vec4 vObjColor = vec4( vColor, 1.0 );
	vec3 kCool;
	vec3 kWarm;
	vec4 textureColor;
	vec3 diffuse = vec3(0.0);
	float specular = 0.0;
	float e = 1.5;

	vec3 R = reflect(-L, N);

	kCool = vec3( 0.0, 0.0, b) + (alpha * vObjColor.rgb);
	kWarm = vec3( y, y, 0.0 ) + (beta * vObjColor.rgb);
	
	float widthScalar = clamp(dot(normalize(-P),N),0.5,1.0);

    if(any(lessThan(vBC, vec3(0.075)*widthScalar)))
    {
		color = vec4(vEdgeColor, 1.0);
	}
	else{
		// Implementing Gooch Shading:
		//		Formula: I = (( 1 - (L.N))/2) * kCool +
		//					  (1 - (1 - (L.N))/2) * kWarm
		float fDotCalc = clamp((1.0 - dot(N,L)) / 2.0, 0.0, 1.0);
		diffuse = (fDotCalc * kCool) + ((1 - fDotCalc) * kWarm);
		specular = pow(max( dot(R,normalize(-P)), 0.0 ), e);

		// Diffuse
		float kd = 0.5;
		//diffuse = kd*vObjColor.rgb*max( 0.0, dot( N, normalize(L - P)));
		//color = vec4( vObjColor.rgb, 1.0 );
		// Specular
		color = vec4(clamp(diffuse + vObjColor.xyz*specular, 0.0, 1.0), 1.0);
		//color = vec4(N,1.0);//vObjColor;
	}

}

