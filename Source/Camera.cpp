#include "Camera.h"

/***********\
* Defines *
\***********/
#define PI					3.14159265f
#define VERT_LOWER_LIMIT	1.f - FLT_EPSILON
#define VERT_UPPER_LIMIT	180.f - VERT_LOWER_LIMIT
#define HORIZONTAL_LIMIT	360.f
#define ZOOM_MIN			0.05f
#define FOV_Y				60.f
#define Z_CLOSE				0.01f
#define Z_FAR				10000.f
#define START_RADIUS		2.5f
#define ORTHO_LR_ZOOM_RATE	2.0f
#define ORTHO_TB_ZOOM_RATE	2.0f
#define ORTHO_DEFAULT_RANGE	7.0f

// Vector Indexing
#define I_THETA				0		// Spherical
#define I_PHI				1
#define I_RADIUS			2
#define I_X					0		// Cartesian
#define I_Y					1
#define I_Z					2

// Shortened Indexing
#define PHI					m_vSPos[I_PHI]
#define THETA				m_vSPos[I_THETA]
#define RADIUS				m_vSPos[I_RADIUS]
#define X_LOOK				m_vSPos[I_X]
#define Y_LOOK				m_vSPos[I_Y]
#define Z_LOOK				m_vSPos[I_Z]
#define ORTHO_LEFT			m_fOrthoCoords[0]
#define ORTHO_RIGHT			m_fOrthoCoords[1]
#define ORTHO_BOTTOM		m_fOrthoCoords[2]
#define ORTHO_TOP			m_fOrthoCoords[3]
#define ORTHO_Z_NEAR		m_fOrthoCoords[4]
#define ORTHO_Z_FAR			m_fOrthoCoords[5]

// Constructor
Camera::Camera( int iHeight, int iWidth )
{
	m_b2DCamera = false;

		//m_vSPos = vec3( 90.f, 8.5308f, START_RADIUS );	// (Theta, Phi, Radius)
	m_vSPos = vec3( 180.f, 90.f, START_RADIUS );	// (Theta, Phi, Radius) facing down Z-Axis.
	m_vWorldLookAt = vec3( 0.f, 0.f, 0.f );		// (X, Y, Z)
	updateHxW( iHeight, iWidth );

	// Default Orthogonal dimensions
	ORTHO_LEFT		= -ORTHO_DEFAULT_RANGE;
	ORTHO_RIGHT		= ORTHO_DEFAULT_RANGE;
	ORTHO_BOTTOM	= -ORTHO_DEFAULT_RANGE;
	ORTHO_TOP		= ORTHO_DEFAULT_RANGE;
	ORTHO_Z_NEAR	= -5.0f;
	ORTHO_Z_FAR		= 5.0f;

	m_bUpdated = true;
}

// Copy Constructor
Camera::Camera( Camera& pCopy )
{
	m_vSPos = pCopy.m_vSPos;
	m_vWorldLookAt = pCopy.m_vWorldLookAt;
	m_bUpdated = pCopy.m_bUpdated;
	m_b2DCamera = pCopy.m_b2DCamera;

	ORTHO_LEFT		= pCopy.ORTHO_LEFT;
	ORTHO_RIGHT		= pCopy.ORTHO_RIGHT;
	ORTHO_BOTTOM	= pCopy.ORTHO_BOTTOM;
	ORTHO_TOP		= pCopy.ORTHO_TOP;
	ORTHO_Z_NEAR	= pCopy.ORTHO_Z_NEAR;
	ORTHO_Z_FAR		= pCopy.ORTHO_Z_FAR;
}

// Destructor
Camera::~Camera()
{
	// Nothing to Destruct
}

// Generates toCamera Matrix and updates Uniform in ShaderManager.
mat4 Camera::getToCameraMat()
{
	vec3 vCamCPos = getCartesianPos();
	return lookAt( vCamCPos, m_vWorldLookAt, vec3( 0.f, 1.f, 0.f ) );
}

// Generates toCamera Matrix and updates Uniform in ShaderManager.
mat4 Camera::getPerspectiveMat()
{
	// vec3 vCamCPos = getCartesianPos();

	if ( !m_b2DCamera )
		return perspective( FOV_Y * PI / 180.f, m_fAspectRatio, Z_CLOSE, Z_FAR );
	else
	{
		float target_Height = (ORTHO_TOP - ORTHO_BOTTOM);
		float target_Width = (ORTHO_RIGHT - ORTHO_LEFT);
		float P = (ORTHO_RIGHT - ORTHO_LEFT) / (ORTHO_TOP - ORTHO_BOTTOM);

		return ortho( m_fAspectRatio * ORTHO_LEFT, m_fAspectRatio * ORTHO_RIGHT,
						m_fAspectRatio * ORTHO_TOP, m_fAspectRatio* ORTHO_BOTTOM,
						ORTHO_Z_NEAR, ORTHO_Z_FAR );

		//return ortho( ORTHO_LEFT, ORTHO_RIGHT, ORTHO_BOTTOM, ORTHO_TOP, ORTHO_Z_NEAR, ORTHO_Z_FAR );
	}
}

// fetches the World Position of the Camera
vec3 Camera::getCameraWorldPos()
{
	return getCartesianPos();
}

// Returns The Look At Vector for the Camera.
const vec3 Camera::getLookAt()
{
	return m_vWorldLookAt - getCartesianPos();
}

// Returns the Current Camera Position in Cartesian Coordinates
vec3 Camera::getCartesianPos()
{
	float fPhi_Rads = PHI * PI / 180.f;
	float fTheta_Rads = THETA * PI / 180.f;
	vec3 vReturn;

	vReturn.z = RADIUS * sin( fPhi_Rads );	// Z = r·sinϕ
	vReturn.x = vReturn.z * sin( fTheta_Rads );		// use Z for X = r·sinϕ·sinθ
	vReturn.x = abs( vReturn.x ) < FLT_EPSILON ? 0.f : vReturn.x;
	vReturn.z *= cos( fTheta_Rads );			// Finish Z: Z = r·sinϕ·cosθ
	vReturn.z = abs( vReturn.z ) < FLT_EPSILON ? 0.f : vReturn.z;
	vReturn.y = RADIUS * cos( fPhi_Rads );	// Y: r·cosϕ
	vReturn.y = abs( vReturn.y ) < FLT_EPSILON ? 0.f : vReturn.y;

	mat4 mLookAtTranslation = translate( mat4( 1.f ), m_vWorldLookAt );
	vec4 mTranslatedPos = mLookAtTranslation * vec4( vReturn, 1.f );
	vReturn = vec3( mTranslatedPos );

	return vReturn;
}

// Handle logic for changing window size.
void Camera::updateHxW( int iHeight, int iWidth )
{
	m_iHeight = iHeight;
	m_iWidth = iWidth;

	m_fAspectRatio = (float)m_iWidth / (float)m_iHeight;
}

/// Camera Manipulation Functions
// Rotatable 360 degrees around.  Laps if it goes over that limit.
void Camera::orbit( vec2 pDelta )
{
	if ( !m_b2DCamera )
	{
		THETA += pDelta.x;
		THETA = THETA > HORIZONTAL_LIMIT ? THETA - HORIZONTAL_LIMIT : THETA;
		THETA = THETA < 0.f ? THETA + HORIZONTAL_LIMIT : THETA;

		PHI += pDelta.y;
		PHI = PHI < VERT_LOWER_LIMIT ? VERT_LOWER_LIMIT : PHI;
		PHI = PHI > VERT_UPPER_LIMIT ? VERT_UPPER_LIMIT : PHI;

		//cout << "CAMERA ORBIT: {" << THETA << ", " << PHI << "}" << endl;
	}
	else // Pan Camera
	{
		m_vWorldLookAt += vec3( -0.05f * pDelta.x, 0.05f * pDelta.y, 0.0f );
		//cout << "LOOK AT: {" << m_vWorldLookAt.x << ", " << m_vWorldLookAt.y << ", " << m_vWorldLookAt.z << "}\n";
	}
}

// Zoom in and out by a given Delta
void Camera::zoom( float fDelta )
{
	if ( !m_b2DCamera )
	{
		RADIUS -= fDelta;
		RADIUS = RADIUS < ZOOM_MIN ? ZOOM_MIN : RADIUS;
		RADIUS = RADIUS > ZOOM_MAX ? ZOOM_MAX : RADIUS;

		//cout << "CAMERA ZOOM: {" << RADIUS << "}" << endl;
	}
	else  // 2D Camera: Shrink orthogonal frustrum.
	{
		// Left/Right
		float fLRZoom = fDelta * ORTHO_LR_ZOOM_RATE;
		ORTHO_LEFT += fLRZoom;
		ORTHO_LEFT = ORTHO_LEFT > 0 ? 0 : ORTHO_LEFT;

		ORTHO_RIGHT -= fLRZoom;
		ORTHO_RIGHT = ORTHO_RIGHT < 0 ? 0 : ORTHO_RIGHT;

		// Top/Bottom
		float fTBZoom = fDelta * ORTHO_TB_ZOOM_RATE;
		ORTHO_TOP -= fTBZoom;
		ORTHO_TOP = ORTHO_TOP < 0 ? 0 : ORTHO_TOP;

		ORTHO_BOTTOM += fTBZoom;
		ORTHO_BOTTOM = ORTHO_BOTTOM > 0 ? 0 : ORTHO_BOTTOM;

		/*
		cout << "CAMERA ZOOM: { L=" << ORTHO_LEFT << "; R=" << ORTHO_RIGHT
			   << "; T=" << ORTHO_TOP << "; B=" << ORTHO_BOTTOM
			   << "; Near=" << ORTHO_Z_NEAR << "; Far=" << ORTHO_Z_FAR << " }\n" ;//*/
	}
}

// Ray Casting:
//	Tutorial from: antongerdelan.net/opengl/raycasting.html
vec3 Camera::getRay( float fX, float fY )
{
	// Convert to 3D Normalized Device Coordinates
	float fMouseX = (2.0f * fX / (float)m_iWidth) - 1.0f;
	float fMouseY = (-2.0f * fY / (float)m_iHeight) + 1.0f;

	// Create 4D Homogeneous Clip Coordinates
	vec4 ray_clip = vec4( fMouseX, fMouseY, -1.0f, 1.0f );

	// 4D Eye (Camera) Coordinates - Inverse Projection Matrix
	vec4 ray_eye = inverse( getPerspectiveMat() ) * ray_clip;

	// Un-project the x,y part -> manually set z,w
	ray_eye = vec4( ray_eye.x, ray_eye.y, -1.0, 0.0 );

	// 4D World Coordinates
	return normalize( vec3( inverse( getToCameraMat() ) * ray_eye ) );
}
