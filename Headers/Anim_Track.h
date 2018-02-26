#pragma once
#include "stdafx.h"
#include "MeshManager.h"
#include "TextureManager.h"
#include <ctime>

class Anim_Track 
{
public:
	Anim_Track( long lID, const string* sContourFile, 
				const string* sMeshFile = nullptr, const string* sTextureFile = nullptr,  bool bOpen = false );
	Anim_Track( const Anim_Track& pRHS );
	Anim_Track& operator=( const Anim_Track& pRHS );
	~Anim_Track();

	// Public Methods
	void draw( );
	void loadAnimTrack( const string& pContourFile );
	void animate();
	vec3 getPosition() { return getPosition( m_fCurrDist ); }
	mat4 getFreNetFrames() { return m_m4CurrentFrenetFrame; }

private:
	void initializeTrack();
	// Vector of Points on Curve
	vector< vec3 > m_vKeyFrames, m_vTrackFrames[2];
	string m_sContourFile, m_sMeshFile, m_sTextureFile;

	//  Vertex Array and Buffer
	GLuint m_iVertexArray, m_iVertexBuffer;

	// Private Variables
	float m_fCurveLength;
	long m_lID;
	bool m_bOpenCurve;
	float m_fCurrDist, m_fDistToFreeFall;
	float m_fCurrHeight, m_fMaxHeight, m_fMinHeight;
	mat4 m_m4CurrentFrenetFrame;
	vec3 m_vDecelStartPosition;
	enum eCurrVelocityState
	{
		LIFTING_STATE = 0,
		GRAVITY_FREE_FALL,
		DECELERATION
	} m_eCurrentState;
	eCurrVelocityState getState( float fDist );

	// Rendering Members
	Mesh* m_pMesh;
	Texture* m_pTexture;

	// Private Functions
	void smoothCurve();
	void preProcessCurve();
	void generateTrackFrames();
	vec3 getPosition( float sDist );
	vec3 getTangent( float fDist );
	float getVelocity( float fDist );
	float getHorizonAngle() { return dot( getTangent( m_fCurrDist ), vec3( 1.0f, 0.0f, 1.0f ) ); }

	void computeFreNetFrames();
	vec3 computeBiNormal( const vec3& vTangent );
	vec3 computeNormal( float fDist );
	vec3 computeNormal() { return computeNormal( m_fCurrDist ); }
	vec3 getCentripetalAcce( float fDist );
	vec3 getCentripetalAcce() { return getCentripetalAcce( m_fCurrDist ); }

	inline float wrap( float s )
	{
		if ( s > m_fCurveLength )
			s -= m_fCurveLength;
		if ( s < 0 )
			s += m_fCurveLength;
		return s;
	}

	inline void evalHeight( float fHeight, float fDist )
	{
		if ( fHeight > m_fMaxHeight )
		{
			m_fMaxHeight = fHeight;
			m_fDistToFreeFall = fDist;
		}
		m_fMinHeight = (fHeight < m_fMinHeight) ? fHeight : m_fMinHeight;
	}
	
};

