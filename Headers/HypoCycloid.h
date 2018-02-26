#pragma once
#include "Object.h"
class HypoCycloid
{
public:
	HypoCycloid();
	~HypoCycloid();

	void draw( );
	string getType() { return "Hypocycloid"; }
	string getDebugOutput() { return "No Debug Output yet.\n"; }

	// Manipulation functions
	void modifySmallR( float fDelta );
	void modifyBigR( float fDelta );
	void modifyNumCycles( int iDelta );
	void cycleAnimationModes()
	{
		if( m_bAnimating )
			m_bPaused = !m_bPaused;
	}
	void toggleAnimation()
	{
		m_bAnimating = !m_bAnimating;
		m_bPaused = false;
		resetAnimation();
	}

private:
	float m_fSmallR, m_fBigR;
	int m_iNumCycles;
	vec3 m_vCurrSmallPos;
	float m_fCurrSample;
	GLuint m_iVertexArray, m_iVertexBuffer;

	// Animation Functionality
	bool m_bAnimating;
	bool m_bPaused;

	vector<vec3> m_vPattern, m_vLargeCircle, m_vSmallCircle, m_vSmallLine;

	// Calculations functions.
	void calculateBigCircle();
	void calculateSmallCircle();
	void calculateHypoCycloid();
	void resetAnimation();
};

