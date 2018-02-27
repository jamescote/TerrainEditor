#pragma once
#include "stdafx.h"

#define U_INC		0.001f
#define U_INC_MIN	0.001f
#define U_INC_MAX	0.05f

// NURBS Class
// Positions a NURBS in the world at some spherical coordinates 
class NURBS
{
public:
	NURBS( );			// Default Constructor
	NURBS( NURBS& pCopy );					// Overloaded Copy Constructor.
	~NURBS();

	// Updating Functions:
	void modifyU( char cDirection ) 
	{ 
		m_fUInc += (float)(cDirection) * U_INC;

		if ( m_fUInc > U_INC_MAX )
			m_fUInc = U_INC_MAX;
		if ( m_fUInc < U_INC_MIN )
			m_fUInc = U_INC_MIN;

		m_bUpdateNeeded = true;
	}
	void modifyOrder( int iDirection ) 
	{ 
		m_iOrder += iDirection;
		if ( m_iOrder < 2 )
			m_iOrder = 2;
		cout << "New Order: " << m_iOrder << endl;
		m_bUpdateNeeded = true;
	}
	void moveTarget( vec3 vPos );
	void selectAdd( vec3 vPosition );
	void releaseSelection() { m_pTargetPoint = nullptr; }
	void undoAdd();
	void advanceNurbMan();
	void retractNurbMan();
	bool modifyWeight( vec3 vPos, float fVal );
	void toggleDrawAffine() { m_bDrawAffineOnly = !m_bDrawAffineOnly; }

	// Evaluation Functions
	void drawNURBS();
	void generateCurve();

private:
	vector< vec3 > m_vControlPoints;
	vector< float > m_vWeights, m_vKnotSequence;
	vector< vec3 > m_vAffineCurve, m_vCurve;
	vector< vec3 > m_vGeometrix;
	vec3 m_vCurrPos;
	vec3* m_pTargetPoint;
	GLuint m_iVertexArray, m_iVertexBuffer;
	float m_fU;
	int m_iCurrDelta;
	float m_fUInc;
	unsigned int m_iOrder;
	bool m_bUpdateNeeded;
	bool m_bDrawAffineOnly;

	// NURBS calculation function
	unsigned int delta( float fCurrU );
	void calculateKnots();
	void generateNurbMan();
	void clearNurbMan();
	vec3 getPointOnCurve( float fCurrU, unsigned int iDelta );
	float getWeightDenom( float fCurrU, unsigned int iDelta );
	vec3 getWeightNumer( float fCurrU, unsigned int iDelta );
};

