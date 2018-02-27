#include "NURBS.h"
#include "ShaderManager.h"

/***********\
 * Defines *
\***********/
#define POINT_WIDTH 0.05

// Constructor
NURBS::NURBS( )
{
	m_fUInc = 0.001f;
	m_fU = -1.0f;
	m_iOrder = 2;
	m_iCurrDelta = m_iOrder;
	m_pTargetPoint = nullptr;
	m_bUpdateNeeded = false;
	m_bDrawAffineOnly = true;

	glGenVertexArrays( 1, &m_iVertexArray );

	m_iVertexBuffer = ShaderManager::getInstance()->genVertexBuffer( m_iVertexArray, 0, 3, m_vCurve.data(), sizeof( vec3 ) * m_vCurve.size(), GL_DYNAMIC_DRAW );
}

// Copy Constructor
NURBS::NURBS( NURBS& pCopy )
{
	m_vControlPoints	= pCopy.m_vControlPoints;
	m_vAffineCurve		= pCopy.m_vAffineCurve;
	m_vCurve			= pCopy.m_vCurve;
	m_fU				= pCopy.m_fU;
	m_iOrder			= pCopy.m_iOrder;
	m_pTargetPoint		= nullptr;

	// Generate new OpenGL handles.
	glGenVertexArrays( 1, &m_iVertexArray );
	m_iVertexBuffer = ShaderManager::getInstance()->genVertexBuffer( m_iVertexArray, 0, 3, m_vCurve.data(), sizeof( vec3 ) * m_vCurve.size(), GL_DYNAMIC_DRAW );
}

// Destructor
NURBS::~NURBS()
{
	glDeleteBuffers( 1, &m_iVertexBuffer );
	glDeleteVertexArrays( 1, &m_iVertexArray );
}

void NURBS::drawNURBS()
{
	generateCurve();
	// Init OpenGl for drawing HypoCycloid
	glBindVertexArray( m_iVertexArray );
	glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::HYPO_SHDR ) );
	vec3 WHITE = vec3( 1.0f );
	vec3 SELECTED = vec3( 0.0f, 0.44f, 0.87f );
	vec3 LINE_COLOUR = vec3( 0.58f, 0.51f, 0.79f );
	vec3 CURVE_COLOUR = vec3( 0.78f, 0.61f, 0.43f );
	vec3 GEOMETRIX_COLOUR = vec3( 0.41f, 0.8f, 0.94f );
	vec3 GEOMETRIX_POINT_COLOUR = vec3( 0.67f, 0.83f, 0.45f );
	vec3 NURBMAN_COLOUR = vec3( 0.77f, 0.12f, 0.23f );
	vec3 NON_AFFINE_COLOUR = vec3( 1.0f, 0.49f, 0.04f );

	ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::HYPO_SHDR, "vColor", &SELECTED );

	// Draw Target Point being manipulated.
	if ( nullptr != m_pTargetPoint )
	{
		glBufferData( GL_ARRAY_BUFFER, sizeof( vec3 ), m_pTargetPoint, GL_DYNAMIC_DRAW );

		glPointSize( 10.0f );
		glDrawArrays( GL_POINTS, 0, 1 );
		glPointSize( 1.0f );
	}

	// Bind and set up Control Point Buffer Data.
	glBindBuffer( GL_ARRAY_BUFFER, m_iVertexBuffer );
	glBufferData( GL_ARRAY_BUFFER, m_vControlPoints.size() * sizeof( glm::vec3 ), m_vControlPoints.data(), GL_DYNAMIC_DRAW );

	// Draw Control Points
	ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::HYPO_SHDR, "vColor", &WHITE );
	for ( unsigned int i = 0; i < m_vControlPoints.size(); ++i )
	{
		glPointSize( 5.0f * m_vWeights[i] );
		glDrawArrays( GL_POINTS, i, 1 );
		glPointSize( 1.0f );
	}

	// Draw Convex Hull
	ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::HYPO_SHDR, "vColor", &LINE_COLOUR );
	glDrawArrays( GL_LINE_STRIP, 0, m_vControlPoints.size() );

	// Draw Curve
	ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::HYPO_SHDR, "vColor", &CURVE_COLOUR );
	glBufferData( GL_ARRAY_BUFFER, m_vAffineCurve.size() * sizeof( vec3 ), m_vAffineCurve.data(), GL_DYNAMIC_DRAW );
	glDrawArrays( GL_LINE_STRIP, 0, m_vAffineCurve.size() );

	// Draw Non Affine Curve
	if ( !m_bDrawAffineOnly )
	{
		ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::HYPO_SHDR, "vColor", &NON_AFFINE_COLOUR );
		glBufferData( GL_ARRAY_BUFFER, m_vCurve.size() * sizeof( vec3 ), m_vCurve.data(), GL_DYNAMIC_DRAW );
		glDrawArrays( GL_LINE_STRIP, 0, m_vCurve.size() );
	}

	// Draw NurbMan
	if ( 0.0f <= m_fU && m_fU < 1.0f && m_vControlPoints.size() >= m_iOrder )
	{
		ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::HYPO_SHDR, "vColor", &NURBMAN_COLOUR );
		glBufferData( GL_ARRAY_BUFFER, sizeof( vec3 ), &m_vCurrPos, GL_DYNAMIC_DRAW );
		glPointSize( 10.0f );
		glDrawArrays( GL_POINTS, 0, 1 );

		ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::HYPO_SHDR, "vColor", &GEOMETRIX_COLOUR );
		glBufferData( GL_ARRAY_BUFFER, m_vGeometrix.size() * sizeof( vec3 ), m_vGeometrix.data(), GL_DYNAMIC_DRAW );
		glDrawArrays( GL_LINES, 0, m_vGeometrix.size() );
		ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::HYPO_SHDR, "vColor", &GEOMETRIX_POINT_COLOUR );
		glDrawArrays( GL_POINTS, 0, m_vGeometrix.size() );
		glPointSize( 1.0f );
	}
	else
		clearNurbMan();


	// Clean Up.
	glUseProgram( 0 );
	glBindVertexArray( 0 );
}

void NURBS::generateCurve()
{
	// Local Variables
	unsigned int iDelta = m_iOrder - 1;
	float fCurrU;

	if ( m_bUpdateNeeded )
	{
		m_vAffineCurve.clear();
		m_vCurve.clear();
		clearNurbMan();

		if ( m_vControlPoints.size() >= m_iOrder )
		{
			// Compute a Standard Knot Sequence for Curve
			calculateKnots();
			fCurrU = m_vKnotSequence[ m_iOrder - 1 ];

			// Compute Curve
			while ( fCurrU < m_vKnotSequence[ m_vControlPoints.size() ] )
			{
				while ( fCurrU > m_vKnotSequence[ iDelta + 1 ] )
					iDelta++;

				vec3 vNumerator = getWeightNumer( fCurrU, iDelta );
				float vDenom = getWeightDenom( fCurrU, iDelta );
				m_vAffineCurve.push_back( vNumerator / vDenom );
				m_vCurve.push_back( vNumerator );

				fCurrU += m_fUInc;
			}

			if( -1.0f != m_fU )
				generateNurbMan();
		}

		// Flag that no update is needed to avoid excessive computation
		m_bUpdateNeeded = false;
	}
}

// Function for Calculating the Weighted Denominator
float NURBS::getWeightDenom( float fCurrU, unsigned int iDelta )
{
	vector< float > c;
	for ( unsigned int i = 0; i < m_iOrder; ++i )
		c.push_back( m_vWeights[ iDelta - i ] );
	for ( unsigned int r = m_iOrder; r >= 2; r -= 1 )
	{
		int i = iDelta;
		for ( unsigned int s = 0; s <= r - 2; ++s )
		{
			float fOmega = (fCurrU - m_vKnotSequence[ i ]) / (m_vKnotSequence[ i + r - 1 ] - m_vKnotSequence[ i ]);
			c[ s ] = fOmega* c[ s ] + (1 - fOmega)* c[ s + 1 ];
			i--;
		}
	}

	return c.front();
}

// Function for Calculating the Weighted Numerator
vec3 NURBS::getWeightNumer( float fCurrU, unsigned int iDelta )
{
	vector< vec3 > c;
	for ( unsigned int i = 0; i < m_iOrder; ++i )
		c.push_back( m_vControlPoints[ iDelta - i ] * m_vWeights[ iDelta - i ] );
	for ( unsigned int r = m_iOrder; r >= 2; r -= 1 )
	{
		int i = iDelta;
		for ( unsigned int s = 0; s <= r - 2; ++s )
		{
			float fOmega = (fCurrU - m_vKnotSequence[ i ]) / (m_vKnotSequence[ i + r - 1 ] - m_vKnotSequence[ i ]);
			c[ s ] = fOmega* c[ s ] + (1 - fOmega)* c[ s + 1 ];
			i--;
		}
	}

	return c.front();
}

// Moves targeted point to new position
void NURBS::moveTarget( vec3 vPos )
{
	if ( nullptr != m_pTargetPoint )
	{
		(*m_pTargetPoint) = vPos;
		m_bUpdateNeeded = true;
	}
}
// Advances NurbMan further along the curve.
void NURBS::advanceNurbMan()
{
	if ( m_vControlPoints.size() >= m_iOrder )
	{
		if ( m_vKnotSequence.empty() )
			calculateKnots();

		if ( -1.0f == m_fU )
			m_fU = m_vKnotSequence[ m_iOrder - 1 ];
		else
			m_fU += m_fUInc;

		if ( m_fU >= m_vKnotSequence[ m_vControlPoints.size() ] )
			m_fU = -1.0f;

		if ( -1.0f != m_fU )
			generateNurbMan();
	}
}

// Moves Nurbman back across the curve.
void NURBS::retractNurbMan()
{
	if ( m_vControlPoints.size() >= m_iOrder )
	{
		if ( m_vKnotSequence.empty() )
			calculateKnots();

		if ( -1.0f == m_fU )
			m_fU = m_vKnotSequence[ m_vControlPoints.size() ];
		else
			m_fU -= m_fUInc;

		if ( m_fU <= m_vKnotSequence[ m_iOrder - 1 ] )
			m_fU = -1.0f;

		if ( -1.0f != m_fU )
			generateNurbMan();
	}
}

// Calculates the current point specified internally as the
// position for NurbMan. Also calculates the geometric calculation.
void NURBS::generateNurbMan()
{
	unsigned int iDelta = delta( m_fU );
	vector< vec3 > c;
	vector< float > d;
	m_vGeometrix.clear(); // Clear current values

	if ( iDelta < m_vKnotSequence.size() )
	{
		for ( unsigned int i = 0; i < m_iOrder; ++i )
		{
			c.push_back( m_vControlPoints[ iDelta - i ] * m_vWeights[ iDelta - i ] );
			d.push_back( m_vWeights[ iDelta - i ] );
		}
		for ( unsigned int r = m_iOrder; r >= 2; r -= 1 )
		{
			int i = iDelta;
			for ( unsigned int s = 0; s <= r - 2; ++s )
			{
				float fOmega = (m_fU - m_vKnotSequence[ i ]) / (m_vKnotSequence[ i + r - 1 ] - m_vKnotSequence[ i ]);
				// Store Geometry
				m_vGeometrix.push_back( c[ s ] / d[ s ] );
				m_vGeometrix.push_back( c[ s + 1 ] / d[ s + 1 ] );
				c[ s ] = fOmega* c[ s ] + (1 - fOmega)* c[ s + 1 ];
				d[ s ] = fOmega* d[ s ] + (1 - fOmega)* d[ s + 1 ];
				i--;
			}
		}
		// Store NurbMan's Position.
		m_vCurrPos = c[ 0 ] / d[ 0 ];
	}
}

// Clears NurbMan to stop rendering him
void NURBS::clearNurbMan()
{
	m_vCurrPos = vec3( 0.0f );
	m_vGeometrix.clear();
}

// Selects a nearby control point or adds a new control point at
//	the specified position.
void NURBS::selectAdd( vec3 vPosition )
{
	// Search for an established control point.
	for ( vector< vec3 >::iterator iter = m_vControlPoints.begin();
		 iter != m_vControlPoints.end();
		 ++iter )
	{
		// Grab that control point
		if ( distance( (*iter), vPosition ) <= POINT_WIDTH )
			m_pTargetPoint = &(*iter);
	}

	// if no control point found, add a new one and grab that one.
	if ( nullptr == m_pTargetPoint )
	{
		m_vControlPoints.push_back( vPosition );
		m_vWeights.push_back( 1.0f );
		m_pTargetPoint = &m_vControlPoints.back();	
		m_bUpdateNeeded = true;
	}
}

bool NURBS::modifyWeight( vec3 vPos, float fVal )
{
	unsigned int i = 0;
	// Search for an established control point.
	for ( i; i < m_vControlPoints.size(); ++i )
	{
		// Grab that control point
		if ( distance( m_vControlPoints[ i ], vPos ) <= POINT_WIDTH )
			break;
	}

	if ( i < m_vControlPoints.size() )
	{
		m_vWeights[ i ] += fVal;

		if ( m_vWeights[ i ] < 0.0f )
			m_vWeights[ i ] = 0.0f;

		m_bUpdateNeeded = true;
		return true;
	}

	return false;
}

// Removes the last added Control Point
void NURBS::undoAdd()
{
	if ( !m_vControlPoints.empty() )
	{
		m_vControlPoints.pop_back();
		m_vWeights.pop_back();
		m_bUpdateNeeded = true;
	}
}

/*************************************************************\
 * NURBS Calculation Methods                                 *
\*************************************************************/

// Finds the delta index to use for the NURBS calculation
//	Returns -1 if no delta found.
unsigned int NURBS::delta( float fCurrU )
{
	unsigned int m = m_vControlPoints.size() - 1;

	for ( unsigned int i = 0; i < (m + m_iOrder - 1); ++i )
		if ( fCurrU >= m_vKnotSequence[ i ] && fCurrU < m_vKnotSequence[ i + 1 ] )
			return i;

	return -1;
	
}

// Generates a Uniform Knot Sequence. 
void NURBS::calculateKnots()
{
	// Determine spacing of Knot Values ( 1 / m-1 )
	float fUniformStep = 1.0f / (float) (m_vControlPoints.size() - m_iOrder + 1);

	// Start at 0.0f
	m_vKnotSequence.assign( 1, 0.0f );

	// Place First values in back.
	for ( unsigned int i = 0; i < m_iOrder - 1; ++i )
		m_vKnotSequence.push_back( 0.0f );

	// Incrementally add knots until 1.0f is reached.
	while ( m_vKnotSequence.back() < 1.0f )
		m_vKnotSequence.push_back( m_vKnotSequence.back() + fUniformStep );

	// Place Last Values into sequence.
	for ( unsigned int i = 0; i < m_iOrder - 1; ++i )
		m_vKnotSequence.push_back( 1.0f );
}