/*********************************************************\
 * NAME: HypoCycloid class
 * Written By: James Cote
 * Description: Implements logic for the HypoCycloid algorithm taught in Modelling (CPSC 587)
 *					Modifiable values:
 *						- Big Radius: Radius of the big circle
 *					    - Small Radius: Radius of the Small Circle
 *						- Number of Cycles: Number of Rotations for Small Circle to do about the Big Circle
\*********************************************************/
#include "HypoCycloid.h"
#include "ShaderManager.h"

/***********\
 * DEFINES *
\***********/
#define PI				3.14159265f
#define Z_POS			-1.0f
#define SAMPLE_RATE		(float)(0.1f / PI)
#define SAMPLE_LIMIT	2.0f * PI

/*************\
 * CONSTANTS *
\*************/
const vec3 BIG_CIRCLE_COLOR = vec3( 0.0f, 0.0f, 1.0f );
const vec3 SMALL_CIRCLE_COLOR = vec3( 0.0f, 1.0f, 0.0f );
const vec3 HYPO_COLOR = vec3( 1.0f, 0.0f, 0.0f );

// Default Constructor
HypoCycloid::HypoCycloid()
{
	// Initialize Values.
	glGenVertexArrays( 1, &m_iVertexArray );
	m_fBigR = 3.0f;
	m_fSmallR = 1.0f;
	m_iNumCycles = 1;
	m_fCurrSample = 0.0f;
	m_bAnimating = true;
	m_bPaused = false;

	// Generate Vertex Buffer
	m_iVertexBuffer = ShaderManager::getInstance()->genVertexBuffer( m_iVertexArray, 0, 3, nullptr, 0, GL_DYNAMIC_DRAW );

	// Do initial data calculations
	calculateHypoCycloid();
	calculateBigCircle();
	calculateSmallCircle();
}

// Destructor
HypoCycloid::~HypoCycloid()
{
	// Clean Up buffers.
	glDeleteBuffers( 1, &m_iVertexBuffer );
	glDeleteVertexArrays( 1, &m_iVertexArray );
}

void HypoCycloid::draw( )
{
	if ( m_bAnimating && !m_bPaused )
	{
		m_fCurrSample += SAMPLE_RATE;
		calculateHypoCycloid();
		calculateSmallCircle();
	}

	// Init OpenGl for drawing HypoCycloid
	glBindVertexArray( m_iVertexArray );
	glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::HYPO_SHDR ) );
	vec3 vHypoPos;
	vec3 WHITE = vec3( 1.0f );

	// Draw Big Circle:
	ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::HYPO_SHDR, "vColor", &BIG_CIRCLE_COLOR );
	glBindBuffer( GL_ARRAY_BUFFER, m_iVertexBuffer );
	glBufferData( GL_ARRAY_BUFFER, m_vLargeCircle.size() * sizeof( glm::vec3 ), m_vLargeCircle.data(), GL_DYNAMIC_DRAW );

	glDrawArrays( GL_LINE_LOOP, 0, m_vLargeCircle.size() );

	// Draw Small Circle
	ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::HYPO_SHDR, "vColor", &SMALL_CIRCLE_COLOR );
	glBindBuffer( GL_ARRAY_BUFFER, m_iVertexBuffer );
	glBufferData( GL_ARRAY_BUFFER, m_vSmallCircle.size() * sizeof( glm::vec3 ), m_vSmallCircle.data(), GL_DYNAMIC_DRAW );

	glDrawArrays( GL_LINE_LOOP, 0, m_vSmallCircle.size() );

	// Draw Hypocycloid
	ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::HYPO_SHDR, "vColor", &HYPO_COLOR );
	glBindBuffer( GL_ARRAY_BUFFER, m_iVertexBuffer );
	glBufferData( GL_ARRAY_BUFFER, m_vPattern.size() * sizeof( glm::vec3 ), m_vPattern.data(), GL_DYNAMIC_DRAW );

	glDrawArrays( GL_LINE_STRIP, 0, m_vPattern.size() );

	// Draw Line
	ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::HYPO_SHDR, "vColor", &WHITE );
	glBindBuffer( GL_ARRAY_BUFFER, m_iVertexBuffer );
	glBufferData( GL_ARRAY_BUFFER, m_vSmallLine.size() * sizeof( glm::vec3 ), m_vSmallLine.data(), GL_DYNAMIC_DRAW );

	glDrawArrays( GL_LINES, 0, m_vSmallLine.size() );

	glPointSize( 5.0f );
	glDrawArrays( GL_POINTS, 0, m_vSmallLine.size() );
	glPointSize( 1.0f );

	// Clean Up.
	glUseProgram( 0 );
	glBindVertexArray( 0 );
}

/************************************************************\
 * Modifier functions                                       *
\************************************************************/

void HypoCycloid::resetAnimation()
{
	m_fCurrSample = 0.0f;
	m_vPattern.clear();
	calculateHypoCycloid();
	calculateSmallCircle();
}

// Modifies the smaller radius by some delta.
void HypoCycloid::modifySmallR( float fDelta )
{
	m_fSmallR += fDelta;
	resetAnimation();
	calculateHypoCycloid();
	calculateSmallCircle();
}

// Modifies the larger radius by some delta.
void HypoCycloid::modifyBigR( float fDelta )
{
	// Change Big Circle, then Recalculate Hypocycloid, big circle and small circle.
	m_fBigR += fDelta;
	resetAnimation();
	calculateHypoCycloid();
	calculateBigCircle();
	calculateSmallCircle();
}

// Modifies the number of cycles by some delta integer. Min = 1.
void HypoCycloid::modifyNumCycles( int iDelta )
{
	// Modify Num Circles, min = 1. Then recalculate HypoCycloid.
	m_iNumCycles += iDelta;
	m_iNumCycles = m_iNumCycles < 1 ? 1 : m_iNumCycles;
	resetAnimation();
	calculateHypoCycloid();
}

/*********************************************************************************************\
 * Calculation Functions																	 *
\*********************************************************************************************/

// Function to recalculate the Big circle.
void HypoCycloid::calculateBigCircle()
{
	// Local Variables
	float fCurrSample = 0.0f;

	// Clear Previous Calculations
	m_vLargeCircle.clear();

	// Loop Through and calculate new Circle
	while ( fCurrSample <= SAMPLE_LIMIT )
	{
		// Add New Position on Circle
		m_vLargeCircle.push_back( vec3( m_fBigR * cos( fCurrSample ), m_fBigR * sin( fCurrSample ), Z_POS ) );

		// Move to next sample
		fCurrSample += SAMPLE_RATE;
	}
}

// Function to recalculate the Small Circle, translated to the small circle position.
void HypoCycloid::calculateSmallCircle()
{
	// Local Variables
	mat4 mTranslateMat4 = glm::translate( m_vCurrSmallPos );	// Translation matrix to current Circle Position.
	float fCurrSample = 0.0f;

	// Clear previous Calculations
	m_vSmallCircle.clear();

	// Iterate through and Generate Circle
	while ( fCurrSample <= SAMPLE_LIMIT )
	{
		// Calculate Position of circle and translate to circle position. Save in Vector.
		m_vSmallCircle.push_back( vec3( mTranslateMat4 * vec4( m_fSmallR * cos( fCurrSample ), m_fSmallR * sin( fCurrSample ), Z_POS, 1.0f ) ) );

		// Increment.
		fCurrSample += SAMPLE_RATE;
	}

}

// Function to calculate the HypoCycloid.
void HypoCycloid::calculateHypoCycloid()
{
	// Local Variables
	float fCurrSample = 0.0f;
	vec2 vHypoPos;
	float fBigMinusSmall = m_fBigR - m_fSmallR;
	float fBigMinusSmallOverSmall = fBigMinusSmall / m_fSmallR;

	// Reset Previous Calculations
	m_vSmallLine.clear();

	if ( !m_bAnimating )
	{
		m_vPattern.clear();

		// Iterate through Number of Cycles times and calculate the HypoCycloid.
		while ( fCurrSample <= (SAMPLE_LIMIT * m_iNumCycles) )
		{
			// Cycloid Position at current step.
			float fBMSSBySample = fCurrSample * fBigMinusSmallOverSmall;
			vHypoPos.x = ((fBigMinusSmall) * cos( fCurrSample )) + (m_fSmallR * (cos( fBMSSBySample )));
			vHypoPos.y = ((fBigMinusSmall) * sin( fCurrSample )) - (m_fSmallR * (sin( fBMSSBySample )));

			// Store in Pattern Vector
			m_vPattern.push_back( vec3( vHypoPos, Z_POS ) );

			// Increment step
			fCurrSample += SAMPLE_RATE;
		}

		// Calculate base HypoCycloid position for drawing line in small circle
		vHypoPos.x = ((fBigMinusSmall) * cos( m_fCurrSample )) + (m_fSmallR * (cos( m_fCurrSample * fBigMinusSmallOverSmall )));
		vHypoPos.y = ((fBigMinusSmall) * sin( m_fCurrSample )) - (m_fSmallR * (sin( m_fCurrSample * fBigMinusSmallOverSmall )));

		// Calculate current small circle position (would be updated if animating)
		m_vCurrSmallPos = vec3( fBigMinusSmall * cos( 0.0f ), fBigMinusSmall * sin( 0.0f ), Z_POS );

		// Store line information.
		m_vSmallLine.push_back( vec3( vHypoPos, Z_POS ) );
		m_vSmallLine.push_back( m_vCurrSmallPos );
	}
	else if( !m_bPaused ) // Calculate animating
	{
		if ( m_fCurrSample > (SAMPLE_LIMIT * m_iNumCycles) )
		{
			m_fCurrSample = 0.0f;
			m_vPattern.clear();
		}

		// Cycloid Position at current step.
		float fBMSSBySample = m_fCurrSample * fBigMinusSmallOverSmall;
		vHypoPos.x = ((fBigMinusSmall) * cos( m_fCurrSample )) + (m_fSmallR * (cos( fBMSSBySample )));
		vHypoPos.y = ((fBigMinusSmall) * sin( m_fCurrSample )) - (m_fSmallR * (sin( fBMSSBySample )));

		// Store in Pattern Vector
		m_vPattern.push_back( vec3( vHypoPos, Z_POS ) );

		// Calculate current small circle position (would be updated if animating)
		m_vCurrSmallPos = vec3( fBigMinusSmall * cos( m_fCurrSample ), fBigMinusSmall * sin( m_fCurrSample ), 0 );

		// Store line information.
		m_vSmallLine.push_back( vec3( vHypoPos, Z_POS ) );
		m_vSmallLine.push_back( vec3( m_vCurrSmallPos.x, m_vCurrSmallPos.y, Z_POS ) );
	}
}
