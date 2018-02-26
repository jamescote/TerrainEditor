#include "Anim_Track.h"
#include "ShaderManager.h"
#include <sstream>

/***********\
 * DEFINES *
\***********/
#define DEFAULT_DELTA 0.05f
#define SUBDIVISION_COUNT 8
#define DELTA_S 1.f
#define DELTA_ST 0.01f
#define DELTA_T 0.01f
#define MIN_VELOCITY 0.25f
#define SPEED_CONST 5.0f
#define TRACK_WIDTH 0.05f
#define DECEL_BEGIN_ANGLE -5.0f
#define DECEL_THRESHOLD(fD) (fD * 0.9f)
#define GRAVITY 9.81f
#define POSITION_OFFSET 0.025f
#define LEFT_TRACK m_vTrackFrames[0]
#define RIGHT_TRACK m_vTrackFrames[1]

// Default Constructor
Anim_Track::Anim_Track( long lID, const string* sContourFile, const string* pMeshLoc, 
						const string* pTexLoc, bool bOpen )
{
	// Store generation values
	m_lID = lID;
	m_sContourFile = *sContourFile;
	m_sMeshFile = *pMeshLoc;
	m_sTextureFile = *pTexLoc;
	m_bOpenCurve = bOpen;
	
	initializeTrack();
}

// Copy Constructor Overload
Anim_Track::Anim_Track( const Anim_Track& pRHS )
{
	// Utilize operator= overload.
	(*this) = pRHS;
}

// Assignment operator overload.
Anim_Track& Anim_Track::operator=( const Anim_Track& pRHS )
{
	// Copy over initialization parameters
	this->m_lID = pRHS.m_lID;
	this->m_sContourFile = pRHS.m_sContourFile;
	this->m_sMeshFile = pRHS.m_sMeshFile;
	this->m_sTextureFile = pRHS.m_sTextureFile;
	this->m_bOpenCurve = pRHS.m_bOpenCurve;

	// Initialize Track based on same input files.
	initializeTrack();

	return (*this);
}

void Anim_Track::initializeTrack()
{
	loadAnimTrack( m_sContourFile );
	assert( m_vKeyFrames.size() != 0 );
	m_fCurrDist = 0.f;
	m_fCurrHeight = m_vKeyFrames.front().y;
	m_eCurrentState = LIFTING_STATE;

	// Create Vertex Array for Rendering Track
	glGenVertexArrays( 1, &m_iVertexArray );

	// Load a Mesh if mesh is specified.
	if ( !m_sMeshFile.empty() )
	{
		m_pMesh = MeshManager::getInstance()->loadMesh( m_sMeshFile, m_lID );

		if ( nullptr != m_pMesh )
			m_pMesh->initMesh();
	}
	else m_pMesh = nullptr;

	if ( !m_sTextureFile.empty() )
		m_pTexture = TextureManager::getInstance()->loadTexture( m_sTextureFile, m_lID );
	else m_pTexture = nullptr;

	// Generate Vertex buffer for curve.
	m_iVertexBuffer = ShaderManager::getInstance()->genVertexBuffer( m_iVertexArray, 0, 3, m_vKeyFrames.data(), m_vKeyFrames.size() * sizeof( vec3 ), GL_STATIC_DRAW );
}

// Destructor
Anim_Track::~Anim_Track()
{
	if ( nullptr != m_pMesh )
		MeshManager::getInstance()->unloadMesh( m_pMesh->getFileName(), m_lID );

	if ( nullptr != m_pTexture )
		TextureManager::getInstance()->unloadTexture( m_pTexture->getFileName(), m_lID );

	glDeleteBuffers( 1, &m_iVertexBuffer );
	glDeleteVertexArrays( 1, &m_iVertexArray );
}

// loadAnimTrack
// Takes in a set of keyframes loaded from a file.  Temporary function, may switch to take in a file to load the keyframes from.
// Function taken from provided Vec3f_FileIO.h:
// i) Reads in file (e.g. text file (.txt) or contour file (.con)).
// ii) Parses it line by line.
// iii) Extracts the first three floating point values (ignoring all other values)
// iv) Stores these values as a Vec3f in a vector container
void Anim_Track::loadAnimTrack( const string& pContourFile )
{
	using std::string;
	using std::stringstream;
	using std::istream_iterator;

	// Set up input stream
	std::ifstream file( pContourFile );

	// Error opening file.
	if ( !file )
	{
		cout << "Could not open " << pContourFile << " for Animation Track." << endl;
		return;
	}

	// Local Variables.
	vec3 v;
	string line;
	size_t index;
	stringstream ss( std::ios_base::in );
	size_t lineNum = 0;

	// Clear keyframes
	m_vKeyFrames.clear();

	// Get line by line
	while ( getline( file, line ) )
	{
		++lineNum;

		// remove comments	
		index = line.find_first_of( "#" );
		if ( index != string::npos )
		{
			line.erase( index, string::npos );
		}

		// removes leading/tailing junk
		line.erase( 0, line.find_first_not_of( " \t\r\n\v\f" ) );
		index = line.find_last_not_of( " \t\r\n\v\f" ) + 1;
		if ( index != string::npos )
		{
			line.erase( index, string::npos );
		}

		if ( line.empty() )
		{
			continue; // empty or commented out line
		}

		// read line into string stream and clear any flags
		ss.str( line );
		ss.clear();

		// store the position into the vec3f
		ss >> v.x;
		ss >> v.y;
		ss >> v.z;

		// No issues, store the position
		if ( ss.good() )
			m_vKeyFrames.push_back( v );
	}
	file.close();

	// Smooth out curve.
	for ( int i = 0; i < SUBDIVISION_COUNT; ++i )
		smoothCurve();

	// PreProcess for easier position access and evaluation
	preProcessCurve();

	// Calculate Total Curve Length
	m_fCurveLength = ((float)(m_vKeyFrames.size() - 1) * DELTA_ST)
		+ length( m_vKeyFrames.front() - m_vKeyFrames.back() );

	// Generate Tracks for the Car.
	generateTrackFrames();
}

// render
// Renders the Animation track.  Can be turned on or off depending on whether it is necessary to render the track.
// TODO: want to extend to be able to repeat a model along the track such as a rollercoaster track.
void Anim_Track::draw( )
{
	// Get current bindings as a restore point
	GLint iCurrProgBinding = 0, iCurrVABinding = 0;
	glGetIntegerv( GL_VERTEX_ARRAY_BINDING, &iCurrVABinding );
	glGetIntegerv( GL_CURRENT_PROGRAM, &iCurrProgBinding );

	glBindVertexArray( m_iVertexArray );


	// Upload new Indices
#ifdef DEBUG	// For Debugging
	vec3 vColor( 0.5 );
	glPointSize( 10.f );
	glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::WORLD_SHDR ) );
	vector< vec3 > vTemps = { vec3( 0.0f ), getPosition() };
	glBindBuffer( GL_ARRAY_BUFFER, m_iVertexBuffer );
	//glBufferData( GL_ARRAY_BUFFER, vTemps.size() * sizeof( vec3 ), vTemps.data(), GL_DYNAMIC_DRAW );
	//glDrawArrays( GL_POINTS, 0, 2 );
	//glDrawArrays( GL_LINES, 0, 2 );

	vTemps.clear();
	vTemps.push_back( vec3( m_m4CurrentFrenetFrame[ 3 ] ) );
	vTemps.push_back( vec3( m_m4CurrentFrenetFrame[ 3 ] ) + getCentripetalAcce() );
	glBufferData( GL_ARRAY_BUFFER, vTemps.size() * sizeof( vec3 ), vTemps.data(), GL_DYNAMIC_DRAW );
	ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::WORLD_SHDR, "vColor", &vColor );
	//glDrawArrays( GL_POINTS, 0, vTemps.size() );
	//glDrawArrays( GL_LINES, 0, vTemps.size() );

	// Prev Pos
	//vTemps.clear();
	//vTemps.push_back( vec3( m_m4CurrentFrenetFrame[ 3 ] ) );
	//vTemps.push_back( vec3( getPosition( m_fCurrDist - (100.f * DELTA_S) ) ) );
	//glBufferData( GL_ARRAY_BUFFER, vTemps.size() * sizeof( vec3 ), vTemps.data(), GL_DYNAMIC_DRAW );
	//ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::WORLD_SHDR, "vColor", &vec3( 0.25f ) );
	//glDrawArrays( GL_POINTS, 0, vTemps.size() );
	//glDrawArrays( GL_LINES, 0, vTemps.size() );

	// Next Pos
	//vTemps.clear();
	//vTemps.push_back( vec3( m_m4CurrentFrenetFrame[ 3 ] ) );
	//vTemps.push_back( vec3( getPosition( m_fCurrDist + (100.f * DELTA_S) ) ) );
	//glBufferData( GL_ARRAY_BUFFER, vTemps.size() * sizeof( vec3 ), vTemps.data(), GL_DYNAMIC_DRAW );
	//ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::WORLD_SHDR, "vColor", &vec3( 0.75f ) );
	//glDrawArrays( GL_POINTS, 0, vTemps.size() );
	//glDrawArrays( GL_LINES, 0, vTemps.size() );

	mat3 mTemp = mat3( 1.0f );

	for ( int i = 0; i < 3; ++i )
	{
		// Debug FreNet Frame Lines:
		//	Red = BiNormal
		//	Green = Normal
		//	Blue = Tangent
		vTemps.clear();
		vTemps.push_back( vec3( m_m4CurrentFrenetFrame[ 3 ] ) );
		vTemps.push_back( vec3( m_m4CurrentFrenetFrame[ 3 ] + (m_m4CurrentFrenetFrame[ i ] * 0.5) ) );
		glBufferData( GL_ARRAY_BUFFER, vTemps.size() * sizeof( vec3 ), vTemps.data(), GL_DYNAMIC_DRAW );
		ShaderManager::getInstance()->setUniformVec3( ShaderManager::eShaderType::WORLD_SHDR, "vColor", &mTemp[ i ] );
		glDrawArrays( GL_POINTS, 0, vTemps.size() );
		glDrawArrays( GL_LINES, 0, vTemps.size() );
	}
	glPointSize( 1.0f );
#endif

	// Draw Main Curve
	glUseProgram( ShaderManager::getInstance()->getProgram( ShaderManager::eShaderType::RC_TRACK_SHDR ) );
	glBindBuffer( GL_ARRAY_BUFFER, m_iVertexBuffer );
	glBufferData( GL_ARRAY_BUFFER, m_vKeyFrames.size() * sizeof( vec3 ), m_vKeyFrames.data(), GL_DYNAMIC_DRAW );
	glDrawArrays( GL_LINE_LOOP, 0, m_vKeyFrames.size() );

	// Draw Track
	for ( int i = 0; i < 2; ++i )
	{
		glBindBuffer( GL_ARRAY_BUFFER, m_iVertexBuffer );
		glBufferData( GL_ARRAY_BUFFER, m_vTrackFrames[i].size() * sizeof( vec3 ), m_vTrackFrames[ i ].data(), GL_DYNAMIC_DRAW );
		glDrawArrays( GL_LINE_LOOP, 0, m_vTrackFrames[ i ].size() );
	}

	glUseProgram( iCurrProgBinding );
	glBindVertexArray( iCurrVABinding );
}

// animate
//		Called from container object to update the Animation.
//		Updates the Current Distance by Velocity and evaluates any stat changes.
void Anim_Track::animate()
{
	m_fCurrDist = wrap( m_fCurrDist + (getVelocity( m_fCurrDist ) * DELTA_T) );
	m_fCurrHeight = getPosition().y;
	computeFreNetFrames();
	m_eCurrentState = getState( m_fCurrDist );
}

/*******************************************************************************\
* Private Functions  														   *
\*******************************************************************************/

// smoothCurve: pre-processes the curve by splitting the curve into several smaller segments
//		Each of which have, at most, size delta which is evaluated as the smallest distance between
//		two points on the initial set of key frames.
void Anim_Track::smoothCurve()
{
	// Local Variables
	vector< vec3 > vNewCurve, vShiftedCurve;		// New working curve
	vNewCurve.reserve( m_vKeyFrames.size() << 1 );	// Resize to 2n
	vector< vec3 >::const_iterator v_iCurrPoint, v_iNextPoint;	// Point iterators.
	v_iNextPoint = m_vKeyFrames.begin();
	v_iCurrPoint = v_iNextPoint++;
	vNewCurve.push_back( (*v_iCurrPoint) );

	// Steps:
	//	add new points between every 2 points on the curve
	//  move each point halfway forward between each set of points

	// 1st pass
	while ( v_iNextPoint != m_vKeyFrames.end() )
	{
		vNewCurve.push_back( ((*v_iNextPoint) + (*v_iCurrPoint)) * 0.5f );
		vNewCurve.push_back( (*v_iNextPoint) );
		v_iCurrPoint = v_iNextPoint++;
	}

	// Add one last point between the first and last points of the curve (if it's a closed curve)
	if ( !m_bOpenCurve )
		vNewCurve.push_back( (m_vKeyFrames.front() + (*v_iCurrPoint)) * 0.5f );
	else
		vShiftedCurve.push_back( vNewCurve.front() );

	// 2nd pass: shift everything over.
	v_iNextPoint = vNewCurve.begin();
	v_iCurrPoint = v_iNextPoint++;

	while ( v_iNextPoint != vNewCurve.end() )
	{
		vShiftedCurve.push_back( ((*v_iCurrPoint) + (*v_iNextPoint)) * 0.5f );
		v_iCurrPoint = v_iNextPoint++;
	}

	// Shift point over between first and last points on a closed curve.
	if ( !m_bOpenCurve )
		vShiftedCurve.push_back( ((*v_iCurrPoint) + vNewCurve.front()) * 0.5f);

	// Finally, swap new set to stored set
	m_vKeyFrames.swap( vShiftedCurve );
}

// preProcessCurve: creates a new curve with a uniform distance of Delta_S between points.
void Anim_Track::preProcessCurve( )
{
	// Locals
	vector< vec3 > vNewCurve;
	vector< vec3 >::iterator v_iNextPos;
	vec3 vCurrPos;
	v_iNextPos = m_vKeyFrames.begin() + 1;
	vCurrPos = m_vKeyFrames.front();
	m_fMaxHeight = vCurrPos.y;
	m_fMinHeight = vCurrPos.y;
	m_fDistToFreeFall = 0.0f;

	// Push on current point, evaluate from here
	vNewCurve.push_back( vCurrPos );
	
	// Go through each point to set a set distance between each point.
	while ( v_iNextPos != m_vKeyFrames.end() )
	{
		// Distance too small? evaluate next point and try again
		if ( length( *v_iNextPos - vCurrPos ) < DELTA_ST )
			v_iNextPos++;
		else         
		{
			// If the distance is too long, truncate it to the set size and add new point to list
			vNewCurve.push_back( vCurrPos + ((*v_iNextPos - vCurrPos) * DELTA_ST ) );
			vCurrPos = vNewCurve.back();
			
			// Evaluate min and max height as curve is pieced through
			evalHeight( vNewCurve.back().y, vNewCurve.size() * DELTA_ST );
		}
	}

	// Evaluate edge case (between current end and the beginning)
	v_iNextPos = vNewCurve.begin();

	// Stop when the last segment is smaller than the delta.
	while ( length( vNewCurve.front() - vNewCurve.back() ) >= DELTA_ST )
	{
		// Same algorithm
		if ( length( *v_iNextPos - vNewCurve.back() ) < DELTA_ST )
			v_iNextPos++;
		else
		{
			vNewCurve.push_back( vNewCurve.back() + ((*v_iNextPos - vNewCurve.back()) * DELTA_ST) );
			evalHeight( vNewCurve.back().y, vNewCurve.size() * DELTA_ST );
		}
	}

	// Store the new curve.
	m_vKeyFrames.swap( vNewCurve );
}

// After everything has been pre-processed, evaluate the tracks using the binormal at each position on the curve.
void Anim_Track::generateTrackFrames()
{
	// Locals
	vector< vec3 >::iterator v_iCurrPos, v_iNextPos;
	vec3 vBiNormal, vCurrPos;
	LEFT_TRACK.reserve( m_vKeyFrames.size() );		// Reserve necessary space for the tracks
	RIGHT_TRACK.reserve( m_vKeyFrames.size() );		// Reserve necessary space for the tracks
	v_iNextPos = m_vKeyFrames.begin();
	v_iCurrPos = v_iNextPos++;
	float fCurrDist = 0.f;
	float fHalfWidth = TRACK_WIDTH * 0.5f;

	// Track lengths will be the same so only one loop needed
	while ( v_iNextPos != m_vKeyFrames.end() )
	{
		// get binormal and position then create the track points
		vBiNormal = normalize( cross( getTangent( fCurrDist ), computeNormal( fCurrDist ) ) );
		vCurrPos = getPosition(fCurrDist);
		LEFT_TRACK.push_back( vCurrPos + (vBiNormal * -fHalfWidth ) );
		RIGHT_TRACK.push_back( vCurrPos + (vBiNormal * fHalfWidth ) );

		// Iterate through
		fCurrDist += DELTA_ST;
		v_iCurrPos = v_iNextPos++;
	}

	// Last Point
	vBiNormal = normalize( cross( getTangent( fCurrDist ), computeNormal( fCurrDist ) ) );
	vCurrPos = getPosition( fCurrDist );
	LEFT_TRACK.push_back( vCurrPos + (vBiNormal * -fHalfWidth) );
	RIGHT_TRACK.push_back( vCurrPos + (vBiNormal * fHalfWidth) );
}

// Returns an Arc-Length distance at a certain point along the curve.
vec3 Anim_Track::getPosition( float sDist )
{
	// Locals
	vec3 fReturnDist = vec3( 0.f );

	// Evaluate the index from the current distance on the track
	sDist = wrap( sDist );
	float fIndexedDist = (float)sDist / (float)DELTA_ST;
	int iIndex = (unsigned int) (fIndexedDist >= m_vKeyFrames.size() ? 0 : fIndexedDist);
	int iNextIndex = (unsigned int) (iIndex + 1) < m_vKeyFrames.size() ? (iIndex + 1) : 0;

	// Get interpolated distance
	fReturnDist = m_vKeyFrames[ iIndex ] + 
		((m_vKeyFrames[ iNextIndex ] - m_vKeyFrames[ iIndex ]) * (fIndexedDist - (float)iIndex));

	return fReturnDist;
}

// Evaluate the current roller coaster state.
Anim_Track::eCurrVelocityState Anim_Track::getState( float fDist )
{
	vec3 vPosAtDist = getPosition( fDist );
	eCurrVelocityState eReturnState = GRAVITY_FREE_FALL;

	if ( fDist >= DECEL_THRESHOLD( m_fCurveLength ) )	// In the Deceleration zone
		eReturnState = DECELERATION;
	else if ( fDist >= m_fDistToFreeFall )				// In the Free Fall zone
		eReturnState = GRAVITY_FREE_FALL;
	else												// Anything else is the Lifting Stage
		eReturnState = LIFTING_STATE;

	return eReturnState;
}

// Returns the Velocity at a given position.
float Anim_Track::getVelocity( float fDist )
{
	// Locals
	float fReturnVelocity = MIN_VELOCITY;
	float fDecel_Start_Velocity;
	float fDecelDistance;

	// Evaluate Velocity based on current state.
	switch ( getState( fDist ) )
	{
		case GRAVITY_FREE_FALL:	// Free Fall, just use Gravity
			fReturnVelocity = sqrt( (2.0f * GRAVITY) * (m_fMaxHeight - m_fCurrHeight) );
			break;
		case DECELERATION:	// Deceleration Stage -> coming to a stop.
			fDecel_Start_Velocity = 
				sqrt( (2.0f * GRAVITY) * (m_fMaxHeight - getPosition( DECEL_THRESHOLD( m_fCurveLength ) ).y ) );
			fDecelDistance = m_fCurveLength - DECEL_THRESHOLD( m_fCurveLength );
			fReturnVelocity = fDecel_Start_Velocity * ((m_fCurveLength - m_fCurrDist) / fDecelDistance);
			break;
	}

	// Truncate the Velocity to a Minimum.
	fReturnVelocity = fReturnVelocity < MIN_VELOCITY ? MIN_VELOCITY : fReturnVelocity;

	return fReturnVelocity;
}

// Computes FreNet Frames for the Containing Object to position its model with
void Anim_Track::computeFreNetFrames()
{
	vec3 vTangent = getTangent( m_fCurrDist );
	vec3 vBiNormal = computeBiNormal( vTangent );
	vec3 vNormal = normalize( cross( vBiNormal, vTangent ) );

	// mat4 is column wise
	m_m4CurrentFrenetFrame[ 0 ] = vec4( vBiNormal, 0.0f );		// Column 1: Binormal
	m_m4CurrentFrenetFrame[ 1 ] = vec4( vNormal, 0.0f );		// Column 2: Normal
	m_m4CurrentFrenetFrame[ 2 ] = vec4( vTangent, 0.0f );		// Column 3: Tangent
	m_m4CurrentFrenetFrame[ 3 ] = vec4( getPosition() + (vNormal * POSITION_OFFSET), 1.0f );	// Column 4: Translation to Position

}

// Returns the Centripetal Acceleration based on the curvature at the given position on the animation track
vec3 Anim_Track::getCentripetalAcce( float fDist )
{
	// Calculate the Centripetal Acceleration as: **THIS WASN'T WORKING**
	//		(1 / (x^2 + c^2)) * (p'[i+1] - 2p[i] + p'[i-1])
	vec3 vReturnValue;/* = (getPosition( fDist + DELTA_S ) - 2.0f * getPosition( fDist ) + getPosition( fDist - DELTA_S )) / (DELTA_S * DELTA_S);
	float fX = length( vReturnValue ) * 0.5f;												// 1/2 * (p'[i+1] - 2p[i] + p'[i-1])
	float fC = length( getPosition( fDist + DELTA_S ) - getPosition( fDist - DELTA_S ) )	// 1/2 * (p'[i+1] - p'[i-1])
		* 0.5f;
	vReturnValue *= 2.0f*fX / ((fX*fX) + (fC*fC));*/

	vReturnValue = getVelocity( fDist ) * getVelocity( fDist ) *
		(getPosition( fDist + DELTA_S ) - 2.0f * getPosition( fDist ) + getPosition( fDist - DELTA_S )) / (DELTA_S * DELTA_S) +
		vec3( 0.0f, -GRAVITY, 0.0f ) + ((getVelocity( fDist + getVelocity( fDist ) * DELTA_T ) - getVelocity( fDist )) / DELTA_T) *
		getTangent( fDist );  /* Recreation of Code from handout -> This Works, Previous one did not */

	return vReturnValue;
}

// Computes the BiNormal given a Tangent
vec3 Anim_Track::computeBiNormal( const vec3& vTangent )
{
	return normalize( cross( vTangent, computeNormal() ) );
}

// Computes the Normal at a given point on the track.
vec3 Anim_Track::computeNormal( float fDist )
{
	return normalize( getCentripetalAcce( fDist ) - vec3( 0.0f, -GRAVITY, 0.0f ) );
}

// Return estimated tangent at point fDist.
vec3 Anim_Track::getTangent( float fDist )
{
	return normalize( (getPosition( fDist + DELTA_S ) - getPosition( fDist )) / DELTA_S );
}