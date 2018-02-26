#include "Mesh.h"

Mesh::Mesh( const string &sFileName )
{
	m_sFileName = sFileName;
}

// Delete any buffers that we initialized
Mesh::~Mesh()
{
	glDeleteBuffers( 1, &m_iNormalBuffer );
	glDeleteBuffers( 1, &m_iIndicesBuffer );
	glDeleteBuffers( 1, &m_iVertexBuffer );
	glDeleteVertexArrays( 1, &m_iVertexArray );

	if (nullptr != m_pMesh)
		delete m_pMesh;
}

// Load the Mesh from a given file name
//  Result: Stores the mesh variables into a set of vertices
bool Mesh::genMesh( const string& sFileName )
{
	m_pMesh = trimesh::TriMesh::read(sFileName);
	bool bReturnValue = (nullptr != m_pMesh);

	if ( bReturnValue )
	{
		m_pMesh->need_bbox();
		m_pMesh->need_faces();
		m_pMesh->need_normals();
		m_pMesh->need_bsphere();

		for ( unsigned int i = 0; i < m_pMesh->faces.size(); i++ )
		{
			m_pIndices.push_back( m_pMesh->faces[ i ][ 0 ] );
			m_pIndices.push_back( m_pMesh->faces[ i ][ 1 ] );
			m_pIndices.push_back( m_pMesh->faces[ i ][ 2 ] );
		}
	}

	return bReturnValue;
}

void Mesh::initMesh()
{
	glGenVertexArrays( 1, &m_iVertexArray );

	m_iVertexBuffer = ShaderManager::getInstance()->genVertexBuffer( 
		m_iVertexArray, 
		0, 3, 
		m_pMesh->vertices.data(), 
		m_pMesh->vertices.size() * sizeof( trimesh::point ), 
		GL_STATIC_DRAW );
	m_iNormalBuffer = ShaderManager::getInstance()->genVertexBuffer( 
		m_iVertexArray, 
		1, 3, 
		m_pMesh->normals.data(), 
		m_pMesh->normals.size() * sizeof( trimesh::vec ), 
		GL_STATIC_DRAW );

	m_iIndicesBuffer = ShaderManager::getInstance()->genIndicesBuffer( 
		m_iVertexArray, 
		m_pIndices.data(), 
		m_pIndices.size() * sizeof( unsigned int ), 
		GL_STATIC_DRAW );

}

// draws the Mesh by setting up the Shader
void Mesh::drawMesh( )
{
	ShaderManager* pShdrMngr = ShaderManager::getInstance();

	glBindVertexArray( m_iVertexArray );
	glUseProgram( pShdrMngr->getProgram( ShaderManager::eShaderType::MESH_SHDR ) );

	glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, m_iIndicesBuffer );
	glDrawElements( GL_TRIANGLES, m_pMesh->faces.size() * 3, GL_UNSIGNED_INT, 0 );

	glUseProgram(0);
	glBindVertexArray( 0 );
}