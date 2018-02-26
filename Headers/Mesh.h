#pragma once
#include "stdafx.h"
#include "TriMesh.h"
#include "ShaderManager.h"

//////////////////////////////////////////////////////////////////
// Name: Mesh.h
// Class: Container for TriMesh Objects as well as buffers for normals
//			and indices.
// Written by: James Cote
//////////////////////////////////
class Mesh 
{
public:
	void initMesh( );	// Binds mesh information to a specified Vertex Array.

	// Get Information for drawing elements.
	void bindIndicesBuffer() const { glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, m_iIndicesBuffer ); }
	GLuint getCount() { return m_pMesh->faces.size() * 3; }

	// Gets the file name, only the MeshManager can set this variable.
	const string& getFileName() { return m_sFileName; }

	void drawMesh( );

private:
	// Private Constructor and Copy Constructor to restrict usage to Object_Factory
	Mesh( const string &sFileName );
	Mesh( const Mesh& pCopy );
	Mesh& operator= ( const Mesh& pRHS );
	~Mesh();

	// Private Methods
	bool genMesh( const string& sFileName );

	// Indices for Faces of Mesh and Additional Buffer Addresses on the GPU for
	//	Indices and Normals
	vector<unsigned int> m_pIndices;
	GLuint m_iNormalBuffer, m_iIndicesBuffer, m_iVertexBuffer;
	GLuint m_iVertexArray;
	string m_sFileName;
	ShaderManager* m_pShdrMngr;

	// Mesh Object contains vertices, normals and indices of Mesh
	trimesh::TriMesh* m_pMesh;

	// Friend Class: Object_Factory to create Meshes.
	friend class MeshManager;
};

