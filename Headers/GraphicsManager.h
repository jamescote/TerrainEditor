#pragma once

/* INCLUDES */
#include "stdafx.h"
#include "Camera.h"
#include "HypoCycloid.h"
#include "NURBS.h"
#include "Terrain.h"

// Forward Declarations
class ShaderManager;
class EnvironmentManager;

// Class: Graphics Manager
// Purpose: Acts as the Sinew between all moving parts that are required for drawing
//			with openGL.
// TODO: Set-up a Manager for Geometry (Geometry may be expanded on later).
// Written by: James Cot�
class GraphicsManager
{
public:
	static GraphicsManager* getInstance(GLFWwindow *rWindow);
	~GraphicsManager();

	// Graphics Application
	bool initializeGraphics( string sFileName );
	bool renderGraphics();

	/// HxW Settings
	void resizedWindow( int iHeight, int iWidth ) { m_pCamera->updateHxW( iHeight, iWidth ); }

	// Camera Functions
	//void snapShotCamera();
	void rotateCamera(vec2 pDelta);
	void zoomCamera(float fDelta);
	void moveCamera( vec3 pDirection );

	// Helper Enum for RGB Values
	enum eRGB
	{
		RED = 0,
		GREEN,
		BLUE,
		RGB_MAX
	};

	// Shader Manipulation functions
	void setRGBVal( eRGB eType, float fVal);
	void setBeta(float fVal);
	void setAlpha(float fVal);
	void setBVal(float fVal);
	void setYVal(float fVal);
	void togGooch();
	void togToon();
	void togSpec();
	void setShine( float fVal );
	void setR( float fVal );

	// HypoCycloid Manipulation Functions
	void modifyHypoBigR( float fDelta ) { m_pHypoCyc->modifyBigR( fDelta ); }
	void modifyHypoSmallR( float fDelta ) { m_pHypoCyc->modifySmallR( fDelta ); }
	void modifyHypoN( int iDelta ) { m_pHypoCyc->modifyNumCycles( iDelta ); }
	void modifyHypoAnim() { m_pHypoCyc->cycleAnimationModes(); }
	void toggleHypoAnim() { m_pHypoCyc->toggleAnimation(); }

	// NURBS Manipulation Functions
	void controlPoint( float fX, float fY );
	void dropControlPoint() { m_pNurbs->releaseSelection(); }
	void moveControlPoint( float fX, float fY );
	void undoControlPoint() { m_pNurbs->undoAdd(); }
	void modifyOrder( int iDirection ) { m_pNurbs->modifyOrder( iDirection ); }
	void advanceNurbMan() { m_pNurbs->advanceNurbMan(); }
	void retractNurbMan() { m_pNurbs->retractNurbMan(); }
	void incrementU() { m_pNurbs->modifyU( 1 ); }
	void decrementU() { m_pNurbs->modifyU( -1 ); }
	void modifyWeight( float fX, float fY, float fVal );
	void toggleAffine() { m_pNurbs->toggleDrawAffine(); }
	void resetWeights() { m_pNurbs->resetWeights(); }
	void resetCurve() { m_pNurbs->resetCurve(); }

	// Terrain Manipulation
	void selectFace(float fX, float fY);
	void swapTerrains() { m_iCurrTerrain = (m_iCurrTerrain + 1) % 2; }

private:
	// For Singleton Implementation
	GraphicsManager(GLFWwindow* rWindow);
	GraphicsManager(const GraphicsManager* pCopy);
	static GraphicsManager* m_pInstance;

	// Window Reference
	GLFWwindow* m_pWindow;

	// Axis Buffer/Array Containers
	GLuint m_pVertexArray;
	GLuint m_pVertexBuffer;

	// Camera
	Camera* m_pCamera;
	vec3 getIntersection( float fX, float fY );

	// Render Functions
	void RenderScene();
	void renderAxis();

	// Manages Shaders for all assignments
	ShaderManager* m_pShaderMngr;
	EnvironmentManager* m_pEnvMngr;
	HypoCycloid* m_pHypoCyc;	// Modelling Design
	NURBS* m_pNurbs;

	Terrain* m_pTerrain[2];
	unsigned int m_iCurrTerrain;
};
