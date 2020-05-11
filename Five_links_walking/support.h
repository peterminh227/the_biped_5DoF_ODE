#include <ode/ode.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cmath>
#include <armadillo> // algebra library
#include <drawstuff/drawstuff.h>
using namespace libxl;
using namespace arma;
#ifdef MSVC
#pragma warning(disable:4244 4305)  // for VC++, no precision loss complaints
#endif
#ifdef dDOUBLE
#define dsDrawBox dsDrawBoxD
#define dsDrawSphere dsDrawSphereD
#define dsDrawCylinder dsDrawCylinderD
#define dsDrawCapsule dsDrawCapsuleD
#endif
// define parameter
#define PI					3.14159		// PI
#define DT					0.001	// sampling time
#define GKP					1000000		// eslaticity 1e6
#define GKD					7848.6		// viscocity -> create the critical damped
#define walking_time		10
#define excel 				0// excel record
#define FSM_enable			1
#define en_derivative		1  //
#define en_VBLA				1
#define use_CoM_Torso		0 // 
#define enable_zd			1
#define enable_sw			0
#define enable_sw_JW		1
#define Minh_enable			1
#define Minh_OSC			0
#define T_gap				5
#define t_start				1
#define v_c					0.6
#define v_d					0.85
#define exp					2.718
#define T_gap2				5
#define t_start2			14
#define v_c2				0.65
#define v_d2				0.5
// Information for Robot model
//---------------------------------
//　Body setting
//---------------------------------
// weight
#define Torso_Mass	12.5
#define Thigh_Mass	0.7
#define Shin_Mass	0.7
// using the point foot
#define Foot_Mass	0.1 
// Define box
#define Torso_L		0.1
#define Torso_W		0.1
#define Torso_H		0.42
// Define capsule
#define Thigh_R		0.02
#define Thigh_H		0.19
// Define capsule
#define Shin_R		0.015
#define Shin_H		0.155

// --> using the point foot
#define Foot_R		0.015
#define Foot_H		0.02
// --> define moment of inertia
#define Ib_yy       0.2254 // moi of body
#define If_yy		0.004472
#define Is_yy		0.004472
//

//
int sign(dReal value);
double min(dReal x, dReal y);
double max(dReal x, dReal y);
double getRotateAngle(dVector3 vec1, dVector3 vec2);
double deg2rad(dReal angle);
//void auto_MechEnergy(mat q, mat dq, double& KE_tot, double& PE_tot);
mat get_D(mat q_f);
mat get_jacobian1(mat q_f);
mat get_jacobian2(mat q_f);

