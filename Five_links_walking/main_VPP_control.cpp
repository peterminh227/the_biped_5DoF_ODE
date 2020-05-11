#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "libxl.h"
#include <iomanip>
#include <locale>
#include <sstream>
#include <string> // this should be already included in <sstream>
#include <armadillo> // algebra library
#include <cmath>
#include "support.h"
//---------------------------------
//　Physic world
//---------------------------------
// define the paramter for terrain ground
double gr_W = 0.2, gr_L = 0.2, gr_H[50], height_max = 0.015;
int maker = 0; // for flight phase
double num_box = 20, distance = 1;
// Control parameter
double ks = 10000;
double Lmax = 0.38;
double L0 = Lmax*0.97;
double L0_sw = L0*0.85;
double alpha0	=			1.84; // if use non-VBLA %1.85 OK 1.8326
double cs	=				10;
double cv	=				1;
double mu_fic	=			1;
double offset	=			-PI/80.0;// -PI/40
//double offset = deg2rad(-2.2); // change value
double c_swing		=		10;
double VBLA_coff	=		0.5;
double KE_tot = 0, PE_tot = 0;
double delta_E, delta_Ts;
double v_VBLA;
static double D_n,step_length = 0, v_bar; // capture region
static double pre_t = 0;
static dJointFeedback *feedback1 = new dJointFeedback;
static dJointFeedback *feedback2 = new dJointFeedback;

FILE *fp;
//configuration of physic world
static dWorldID world;
static dSpaceID space;
static dBodyID body[20];
static dJointID joint[20];
static dJointGroupID contactgroup;
static dGeomID ground;
static dSpaceID geom_group,right,left;
static dGeomID box[20];
static dGeomID ground_box[50];
double ground_position[50][3];
int leg_state[3];
// 
static double t=0;
// set of body
double mass;
double mass_part[4];
dQuaternion qrot1,qrot2,qrot3,qrot4,qrot5,qrot6,qrot10,qrot11;
double theta1,theta2,theta3,theta4,theta5,theta6;
double posi0[3],posi1[3],posi2[3],posi3[3],posi4[3],posi5[3],posi6[3];

//sensor
int rf=0,lf=0;
double tpos[3]={0,0,0};
double aja[7];
double ja[6], jar[6];
struct legstruct{
	mat F; // force(2,1)
	mat Tor; // Tor(7,1)
	mat J; //J(2,7)
	mat J_m;
	dReal zeta;
	dReal alpha;
	dReal alpha_dot;
	dReal length;
	dReal theta;
	dReal beta;
	dReal dL;
	dReal dth;
	int state;
	dReal Fs;
	dReal Ft;
	dReal Tau_ff;
	dReal Tau_fb;
	int contact;
	dReal phi_knee_dot; // for swing leg control
	dReal phi_hip_dot; //  for swing leg control 
	dReal Tau_hip;
	dReal Tau_knee;
	mat foot; // foot position
	dReal GRF;
	dReal Fn;
	dReal Fta;
	double Fx; // ground reaction force
	double Fz;
	double Fy;
};

// testing
dReal q1L_pre = 0, q31L_pre = 0, q32L_pre = 0, q41L_pre = 0, q42L_pre = 0;
//----------------------------------------------------------------------------------

double f(double k){if(k<=0.0)return(0.0);else return(k);};
double h(double k){if(k<=0.0)return(0.0);else return(1.0);};
double dir(double k){if(k>=0.0)return(1.0);else return(-1.);};
//----------------------------------------------------
static Book* book = xlCreateBook(); // xlCreateXMLBook() for xlsx
static Sheet* sheet = book->addSheet(L"Sheet1");
static double count = 1;
//----------------------------------------------------
// Camera position
//----------------------------------------------------
static void start()
{
  static float xyz[3] = {0.500f,-2.00f,0.900f};
  static float hpr[3] = {90.0000f,-0.000f,0.0000f};
  
  dsSetViewpoint (xyz,hpr);
}
//----------------------------------------------------
// collision
//----------------------------------------------------
static void nearCallback (void *data, dGeomID o1, dGeomID o2)
{
	
  int i,n,j;
  int o1i, o2i;
  o1i=(int)o1;
  o2i=(int)o2;
  // collision avoidance among body parts
  if(o1i==(int)geom_group || o1i == (int)left || o1i == (int)right){
  	  if( o2i==(int)geom_group || o2i == (int)left || o2i == (int)right)return;
  }
  // only collision with the ground
  //for (int j =0; j <50; j++)
  //{
  //int g1 = (o1 == ground) || (o1 == ground_box[j]);
  //int g2 = (o2 == ground) || (o2 == ground_box[j]);
  //if (!(g1^g2)) return;
  //}

  const int N = 100;
  dContact contact[N];
  n = dCollide (o1,o2,N,&contact[0].geom,sizeof(dContact));
  if (n > 0) {
    for (i=0; i<n; i++) 
		{
		//　Create contact when collision
		contact[i].surface.mode = dContactSoftERP | dContactSoftCFM | dContactApprox1;
		contact[i].surface.mu = 1;
		contact[i].surface.soft_erp = DT*GKP/(DT*GKP+GKD);
		contact[i].surface.soft_cfm = 1/(DT*GKP+GKD);
		dJointID c = dJointCreateContact (world,contactgroup,&contact[i]);
		dJointAttach(c,dGeomGetBody(contact[i].geom.g1),dGeomGetBody(contact[i].geom.g2));
		if ((o1i==(int)right && o2i==(int)ground) || (o2i==(int)right && o1i==(int)ground))
		{dJointSetFeedback(c,feedback1);}
		if ((o1i==(int)left && o2i==(int)ground) || (o2i==(int)left && o1i==(int)ground))
		{dJointSetFeedback(c,feedback2);}
		}
  }
   
  //　contact information of right foot -- adding the ground box
  for (int i = 0; i<num_box;i++)
  {
  if((o1i==(int)right && o2i==(int)ground) || (o1i==(int)right && o2i==(int)ground_box[i]) )
	  rf = n;		
  else if((o2i==(int)right && o1i==(int)ground) || (o2i==(int)right && o1i==(int)ground_box[i]))
	  rf = n;		
  //　contact information of left foot
  if ((o1i==(int)left && o2i==(int)ground) || (o1i==(int)left && o2i==(int)ground_box[i]))
	  lf = n;		
  else if((o2i==(int)left && o1i==(int)ground) || (o2i==(int)left && o1i==(int)ground_box[i]))
	  lf = n;	
  }
}
//----------------------------------------------------
//　Initial setting of robot
//----------------------------------------------------
static void create_model(void){
  //-------------------
  //　Initial setting of physical world
  //-------------------
  world = dWorldCreate();
  space = dHashSpaceCreate(0);
  contactgroup = dJointGroupCreate(0);
  dWorldSetGravity(world,0,0,-9.81);
  dWorldSetCFM(world,1e-5);
  dWorldSetERP(world,0.9);
  ground = dCreatePlane(space,0,0,1,0);
  int i, enable_initilize = 1; 
  dReal height = 0.36;
  dMass m; 
  // create initial angle
  dReal rad_rthigh = PI/6, rad_rshin = 0, rad_lthigh = -PI/10, rad_lshin = -PI/6; // inverse this angle for rotation because in XZ plane
  dReal rad_torso = 0;
  dMatrix3 rot_torso,rot_rthigh, rot_lthigh, rot_rshin, rot_lshin, rot_rfoot, rot_lfoot;
  dRFromAxisAndAngle(rot_rthigh,0,1,0,-rad_rthigh);
  dRFromAxisAndAngle(rot_lthigh,0,1,0,-rad_lthigh);
  dRFromAxisAndAngle(rot_rshin,0,1,0,-rad_rshin);
  dRFromAxisAndAngle(rot_lshin,0,1,0,-rad_lshin);
  dRFromAxisAndAngle(rot_rfoot,0,1,0,-rad_rshin);
  dRFromAxisAndAngle(rot_lfoot,0,1,0,-rad_lshin);
  dRFromAxisAndAngle(rot_torso,0,1,0,-rad_torso);
   
  //-------------------
  //　setting mass
  //-------------------
  mass_part[0]=Torso_Mass;
  mass_part[1]=Thigh_Mass;	
  mass_part[2]=Shin_Mass;
  mass_part[3]=Foot_Mass;
  // Setting of initial position of each object
  if (enable_initilize == 0) // start from the rest
  {
  posi0[0]=0; 
  posi0[1]=0; 
  posi0[2]=height + 0.5*Torso_H;
  // right thigh
  posi1[0]=posi0[0]; 
  posi1[1]=0; 
  posi1[2]=posi0[2]-0.5*Thigh_H - 0.5*Torso_H; 
  // right shin
  posi2[0]=posi1[0];
  posi2[1]=0; 
  posi2[2]=posi1[2]-0.5*Thigh_H-0.5*Shin_H; 
  // right foot
  posi3[0]=posi2[0];
  posi3[1]=0; 
  posi3[2]=posi2[2]-0.5*Shin_H-0.5*Foot_H; 
  // left thigh
  posi4[0]=posi0[0]; 
  posi4[1]=0; 
  posi4[2]=posi0[2]-0.5*Thigh_H -0.5*Torso_H; 
  // left shin
  posi5[0]=posi4[0];
  posi5[1]=0;
  posi5[2]=posi4[2]-0.5*Thigh_H-0.5*Shin_H; 
  // left foot
  posi6[0]=posi5[0];
  posi6[1]=0; 
  posi6[2]=posi5[2]-0.5*Shin_H-0.5*Foot_H; 
  // Object configuration
  body[0] = dBodyCreate (world);
  dBodySetPosition (body[0],posi0[0],posi0[1],posi0[2]);
  //dBodySetLinearVel(body[0],1,0,0);
  dMassSetBox (&m,1,Torso_L,Torso_W,Torso_H);
  mass = mass_part[0];
  dMassAdjust(&m,mass);
  dBodySetMass(body[0], &m);
  box[0] = dCreateBox(0,Torso_L,Torso_W,Torso_H);
  dGeomSetBody (box[0],body[0]);
  // capsule 1: right thigh
  body[1] = dBodyCreate (world);
  dBodySetPosition (body[1],posi1[0],posi1[1],posi1[2]);
  dMassSetCapsuleTotal(&m,mass,3,Thigh_R,Thigh_H);
  mass=mass_part[1];
  dMassAdjust(&m,mass);
  dBodySetMass(body[1],&m);
  box[1] = dCreateCapsule (0,Thigh_R,Thigh_H);
  dGeomSetBody (box[1],body[1]);
  // capsule 2: right shin
  body[2] = dBodyCreate (world);
  dBodySetPosition (body[2],posi2[0],posi2[1],posi2[2]);
  dMassSetCapsuleTotal(&m,mass,3,Shin_R,Shin_H);
  mass=mass_part[2];
  dMassAdjust (&m,mass);
  dBodySetMass (body[2],&m);
  box[2] = dCreateCapsule (0,Shin_R,Shin_H);
  dGeomSetBody (box[2],body[2]);
  // box 3: right foot
  body[3] = dBodyCreate (world);
  dBodySetPosition (body[3],posi3[0],posi3[1],posi3[2]);
  dMassSetCapsuleTotal (&m,mass,3,Foot_R,Foot_H);
  mass=mass_part[3];
  dMassAdjust (&m,mass);
  dBodySetMass (body[3],&m);
  box[3] = dCreateCapsule (0,Foot_R,Foot_H);
  dGeomSetBody (box[3],body[3]);
  // Box 4: left thigh
  body[4] = dBodyCreate (world);
  dBodySetPosition (body[4],posi4[0],posi4[1],posi4[2]);
  dMassSetCapsuleTotal (&m,mass,3,Thigh_R,Thigh_H);
  mass=mass_part[1];
  dMassAdjust (&m,mass);
  dBodySetMass (body[4],&m);
  box[4] = dCreateCapsule (0,Thigh_R,Thigh_H);
  dGeomSetBody (box[4],body[4]);
  // box 5:left shin
  body[5] = dBodyCreate (world);
  dBodySetPosition (body[5],posi5[0],posi5[1],posi5[2]);
  dMassSetCapsuleTotal (&m,mass,3,Shin_R,Shin_H);
  mass=mass_part[2];
  dMassAdjust (&m,mass);
  dBodySetMass (body[5],&m);
  box[5] = dCreateCapsule(0,Shin_R,Shin_H);
  dGeomSetBody (box[5],body[5]);
  // box 6: left foot
  body[6] = dBodyCreate (world);
  dBodySetPosition (body[6],posi6[0],posi6[1],posi6[2]);
  dMassSetCapsuleTotal (&m,mass,3,Foot_R,Foot_H);
  mass=mass_part[3];
  dMassAdjust (&m,mass);
  dBodySetMass (body[6],&m);
  box[6] = dCreateCapsule (0,Foot_R,Foot_H);
  dGeomSetBody (box[6],body[6]);

  //　Joint setting --> motor position
  //Joint 0: right hip joint
  joint[0] = dJointCreateHinge (world,0);
  dJointAttach (joint[0],body[0],body[1]);
  const dReal *a1 = dBodyGetPosition (body[1]);
  dJointSetHingeAnchor (joint[0],a1[0],a1[1],a1[2]+0.5*Thigh_H);
  dJointSetHingeAxis (joint[0],0,1,0);
  dJointSetHingeParam (joint[0],dParamLoStop,-dInfinity);
  dJointSetHingeParam (joint[0],dParamHiStop,dInfinity);
  dJointSetHingeParam (joint[0],dParamFudgeFactor,0.1);
  // Joint 1: Joint right Knee joint
  joint[1] = dJointCreateHinge (world,0);
  dJointAttach (joint[1],body[1],body[2]);
  const dReal *a2 = dBodyGetPosition (body[2]);
  dJointSetHingeAnchor (joint[1],a2[0],a2[1],a2[2]+0.5*Shin_H);
  dJointSetHingeAxis (joint[1],0,1,0);
  dJointSetHingeParam (joint[1],dParamLoStop,-dInfinity);
  dJointSetHingeParam (joint[1],dParamHiStop,dInfinity);
  dJointSetHingeParam (joint[1],dParamFudgeFactor,0.1);
  // Joint 2 right ankle joint
  joint[2] = dJointCreateFixed (world,0);
  dJointAttach (joint[2],body[2],body[3]);
  dJointSetFixed(joint[2]);  

  // for configure the foot rotation
  //const dReal *a3 = dBodyGetPosition (body[3]);
  //dJointSetHingeAnchor (joint[2],a2[0],a3[1],a3[2]+0.5*Foot_H);
  //dJointSetHingeAxis (joint[2],0,1,0);
  //dJointSetHingeParam (joint[2],dParamLoStop,-dInfinity);
  //dJointSetHingeParam (joint[2],dParamHiStop,dInfinity);
  //dJointSetHingeParam (joint[2],dParamFudgeFactor,0.1);

  // Joint hip left joint
  joint[3] = dJointCreateHinge (world,0);
  dJointAttach (joint[3],body[0],body[4]);
  const dReal *b1 = dBodyGetPosition (body[4]);
  dJointSetHingeAnchor (joint[3],b1[0],b1[1],b1[2]+0.5*Thigh_H);
  dJointSetHingeAxis (joint[3],0,1,0);
  dJointSetHingeParam (joint[3],dParamLoStop,-dInfinity);
  dJointSetHingeParam (joint[3],dParamHiStop,dInfinity);
  dJointSetHingeParam (joint[3],dParamFudgeFactor,0.1);
  // Joint 4: left knee joint
  joint[4] = dJointCreateHinge (world,0);
  dJointAttach (joint[4],body[4],body[5]);
  const dReal *b2 = dBodyGetPosition (body[5]);
  dJointSetHingeAnchor (joint[4],b2[0],b2[1],b2[2]+0.5*Shin_H);
  dJointSetHingeAxis (joint[4],0,1,0);	
  dJointSetHingeParam (joint[4],dParamLoStop,-dInfinity);
  dJointSetHingeParam (joint[4],dParamHiStop,dInfinity);
  dJointSetHingeParam (joint[4],dParamFudgeFactor,0.1);
  // Joint 5: ankle left joint
  joint[5] = dJointCreateFixed(world,0);
  dJointAttach (joint[5],body[5],body[6]);
  dJointSetFixed(joint[5]);
  }
  else //**************************************************************************************//
  {
  posi0[0]= 0.5*Torso_H*sin(-rad_torso); 
  posi0[1]= 0; 
  posi0[2]= height + 0.5*Torso_H*cos(rad_torso);

  // right thigh
  posi1[0]= 0  + 0.5*Thigh_H*sin(rad_rthigh); 
  posi1[1]=0; 
  posi1[2]=posi0[2] - 0.5*Torso_H*cos(rad_torso) - 0.5*Thigh_H*cos(rad_rthigh); 
  // right shin
  posi2[0]=posi1[0] + 0.5*Thigh_H*sin(rad_rthigh) + 0.5*Shin_H*sin(rad_rshin);
  posi2[1]=0; 
  posi2[2]=posi1[2] - 0.5*Thigh_H*cos(rad_rthigh) - 0.5*Shin_H*cos(rad_rshin); 
  // right foot
  posi3[0]=posi2[0] + 0.5*Shin_H*sin(rad_rshin) + 0.5*Foot_H*sin(rad_rshin);
  posi3[1]=0; 
  posi3[2]=posi2[2] - 0.5*Shin_H*cos(rad_rshin) - 0.5*Foot_H*cos(rad_rshin); 
  // left thigh
  posi4[0]=0 + 0.5*Thigh_H*sin(rad_lthigh);  
  posi4[1]=0; 
  posi4[2]=posi0[2] - 0.5*Torso_H - 0.5*Thigh_H*cos(rad_lthigh); 
  // left shin
  posi5[0]=posi4[0] +  0.5*Thigh_H*sin(rad_lthigh) + 0.5*Shin_H*sin(rad_lshin);
  posi5[1]=0;
  posi5[2]=posi4[2] -  0.5*Thigh_H*cos(rad_lthigh) - 0.5*Shin_H*cos(rad_lshin); 
  // left foot
  posi6[0]=posi5[0] + 0.5*Shin_H*sin(rad_lshin) + 0.5*Foot_H*sin(rad_lshin);
  posi6[1]=0; 
  posi6[2]=posi5[2] - 0.5*Shin_H*cos(rad_lshin) - 0.5*Foot_H*cos(rad_lshin); 
  // Object configuration
  body[0] = dBodyCreate (world);
  dBodySetPosition (body[0],posi0[0],posi0[1],posi0[2]);
  dBodySetRotation(body[0],rot_torso);
  dBodySetLinearVel(body[0],0.8,0,0);
  dMassSetBox (&m,1,Torso_L,Torso_W,Torso_H); 
  // set moment of inertia 
  //void dMassSetParameters (dMass *, dReal themass,
  //                       dReal cgx, dReal cgy, dReal cgz,
  //                       dReal I11, dReal I22, dReal I33,
  //                       dReal I12, dReal I13, dReal I23);

  mass = mass_part[0];
  //dMassAdjust(&m,mass);
  dMassSetParameters(&m,mass,0,0,0,Ib_yy,Ib_yy,Ib_yy,0,0,0);
  dBodySetMass(body[0], &m);
  box[0] = dCreateBox(0,Torso_L,Torso_W,Torso_H);
  dGeomSetBody (box[0],body[0]);
  //
  body[10] = dBodyCreate (world);
  dBodySetPosition (body[10],2,0,0.5);
  dMassSetBox (&m,1,2,2,0.5);
  mass = 10;
  dMassAdjust(&m,mass);
  dBodySetMass(body[10], &m);
  box[10] = dCreateBox(0,2,2,0.5);
  dGeomSetBody (box[10],body[10]);
  // capsule 1: right thigh
  body[1] = dBodyCreate (world);
  dBodySetPosition(body[1],posi1[0],posi1[1],posi1[2]);
  dBodySetRotation(body[1],rot_rthigh);
  dMassSetCapsuleTotal(&m,mass,3,Thigh_R,Thigh_H);
  mass=mass_part[1];
  //dMassAdjust(&m,mass);
  dMassSetParameters(&m,mass,0,0,0,If_yy,If_yy,If_yy,0,0,0);
  dBodySetMass(body[1],&m);
  box[1] = dCreateCapsule (0,Thigh_R,Thigh_H);
  dGeomSetBody (box[1],body[1]);
  // capsule 2: right shin
  body[2] = dBodyCreate (world);
  dBodySetPosition(body[2],posi2[0],posi2[1],posi2[2]);
  dBodySetRotation(body[2],rot_rshin);
  dMassSetCapsuleTotal(&m,mass,3,Shin_R,Shin_H);
  mass=mass_part[2];
  //dMassAdjust (&m,mass);
  dMassSetParameters(&m,mass,0,0,0,Is_yy,Is_yy,Is_yy,0,0,0);
  dBodySetMass (body[2],&m);
  box[2] = dCreateCapsule (0,Shin_R,Shin_H);
  dGeomSetBody (box[2],body[2]);
  // box 3: right foot
  body[3] = dBodyCreate (world);
  dBodySetPosition(body[3],posi3[0],posi3[1],posi3[2]);
  dBodySetRotation(body[3],rot_rfoot);
  dMassSetCapsuleTotal (&m,mass,3,Foot_R,Foot_H);
  mass=mass_part[3];
  //dMassAdjust (&m,mass);
  dMassSetParameters(&m,mass,0,0,0,If_yy,If_yy,If_yy,0,0,0);
  dBodySetMass (body[3],&m);
  box[3] = dCreateCapsule (0,Foot_R,Foot_H);
  dGeomSetBody (box[3],body[3]);
  // Box 4: left thigh
  body[4] = dBodyCreate (world);
  dBodySetPosition(body[4],posi4[0],posi4[1],posi4[2]);
  dBodySetRotation(body[4],rot_lthigh);
  dMassSetCapsuleTotal (&m,mass,3,Thigh_R,Thigh_H);
  mass=mass_part[1];
  dMassAdjust (&m,mass);
  dBodySetMass (body[4],&m);
  box[4] = dCreateCapsule (0,Thigh_R,Thigh_H);
  dGeomSetBody (box[4],body[4]);
  // box 5:left shin
  body[5] = dBodyCreate (world);
  dBodySetPosition(body[5],posi5[0],posi5[1],posi5[2]);
  dBodySetRotation(body[5],rot_lshin);
  dMassSetCapsuleTotal (&m,mass,3,Shin_R,Shin_H);
  mass=mass_part[2];
  //dMassAdjust (&m,mass);
  dMassSetParameters(&m,mass,0,0,0,Is_yy,Is_yy,Is_yy,0,0,0);
  dBodySetMass (body[5],&m);
  box[5] = dCreateCapsule(0,Shin_R,Shin_H);
  dGeomSetBody (box[5],body[5]);
  // box 6: left foot
  body[6] = dBodyCreate (world);
  dBodySetPosition(body[6],posi6[0],posi6[1],posi6[2]);
  dBodySetRotation(body[6],rot_lfoot);
  dMassSetCapsuleTotal (&m,mass,3,Foot_R,Foot_H);
  mass=mass_part[3];
  dMassAdjust (&m,mass);
  dBodySetMass (body[6],&m);
  box[6] = dCreateCapsule (0,Foot_R,Foot_H);
  dGeomSetBody (box[6],body[6]);
  // create ground by box:
  
  //　Joint setting --> motor position
  //Joint 0: right hip joint
  joint[0] = dJointCreateHinge (world,0);
  dJointAttach (joint[0],body[0],body[1]);
  const dReal *a1 = dBodyGetPosition (body[0]);
  dJointSetHingeAnchor (joint[0],a1[0],a1[1],a1[2]-0.5*Torso_H);
  dJointSetHingeAxis (joint[0],0,1,0);
  dJointSetHingeParam (joint[0],dParamLoStop,-PI);
  dJointSetHingeParam (joint[0],dParamHiStop,PI);
  dJointSetHingeParam (joint[0],dParamFudgeFactor,0.1);
  // Joint 1: Joint right Knee joint
  joint[1] = dJointCreateHinge (world,0);
  dJointAttach (joint[1],body[1],body[2]);
  dVector3  a2;
  dBodyGetRelPointPos(body[1],0,0,-Thigh_H/2,a2);
  dJointSetHingeAnchor (joint[1],a2[0],a2[1],a2[2]);
  dJointSetHingeAxis (joint[1],0,1,0);
  dJointSetHingeParam (joint[1],dParamLoStop,-PI+PI/6);
  dJointSetHingeParam (joint[1],dParamHiStop,PI/6);
  dJointSetHingeParam (joint[1],dParamFudgeFactor,0.1);
  // Joint 2 right ankle joint
  joint[2] = dJointCreateFixed (world,0);
  dJointAttach (joint[2],body[2],body[3]);
  dJointSetFixed(joint[2]);  

  // for configure the foot rotation
  //const dReal *a3 = dBodyGetPosition (body[3]);
  //dJointSetHingeAnchor (joint[2],a2[0],a3[1],a3[2]+0.5*Foot_H);
  //dJointSetHingeAxis (joint[2],0,1,0);
  //dJointSetHingeParam (joint[2],dParamLoStop,-dInfinity);
  //dJointSetHingeParam (joint[2],dParamHiStop,dInfinity);
  //dJointSetHingeParam (joint[2],dParamFudgeFactor,0.1);

  // Joint hip left joint
  joint[3] = dJointCreateHinge (world,0);
  dJointAttach (joint[3],body[0],body[4]);
  dVector3 b1;
  dBodyGetRelPointPos(body[4],0,0,+Thigh_H/2,b1);
  dJointSetHingeAnchor (joint[3],b1[0],b1[1],b1[2]);
  dJointSetHingeAxis (joint[3],0,1,0);
  dJointSetHingeParam (joint[3],dParamLoStop,-PI);
  dJointSetHingeParam (joint[3],dParamHiStop,PI);
  dJointSetHingeParam (joint[3],dParamFudgeFactor,0.1);
  // Joint 4: left knee joint
  joint[4] = dJointCreateHinge (world,0);
  dJointAttach (joint[4],body[4],body[5]);
  dVector3 b2;
  dBodyGetRelPointPos(body[5],0,0,Shin_H/2,b2);
  dJointSetHingeAnchor (joint[4],b2[0],b2[1],b2[2]);
  dJointSetHingeAxis (joint[4],0,1,0);	
  dJointSetHingeParam (joint[4],dParamLoStop,-PI+PI/12);
  dJointSetHingeParam (joint[4],dParamHiStop,PI/12);
  dJointSetHingeParam (joint[4],dParamFudgeFactor,0.1);
  // Joint 5: ankle left joint
  joint[5] = dJointCreateFixed(world,0);
  dJointAttach (joint[5],body[5],body[6]);
  dJointSetFixed(joint[5]);
  }
  // for configure the foot rotation
  //const dReal *b3 = dBodyGetPosition (body[6]);
  //dJointSetHingeAnchor (joint[5],b2[0],b3[1],b3[2]+0.5*Foot_H);
  //dJointSetHingeAxis (joint[5],0,1,0);
  //dJointSetHingeParam (joint[5],dParamLoStop,-dInfinity);
  //dJointSetHingeParam (joint[5],dParamHiStop,dInfinity);
  //dJointSetHingeParam (joint[5],dParamFudgeFactor,0.1);
  //　Orther, intialization parameter
  
  rf=0; lf=0; // sensor to detect wherether leg touch the ground
  tpos[0]=0; 
  tpos[1]=0; 
  tpos[2]=0;
  t=0;
    
  //dBodySetRotation(body[5],R3);
	

  // geometry grouping
  geom_group = dSimpleSpaceCreate(space);  
  dSpaceSetCleanup(geom_group,0);
  dSpaceAdd(geom_group,box[0]);
  dSpaceAdd(geom_group,box[1]);
  dSpaceAdd(geom_group,box[4]);
  dSpaceAdd(geom_group,box[2]);
  dSpaceAdd(geom_group,box[5]);
   //
 
  //--> space ID for foot right and left
  right = dSimpleSpaceCreate(space);  
  dSpaceSetCleanup(right,0);
  dSpaceAdd(right,box[3]);
  left = dSimpleSpaceCreate(space);  
  dSpaceSetCleanup(left,0);
  dSpaceAdd(left,box[6]);
  // environment
  for(int i= 0; i<num_box; i++)
	{
		if ((i % 2) == 1)
		gr_H[i] = height_max ; 
		else gr_H[i] = 0;
		ground_box[i] = dCreateBox(space,gr_L,gr_W,gr_H[i]);
		if (i == 0) {ground_position[i][0] = distance;}
		else {ground_position[i][0] = ground_position[i-1][0] +  gr_L;}
		dGeomSetPosition(ground_box[i],ground_position[i][0],0,0);
	}
  
}

//----------------------------------------------------
// Destroy simulation
//----------------------------------------------------
static void destroy(void){
	int i;
	dJointGroupDestroy (contactgroup);
	for(i=0;i<7;i++){
	  dGeomDestroy (box[i]);
	}
	dSpaceDestroy(space);
	dWorldDestroy(world);
}

// wrap code
template <class T>
inline void wrapTo2PiInPlace(T &a)
{
	bool was_neg = a<0;
	a = fmod(a, static_cast<T>(2.0*PI));
	if (was_neg) a+=static_cast<T>(2.0*PI);
}
template <class T>
inline T wrapTo2Pi(T a)
{
	wrapTo2PiInPlace(a);
	return a;
}
template <class T>
inline T wrapToPi(T a)
{
	return wrapTo2Pi(a+ static_cast<T>(PI)) - static_cast<T>(PI);
}

//----------------------------------------------------
// Simulation loop
//----------------------------------------------------
static void simLoop (int pause)
{	
	//
	// Dynamics calculation and 1 step elapsed
	rf=0;lf=0;  
    dSpaceCollide (space,0,&nearCallback);	// Computation of Collision Detect
    dWorldStep (world,DT);					// Sampling time  
    dJointGroupEmpty(contactgroup);		// remove all contact joints
	int i;
	legstruct leg[3]; // leg[1]: right; leg[2]: left
	dReal L1, L3, L4; // L1: Torso length, L3: thigh, L4: shin
	dVector3 hpos, rthpos, lthpos, rshpos, lshpos, r_f, l_f, nominal,tangital; 
	dVector3 hip_tor, hip_rknee, hip_lknee, rknee_rshin,lknee_lshin; // define vector 
	dVector3 rfoot_hip, lfoot_hip, rfoot_CoM, lfoot_CoM;
	dReal q31L, q41L, q32L, q42L, q1L; // the relative angle 
	dReal q1 = 0, q2 = 0, q3 = 0, q4 = 0,q7 = 0; // q1 q2: right leg, q3 q4: left leg- JW coordinate 
	dReal dq31L, dq41L, dq32L, dq42L, dq1L, dxCoM, dzCoM; // velocity 
	dReal xCoM, zCoM, phi_tilde;
	// Controller varibles
	dReal val_swing, range, SW_L_d; //VBLA
	mat q_f(7,1), dq_f(7,1), R_dX(2,1), L_dX(2,1), V(2,1), G(2,1), O(2,1);
	dReal alpha_d, length; // using for VBLA
	
	//*****************************************************************************************************
	// Sensor capture
	// Joint angle
	ja[0]=dJointGetHingeAngle(joint[0]);
	ja[1]=dJointGetHingeAngle(joint[1]);
	//ja[2]=dJointGetHingeAngle (joint[2]);	right_ankle
	ja[3]=dJointGetHingeAngle(joint[3]);
	ja[4]=dJointGetHingeAngle(joint[4]);	
	//ja[5]=dJointGetHingeAngle (joint[5]);	left_ankle
	// Joint angular velocity
	jar[0]=dJointGetHingeAngleRate(joint[0]);
	jar[1]=dJointGetHingeAngleRate(joint[1]);
	//jar[2]=dJointGetHingeAngleRate (joint[2]);
	jar[3]=dJointGetHingeAngleRate(joint[3]);
	jar[4]=dJointGetHingeAngleRate(joint[4]);	
	//jar[5]=dJointGetHingeAngleRate (joint[5]);
	// absolute position of each object
	const dReal *pos0= dBodyGetPosition(body[0]); //torso
	const dReal *pos1= dBodyGetPosition(body[1]);
	const dReal *pos2= dBodyGetPosition(body[2]);
	const dReal *pos3= dBodyGetPosition(body[3]);// foot position --> center of capsule, not touch point
	const dReal *pos4= dBodyGetPosition(body[4]);
	const dReal *pos5= dBodyGetPosition(body[5]);
	const dReal *pos6= dBodyGetPosition(body[6]);// foot position --> center of capsule, not touch point
	//
	const dReal *pos0rate= dBodyGetLinearVel(body[0]);
	const dReal *pos1rate= dBodyGetLinearVel(body[1]);
	const dReal *pos2rate= dBodyGetLinearVel(body[2]);
	const dReal *pos3rate= dBodyGetLinearVel(body[3]);//
	const dReal *pos4rate= dBodyGetLinearVel(body[4]);
	const dReal *pos5rate= dBodyGetLinearVel(body[5]);
	const dReal *pos6rate= dBodyGetLinearVel(body[6]);// 
	//
	const dReal *angular = dBodyGetAngularVel(body[0]);
	dVector3 CoMpos, dCoM;
	dReal rh;
	for ( int  i=0; i<3; i++) 
	{
	CoMpos[i] = (pos0[i]*Torso_Mass + pos1[i]*Thigh_Mass + pos2[i]*Shin_Mass + pos3[i]*Foot_Mass + pos4[i]*Thigh_Mass + pos5[i]*Shin_Mass + pos6[i]*Foot_Mass)/(Torso_Mass +2*(Thigh_Mass+Shin_Mass+Foot_Mass));
	dCoM[i] = (pos0rate[i]*Torso_Mass + pos1rate[i]*Thigh_Mass + pos2rate[i]*Shin_Mass + pos3rate[i]*Foot_Mass + pos4rate[i]*Thigh_Mass + pos5rate[i]*Shin_Mass + pos6rate[i]*Foot_Mass)/(Torso_Mass +2*(Thigh_Mass+Shin_Mass+Foot_Mass));
	}
	// take the position
	//*****************************************************************************************************//
	
	//xCoM = CoMpos[0];
	//zCoM = CoMpos[2];
	if (use_CoM_Torso == 1)
	{
	xCoM = pos0[0];
	zCoM = pos0[2];
	
	} else {xCoM = CoMpos[0]; zCoM = CoMpos[2];}
	// Absolute coordinate position of the center of mass position
	tpos[0]=pos0[0];
	tpos[1]=pos0[1];
	tpos[2]=pos0[2];
	nominal[0]=0; tangital[0] = 1;
	nominal[1]=0; tangital[1] = 0; 
	nominal[2]=1; tangital[2] = 0; 
	dBodyGetRelPointPos(body[0],0,0,-Torso_H/2,hpos);
	dBodyGetRelPointPos(body[1],0,0, -Thigh_H/2 , rthpos);
	dBodyGetRelPointPos(body[2],0,0,-Shin_H/2 , rshpos);
	dBodyGetRelPointPos(body[3],0,0,-Foot_H/2 - Foot_R, r_f);//
	dBodyGetRelPointPos(body[4],0,0,-Thigh_H/2 , lthpos);
	dBodyGetRelPointPos(body[5],0,0,-Shin_H/2 , lshpos);
	dBodyGetRelPointPos(body[6],0,0,-Foot_H/2 - Foot_R, l_f);//
	// vector calculation
	for ( int  i=0; i<3; i++) {
	hip_tor[i] = tpos[i] - hpos[i];
	hip_rknee[i] = rthpos[i] - hpos[i];
	hip_lknee[i] = lthpos[i] - hpos[i];
	rknee_rshin[i] = rshpos[i] - rthpos[i];
	lknee_lshin[i] = lshpos[i] - lthpos[i];
	rfoot_hip[i] = hpos[i] - r_f[i];
	lfoot_hip[i] = hpos[i] - l_f[i];
	rfoot_CoM[i] = CoMpos[i] - r_f[i];
	lfoot_CoM[i] = CoMpos[i] - l_f[i];
	}
	//
	rh = sqrt(pow(xCoM-hpos[0],2)+pow(zCoM-hpos[2],2));
	//
	q1L = getRotateAngle(nominal,hip_tor);q1L = wrapToPi(q1L); 
	q31L = getRotateAngle(hip_tor,hip_rknee) ; q32L = getRotateAngle(hip_tor,hip_lknee);
	q41L = getRotateAngle(hip_rknee,rknee_rshin); q41L = wrapToPi(q41L);
	q42L = getRotateAngle(hip_lknee,lknee_lshin); q42L = wrapToPi(q42L);
	//--> Minh coordinate 
	q1 = q31L-PI; q2 = q41L; q3 = q32L-PI; q4 = q42L; q7 = q1L; // JW coordinate
	leg[1].zeta = wrapToPi(getRotateAngle(rfoot_CoM,rfoot_hip)); leg[2].zeta = wrapToPi(getRotateAngle(lfoot_CoM, lfoot_hip));
	leg[1].alpha = wrapToPi(getRotateAngle(tangital, rfoot_hip)); leg[2].alpha = wrapToPi(getRotateAngle(tangital, lfoot_hip));
	dq1L = -angular[1]; // xz plane
	dq31L = jar[0]; dq32L = jar[3]; dq41L = jar[1]; dq42L = jar[4]; 
	dxCoM = dCoM[0]; dzCoM = dCoM[2];
	dq_f(0,0) = dq31L; dq_f(1,0) = dq41L; dq_f(2,0) = dq32L; dq_f(3,0) = dq42L; dq_f(4,0) = pos0rate[0]; dq_f(5,0) = pos0rate[2]; dq_f(6,0) = dq1L; //JW
	q_f(0,0) = q1; q_f(1,0) = q2; q_f(2,0) = q3; q_f(3,0) = q4; q_f(4,0) = xCoM; q_f(5,0) = zCoM; q_f(6,0) = q7; //JW
	// for swing leg control -  Peterminh
	leg[1].alpha_dot  = dq31L+dq41L/2+dq1L; leg[2].alpha_dot = dq32L+dq42L/2 + dq1L;
	leg[1].phi_knee_dot = dq41L; leg[2].phi_knee_dot = dq42L;
	//printf("value at source:%f \n",leg[1].phi_knee_dot);
	//printf("value at source:%f \n",leg[2].phi_knee_dot); 
	// testing angle:
	dReal t_dq1L, t_dq31L, t_dq32L, t_dq41L, t_dq42L;
	t_dq1L = (q1L - q1L_pre)/DT; t_dq31L = (q31L - q31L_pre)/DT; t_dq32L = (q32L - q32L_pre)/DT;t_dq41L = (q41L - q41L_pre)/DT;t_dq42L = (q42L - q42L_pre)/DT;
	q1L_pre = q1L; q31L_pre = q31L; q32L_pre  = q32L; q41L_pre = q41L; q42L_pre = q42L; 
	//absolute angle of each object
	const dReal *qua0= dBodyGetQuaternion(body[0]);
	const dReal *qua1= dBodyGetQuaternion(body[1]);
	const dReal *qua2= dBodyGetQuaternion(body[2]);
	const dReal *qua3= dBodyGetQuaternion(body[3]); // 
	const dReal *qua4= dBodyGetQuaternion(body[4]);
	const dReal *qua5= dBodyGetQuaternion(body[5]);
	const dReal *qua6= dBodyGetQuaternion(body[6]); // left foot
	// absolute angle in XZ plane
	aja[0]=fabs(2*acos(qua0[0]))*dir(qua0[2]);
	aja[1]=fabs(2*acos(qua1[0]))*dir(qua1[2]);
	aja[2]=fabs(2*acos(qua2[0]))*dir(qua2[2]);
	aja[3]=fabs(2*acos(qua3[0]))*dir(qua3[2]);
	aja[4]=fabs(2*acos(qua4[0]))*dir(qua4[2]);
	aja[5]=fabs(2*acos(qua5[0]))*dir(qua5[2]);
	aja[6]=fabs(2*acos(qua6[0]))*dir(qua6[2]);
	//
	leg[1].foot = zeros(2,1);
	leg[1].foot(0,0) = r_f[0]; 	leg[1].foot(1,0) = r_f[2];
	leg[2].foot = zeros(2,1);
	leg[2].foot(0,0) = l_f[0]; leg[2].foot(1,0) = l_f[2];
	// energy of model
	//auto_MechEnergy(q_f,dq_f,KE_tot,PE_tot);
	// VPP-JW controller
	L1 = Torso_H; L3 = Thigh_H; L4 = Shin_H +Foot_H + Foot_R;
	//
	leg[1].length = pow((L4*L4 + L3*L3 + 2*L4*L3*cos(q2)),0.5);  
	leg[1].theta = q1 + (PI)/2 + atan2(- L3 - L4*cos(q2), L4*sin(q2));
	//
	leg[2].length = pow((L4*L4 + L3*L3 + 2*L4*L3*cos(q4)),0.5);
	leg[2].theta = q3 + (PI)/2 + atan2(- L3 - L4*cos(q4), L4*sin(q4));
	// initial state for matrix
	leg[1].J = zeros(2,7); leg[2].J = zeros(2,7);
	leg[1].J_m = zeros(2,7); leg[2].J_m = zeros(2,7);
	leg[1].F = zeros(2,1); leg[2].F = zeros(2,1);
	leg[1].Tor = zeros(7,1); leg[2].Tor = zeros(7,1);
	// jacobian for each leg polar jacobian
	leg[1].J(0,1)=-(L4*L3*sin(q2))/pow((L4*L4 + L3*L3 + 2*L4*L3*cos(q2)),0.5);
	leg[1].J(1,0)=1;
	leg[1].J(1,1)=(L4*(L4 + L3*cos(q2)))/(L4*L4 + L3*L3 + 2*L4*L3*cos(q2));
	//
	leg[2].J(0,3)=-(L4*L3*sin(q4))/pow((L4*L4 + L3*L3 + 2*L4*L3*cos(q4)),0.5);
	leg[2].J(1,2)=1;
	leg[2].J(1,3)=(L4*(L4 + L3*cos(q4)))/(L4*L4 + L3*L3 + 2*L4*L3*cos(q4));
	//-> Minh's coordinate
	//leg[1].J_m(0,0)=L3*cos(q31L + q1L) + L4*cos(q31L + q41L + q1L);
	//leg[1].J_m(0,1)=L4*cos(q31L + q41L + q1L);
	//leg[1].J_m(0,6)=L3*cos(q31L + q1L) + L4*cos(q31L + q41L + q1L);
	//leg[1].J_m(1,0)=L3*sin(q31L + q1L) + L4*sin(q31L + q41L + q1L);
	//leg[1].J_m(1,1)=L4*sin(q31L + q41L + q1L);
	//leg[1].J_m(1,6)=L3*sin(q31L + q1L) + L4*sin(q31L + q41L + q1L);
	leg[1].J_m(0,0)=- L3*cos(q7 + q1) - L4*cos(q7 + q2 + q1);
	leg[1].J_m(0,1)=-L4*cos(q7 + q2 + q1);
	leg[1].J_m(0,6)=- L3*cos(q7 + q1) - L4*cos(q7 + q2 + q1);
	//
	leg[1].J_m(1,0)=- L3*sin(q7 + q1) - L4*sin(q7 + q2 + q1);
	leg[1].J_m(1,1)=-L4*sin(q7 + q2 + q1);
	leg[1].J_m(1,6)=- L3*sin(q7 + q1) - L4*sin(q7 + q2 + q1);

	//
	//leg[2].J_m(0,2)=L3*cos(q32L + q1L) + L4*cos(q32L + q42L + q1L);
	//leg[2].J_m(0,3)=L4*cos(q32L + q42L + q1L);
	//leg[2].J_m(0,6)=L3*cos(q32L + q1L) + L4*cos(q32L + q42L + q1L);
	//leg[2].J_m(1,2)=L3*sin(q32L + q1L) + L4*sin(q32L + q42L + q1L);
	//leg[2].J_m(1,3)=L4*sin(q32L + q42L + q1L);
	//leg[2].J_m(1,6)=L3*sin(q32L + q1L) + L4*sin(q32L + q42L + q1L);

	 leg[2].J_m(0,2)=- L3*cos(q7 + q3) - L4*cos(q7 + q4+ q3);
	 leg[2].J_m(0,3)=-L4*cos(q7 + q4 + q3);
 	 leg[2].J_m(0,6)=- L3*cos(q7 + q3) - L4*cos(q7 + q4 + q3);
     leg[2].J_m(1,2)=- L3*sin(q7 + q3) - L4*sin(q7 + q4 + q3);
	 leg[2].J_m(1,3)=-L4*sin(q7 + q4 + q3);
     leg[2].J_m(1,6)=- L3*sin(q7 + q3) - L4*sin(q7 + q4 + q3);
	//
	R_dX = leg[1].J * dq_f; 
	leg[1].dL = R_dX(0,0);
	leg[1].dth = R_dX(1,0);
	//
	L_dX = leg[2].J * dq_f;
	leg[2].dL = L_dX(0,0);
	leg[2].dth = L_dX(1,0);
	// calculate the control angle of trunk
	if (en_derivative == 0) // no need trunk vel
	{
		phi_tilde = -wrapToPi(cs*wrapToPi(q7-offset));
	}
	else { // using trunk velocity
		phi_tilde = -wrapToPi(cs*wrapToPi(q7 - offset)) - cv*dq1L;
	}
	leg[1].beta = leg[1].zeta + phi_tilde;
	leg[2].beta = leg[2].zeta + phi_tilde;
	// limitation of the control angle
	dReal minangle  = PI/2 - atan(mu_fic), maxangle = PI/2 + atan(mu_fic);
	dReal minbeta = deg2rad(-89), maxbeta = deg2rad(89);
	
	for (i=1;i<3;i++)
	{
		if (leg[i].beta > min(leg[i].alpha - minangle,maxbeta))
			{
				leg[i].beta = min(leg[i].alpha - minangle,maxbeta);
			}
		if (leg[i].beta < max(leg[i].alpha - maxangle,minbeta)) 
			{
				leg[i].beta = max(leg[i].alpha - maxangle,minbeta);
			}
	}
	// VBLA controller
	if (t < t_start) {v_VBLA = v_c;}
	if ((t>=t_start) && (t<t_start+T_gap))
	{
		v_VBLA = t*(v_d-v_c)/T_gap + ((T_gap+1)*v_c - v_d)/T_gap;
		//offset = -PI/80 - (v_VBLA -dq_f(4,0));
	}
	if (t>t_start+T_gap)
	{
	v_VBLA = dq_f(4,0);
	//offset = -PI/80 - (v_VBLA -dq_f(4,0));
	}
	v_VBLA = dq_f(4,0);
	if (en_VBLA == 1)
	{	
		if (rf ==0)
		{
			length = leg[1].length ;
		}
		else
		{
			length = leg[2].length;
		}
		V(0,0) = dq_f(4,0) - (v_VBLA - dq_f(4,0))/VBLA_coff; V(1,0) = dq_f(5,0);
		V= V*1/sqrt(9.81*length);
		G(0,0) = 0; G(1,0) = -1;
		O = V*VBLA_coff + G*(1-VBLA_coff);
		//alpha_d = PI + atan2(O(1,0),O(0,0));
		alpha_d = atan2(-O(1,0),-O(0,0));
	}
	else
	{
		alpha_d = alpha0;
	}
	// Finite State Machine - FSM
	int  isDS, isSwCond1,isSwCond2,isSwCond3,isSwCond4, isStCond1;
	int  sub_i;
	int  trace = 1;
	// pre-define ?.? the problem in here haha
	leg[1].state = leg_state[1]; leg[2].state = leg_state[2];
	// Finite state machine:
	// pre-define
	if (rf > 0) {leg[1].contact = 1;} else {leg[1].contact = 0;}
	if (lf > 0) {leg[2].contact = 1;} else {leg[2].contact = 0;}
	if (FSM_enable == 1 )
	{
	if ((leg[1].state == 1) && (leg[2].state) == 1) {isDS = 1;} else {isDS = 0;}
	for (int i = 1; i <3; i++)
		{	
			if (leg[i].state == 1) // stance phase
				{
					if (i == 1) {sub_i = 2;} else {sub_i = 1;} // choose oppsite leg
					if ((leg[sub_i].alpha - PI/2.0) < deg2rad(15)) {isSwCond1 = 1;} else {isSwCond1 = 0;} 
					if (leg[sub_i].length < leg[i].length) {isSwCond2 = 1;} else {isSwCond2 = 0;}
					if (leg[sub_i].alpha > leg[i].alpha) {isSwCond3 = 1;} else {isSwCond3 = 0;}
					if (leg[i].contact  == 0) {isSwCond4 = 1;} else {isSwCond4 = 0;}
					if ((isSwCond1 == 1) && (isSwCond2 == 1)&& (isSwCond3 == 1) && (isSwCond4 == 1)) {leg[i].state = 0;leg_state[i] = 0;}   
				}
			else 
				{	
					if (leg[i].contact == 1) {isStCond1 = 1;} else {isStCond1 = 0;}
					if (isStCond1 == 1) {leg[i].state = 1; maker = i;leg_state[i] = 1;trace = i;}
				}
		}
	}
	//
		dReal DelL0, DelPsi = PI/4.0, VarPsi ;
		dReal L_sw, dL_sw_dpsi, dL_sw;
		DelL0 = L0 - L0_sw;
	for (int i = 1;i < 3;i++)
	{
		if (leg[i].state == 1)
		{
		//case 1 : // stance
			leg[i].Fs = ks * (L0 - leg[i].length) - ks/100 * leg[i].dL;
			if (leg[i].Fs < 0) {leg[i].Fs = 0;}
			leg[i].Ft = leg[i].Fs * tan(leg[i].beta);
			leg[i].Tau_ff = leg[i].Ft * leg[i].length;
			
			// zero dynamics
			
			dReal zd_Kp = 50, zd_Kd = 5, epsilon = 0.1, v; 
			
			if (enable_zd ==1)
			{
				// add desired velocity
				double eta_d,v_d_ZD = 1, k_ve = 0.1;
				eta_d = k_ve*(v_d_ZD - dxCoM)/v_d;
				v = -zd_Kp*(q7+1*PI/180.0)/pow(epsilon,2) - zd_Kd*dq1L/epsilon;
				printf("v=%f\n",v);
				leg[i].Tau_ff = leg[i].length* (v+leg[i].Ft*rh*sin(leg[i].theta))/(leg[i].length+rh*cos(leg[i].theta));
			}
			leg[i].F(0,0) = leg[i].Fs; leg[i].F(1,0) = -leg[i].Tau_ff;
			// store the GRF value:
			if (i==1)
			{
			leg[i].Fx = feedback1->f1[0];
			leg[i].Fz = feedback1->f1[2];
			leg[i].Fy = feedback1->f1[1];
			}
			if (i==2)
			{
			leg[i].Fx = feedback2->f1[0];
			leg[i].Fz = feedback2->f1[2];
			leg[i].Fy = feedback2->f1[1];
			}
		}
		else {
		//case 0 : // swing
			
		
		// prevent flight phase
				VarPsi = leg[i].theta + offset;
				L_sw = L0_sw + DelL0/2 - DelL0/2*cos(VarPsi*(2*PI)/DelPsi);
				dL_sw_dpsi = + DelL0/2*sin(VarPsi*(2*PI)/DelPsi)*(2*PI)/DelPsi;
				dL_sw =  dL_sw_dpsi*leg[i].dth;
				if ((VarPsi > DelPsi/2.0) || (VarPsi < -DelPsi/2.0))
				{
                L_sw = L0;
                dL_sw = 0;
				}
			
			leg[i].Fs = ks*(L_sw-leg[i].length) + ks/100*(dL_sw-leg[i].dL);
			leg[i].Tau_ff = 0;
			leg[i].Tau_fb = -c_swing*(wrapToPi(alpha_d - leg[i].alpha)) + 1*leg[i].dth;
			leg[i].Ft = leg[i].Tau_fb/leg[i].length;
			leg[i].F(0,0) = leg[i].Fs; leg[i].F(1,0) = -(leg[i].Tau_ff+leg[i].Tau_fb);
			// reset beta
			leg[i].beta = 0;
			leg[i].Fx = 0;
			leg[i].Fz = 0;
		}

	}


	// ->XCoM calculation Extrapolated CoM
	dReal XCoM, length_LIP, base_x = 0, base_z = 0,xi = 0;
	dReal delta_leg1_CoM = 0, delta_leg2_CoM = 0;
	if ((leg[1].state ==0) && (leg[2].state==0))
	{xi =0;}
	else
	{
	xi = xCoM + dxCoM*pow(zCoM/9.81,0.5);
	}
	for (int i=1;i<3;i++)
	{
		int sub_i = 1;
		if (i==1) {sub_i = 2;} else {sub_i = 1;}
    	if ((leg[i].state ==1) && (leg[sub_i].state == 1))
			{
				
			
				if (r_f[0] > l_f[0])
				{
					D_n = xi - r_f[0];
				}
				else {D_n = xi - l_f[0];}
			}
		else
		{
			if (leg[i].state ==1)
			{
				if (i==1) {D_n = xi - r_f[0];}
				if (i==2) {D_n = xi - l_f[0];}
			}
		}
	}
	if ((leg[1].state ==1) && (leg[2].state ==1))
	{
		if (t-pre_t >= 0.15) //delta_time_strike min = 0.15
		{
		step_length = abs(r_f[0] - l_f[0]);
		delta_Ts = t - pre_t;
		v_bar = step_length/(t-pre_t);
		//printf("pre_t:%f \n",t-pre_t);
		pre_t = t;
		}
	}
	// Torque calculation
	//***************** --> using polar calculation
	if (Minh_enable ==0)
	{
	leg[1].Tor = trans(leg[1].J) *  leg[1].F;
	leg[2].Tor = trans(leg[2].J) *  leg[2].F;
	leg[1].Tau_hip = leg[1].Tor(0,0);
	leg[1].Tau_knee = leg[1].Tor(1,0);
	leg[2].Tau_hip = leg[2].Tor(2,0);
	leg[2].Tau_knee = leg[2].Tor(3,0);
	}
	else
	{
			for(i=1;i<3;i++)
			{
				if (leg[i].state ==0) // swing
				{
					leg[i].beta = atan2(leg[i].Ft,leg[i].Fs);
				}
				leg[i].GRF = sqrt(leg[i].Ft*leg[i].Ft + leg[i].Fs*leg[i].Fs);
				if (leg[i].state ==1)
				{
				leg[i].Fn = abs(sin(leg[i].alpha - leg[i].beta)*leg[i].GRF);
				}
				else
				{
				leg[i].Fn = sin(leg[i].alpha - leg[i].beta)*leg[i].GRF;
				}
				if (leg[i].state ==0)
				{
				leg[i].Fta = leg[i].Fn/tan(leg[i].alpha-leg[i].beta);
				}
				else
				{
					if ((leg[i].alpha - leg[i].beta) > PI/2.0)
					{
						leg[i].Fta = leg[i].Fn/tan(leg[i].alpha-leg[i].beta);
					}
					else
					{
						leg[i].Fta = leg[i].Fn/tan(leg[i].alpha-leg[i].beta);
					}
				}
				leg[i].F(0,0) = leg[i].Fta;
				leg[i].F(1,0) = leg[i].Fn;
			}
			if (Minh_OSC == 0)
			{
			leg[1].Tor = trans(leg[1].J_m) *  leg[1].F;
			leg[2].Tor = trans(leg[2].J_m) *  leg[2].F;
			}
			else // lets implementation OSC in here
			{
				mat D_f(7,7), jac_c1(2,7), jac_c2(2,7),eye77(7,7), P(7,7), M_c(7,7), N(7,7), J_T_fly(4,7); // jac_c1: for right leg, jac_c2: for leftleg
		        mat J(4,7), B_f(7,7); // jac_c = c1 + c2 by column
				mat coff_matrix(7,7), jac_c;
				
				eye77 = eye(7,7);
				J = join_cols(leg[1].J_m,leg[2].J_m);

				D_f = get_D(q_f); jac_c1 = get_jacobian1(q_f), jac_c2 = get_jacobian2(q_f);

				B_f = zeros(7,7); B_f(0,0) = 1; B_f(1,1) = 1; B_f(2,2) = 1; B_f(3,3) = 1;	
				if ((leg[1].state ==1) && (leg[2].state ==1))
				{jac_c = join_cols(jac_c1,jac_c2);}
				else
				{
					if (leg[1].state ==1) {jac_c = jac_c1;}
					if (leg[2].state ==1) {jac_c = jac_c2;}
				}
				// OSC formulation
				if ((leg[1].state ==1) || (leg[2].state ==1))
				{
				P = eye77 - pinv(jac_c)*jac_c;
				M_c = P*D_f + eye77 - P;	
				J_T_fly = inv(J*inv(M_c)*P*trans(J))*J*inv(M_c)*P;
				N = eye77 - trans(J)*J_T_fly;
				coff_matrix = eye77 - N*pinv((eye77 - B_f)*N);
				if ((leg[1].state ==1) && (leg[2].state ==1)) {coff_matrix = eye77;}
				leg[1].Tor = coff_matrix*trans(J)*join_cols(leg[1].F,leg[2].F);
				leg[2].Tor = leg[1].Tor;
				}
						
			}
			leg[1].Tau_hip = leg[1].Tor(0,0);
			leg[1].Tau_knee = leg[1].Tor(1,0);
			leg[2].Tau_hip = leg[2].Tor(2,0);
			leg[2].Tau_knee = leg[2].Tor(3,0);
	}
	if ((leg[1].state ==0)&& (leg[2].state==0))
	{
			leg[1].Tau_hip = 0;
			leg[1].Tau_knee = 0;
			leg[2].Tau_hip = 0;
			leg[2].Tau_knee = 0;
	}
	//

	dReal Tor_max = 150;
	for (int j = 1; j<3; ++j)
	{
		for (int i =0;i<7;++i)
			{
				if (leg[j].Tor(i,0) >= Tor_max)
					{
						leg[j].Tor(i,0) = Tor_max;
					}
				if (leg[j].Tor(i,0) <= -Tor_max)
					{
						leg[j].Tor(i,0) = -Tor_max;
					}
			}
	}
	//
		for (int j = 1; j<3; ++j)
	{
		if (leg[j].Tau_hip > Tor_max) {leg[j].Tau_hip = Tor_max;}
		if (leg[j].Tau_hip < -Tor_max) {leg[j].Tau_hip = -Tor_max;}
		if (leg[j].Tau_knee > Tor_max) {leg[j].Tau_knee = Tor_max;}
		if (leg[j].Tau_knee < -Tor_max) {leg[j].Tau_knee = -Tor_max;}
	}
	// Add torque to joint
	dJointAddHingeTorque(joint[0],leg[1].Tau_hip);
	//dJointSetHingeParam(joint[0],dParamFMax,400);
	dJointAddHingeTorque(joint[1],leg[1].Tau_knee);
	//dJointSetHingeParam(joint[1],dParamFMax,400);
	dJointAddHingeTorque(joint[3],leg[2].Tau_hip);
	//dJointSetHingeParam(joint[3],dParamFMax,400);
	dJointAddHingeTorque(joint[4],leg[2].Tau_knee);
	//dJointSetHingeParam(joint[4],dParamFMax,400);
	

	//waist
	dsSetColor (1,1,0); // yellow
	dReal sides0[3] = {Torso_L,Torso_W,Torso_H};
	dsDrawBox (dBodyGetPosition(body[0]),dBodyGetRotation(body[0]),sides0);

	//thigh
	dReal sides1[2] = {Thigh_R,Thigh_H};
	dsSetColor (1,1,1);	dsDrawCapsule (dBodyGetPosition(body[1]),dBodyGetRotation(body[1]),Thigh_H,Thigh_R);
	dsSetColor (1,1,1);	dsDrawCapsule (dBodyGetPosition(body[4]),dBodyGetRotation(body[4]),Thigh_H,Thigh_R);
	//Shin
	dReal sides2[3] = {Shin_R,Shin_H};
	dsSetColor (0,0,1);	dsDrawCapsule (dBodyGetPosition(body[2]),dBodyGetRotation(body[2]),Shin_H,Shin_R);
	dsSetColor (0,0,1);	dsDrawCapsule (dBodyGetPosition(body[5]),dBodyGetRotation(body[5]),Shin_H,Shin_R);
	// Foot
	dReal sides3[3] = {Foot_R,Foot_H};
	dsSetColor (1,0,0);	dsDrawCapsule (dBodyGetPosition(body[3]),dBodyGetRotation(body[3]),Foot_H,Foot_R);
	dsSetColor (0,1,0);	dsDrawCapsule (dBodyGetPosition(body[6]),dBodyGetRotation(body[6]),Foot_H,Foot_R);
	// data store process
	//title
	if (excel == 1)
	{
		if ((int(t*1000) % 5) == 0)
		{
	count +=1;
	sheet->writeStr(1, 0,L"time");
	sheet->writeStr(1, 1,L"Hip_angle");sheet->writeStr(1,2,L"R_hip_angle");sheet->writeStr(1,3,L"L_hip_angle");
	sheet->writeStr(1, 4,L"R_knee_angle");sheet->writeStr(1,5,L"L_knee_angle");
	sheet->writeStr(1, 6,L"R_hip_vel");sheet->writeStr(1,7,L"L_hip_vel");sheet->writeStr(1,8,L"R_knee_vel");
	sheet->writeStr(1, 9,L"L_knee_vel");sheet->writeStr(1,10,L"hip_vel");
	sheet->writeStr(1, 15,L"rf"); sheet->writeStr(1, 16,L"lf");
	sheet->writeStr(1, 11,L"Tor_R_Th");sheet->writeStr(1, 12,L"Tor_R_Kn"); sheet->writeStr(1, 13,L"Tor_L_Th"); sheet->writeStr(1, 14,L"Tor_L_Kn");
	sheet->writeStr(1, 17,L"hpos[2]"); sheet->writeStr(1, 18,L"leg[1].length"); sheet->writeStr(1, 19,L"leg[1].beta"); 
	sheet->writeStr(1, 20,L"leg[1].alpha"); sheet->writeStr(1, 21,L"leg[1].Fs"); sheet->writeStr(1, 22,L"leg[1].Tau");
	sheet->writeStr(1, 23,L"leg[2].length") ; sheet->writeStr(1, 24,L"leg[2].beta"); 
	sheet->writeStr(1, 25,L"leg[2].alpha"); sheet->writeStr(1, 26,L"leg[2].Fs"); sheet->writeStr(1, 27,L"leg[2].Tau");
	sheet->writeStr(1, 28,L"KE"); sheet->writeStr(1, 29,L"PE"); sheet->writeStr(1, 30,L"z");
	sheet->writeStr(1, 31,L"dz");
	sheet->writeStr(1,32,L"D_n");
	sheet->writeStr(1,33,L"Forward_Vel"); 
	sheet->writeStr(1,34,L"leg[1].state"); 
	sheet->writeStr(1,35,L"leg[2].state"); 
    sheet->writeStr(1,36,L"step_length"); 
	sheet->writeStr(1,37,L"v_bar"); 
	sheet->writeStr(1,38,L"Step_time"); 
	sheet->writeStr(1,39,L"leg[1].Fx"); 
	sheet->writeStr(1,40,L"leg[1].Fz"); 
	sheet->writeStr(1,41,L"leg[2].Fx"); 
	sheet->writeStr(1,42,L"leg[2].Fz"); 
		
	// write data to excel file
	sheet->writeNum(count, 0 ,t);
	
	sheet->writeNum(count, 1,q1L*180/PI); 
	sheet->writeNum(count, 2,q31L*180/PI);
	sheet->writeNum(count, 3,q32L*180/PI); 
	sheet->writeNum(count, 4,q41L*180/PI); 
	sheet->writeNum(count, 5,q42L*180/PI); 
	//
	sheet->writeNum(count, 6,dq31L); 
	sheet->writeNum(count, 7,dq32L); 
	sheet->writeNum(count, 8,dq41L); 
	sheet->writeNum(count, 9,dq42L);
	sheet->writeNum(count, 10,dq1L);
	//
	sheet->writeNum(count, 11,leg[1].Tau_hip); 
	sheet->writeNum(count, 12,leg[1].Tau_knee); 
	sheet->writeNum(count, 13,leg[2].Tau_hip); 
	sheet->writeNum(count, 14,leg[2].Tau_knee); 
	//
	sheet->writeNum(count, 15,rf); 
	sheet->writeNum(count, 16,lf);
	//
	sheet->writeNum(count, 17,hpos[2]); 
	//
	sheet->writeNum(count, 18,leg[1].length);
	sheet->writeNum(count, 19,leg[1].beta*180/PI);
	sheet->writeNum(count, 20,leg[1].alpha*180/PI);
	sheet->writeNum(count, 21,leg[1].F(0,0));
	sheet->writeNum(count, 22,leg[1].F(1,0));
	//
	sheet->writeNum(count, 23,leg[2].length); 
	sheet->writeNum(count, 24,leg[2].beta*180/PI);
	sheet->writeNum(count, 25,leg[2].alpha*180/PI);
	sheet->writeNum(count, 26,leg[2].F(0,0));
	sheet->writeNum(count, 27,leg[2].F(1,0));
	sheet->writeNum(count, 28,KE_tot);
	sheet->writeNum(count, 29,PE_tot);
	sheet->writeNum(count, 30,zCoM);
	sheet->writeNum(count, 31,dzCoM);
	//
	sheet->writeNum(count, 32,D_n);
	
	//
	sheet->writeNum(count, 33,dxCoM);
	sheet->writeNum(count, 34,leg[1].state);
	sheet->writeNum(count, 35,leg[2].state);
	sheet->writeNum(count, 36,step_length);
	sheet->writeNum(count, 37,v_bar);
	sheet->writeNum(count,38,delta_Ts);
	//
	sheet->writeNum(count, 39,leg[1].Fx);
	sheet->writeNum(count, 40,leg[1].Fz);
	sheet->writeNum(count, 41,leg[2].Fx);
	sheet->writeNum(count, 42,leg[2].Fz);
	sheet->writeNum(count, 43,leg[1].Fy);
	sheet->writeNum(count, 44,leg[2].Fy);

	book->save(L"GRF.xls"); // change name
		}
	}

	// update state
	//q31L_pre = q31L; q41L_pre = q41L; q32L_pre = q32L; q42L_pre = q42L; q1L_pre = q1L;
	//xCoM_pre = xCoM; zCoM_pre = zCoM;
	
	
	// Incline
	for (int i =0; i < num_box;i++)
	{
		if ((i % 2) == 0)
		dsSetColor(250,250 ,210);
		else {dsSetColor(0,1,0);}
		dVector3 ss;
		dGeomBoxGetLengths(ground_box[i],ss);
		dsDrawBox (dGeomGetPosition(ground_box[i]),dGeomGetRotation(ground_box[i]),ss); 
	}
	float xyz1[3] = {0.500 + xCoM,-2.00f,0.900f};
    float hpr1[3] = {110 ,0, 0};
    dsSetViewpoint(xyz1,hpr1);
	t += DT;
	if(t>walking_time) dsStop();
	printf("time = %f \n",t);
}
//----------------------------------------------------
// Main loop
//----------------------------------------------------
int main (int argc, char **argv)
{
  leg_state[1] = 0; leg_state[2] = 0;
  int i,j,k;
  	// random parameter
	srand((unsigned)time(NULL));
  // setup pointers to drawstuff callback functions
  dsFunctions fn;
  fn.version = DS_VERSION;
  fn.start = &start;
  fn.step = simLoop;
  fn.command = 0;
  fn.stop = 0;
  fn.path_to_textures = "./drawstuff/textures";	
  dInitODE();
  // start simulation
  create_model();
  dsSimulationLoop(argc,argv,640,480,&fn);
  destroy();
  return 0;
}
