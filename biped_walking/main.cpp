#include <ode/ode.h>
#include <drawstuff/drawstuff.h>
#include <math.h>
#include <stdio.h>
#include <windows.h>
#include <time.h>
#include <stdlib.h>
#ifdef MSVC
#pragma warning(disable:4244 4305)  // for VC++, no precision loss complaints
#endif
#ifdef dDOUBLE
#define dsDrawBox dsDrawBoxD
#define dsDrawSphere dsDrawSphereD
#define dsDrawCylinder dsDrawCylinderD
#endif
//---------------------------------
//　Physic world
//---------------------------------
#define GA_SIM				// GA:0 / SIMULATION:1 
#define PI		3.14159		// PI
#define DT		0.01		// sampling time
#define GKP		30000		// eslaticity
#define GKD		30000		// viscocity

//---------------------------------
//　GA parameter
//---------------------------------
#define SEED	5		// number of species
#define GENE	80		// number of generation
#define POP		40		// number of individuals( mul 4)
#define LEN		10		// Gen parameter +1
#define MU_NUM	2		// the number of mutation
//---------------------------------
//　Control setting
//---------------------------------
#define u0			5.		// Constant into force
#define beta		2.5		// constant
#define tau0		0.12	// time constant
#define tau_dash0	1.4		// 
#define PGAIN		10		// P gain control

//---------------------------------
//　Body setting
//---------------------------------
// weight
#define Torso_Mass	8.5
#define Thigh_Mass	0.9
#define Shin_Mass	0.9
#define Foot_Mass	0.1
// size
#define Torso_R		0.05
#define Torso_L		0.50
#define Thigh_L		0.03
#define Thigh_W		0.03
#define Thigh_H		0.4
#define Shin_L		0.03
#define Shin_W		0.03
#define Shin_H		0.4
#define Foot_L		0.2
#define Foot_W		0.03
#define Foot_H		0.02	
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
static dGeomID ground_box;
// 
double t=0.;
int dflag_key2=0;
//setting of control
double param[LEN];
double tr[6], da[6];
double a[9];
double y[12], u[12], v[12];
// set of body
double mass=0;
double mass_part[4];
dQuaternion qrot1,qrot2,qrot3,qrot4,qrot5,qrot6,qrot10,qrot11;
double theta1,theta2,theta3,theta4,theta5,theta6;
double posi0[3],posi1[3],posi2[3],posi3[3],posi4[3],posi5[3],posi6[3];
//sensor
int rf=0,lf=0;
double tpos[3]={0,0,0};
double aja[7];
double ja[6], jar[6];
// setting of GA
int ipop=0, igene=0, iseed=0;
double *p[POP];
double ss[POP][LEN];
double *st;
int a2;
int a2_;
//----------------------------------------------------------------------------------
double cpg2ga[LEN-1]={-0.26,-0.21,1.00,0.69,0.61,0.80,-0.90,0.25,-0.75};	// 4.31
double feed(int i);
double f(double k){if(k<=0.0)return(0.0);else return(k);};
double h(double k){if(k<=0.0)return(0.0);else return(1.0);};
double dir(double k){if(k>=0.0)return(1.0);else return(-1.);};
//----------------------------------------------------
// Decoding GA 
//----------------------------------------------------
void parameters(double pa[]){
	int k;
	k=0;			

	a[0]=5*(pa[k]);				k++;
	a[1]=5*(pa[k]);				k++;
	a[2]=5*(pa[k]);				k++;
	a[3]=5*(pa[k]);				k++;
	a[4]=5*(pa[k]);				k++;
	a[5]=5*(pa[k]);				k++;
	a[6]=5*(pa[k]);				k++;
	a[7]=5*(pa[k]);				k++;
	a[8]=5*(pa[k]);				k++;
};

//----------------------------------------------------
// GA
//----------------------------------------------------
	void ga(int aa, int bb, double *gstring[]){

	// random parameter
	srand((unsigned)time(NULL));
	int i=0, j=0, k=0;
	while(1){						// Sort: Changing the order followed by fitness values
		if(gstring[i][0]<gstring[i+1][0]){
			st=gstring[i];
			gstring[i]=gstring[i+1];
			gstring[i+1]=st;
			i=-1;
		}
		i++;
		if(i>=aa-1)break;
	}
	a2=(int)(aa/4);
	for(i=0;i<a2;i++){
		for(j=0;j<bb;j++){
			gstring[i+1*a2][j]=gstring[i][j];			// Best-half is copied to worst-half
			gstring[i+2*a2][j]=gstring[i][j];			// Best-half is copied to worst-half
			gstring[i+3*a2][j]=gstring[i][j];			// Best-half is copied to worst-half
		}
		gstring[i+1*a2][0]=0;							// Half Fitness initialization
		gstring[i+2*a2][0]=0;							// Half Fitness initialization
		gstring[i+3*a2][0]=0;							// Half Fitness initialization
	}

	for(k=1;k<4;k++){
		for(i=k*a2;i<(k+1)*a2;i=i+2){
			j=0;
			// mutation
			for(j=0;j<MU_NUM;j++){
				gstring[i][rand()%(bb-1)+1]=((rand()%201)-100)/100.;
				gstring[i+1][rand()%(bb-1)+1]=((rand()%201)-100)/100.;
			}
		}
	}
}
//----------------------------------------------------
// Camera position
//----------------------------------------------------
static void start()
{
  static float xyz[3] = {1.50f,-3.00f,0.600f};
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
  const int N = 30;
  dContact contact[N];
  n = dCollide (o1,o2,N,&contact[0].geom,sizeof(dContact));
  if (n > 0) {
    for (i=0; i<n; i++) {
		//　Create contact when collision
		contact[i].surface.mode = dContactSoftERP | dContactSoftCFM | dContactApprox1;
		contact[i].surface.mu = 1;
		contact[i].surface.soft_erp = DT*GKP/(DT*GKP+GKD);
		contact[i].surface.soft_cfm = 1/(DT*GKP+GKD);
		dJointID c = dJointCreateContact (world,contactgroup,&contact[i]);
		dJointAttach (c,
			dGeomGetBody(contact[i].geom.g1),
			dGeomGetBody(contact[i].geom.g2));
	}
  }
   
  //　contact information of right foot
  if(o1i==(int)right && o2i==(int)ground)
	  rf=n;		
  else if(o2i==(int)right && o1i==(int)ground)
	  rf=n;		
  //　contact information of left foot
  if(o1i==(int)left && o2i==(int)ground)
	  lf=n;		
  else if(o2i==(int)left && o1i==(int)ground)
	  lf=n;	
}

//----------------------------------------------------
//　Initial setting of robot
//----------------------------------------------------
static void reset(void){
  int i;
  dMass m;

  //-------------------
  //　Initial setting of physical world
  //-------------------
  world = dWorldCreate();
  space = dHashSpaceCreate (0);
  contactgroup = dJointGroupCreate (0);
  dWorldSetGravity (world,0,0,-9.8);
  dWorldSetCFM (world,1e-5);
  dWorldSetERP (world,0.2);
  ground = dCreatePlane (space,0,0,1,0);
  //-------------------
  //　setting mass
  //-------------------
  mass_part[0]=Torso_Mass;
  mass_part[1]=Thigh_Mass;	
  mass_part[2]=Shin_Mass;
  mass_part[3]=Foot_Mass;
  // Setting of initial position of each object
  posi0[0]=0; 
  posi0[1]=0; 
  posi0[2]=0.9+ 0.5*Torso_L;

  posi1[0]=posi0[0]; 
  posi1[1]=0; 
  posi1[2]=posi0[2]-0.5*Thigh_H - 0.5*Torso_L; 

  posi2[0]=posi1[0];
  posi2[1]=0; 
  posi2[2]=posi1[2]-0.5*Thigh_H-0.5*Shin_H; 

  posi3[0]=posi2[0];
  posi3[1]=0; 
  posi3[2]=posi2[2]-0.5*Shin_H-0.5*Foot_H; 

  posi4[0]=posi0[0]; 
  posi4[1]=0; 
  posi4[2]=posi0[2]-0.5*Thigh_H -0.5*Torso_L; 

  posi5[0]=posi4[0];
  posi5[1]=0;
  posi5[2]=posi4[2]-0.5*Thigh_H-0.5*Shin_H; 

  posi6[0]=posi5[0];
  posi6[1]=0; 
  posi6[2]=posi5[2]-0.5*Shin_H-0.5*Foot_H; 
  // Object configuration
  body[0] = dBodyCreate (world);
  dBodySetPosition (body[0],posi0[0],posi0[1],posi0[2]);
  mass=mass_part[0];
  dMassSetZero(&m); //
  //dMassSetParameters(&m,10,0,0,0,0,0,0,0,0,1.05);
  dMassSetCapsuleTotal(&m,mass,3,Torso_R,Torso_L);
  dBodySetMass(body[0], &m);
  box[0] = dCreateCapsule(0,Torso_R,Torso_L);
  dGeomSetBody (box[0],body[0]);
  // Box1: right thigh
  body[1] = dBodyCreate (world);
  dBodySetPosition (body[1],posi1[0],posi1[1],posi1[2]);
  dMassSetBox (&m,1,Thigh_L,Thigh_W,Thigh_H);
  mass=mass_part[1];
  dMassAdjust (&m,mass);
  dBodySetMass (body[1],&m);
  box[1] = dCreateBox (0,Thigh_L,Thigh_W,Thigh_H);
  dGeomSetBody (box[1],body[1]);
  // Box2: right shin
  body[2] = dBodyCreate (world);
  dBodySetPosition (body[2],posi2[0],posi2[1],posi2[2]);
  dMassSetBox (&m,1,Shin_L,Shin_W,Shin_H);
  mass=mass_part[2];
  dMassAdjust (&m,mass);
  dBodySetMass (body[2],&m);
  box[2] = dCreateBox (0,Shin_L,Shin_W,Shin_H);
  dGeomSetBody (box[2],body[2]);
  // box 3: right foot
  body[3] = dBodyCreate (world);
  dBodySetPosition (body[3],posi3[0],posi3[1],posi3[2]);
  dMassSetBox (&m,1,Foot_L,Foot_W,Foot_H);
  mass=mass_part[3];
  dMassAdjust (&m,mass);
  dBodySetMass (body[3],&m);
  box[3] = dCreateBox (0,Foot_L,Foot_W,Foot_H);
  dGeomSetBody (box[3],body[3]);
  // Box 4: left thigh
  body[4] = dBodyCreate (world);
  dBodySetPosition (body[4],posi4[0],posi4[1],posi4[2]);
  dMassSetBox (&m,1,Thigh_L,Thigh_W,Thigh_H);
  mass=mass_part[1];
  dMassAdjust (&m,mass);
  dBodySetMass (body[4],&m);
  box[4] = dCreateBox (0,Thigh_L,Thigh_W,Thigh_H);
  dGeomSetBody (box[4],body[4]);
  // box 5:left shin
  body[5] = dBodyCreate (world);
  dBodySetPosition (body[5],posi5[0],posi5[1],posi5[2]);
  dMassSetBox (&m,1,Shin_L,Shin_W,Shin_H);
  mass=mass_part[2];
  dMassAdjust (&m,mass);
  dBodySetMass (body[5],&m);
  box[5] = dCreateBox (0,Shin_L,Shin_W,Shin_H);
  dGeomSetBody (box[5],body[5]);
  // box 6: left foot
  body[6] = dBodyCreate (world);
  dBodySetPosition (body[6],posi6[0],posi6[1],posi6[2]);
  dMassSetBox (&m,1,Foot_L,Foot_W,Foot_H);
  mass=mass_part[3];
  dMassAdjust (&m,mass);
  dBodySetMass (body[6],&m);
  box[6] = dCreateBox (0,Foot_L,Foot_W,Foot_H);
  dGeomSetBody (box[6],body[6]);
  //　Joint setting
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
  joint[2] = dJointCreateHinge (world,0);
  dJointAttach (joint[2],body[2],body[3]);
  const dReal *a3 = dBodyGetPosition (body[3]);
  dJointSetHingeAnchor (joint[2],a2[0],a3[1],a3[2]+0.5*Foot_H);
  dJointSetHingeAxis (joint[2],0,1,0);
  dJointSetHingeParam (joint[2],dParamLoStop,-dInfinity);
  dJointSetHingeParam (joint[2],dParamHiStop,dInfinity);
  dJointSetHingeParam (joint[2],dParamFudgeFactor,0.1);
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
  joint[5] = dJointCreateHinge (world,0);
  dJointAttach (joint[5],body[5],body[6]);
  const dReal *b3 = dBodyGetPosition (body[6]);
  dJointSetHingeAnchor (joint[5],b2[0],b3[1],b3[2]+0.5*Foot_H);
  dJointSetHingeAxis (joint[5],0,1,0);
  dJointSetHingeParam (joint[5],dParamLoStop,-dInfinity);
  dJointSetHingeParam (joint[5],dParamHiStop,dInfinity);
  dJointSetHingeParam (joint[5],dParamFudgeFactor,0.1);
  // geometry grouping
  geom_group = dSimpleSpaceCreate (space);  
  dSpaceSetCleanup (geom_group,0);
  dSpaceAdd   (geom_group,box[0]);
  dSpaceAdd   (geom_group,box[1]);
  dSpaceAdd   (geom_group,box[4]);
  dSpaceAdd   (geom_group,box[2]);
  dSpaceAdd   (geom_group,box[5]);
  //--> space ID for foot right and left
  right = dSimpleSpaceCreate (space);  
  dSpaceSetCleanup (right,0);
  dSpaceAdd   (right,box[3]);
  left = dSimpleSpaceCreate (space);  
  dSpaceSetCleanup (left,0);
  dSpaceAdd   (left,box[6]);
  //　Orther, intialization parameter
  for(i=0;i<6;i++){
	  tr[i]=0;
	  ja[i]=0;
	  jar[i]=0;
	  da[i]=0;
  }
  for(i=0;i<7;i++){
	  aja[i]=0;
  }
  for(i=0;i<12;i++){
		y[i]=0;
		u[i]=0;
		v[i]=0;
   }
  rf=0; lf=0; // sensor to detect wherether leg touch the ground
  tpos[0]=0; 
  tpos[1]=0; 
  tpos[2]=0;
  t=0;
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
	dSpaceDestroy (space);
	dWorldDestroy (world);
}

//-----------------------------
//　Sensor feedback function
//-----------------------------
double feed(int i){
	double feed0=0.0;
	switch(i){
	//Move the thigh joint rules: their own angle, reserve angle thigh, shine angle = sensor feedback system
	case 0:		feed0=(a[0]*aja[1]+a[1]*aja[4]+a[2]*aja[2]*h(rf)+a[3]*h(lf));		return(feed0);
	case 1:		feed0=(a[0]*aja[1]+a[1]*aja[4]+a[2]*aja[2]*h(rf)+a[3]*h(lf))*(-1);	return(feed0);
	case 2:		feed0=(a[0]*aja[4]+a[1]*aja[1]+a[2]*aja[5]*h(lf)+a[3]*h(rf));		return(feed0);
	case 3:		feed0=(a[0]*aja[4]+a[1]*aja[1]+a[2]*aja[5]*h(lf)+a[3]*h(rf))*(-1);	return(feed0);	
   // In response to the reverse pendulum state of reverse thigh, bend rules
	case 4:		feed0=(a[4]*aja[4])*h(lf);				return(feed0);		
	case 5:		feed0=(a[4]*aja[4])*h(lf)*(-1);			return(feed0);		
	case 6:		feed0=(a[4]*aja[1])*h(rf);				return(feed0);		
	case 7:		feed0=(a[4]*aja[1])*h(rf)*(-1);			return(feed0);		
    // Adding counter torque when the torque not in contact
	case 8:		feed0=(a[5]*aja[3]*h(rf)+a[6]*aja[6]*h(lf)+a[7]*h(rf));		return(feed0);
	case 9:		feed0=(a[5]*aja[3]*h(rf)+a[6]*aja[6]*h(lf)+a[7]*h(rf))*(-1);	return(feed0);
	case 10:	feed0=(a[5]*aja[6]*h(lf)+a[6]*aja[3]*h(rf)+a[7]*h(lf));		return(feed0);
	case 11:	feed0=(a[5]*aja[6]*h(lf)+a[6]*aja[3]*h(rf)+a[7]*h(lf))*(-1);	return(feed0);
	}
	return(feed0);
}
//-----------------------------
//　CPG
//-----------------------------
void cpgp(void){
	int i,j;
	double ww;
	double tau[12]={tau0,tau0,tau0,tau0,tau0,tau0,tau0,tau0,tau0,tau0,tau0,tau0};
	double tau_dash[12]={tau_dash0,tau_dash0,tau_dash0,tau_dash0,tau_dash0,tau_dash0,tau_dash0,tau_dash0,tau_dash0,tau_dash0,tau_dash0,tau_dash0};
	//-->From G.Taga 1991
	double w[12][12]={0.,-2.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
					 -2., 0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 0.,
					 -1., 0., 0.,-2., 0., 0., 0., 0., 0., 0., 0., 0.,
					  0.,-1.,-2., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
					  0., 0., 0., 0., 0.,-2., 0., 0., 0., 0., 0., 0.,
					 -1.,-1., 0., 0.,-2., 0., 0., 0., 0., 0., 0., 0.,
					  0., 0., 0., 0., 0., 0., 0.,-2., 0., 0., 0., 0.,
					  0., 0.,-1.,-1., 0., 0.,-2., 0., 0., 0., 0., 0.,
					 -1., 0., 0., 0., 0., 0., 0., 0., 0.,-2., 0., 0.,
					  0.,-1., 0., 0., 0., 0., 0., 0.,-2., 0., 0., 0.,
					  0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 0.,-2.,
					  0., 0., 0.,-1., 0., 0., 0., 0., 0., 0.,-2., 0.};	
	// CPG equation = non- linear second-order differential equation
	for(i=0;i<12;i++){
		ww=0.0;
		for(j=0;j<12;j++)ww=ww+w[i][j]*f(u[j]);
		u[i]=u[i]+DT*(-u[i]+ww-beta*v[i]+u0+feed(i))/tau[i];
		v[i]=v[i]+DT*(-v[i]+f(u[i]))/tau_dash[i];
	}
	for(i=0;i<12;i++)y[i]=f(u[i]);
   //　It calculates the target angle of each joint from equation
	da[0]=PI/6/5*(y[0]-y[1]);			
	da[1]=PI/2/5*(y[4]-y[5]);	
	if(da[1]>0)da[1]=0;	
	da[2]=PI/6/5*(y[8]-y[9]);	
	da[3]=PI/6/5*(y[2]-y[3]);			
	da[4]=PI/2/5*(y[6]-y[7]);	
	if(da[4]>0)da[4]=0;	
	da[5]=PI/6/5*(y[10]-y[11]);
}
//----------------------------------------------------
// Simulation loop
//----------------------------------------------------
static void simLoop (int pause)
{
	int i;
	// Sensor capture
	// Joint angle
	ja[0]=dJointGetHingeAngle (joint[0]);
	ja[1]=dJointGetHingeAngle (joint[1]);
	ja[2]=dJointGetHingeAngle (joint[2]);	
	ja[3]=dJointGetHingeAngle (joint[3]);
	ja[4]=dJointGetHingeAngle (joint[4]);	
	ja[5]=dJointGetHingeAngle (joint[5]);	
	// Joint angular velocity
	jar[0]=dJointGetHingeAngleRate (joint[0]);
	jar[1]=dJointGetHingeAngleRate (joint[1]);
	jar[2]=dJointGetHingeAngleRate (joint[2]);	
	jar[3]=dJointGetHingeAngleRate (joint[3]);
	jar[4]=dJointGetHingeAngleRate (joint[4]);	
	jar[5]=dJointGetHingeAngleRate (joint[5]);
	// absolute position of each object
	const dReal *pos0= dBodyGetPosition(body[0]);
	const dReal *pos1= dBodyGetPosition(body[1]);
	const dReal *pos2= dBodyGetPosition(body[2]);
	const dReal *pos3= dBodyGetPosition(body[3]);
	const dReal *pos4= dBodyGetPosition(body[4]);
	const dReal *pos5= dBodyGetPosition(body[5]);
	const dReal *pos6= dBodyGetPosition(body[6]);
	//absolute angle of each object
	const dReal *qua0= dBodyGetQuaternion(body[0]);
	const dReal *qua1= dBodyGetQuaternion(body[1]);
	const dReal *qua2= dBodyGetQuaternion(body[2]);
	const dReal *qua3= dBodyGetQuaternion(body[3]);
	const dReal *qua4= dBodyGetQuaternion(body[4]);
	const dReal *qua5= dBodyGetQuaternion(body[5]);
	const dReal *qua6= dBodyGetQuaternion(body[6]);
	// absolute angle in XZ plane
	aja[0]=fabs(2*acos(qua0[0]))*dir(qua0[2]);
	aja[1]=fabs(2*acos(qua1[0]))*dir(qua1[2]);
	aja[2]=fabs(2*acos(qua2[0]))*dir(qua2[2]);
	aja[3]=fabs(2*acos(qua3[0]))*dir(qua3[2]);
	aja[4]=fabs(2*acos(qua4[0]))*dir(qua4[2]);
	aja[5]=fabs(2*acos(qua5[0]))*dir(qua5[2]);
	aja[6]=fabs(2*acos(qua6[0]))*dir(qua6[2]);
	// Absolute coordinate position of the waist position (HIP)S
	tpos[0]=pos0[0];
	tpos[1]=pos0[1];
	tpos[2]=pos0[2];
	//if(t==0){
		//printf("%2.4f	%2.4f	%2.4f	%2.4f\n",t,ja[0],ja[1],ja[2]);
	//}
	//Initialization
	for(i=0;i<6;i++){
		tr[i]=0;
		da[i]=0;
	}
	// CPG calculation
	cpgp();
	// Position control (PI control)
	tr[0]=PGAIN*(da[0]-(ja[0]));	
	tr[1]=PGAIN*(da[1]-ja[1]);	
	tr[2]=PGAIN*(da[2]-ja[2]);	
	tr[3]=PGAIN*(da[3]-(ja[3]));	
	tr[4]=PGAIN*(da[4]-ja[4]);
	tr[5]=PGAIN*(da[5]-ja[5]);

	//Motor command
	for(i=0;i<6;i++){
		dJointSetHingeParam (joint[i],dParamVel,tr[i]);
		dJointSetHingeParam (joint[i],dParamFMax,500.0);
	}
	// Dynamics calculation and 1 stap elapsed
	rf=0;lf=0;  
    dSpaceCollide (space,0,&nearCallback);	// Computation of Collision Detect
    dWorldStep (world,DT);					// Sampling time  
    dJointGroupEmpty (contactgroup);		// remove all contact joints
	t=t+DT;
	/*
	// Drawing
	if(GA_SIM==0){
		// Drawing off
		if(dflag_key2==0){
			keybd_event(VK_CONTROL,0,0,0);
			keybd_event('T',0,0,0);
			keybd_event('T',0, KEYEVENTF_KEYUP, 0);
			keybd_event('S',0, 0, 0);
			keybd_event('S',0, KEYEVENTF_KEYUP, 0);
			keybd_event(VK_CONTROL,0, KEYEVENTF_KEYUP, 0);
			dflag_key2=1;
		}

		//　Falling and exit conditions of individual
		if(tpos[2]<0.7){
			p[ipop][0]=tpos[0];
			dsStop();
		}
		if(tpos[0]<-0.6){
			p[ipop][0]=-1;
			dsStop();
		}

		if(t>6){
			p[ipop][0]=tpos[0];
			dsStop();
		}
	}
	else {
		if(t>20)dsStop();
	}
	*/
	//waist
	dsSetColor (1,0,0);
	dsDrawCapsuleD (dBodyGetPosition(body[0]),dBodyGetRotation(body[0]),Torso_L,Torso_R);
	//thigh
	dReal sides1[3] = {Thigh_L,Thigh_W,Thigh_H};
	dsSetColor (1,1,1);	dsDrawBox (dBodyGetPosition(body[1]),dBodyGetRotation(body[1]),sides1);
	dsSetColor (1,1,1);	dsDrawBox (dBodyGetPosition(body[4]),dBodyGetRotation(body[4]),sides1);
	//Shin
	dReal sides2[3] = {Shin_L,Shin_W,Shin_H};
	dsSetColor (0,0,1);	dsDrawBox (dBodyGetPosition(body[2]),dBodyGetRotation(body[2]),sides2);
	dsSetColor (0,0,1);	dsDrawBox (dBodyGetPosition(body[5]),dBodyGetRotation(body[5]),sides2);
	// Foot
	dReal sides3[3] = {Foot_L,Foot_W,Foot_H};
	dsSetColor (0,1,0);	dsDrawBox (dBodyGetPosition(body[3]),dBodyGetRotation(body[3]),sides3);
	dsSetColor (0,1,0);	dsDrawBox (dBodyGetPosition(body[6]),dBodyGetRotation(body[6]),sides3);
}
//----------------------------------------------------
// Main loop
//----------------------------------------------------
int main (int argc, char **argv)
{
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
  fn.path_to_textures = "E:/ODE_seven_links_biped/biped_walking/ODE_library/drawstuff/textures";	
  dInitODE();
  //search GA
  /*
  if(GA_SIM==0) {
	  
	  //　start of species
		for(iseed=0;iseed<SEED;iseed++){

			for(i=0;i<POP;i++)
				for(j=0;j<LEN;j++)
					ss[i][j]=((rand()%201)-100)/100.;

			for(i=0;i<POP;i++){
				//ss[i][0]=0;
				p[i]=ss[i];
			}

			//start of generation
			for(igene=0;igene<GENE;igene++){
				
				//start of invididual
				for(ipop=0;ipop<POP;ipop++){
					if(ipop==0 && igene!=0){ipop=POP/4;}

						//　variable setting
						reset();

						for(k=0;k<LEN-1;k++){param[k]=p[ipop][k+1];}
						parameters(param);					
						p[ipop][0]=0;
						dsSimulationLoop (argc,argv,640,480,&fn);
						destroy();			
				}
				//　the end of individual
				
				ga(POP,LEN,p);

				printf("	%d	%d	%2.2f	\n",iseed,igene,p[0][0]);
				fp=fopen("gene_data.txt","a");
				if(fp==NULL){
					printf("failure\n");
					exit(-1);
				}	

				// when the fitness is 3.5 or more, kill the species
				if(p[0][0]>3.5){igene=GENE;}
			}
			//　the end of generation
			fp=fopen("seed_gene.txt","a");
				if(fp==NULL){
				printf("failure\n");
				exit(-1);
			}	

			fprintf(fp,"//double cpg2ga[LEN-1]={");
			for(i=1;i<LEN;i++){
				fprintf(fp,"%2.2f",p[0][i]);
				if(i!=LEN-1)fprintf(fp,",");	
			}
			fprintf(fp,"};	// %2.2f\n",p[0][0]);
			fclose(fp);
		}
		//　the end of species
		
	}
  */
  

  // start simulation
  reset();
  parameters(cpg2ga);
  dsSimulationLoop (argc,argv,640,480,&fn);
  destroy();
  return 0;
}
