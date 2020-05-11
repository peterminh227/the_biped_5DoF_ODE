//BTSLIP
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
using namespace libxl;

//#include <PGL/Pgl/PGL.h>

#define NUM 2
#ifdef  dDOUBLE // adaptation for single and double precision
#define dsDrawSphere  dsDrawSphereD
#define dsDrawCapsule dsDrawCapsuleD
#define dsDrawCylinder dsDrawCylinderD
#define PI 3.14159265
#define DT		0.005 //sample time
#define SKP  13000 // spring cofficient
#define SKD  0
#endif
static double THETA[NUM] = { 0.0, 0.0}; // target angle of joints; 
static dWorldID world;
static dSpaceID space;
static dSpaceID legbottom;
static dGeomID  ground;
static dJointGroupID contactgroup;
static double count = 1, Collide_index, count1=0;
static int n; //Collision detection
dsFunctions fn;
static dReal leg_r = 0.02, leg_l = 0.5;
static dReal leng_l1,leng_l2, angle1, angle2, torque1, torque2, ks= 13000, f[2], t=0, off_set = 0;
static dReal zeta1_tilde, zeta2_tilde, beta1_tilde, beta2_tilde, phi_tilde,phi_tilde_pre = 0, vel_phi_tilde;
static dVector3 l1,l2,f1_VPP,f2_VPP, CoM_f1, CoM_f2, trunk_ori, vertical_vec;
dReal timeS = 0.0;
static Book* book = xlCreateBook(); // xlCreateXMLBook() for xlsx
static Sheet* sheet = book->addSheet(L"Sheet1");
static const dReal alpha = 69;
dJointFeedback *feedback_slider0 = new dJointFeedback, *feedback_slider1 = new dJointFeedback;
typedef struct {
    dBodyID body;
    dGeomID geom;
    dReal   radius;
    dReal   length;
    dReal   mass;
} myLink;

myLink leg1[2],leg2[2],torso;

dJointID joint[NUM], jointslider[NUM];


static void nearCallback(void *data, dGeomID o1, dGeomID o2)
{
    const int N = 100;
    dContact contact[N];
	int isGround = ((ground == o1) || (ground == o2));
	dBodyID b1 = dGeomGetBody(o1);
	dBodyID b2 = dGeomGetBody(o2);
	if (b1 && b2 && dAreConnectedExcluding(b1,b2,dJointTypeContact)) return;
	n =  dCollide(o1,o2,N,&contact[0].geom,sizeof(dContact));
	std ::string s;
	std::stringstream out;
	out << n;
	s = out.str();
	Collide_index = ::atof(s.c_str());
	if (isGround)
    {
        for (int i = 0; i < n; i++)
        {
            
			contact[i].surface.mode = dContactSoftERP|dContactSoftCFM|dContactApprox1|dContactBounce;
            contact[i].surface.mu   = dInfinity; // ;
            contact[i].surface.bounce = 0; // (0.0~1.0)
			//contact[i].surface.bounce_vel = 0.0;
			// 0.97-1e-3
            contact[i].surface.soft_erp = 0.999;
            contact[i].surface.soft_cfm = 1e-3;
            dJointID c = dJointCreateContact(world,contactgroup,&contact[i]);
            dJointAttach (c,dGeomGetBody(contact[i].geom.g1),
           dGeomGetBody(contact[i].geom.g2));
						  
	
        }
    }
	// Collide avoid  --> setting CategoryBits and CollideBits

}
void control()    /***  p control; P  ****/
{
    static int step = 0;     // step number; 
    double k1 =  10.0,  fMax  = 100.0; // k1:gain, fMax: max torqeu; 
    printf("\r%6d:",step++);
    for (int j = 1; j <NUM; j++)
    {
        double tmpAngle = dJointGetHingeAngle(joint[j]);  // current joint angle;
        double z = THETA[j] - tmpAngle;  // z= target - current; z: 
        dJointSetHingeParam(joint[j],  dParamVel,  k1*z); // angular velocity; 
        dJointSetHingeParam(joint[j], dParamFMax, fMax); // max torqeu;
    }
}
void command(int cmd)   /***  key control function; ***/
{
    printf("\n control key");
	switch (cmd)
    {
    case 'j':
        THETA[1] += 0.05; //j
        // increases THETA[1] when j key is pressed
        break;
    case 'f':
        THETA[1] -= 0.05;
        break;
    case 'k':
        THETA[2] += 0.05;
        break;
    case 'd':
        THETA[2] -= 0.05;
        break;
   
    }

    // English: limit target angles not to destroy a robot
    // Japanese:
    if (THETA[1] <  - M_PI)    THETA[1] =  - M_PI; // M_PI：
    if (THETA[1] >    M_PI)    THETA[1] =    M_PI;
    if (THETA[2] < -2*M_PI/3)  THETA[2] =  - 2*M_PI/3;
    if (THETA[2] >  2*M_PI/3)  THETA[2] =    2*M_PI/3;
   
}
void controlSlider(dReal target)
{
  static dReal kp   = 13000;                       
  static dReal fmax = 6500;                        

  dReal tmp1  = dJointGetSliderPosition(jointslider[0]);  
  dReal tmp2  = dJointGetSliderPosition(jointslider[1]);  
  //printf("\n control key %2.2f",tmp);
  dReal u1    = -kp * (target - tmp1);               
  dReal u2    = -kp * (target -tmp2);
  dJointSetSliderParam(jointslider[0], dParamVel,  u1);
  dJointSetSliderParam(jointslider[1], dParamVel,  u2);
  dJointSetSliderParam(jointslider[0], dParamFMax, fmax);
  dJointSetSliderParam(jointslider[1], dParamFMax, fmax);
}

int get_state(dReal state[5][3])
{
	dVector3 foot1_p,foot2_p,hip_position,VPP_position,COM_position;
	dBodyGetRelPointPos(leg1[0].body,0,0,-leg_l/2-leg_r,foot1_p);
	dBodyGetRelPointPos(leg2[0].body,0,0,-leg_l/2-leg_r,foot2_p);
	dBodyGetRelPointPos(torso.body,0,0,-0.1,hip_position);
	dBodyGetRelPointPos(torso.body,0,0,0,COM_position);
	dBodyGetRelPointPos(torso.body,0,0,0.1,VPP_position);
	state[0][0] = foot1_p[0];
	state[0][1] = foot1_p[1];
	state[0][2] = foot1_p[2];
	state[1][0] = foot2_p[0];
	state[1][1] = foot2_p[1];
	state[1][2] = foot2_p[2];
	state[2][0] = hip_position[0];
	state[2][1] = hip_position[1];
	state[2][2] = hip_position[2];
	state[3][0] = COM_position[0];
	state[3][1] = COM_position[1];
	state[3][2] = COM_position[2];
	state[4][0] = VPP_position[0];
	state[4][1] = VPP_position[1];
	state[4][2] = VPP_position[2];
	
	return 1;
}
double getRotateAngle(dVector3 vec1, dVector3 vec2)
{
	const double epsilon = 1.0e-6;
	double angle = 0;
	// normalize
	dVector3 norVec1, norVec2;
	norVec1[0] = vec1[0] / sqrt(pow(vec1[0],2)+ pow(vec1[1],2) + pow(vec1[2],2));
	norVec1[1] = vec1[1] / sqrt(pow(vec1[0],2)+ pow(vec1[1],2) + pow(vec1[2],2));
	norVec1[2] = vec1[2] / sqrt(pow(vec1[0],2)+ pow(vec1[1],2) + pow(vec1[2],2));
	norVec2[0] = vec2[0] / sqrt(pow(vec2[0],2)+ pow(vec2[1],2) + pow(vec2[2],2));
	norVec2[1] = vec2[1] / sqrt(pow(vec2[0],2)+ pow(vec2[1],2) + pow(vec2[2],2));
	norVec2[2] = vec2[2] / sqrt(pow(vec2[0],2)+ pow(vec2[1],2) + pow(vec2[2],2));
	// dot product
	double dotProd = (norVec1[0]*norVec2[0] + norVec1[1]*norVec2[1] + norVec1[2]*norVec2[2]);
	if ( abs(dotProd - 1.0) <= epsilon )
	angle = 0;
	else if ( abs(dotProd + 1.0) <= epsilon )
	angle = PI;
	else {
		double cross_x = 0;
		angle = acos(dotProd);
		cross_x = (norVec1[1]*norVec2[2]-norVec1[2]*norVec2[1]);
		if (cross_x < 0) // vec1 rotate clockwise to vec2
		angle = 2 * PI - angle; 

	}
	return angle;
}
double cal_springforce(dReal f[2])
{
	
	dReal tmp2  = dJointGetSliderPosition(jointslider[0]); 
	dReal tmp3 = dJointGetSliderPosition(jointslider[1]);
	const dReal *p1 = dBodyGetPosition(leg1[0].body);
	const dReal *p2 = dBodyGetPosition(leg1[1].body);
	const dReal *p3 = dBodyGetPosition(leg2[0].body);
	const dReal *p4 = dBodyGetPosition(leg2[1].body);
	dReal delta_x1 = sqrt(pow(p2[2] - p1[2],2)+pow(p2[1] - p1[1],2)+pow(p2[0] - p1[0],2)) - 0.5;
	dReal delta_x2 = sqrt(pow(p4[2] - p3[2],2)+pow(p4[1] - p3[1],2)+pow(p4[0] - p3[0],2)) - 0.5;
	f[0] = ks*(-tmp2);
	f[1] = ks*(-tmp3);
	//if (f[0] < 0) f[0] = 0;
	//if (f[1] < 0) f[1] = 0;
	//f[0] = ks*(-delta_x1);
	//f[1] = ks*(-delta_x2);
	return 1;
}
static void simLoop (int pause) {
    const dReal *pos[5],*R[5];
	t = t+ DT;
	dReal state[5][3];
	get_state(state);
	for ( int  i=0; i<3; i++) {
	l1[i] = state[2][i] - state[0][i];
	l2[i] = state[2][i] - state[1][i];
	f1_VPP[i] = state[4][i] - state[0][i];
	f2_VPP[i] = state[4][i] - state[1][i];
	CoM_f1[i] = state[3][i] - state[0][i]; // f1_CoM vector
	CoM_f2[i] = state[3][i] - state[1][i]; // f2_CoM vector
	trunk_ori[i] = state[3][i] - state[2][i]; 
	}
	leng_l1 = sqrt ( pow(l1[0],2)+pow(l1[1],2)+pow(l1[2],2));
	leng_l2 = sqrt ( pow(l2[0],2)+pow(l2[1],2)+pow(l2[2],2));
	//--> calculate the VPP torque f*L*VPPxL/(VPP.L)
	zeta1_tilde = getRotateAngle(CoM_f1,l1);
	zeta2_tilde = getRotateAngle(CoM_f2,l2);
	vertical_vec[0] = 0; vertical_vec[1] = 1; vertical_vec[2] = 0;
	phi_tilde = getRotateAngle(vertical_vec,trunk_ori);

	if (t == 0) vel_phi_tilde = 0;
	if (t > 0 ) vel_phi_tilde = (phi_tilde - phi_tilde_pre)/DT;
	//printf("velocity of trunk = %f \n",vel_phi_tilde);
	angle1 = zeta1_tilde + (-2)*(phi_tilde - off_set) - 1*vel_phi_tilde; 
	angle2 = zeta2_tilde + (-2)*(phi_tilde - off_set) - 1*vel_phi_tilde;
	//angle1 = getRotateAngle(f1_VPP,l1);
	//angle2 = getRotateAngle(f2_VPP,l2);
	
	dContact contact_temp1[20], contact_temp2[20];
	// finite state machine 
	int n1 = dCollide(leg1[0].geom,ground,20,&contact_temp1[0].geom,sizeof(dContact));
	int n2 = dCollide(leg2[0].geom,ground,20,&contact_temp2[0].geom,sizeof(dContact));
	if (state[1][2] < 1e-4) 
	{//n2 =1;
	}
	if (state[0][2] < 1e-4)
	{//n1=1;
	}
	
	if ((n1==0) && (n2==1))
	{	
		
		// leg1 goes under surface -- an airbone
		dGeomDisable(leg1[0].geom);
		//f1=0;
		dGeomEnable(leg2[0].geom);
		dJointAddHingeTorque(joint[0],torque1);
		dJointAddHingeTorque(joint[1],torque2);
		dBodySetPosition(leg1[1].body,0,state[2][1]+ 0.5*0.5*cos(alpha*PI/180),state[2][2]-0.5*0.5*sin(alpha*PI/180));
		dBodySetPosition(leg1[0].body,0,state[2][1]+ 0.75*cos(alpha*PI/180),state[2][2]-0.75*sin(alpha*PI/180));
		cal_springforce(f);
		f[0] = 0;
		dJointAddSliderForce(jointslider[0],f[0]);
		dJointAddSliderForce(jointslider[1],f[1]);
		torque1 = 0;
		torque2 = -f[1]*leng_l2*tan(angle2);
		dJointAddHingeTorque(joint[0],torque1);
		dJointAddHingeTorque(joint[1],torque2);

	}
	if ((n1==1) && (n2==1))
	{
		dGeomEnable(leg1[0].geom);
		dGeomEnable(leg2[0].geom);
		cal_springforce(f);
		dJointAddSliderForce(jointslider[0],f[0]);
		dJointAddSliderForce(jointslider[1],f[1]);
		torque1 = -f[0]*leng_l1*tan(angle1);
		torque2 = -f[1]*leng_l2*tan(angle2);
		dJointAddHingeTorque(joint[0],torque1);
		dJointAddHingeTorque(joint[1],torque2);
	}
	if ((n1==1) && (n2==0))
	{
		// leg1 goes under surface -- an airbone
		dGeomDisable(leg2[0].geom);
		dGeomEnable(leg1[0].geom);
		dBodySetPosition(leg2[1].body,0,state[2][1]+ 0.5*0.5*cos(alpha*PI/180),state[2][2]-0.5*0.5*sin(alpha*PI/180));
		dBodySetPosition(leg2[0].body,0,state[2][1]+ 0.75*cos(alpha*PI/180),state[2][2]-0.75*sin(alpha*PI/180));
		cal_springforce(f);
		f[1] =0;
		dJointAddSliderForce(jointslider[0],f[0]);
		dJointAddSliderForce(jointslider[1],f[1]);
		torque1 = -f[0]*leng_l1*tan(angle1);
		torque2 = 0;
		dJointAddHingeTorque(joint[0],torque1);
		dJointAddHingeTorque(joint[1],torque2);
	}
	if ((n1==0) && (n2==0))
	{
		// leg1 goes under surface -- an airbone
		dGeomEnable(leg2[0].geom);
		dGeomEnable(leg1[0].geom);
	}
	// when it fall down - turn off motor
	if (state[3][2] < state[2][2]) // CoM vertical < Hip vervical 
	{
		dGeomEnable(leg2[0].geom);
		dGeomEnable(leg1[0].geom);
		dJointAddSliderForce(jointslider[0],0);
		dJointAddSliderForce(jointslider[1],0);
		dJointAddHingeTorque(joint[0],0);
		dJointAddHingeTorque(joint[1],0);
		printf("\n Robot fall down at time %f",t);
	
	}
	dSpaceCollide(space,0,&nearCallback);
	dWorldStep(world,DT);
	//add sheet to excel
	//printf("\n export dCollide %d",Collide_index);
	
    //book->release();
    dJointGroupEmpty(contactgroup);
	feedback_slider0 = dJointGetFeedback(jointslider[0]);
	feedback_slider1 = dJointGetFeedback(jointslider[1]);
	//printf("%f2.5",feedback_slider->f2[2]);
	  dsSetColor(0,0,1);
    // draw a torso
    pos[4] = dBodyGetPosition(torso.body);
    R[4]  = dBodyGetRotation(torso.body);
    dsDrawCapsule(pos[4],R[4],torso.length,torso.radius);
	
    // draw a capsule
	dsSetColor(1,0,0); //red --> bottom
    pos[0] = dBodyGetPosition(leg1[0].body);
    R[0]   = dBodyGetRotation(leg1[0].body);
    dsDrawCapsuleD(pos[0],R[0],leg1[0].length,leg1[0].radius);
	//
	dsSetColor(1,1,1); //white --> 
    pos[1] = dBodyGetPosition(leg1[1].body);
    R[1]   = dBodyGetRotation(leg1[1].body);
    dsDrawCapsule(pos[1],R[1],leg1[1].length,leg1[1].radius);
	 // draw a capsule
	  // draw a capsule
	dsSetColor(1,0,0);
    pos[2] = dBodyGetPosition(leg2[0].body);
    R[2]   = dBodyGetRotation(leg2[0].body);
    dsDrawCapsule(pos[2],R[2],leg2[0].length,leg2[0].radius);
	//
	dsSetColor(1,1,1);
    pos[3] = dBodyGetPosition(leg2[1].body);
    R[3]   = dBodyGetRotation(leg2[1].body);
    dsDrawCapsule(pos[3],R[3],leg2[1].length,leg2[1].radius);
	//
	//printf("\n %d:%d",n1,n2);
	//
	double fspring1 = sqrt(pow(feedback_slider0->f1[0],2)+pow(feedback_slider0->f1[1],2)+pow(feedback_slider0->f1[2],2));
	double fspring2 = sqrt(pow(feedback_slider1->f1[0],2)+pow(feedback_slider1->f1[1],2)+pow(feedback_slider1->f1[2],2));
	angle1=angle1*180/PI;
	angle2=angle2*180/PI;
	//title
	sheet->writeStr(1, 0,L"time");
	sheet->writeStr(1, 1,L"leg1 touch");sheet->writeStr(1,2,L"leg2 touch");
	sheet->writeStr(1, 3,L"foot1[x]");sheet->writeStr(1,4,L"foot1[y]");sheet->writeStr(1,5,L"foot1[z]");
	sheet->writeStr(1, 6,L"foot2[x]");sheet->writeStr(1,7,L"foot2[y]");sheet->writeStr(1,8,L"foot2[z]");
	sheet->writeStr(1, 9,L"hip_height");
	sheet->writeStr(1, 10,L"torque1");sheet->writeStr(1,11,L"torque2");
	sheet->writeStr(1, 12,L"f[0]");sheet->writeStr(1,13,L"f[1]");
	sheet->writeStr(1, 14,L"fspring1");sheet->writeStr(1,15,L"fspring2");
	sheet->writeStr(1, 16,L"angle1");sheet->writeStr(1,17,L"angle2");
	sheet->writeStr(1, 18,L"phi_tilde");
	count +=1; 
	sheet->writeNum(count, 0 ,timeS);
	sheet->writeNum(count, 1, n1);
	sheet->writeNum(count, 2, n2);
	sheet->writeNum(count, 3,state[0][0]);
    sheet->writeNum(count, 4,state[0][1]);
	sheet->writeNum(count, 5,state[0][2]);
	sheet->writeNum(count, 6,state[1][0]);
	sheet->writeNum(count, 7,state[1][1]);
	sheet->writeNum(count, 8,state[1][2]);
	sheet->writeNum(count, 9,state[2][2]);
	sheet->writeNum(count,10,torque1);
	sheet->writeNum(count,11, torque2);
	sheet->writeNum(count,12,f[0]);
	sheet->writeNum(count,13,f[1]);
	sheet->writeNum(count,14,fspring1);
	sheet->writeNum(count,15,fspring2);
	sheet->writeNum(count,16,angle1);
	sheet->writeNum(count,17,angle2);
	sheet->writeNum(count,18,phi_tilde);
	book->save(L"minh.xls");
	timeS += DT;
	phi_tilde_pre = phi_tilde;
	if(t>5) dsStop();
	

}

void start() {
    static float xyz[3] = {5,0.0,1.5};
    static float hpr[3] = {-180,0.0,0.0};
    dsSetViewpoint (xyz,hpr);
}
void stop () {
	
}
// Create a ball and a pole
void createBTSLIP() {
    dMass mleg,mtorso;
	
    dReal x[5] = {0.0};
	dReal y[5] = {0.0};
	dReal z[5] = {leg_l/2+leg_r,1.5*leg_l+leg_r,leg_l/2+leg_r,1.5*leg_l+leg_r,2*leg_l+leg_r+0.1};
	//dReal z[5] = {leg_l/2,1.5*leg_l,leg_l/2,1.5*leg_l,2*leg_l+0.1};
	//printf("z0 \t",z[0]);
	// leg1
	for (int j=0; j<2; j++)
	{
    leg1[j].radius = leg_r;
    leg1[j].mass   = 0.5;
	leg1[j].length = leg_l;
    leg1[j].body = dBodyCreate(world);
    dMassSetZero(&mleg);
   	dMassSetCapsuleTotal(&mleg,leg1[j].mass,3,leg1[j].radius,leg1[j].length);
    dBodySetMass(leg1[j].body, &mleg);
    dBodySetPosition(leg1[j].body, x[j], y[j], z[j]);
	leg1[j].geom = dCreateCapsule(space,leg1[j].radius,leg1[j].length);
	dGeomSetBody(leg1[j].geom,leg1[j].body);
	//
	leg2[j].radius = leg_r;
    leg2[j].mass   = 0.5;
	leg2[j].length = leg_l;
    leg2[j].body = dBodyCreate(world);
    dMassSetZero(&mleg);
   	dMassSetCapsuleTotal(&mleg,leg2[j].mass,3,leg2[j].radius,leg2[j].length);
    dBodySetMass(leg2[j].body, &mleg);
	const dReal alpha = 69; // Angle of Attack
    dBodySetPosition(leg2[j].body, x[j+2], y[j+2], z[j+2]);
	
	//
	leg2[j].geom = dCreateCapsule(space,leg2[j].radius,leg2[j].length);
	dGeomSetBody(leg2[j].geom,leg2[j].body);
	}
		
	// torso
	torso.radius = 0.125;
    torso.mass   = 80;
	torso.length = 0.8;
    torso.body = dBodyCreate(world);
    dMassSetZero(&mtorso);
   	dMassSetCylinderTotal(&mtorso,torso.mass,3,torso.radius,torso.length);
    dBodySetMass(torso.body, &mtorso);
    dBodySetPosition(torso.body, x[4], y[4], z[4]);
	dBodySetLinearVel(torso.body,0,1.15,0);
	//dBodySetAngularVel (torso.body, 0.11,0,0);
	torso.geom = dCreateCylinder(space,torso.radius,torso.length);
    dGeomSetBody(torso.geom,torso.body);
	
    // hinge joint
    joint[0] = dJointCreateHinge(world, 0);
    dJointAttach(joint[0], leg1[1].body, torso.body);
	dJointSetHingeAnchor(joint[0], 0, 0, 1.1);
    dJointSetHingeAxis(joint[0], 1, 0, 0);
	// hinge joint 1
	joint[1] = dJointCreateHinge(world, 0);
    dJointAttach(joint[1], leg2[1].body, torso.body);
	dJointSetHingeAnchor(joint[1], 0, 0, 1.1);
    dJointSetHingeAxis(joint[1], 1, 0, 0);

	//
	//	slider joint
	jointslider[0] = dJointCreateSlider(world, 0);
    dJointAttach(jointslider[0], leg1[1].body, leg1[0].body);
	dJointSetSliderAxis(jointslider[0],0,0,1);

	//dJointSetSliderParam (jointslider[0],dParamStopCFM,1/(DT*SKP+SKD));
	//dJointSetSliderParam (jointslider[0],dParamStopERP,(DT*SKP)/(DT*SKP+SKD)); 
	dJointSetSliderParam(jointslider[0], dParamLoStop, -0.25);
	dJointSetSliderParam(jointslider[0], dParamHiStop,  0.25);
	//dJointSetSliderParam(jointslider[0], dParamFMax, 2000);
	//
	//
	jointslider[1] = dJointCreateSlider(world, 0);
    dJointAttach(jointslider[1], leg2[1].body, leg2[0].body);
	dJointSetSliderAxis(jointslider[1],0,0,1);

	//--> try to add spring damp to spring joint
	//dJointSetSliderParam (jointslider[1],dParamStopCFM,1/(DT*SKP+SKD));
	//dJointSetSliderParam (jointslider[1],dParamStopERP,(DT*SKP)/(DT*SKP+SKD)); 
	dJointSetSliderParam(jointslider[1], dParamLoStop, -0.25);
	dJointSetSliderParam(jointslider[1], dParamHiStop,  0.25);
	//dJointSetSliderParam(jointslider[1], dParamFMax, 200);
	//first position
	

}

void  setDrawStuff() {
    fn.version = DS_VERSION;
    fn.start   = &start;
    fn.step    = &simLoop;
    fn.command = &command;
    fn.stop    = NULL;
    fn.path_to_textures = "C:/Lib/ode-0.13/drawstuff/textures";		
}

int main (int argc, char **argv) {
    //count = count +1;
	
	dReal state[5][3];
	//
	setDrawStuff();
    dInitODE();
    world = dWorldCreate();
    space = dHashSpaceCreate(0);
    contactgroup = dJointGroupCreate(0);

    dWorldSetGravity(world,0,0,-9.8);
	dWorldSetERP(world,0.9);          // ERP
	dWorldSetCFM(world,1e-3);          // CFM
	// Create a ground
    ground = dCreatePlane(space,0,0,1,0);
	

    // create an object
    createBTSLIP();
	//leg2 with 69 degree --> to contact with the ground.
	get_state(state);
	// Angle of Attack
	dBodySetPosition(leg2[1].body,0,0.5*0.5*cos(alpha*PI/180),state[2][2]-0.5*0.5*sin(alpha*PI/180));
	dBodySetPosition(leg2[0].body,0,0.75*cos(alpha*PI/180),state[2][2]-0.75*sin(alpha*PI/180));
	//1cm compress leg1++++
	dBodySetPosition(leg1[1].body,0,0,1.5*leg_l+leg_r-0.01);
	dBodySetPosition(leg1[0].body,0,0,0.5*leg_l+leg_r-0.00001);
	dJointSetFeedback(jointslider[0],feedback_slider0);
	dJointSetFeedback(jointslider[1],feedback_slider1);

	dsSimulationLoop (argc,argv,1024,768,&fn);
	
    dWorldDestroy (world);
    dCloseODE();
	

    return 0;
}
