#include <ode/ode.h>
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

// Define parameters
double L_PT_x = 0; 
double L_PT_z = Torso_H/2.0;
// Thigh (Femur)
double L_T = Thigh_H;
// Leg (Tibia)
double L_L = Shin_H+Foot_H+Foot_R;
// Thigh (Femur)
double M_T = Thigh_Mass;
double I_T = If_yy;
double cx_T = 0,  cz_T = L_T/2.0;
// Leg (Tibia)
double M_L = Shin_Mass;
double I_L = Is_yy;
double cx_L = 0,  cz_L = L_L/2.0; 
//Pelvis
double M_P = Torso_Mass;
double I_P = Ib_yy;
//
/*
void auto_MechEnergy(mat q, mat dq, double& KE_tot, double& PE_tot)
{
	double L_PT_x = 0, L_PT_z = 210/1000;
	double L_T = 0.190, L_L = 0.190 + 0.02 + 0.015;
	double M_T = 0.7, I_T = 447.42e-5, cx_T = 0,   cz_T = 0.12577;
	double M_L = 0.7, I_L = 447.42e-5, cx_L = 0,   cz_L = 125.77/1000;
	double M_P = 12.5, I_P = M_P/12*(0.42*0.42+0.2*2);

  KE_tot = (M_L*(pow((dq(4)+dq(0)*(L_T*cos(q(6)+q(0))+cz_L*cos(q(6)+q(1)+q(0))-cx_L*sin(q(6)+q(1)+q(0)))+dq(6)*(L_T*cos(q(6)+q(0))+L_PT_z*cos(q(6))-L_PT_x*sin(q(6))+cz_L*cos(q(6)+q(1)+q(0))-cx_L*sin(q(6)+q(1)+q(0)))+dq(1)*(cz_L*cos(q(6)+q(1)+q(0))-cx_L*sin(q(6)+q(1)+q(0)))), 2)+pow((dq(5)+dq(0)*(L_T*sin(q(6)+q(0))+cx_L*cos(q(6)+q(1)+q(0))+cz_L*sin(q(6)+q(1)+q(0)))+dq(6)*(L_T*sin(q(6)+q(0))+L_PT_x*cos(q(6))+L_PT_z*sin(q(6))+cx_L*cos(q(6)+q(1)+q(0))+cz_L*sin(q(6)+q(1)+q(0)))+dq(1)*(cx_L*cos(q(6)+q(1)+q(0))+cz_L*sin(q(6)+q(1)+q(0)))), 2)))/2.0+(M_L*(pow((dq(4)+dq(2)*(L_T*cos(q(6)+q(2))+cz_L*cos(q(6)+q(3)+q(2))-cx_L*sin(q(6)+q(3)+q(2)))+dq(6)*(L_T*cos(q(6)+q(2))+L_PT_z*cos(q(6))-L_PT_x*sin(q(6))+cz_L*cos(q(6)+q(3)+q(2))-cx_L*sin(q(6)+q(3)+q(2)))+dq(3)*(cz_L*cos(q(6)+q(3)+q(2))-cx_L*sin(q(6)+q(3)+q(2)))), 2)+pow((dq(5)+dq(2)*(L_T*sin(q(6)+q(2))+cx_L*cos(q(6)+q(3)+q(2))+cz_L*sin(q(6)+q(3)+q(2)))+dq(6)*(L_T*sin(q(6)+q(2))+L_PT_x*cos(q(6))+L_PT_z*sin(q(6))+cx_L*cos(q(6)+q(3)+q(2))+cz_L*sin(q(6)+q(3)+q(2)))+dq(3)*(cx_L*cos(q(6)+q(3)+q(2))+cz_L*sin(q(6)+q(3)+q(2)))), 2)))/2.0+(I_T*pow((dq(6)+dq(0)), 2))/2.0+(I_T*pow((dq(6)+dq(2)), 2))/2.0+(I_P*pow(dq(6), 2))/2.0+(M_T*(pow((dq(4)+dq(0)*(cz_T*cos(q(6)+q(0))-cx_T*sin(q(6)+q(0)))+dq(6)*(cz_T*cos(q(6)+q(0))-cx_T*sin(q(6)+q(0))+L_PT_z*cos(q(6))-L_PT_x*sin(q(6)))), 2)+pow((dq(5)+dq(0)*(cx_T*cos(q(6)+q(0))+cz_T*sin(q(6)+q(0)))+dq(6)*(cx_T*cos(q(6)+q(0))+cz_T*sin(q(6)+q(0))+L_PT_x*cos(q(6))+L_PT_z*sin(q(6)))), 2)))/2.0+(M_T*(pow((dq(4)+dq(2)*(cz_T*cos(q(6)+q(2))-cx_T*sin(q(6)+q(2)))+dq(6)*(cz_T*cos(q(6)+q(2))-cx_T*sin(q(6)+q(2))+L_PT_z*cos(q(6))-L_PT_x*sin(q(6)))), 2)+pow((dq(5)+dq(2)*(cx_T*cos(q(6)+q(2))+cz_T*sin(q(6)+q(2)))+dq(6)*(cx_T*cos(q(6)+q(2))+cz_T*sin(q(6)+q(2))+L_PT_x*cos(q(6))+L_PT_z*sin(q(6)))), 2)))/2.0+(M_P*(pow(dq(4), 2)+pow(dq(5), 2)))/2.0+(I_L*pow((dq(6)+dq(1)+dq(0)), 2))/2.0+(I_L*pow((dq(6)+dq(3)+dq(2)), 2))/2.0 ;
  PE_tot = (981*M_L*q(5))/50.0+(981*M_P*q(5))/100.0+(981*M_T*q(5))/50.0-(981*M_L*cz_L*cos(q(6)+q(1)+q(0)))/100.0-(981*M_L*cz_L*cos(q(6)+q(3)+q(2)))/100.0+(981*M_L*cx_L*sin(q(6)+q(1)+q(0)))/100.0+(981*M_L*cx_L*sin(q(6)+q(3)+q(2)))/100.0-(981*L_T*M_L*cos(q(6)+q(0)))/100.0-(981*L_T*M_L*cos(q(6)+q(2)))/100.0-(981*M_T*cz_T*cos(q(6)+q(0)))/100.0-(981*M_T*cz_T*cos(q(6)+q(2)))/100.0+(981*M_T*cx_T*sin(q(6)+q(0)))/100.0+(981*M_T*cx_T*sin(q(6)+q(2)))/100.0-(981*L_PT_z*M_L*cos(q(6)))/50.0-(981*L_PT_z*M_T*cos(q(6)))/50.0+(981*L_PT_x*M_L*sin(q(6)))/50.0+(981*L_PT_x*M_T*sin(q(6)))/50.0 ;
}
*/
mat get_D(mat q_f)
{	
	double q1, q2, q3, q4, q5, q6, q7;
	mat D_f;
	q1 = q_f(0,0); q2 = q_f(1,0); q3 = q_f(2,0); q4 = q_f(3,0); q5 = q_f(4,0); q6 = q_f(5,0); q7 = q_f(6,0);
	D_f = zeros(7,7);
  D_f(0, 0) = I_L+I_T+pow(L_T, 2)*M_L+M_L*pow(cx_L, 2)+M_T*pow(cx_T, 2)+M_L*pow(cz_L, 2)+M_T*pow(cz_T, 2)+2*L_T*M_L*cz_L*cos(q2)-2*L_T*M_L*cx_L*sin(q2) ;
  D_f(0, 1) = I_L+M_L*pow(cx_L, 2)+M_L*pow(cz_L, 2)+L_T*M_L*cz_L*cos(q2)-L_T*M_L*cx_L*sin(q2) ;
  D_f(0, 2) = 0 ;
  D_f(0, 3) = 0 ;
  D_f(0, 4) = (M_L*(2*L_T*cos(q7+q1)+2*cz_L*cos(q7+q2+q1)-2*cx_L*sin(q7+q2+q1)))/2.0+(M_T*(2*cz_T*cos(q7+q1)-2*cx_T*sin(q7+q1)))/2.0 ;
  D_f(0, 5) = (M_L*(2*L_T*sin(q7+q1)+2*cx_L*cos(q7+q2+q1)+2*cz_L*sin(q7+q2+q1)))/2.0+(M_T*(2*cx_T*cos(q7+q1)+2*cz_T*sin(q7+q1)))/2.0 ;
  D_f(0, 6) = I_L+I_T+pow(L_T, 2)*M_L+M_L*pow(cx_L, 2)+M_T*pow(cx_T, 2)+M_L*pow(cz_L, 2)+M_T*pow(cz_T, 2)+L_PT_x*M_L*cx_L*cos(q2+q1)+L_PT_z*M_L*cz_L*cos(q2+q1)-L_PT_z*M_L*cx_L*sin(q2+q1)+L_PT_x*M_L*cz_L*sin(q2+q1)+L_PT_z*L_T*M_L*cos(q1)+L_PT_x*L_T*M_L*sin(q1)+L_PT_x*M_T*cx_T*cos(q1)+2*L_T*M_L*cz_L*cos(q2)+L_PT_z*M_T*cz_T*cos(q1)-2*L_T*M_L*cx_L*sin(q2)-L_PT_z*M_T*cx_T*sin(q1)+L_PT_x*M_T*cz_T*sin(q1) ;
  D_f(1, 0) = I_L+M_L*pow(cx_L, 2)+M_L*pow(cz_L, 2)+L_T*M_L*cz_L*cos(q2)-L_T*M_L*cx_L*sin(q2) ;
  D_f(1, 1) = I_L+M_L*pow(cx_L, 2)+M_L*pow(cz_L, 2) ;
  D_f(1, 2) = 0 ;
  D_f(1, 3) = 0 ;
  D_f(1, 4) = M_L*cz_L*cos(q7+q2+q1)-M_L*cx_L*sin(q7+q2+q1) ;
  D_f(1, 5) = M_L*cx_L*cos(q7+q2+q1)+M_L*cz_L*sin(q7+q2+q1) ;
  D_f(1, 6) = I_L+M_L*pow(cx_L, 2)+M_L*pow(cz_L, 2)+L_PT_x*M_L*cx_L*cos(q2+q1)+L_PT_z*M_L*cz_L*cos(q2+q1)-L_PT_z*M_L*cx_L*sin(q2+q1)+L_PT_x*M_L*cz_L*sin(q2+q1)+L_T*M_L*cz_L*cos(q2)-L_T*M_L*cx_L*sin(q2) ;
  D_f(2, 0) = 0 ;
  D_f(2, 1) = 0 ;
  D_f(2, 2) = I_L+I_T+pow(L_T, 2)*M_L+M_L*pow(cx_L, 2)+M_T*pow(cx_T, 2)+M_L*pow(cz_L, 2)+M_T*pow(cz_T, 2)+2*L_T*M_L*cz_L*cos(q4)-2*L_T*M_L*cx_L*sin(q4) ;
  D_f(2, 3) = I_L+M_L*pow(cx_L, 2)+M_L*pow(cz_L, 2)+L_T*M_L*cz_L*cos(q4)-L_T*M_L*cx_L*sin(q4) ;
  D_f(2, 4) = (M_L*(2*L_T*cos(q7+q3)+2*cz_L*cos(q7+q4+q3)-2*cx_L*sin(q7+q4+q3)))/2.0+(M_T*(2*cz_T*cos(q7+q3)-2*cx_T*sin(q7+q3)))/2.0 ;
  D_f(2, 5) = (M_L*(2*L_T*sin(q7+q3)+2*cx_L*cos(q7+q4+q3)+2*cz_L*sin(q7+q4+q3)))/2.0+(M_T*(2*cx_T*cos(q7+q3)+2*cz_T*sin(q7+q3)))/2.0 ;
  D_f(2, 6) = I_L+I_T+pow(L_T, 2)*M_L+M_L*pow(cx_L, 2)+M_T*pow(cx_T, 2)+M_L*pow(cz_L, 2)+M_T*pow(cz_T, 2)+L_PT_x*M_L*cx_L*cos(q4+q3)+L_PT_z*M_L*cz_L*cos(q4+q3)-L_PT_z*M_L*cx_L*sin(q4+q3)+L_PT_x*M_L*cz_L*sin(q4+q3)+L_PT_z*L_T*M_L*cos(q3)+L_PT_x*L_T*M_L*sin(q3)+L_PT_x*M_T*cx_T*cos(q3)+2*L_T*M_L*cz_L*cos(q4)+L_PT_z*M_T*cz_T*cos(q3)-2*L_T*M_L*cx_L*sin(q4)-L_PT_z*M_T*cx_T*sin(q3)+L_PT_x*M_T*cz_T*sin(q3) ;
  D_f(3, 0) = 0 ;
  D_f(3, 1) = 0 ;
  D_f(3, 2) = I_L+M_L*pow(cx_L, 2)+M_L*pow(cz_L, 2)+L_T*M_L*cz_L*cos(q4)-L_T*M_L*cx_L*sin(q4) ;
  D_f(3, 3) = I_L+M_L*pow(cx_L, 2)+M_L*pow(cz_L, 2) ;
  D_f(3, 4) = M_L*cz_L*cos(q7+q4+q3)-M_L*cx_L*sin(q7+q4+q3) ;
  D_f(3, 5) = M_L*cx_L*cos(q7+q4+q3)+M_L*cz_L*sin(q7+q4+q3) ;
  D_f(3, 6) = I_L+M_L*pow(cx_L, 2)+M_L*pow(cz_L, 2)+L_PT_x*M_L*cx_L*cos(q4+q3)+L_PT_z*M_L*cz_L*cos(q4+q3)-L_PT_z*M_L*cx_L*sin(q4+q3)+L_PT_x*M_L*cz_L*sin(q4+q3)+L_T*M_L*cz_L*cos(q4)-L_T*M_L*cx_L*sin(q4) ;
  D_f(4, 0) = (M_L*(2*L_T*cos(q7+q1)+2*cz_L*cos(q7+q2+q1)-2*cx_L*sin(q7+q2+q1)))/2.0+(M_T*(2*cz_T*cos(q7+q1)-2*cx_T*sin(q7+q1)))/2.0 ;
  D_f(4, 1) = M_L*cz_L*cos(q7+q2+q1)-M_L*cx_L*sin(q7+q2+q1) ;
  D_f(4, 2) = (M_L*(2*L_T*cos(q7+q3)+2*cz_L*cos(q7+q4+q3)-2*cx_L*sin(q7+q4+q3)))/2.0+(M_T*(2*cz_T*cos(q7+q3)-2*cx_T*sin(q7+q3)))/2.0 ;
  D_f(4, 3) = M_L*cz_L*cos(q7+q4+q3)-M_L*cx_L*sin(q7+q4+q3) ;
  D_f(4, 4) = 2*M_L+M_P+2*M_T ;
  D_f(4, 5) = 0 ;
  D_f(4, 6) = M_L*cz_L*cos(q7+q2+q1)+M_L*cz_L*cos(q7+q4+q3)-M_L*cx_L*sin(q7+q2+q1)-M_L*cx_L*sin(q7+q4+q3)+L_T*M_L*cos(q7+q1)+L_T*M_L*cos(q7+q3)+M_T*cz_T*cos(q7+q1)+M_T*cz_T*cos(q7+q3)-M_T*cx_T*sin(q7+q1)-M_T*cx_T*sin(q7+q3)+2*L_PT_z*M_L*cos(q7)+2*L_PT_z*M_T*cos(q7)-2*L_PT_x*M_L*sin(q7)-2*L_PT_x*M_T*sin(q7) ;
  D_f(5, 0) = (M_L*(2*L_T*sin(q7+q1)+2*cx_L*cos(q7+q2+q1)+2*cz_L*sin(q7+q2+q1)))/2.0+(M_T*(2*cx_T*cos(q7+q1)+2*cz_T*sin(q7+q1)))/2.0 ;
  D_f(5, 1) = M_L*cx_L*cos(q7+q2+q1)+M_L*cz_L*sin(q7+q2+q1) ;
  D_f(5, 2) = (M_L*(2*L_T*sin(q7+q3)+2*cx_L*cos(q7+q4+q3)+2*cz_L*sin(q7+q4+q3)))/2.0+(M_T*(2*cx_T*cos(q7+q3)+2*cz_T*sin(q7+q3)))/2.0 ;
  D_f(5, 3) = M_L*cx_L*cos(q7+q4+q3)+M_L*cz_L*sin(q7+q4+q3) ;
  D_f(5, 4) = 0 ;
  D_f(5, 5) = 2*M_L+M_P+2*M_T ;
  D_f(5, 6) = M_L*cx_L*cos(q7+q2+q1)+M_L*cx_L*cos(q7+q4+q3)+M_L*cz_L*sin(q7+q2+q1)+M_L*cz_L*sin(q7+q4+q3)+L_T*M_L*sin(q7+q1)+L_T*M_L*sin(q7+q3)+M_T*cx_T*cos(q7+q1)+M_T*cx_T*cos(q7+q3)+M_T*cz_T*sin(q7+q1)+M_T*cz_T*sin(q7+q3)+2*L_PT_x*M_L*cos(q7)+2*L_PT_x*M_T*cos(q7)+2*L_PT_z*M_L*sin(q7)+2*L_PT_z*M_T*sin(q7) ;
  D_f(6, 0) = I_L+I_T+pow(L_T, 2)*M_L+M_L*pow(cx_L, 2)+M_T*pow(cx_T, 2)+M_L*pow(cz_L, 2)+M_T*pow(cz_T, 2)+L_PT_x*M_L*cx_L*cos(q2+q1)+L_PT_z*M_L*cz_L*cos(q2+q1)-L_PT_z*M_L*cx_L*sin(q2+q1)+L_PT_x*M_L*cz_L*sin(q2+q1)+L_PT_z*L_T*M_L*cos(q1)+L_PT_x*L_T*M_L*sin(q1)+L_PT_x*M_T*cx_T*cos(q1)+2*L_T*M_L*cz_L*cos(q2)+L_PT_z*M_T*cz_T*cos(q1)-2*L_T*M_L*cx_L*sin(q2)-L_PT_z*M_T*cx_T*sin(q1)+L_PT_x*M_T*cz_T*sin(q1) ;
  D_f(6, 1) = I_L+M_L*pow(cx_L, 2)+M_L*pow(cz_L, 2)+L_PT_x*M_L*cx_L*cos(q2+q1)+L_PT_z*M_L*cz_L*cos(q2+q1)-L_PT_z*M_L*cx_L*sin(q2+q1)+L_PT_x*M_L*cz_L*sin(q2+q1)+L_T*M_L*cz_L*cos(q2)-L_T*M_L*cx_L*sin(q2) ;
  D_f(6, 2) = I_L+I_T+pow(L_T, 2)*M_L+M_L*pow(cx_L, 2)+M_T*pow(cx_T, 2)+M_L*pow(cz_L, 2)+M_T*pow(cz_T, 2)+L_PT_x*M_L*cx_L*cos(q4+q3)+L_PT_z*M_L*cz_L*cos(q4+q3)-L_PT_z*M_L*cx_L*sin(q4+q3)+L_PT_x*M_L*cz_L*sin(q4+q3)+L_PT_z*L_T*M_L*cos(q3)+L_PT_x*L_T*M_L*sin(q3)+L_PT_x*M_T*cx_T*cos(q3)+2*L_T*M_L*cz_L*cos(q4)+L_PT_z*M_T*cz_T*cos(q3)-2*L_T*M_L*cx_L*sin(q4)-L_PT_z*M_T*cx_T*sin(q3)+L_PT_x*M_T*cz_T*sin(q3) ;
  D_f(6, 3) = I_L+M_L*pow(cx_L, 2)+M_L*pow(cz_L, 2)+L_PT_x*M_L*cx_L*cos(q4+q3)+L_PT_z*M_L*cz_L*cos(q4+q3)-L_PT_z*M_L*cx_L*sin(q4+q3)+L_PT_x*M_L*cz_L*sin(q4+q3)+L_T*M_L*cz_L*cos(q4)-L_T*M_L*cx_L*sin(q4) ;
  D_f(6, 4) = M_L*cz_L*cos(q7+q2+q1)+M_L*cz_L*cos(q7+q4+q3)-M_L*cx_L*sin(q7+q2+q1)-M_L*cx_L*sin(q7+q4+q3)+L_T*M_L*cos(q7+q1)+L_T*M_L*cos(q7+q3)+M_T*cz_T*cos(q7+q1)+M_T*cz_T*cos(q7+q3)-M_T*cx_T*sin(q7+q1)-M_T*cx_T*sin(q7+q3)+2*L_PT_z*M_L*cos(q7)+2*L_PT_z*M_T*cos(q7)-2*L_PT_x*M_L*sin(q7)-2*L_PT_x*M_T*sin(q7) ;
  D_f(6, 5) = M_L*cx_L*cos(q7+q2+q1)+M_L*cx_L*cos(q7+q4+q3)+M_L*cz_L*sin(q7+q2+q1)+M_L*cz_L*sin(q7+q4+q3)+L_T*M_L*sin(q7+q1)+L_T*M_L*sin(q7+q3)+M_T*cx_T*cos(q7+q1)+M_T*cx_T*cos(q7+q3)+M_T*cz_T*sin(q7+q1)+M_T*cz_T*sin(q7+q3)+2*L_PT_x*M_L*cos(q7)+2*L_PT_x*M_T*cos(q7)+2*L_PT_z*M_L*sin(q7)+2*L_PT_z*M_T*sin(q7) ;
  D_f(6, 6) = 2*I_L+I_P+2*I_T+2*pow(L_PT_x, 2)*M_L+2*pow(L_PT_z, 2)*M_L+2*pow(L_T, 2)*M_L+2*pow(L_PT_x, 2)*M_T+2*pow(L_PT_z, 2)*M_T+2*M_L*pow(cx_L, 2)+2*M_T*pow(cx_T, 2)+2*M_L*pow(cz_L, 2)+2*M_T*pow(cz_T, 2)+2*L_PT_x*M_L*cx_L*cos(q2+q1)+2*L_PT_x*M_L*cx_L*cos(q4+q3)+2*L_PT_z*M_L*cz_L*cos(q2+q1)+2*L_PT_z*M_L*cz_L*cos(q4+q3)-2*L_PT_z*M_L*cx_L*sin(q2+q1)-2*L_PT_z*M_L*cx_L*sin(q4+q3)+2*L_PT_x*M_L*cz_L*sin(q2+q1)+2*L_PT_x*M_L*cz_L*sin(q4+q3)+2*L_PT_z*L_T*M_L*cos(q1)+2*L_PT_z*L_T*M_L*cos(q3)+2*L_PT_x*L_T*M_L*sin(q1)+2*L_PT_x*L_T*M_L*sin(q3)+2*L_PT_x*M_T*cx_T*cos(q1)+2*L_PT_x*M_T*cx_T*cos(q3)+2*L_T*M_L*cz_L*cos(q2)+2*L_T*M_L*cz_L*cos(q4)+2*L_PT_z*M_T*cz_T*cos(q1)+2*L_PT_z*M_T*cz_T*cos(q3)-2*L_T*M_L*cx_L*sin(q2)-2*L_T*M_L*cx_L*sin(q4)-2*L_PT_z*M_T*cx_T*sin(q1)-2*L_PT_z*M_T*cx_T*sin(q3)+2*L_PT_x*M_T*cz_T*sin(q1)+2*L_PT_x*M_T*cz_T*sin(q3) ;
  return D_f;
}
mat get_jacobian1(mat q_f)
{	
	double q1, q2, q3, q4, q5, q6, q7;
	mat jac_c1;
	q1 = q_f(0,0); q2 = q_f(1,0); q3 = q_f(2,0); q4 = q_f(3,0); q5 = q_f(4,0); q6 = q_f(5,0); q7 = q_f(6,0);
	jac_c1 = zeros(2,7);
  jac_c1(0, 0) = L_T*cos(q7+q1)+L_L*cos(q7+q2+q1) ;
  jac_c1(0, 1) = L_L*cos(q7+q2+q1) ;
  jac_c1(0, 2) = 0 ;
  jac_c1(0, 3) = 0 ;
  jac_c1(0, 4) = 1 ;
  jac_c1(0, 5) = 0 ;
  jac_c1(0, 6) = L_T*cos(q7+q1)+L_PT_z*cos(q7)-L_PT_x*sin(q7)+L_L*cos(q7+q2+q1) ;
  jac_c1(1, 0) = L_T*sin(q7+q1)+L_L*sin(q7+q2+q1) ;
  jac_c1(1, 1) = L_L*sin(q7+q2+q1) ;
  jac_c1(1, 2) = 0 ;
  jac_c1(1, 3) = 0 ;
  jac_c1(1, 4) = 0 ;
  jac_c1(1, 5) = 1 ;
  jac_c1(1, 6) = L_T*sin(q7+q1)+L_PT_x*cos(q7)+L_PT_z*sin(q7)+L_L*sin(q7+q2+q1) ;

	return jac_c1;
}
mat get_jacobian2(mat q_f)
{	
	double q1, q2, q3, q4, q5, q6, q7;
	mat jac_c2;
	q1 = q_f(0,0); q2 = q_f(1,0); q3 = q_f(2,0); q4 = q_f(3,0); q5 = q_f(4,0); q6 = q_f(5,0); q7 = q_f(6,0);
	jac_c2 = zeros(2,7);
  jac_c2(0, 0) = 0 ;
  jac_c2(0, 1) = 0 ;
  jac_c2(0, 2) = L_T*cos(q7+q3)+L_L*cos(q7+q4+q3) ;
  jac_c2(0, 3) = L_L*cos(q7+q4+q3) ;
  jac_c2(0, 4) = 1 ;
  jac_c2(0, 5) = 0 ;
  jac_c2(0, 6) = L_T*cos(q7+q3)+L_PT_z*cos(q7)-L_PT_x*sin(q7)+L_L*cos(q7+q4+q3) ;
  jac_c2(1, 0) = 0 ;
  jac_c2(1, 1) = 0 ;
  jac_c2(1, 2) = L_T*sin(q7+q3)+L_L*sin(q7+q4+q3) ;
  jac_c2(1, 3) = L_L*sin(q7+q4+q3) ;
  jac_c2(1, 4) = 0 ;
  jac_c2(1, 5) = 1 ;
  jac_c2(1, 6) = L_T*sin(q7+q3)+L_PT_x*cos(q7)+L_PT_z*sin(q7)+L_L*sin(q7+q4+q3) ;
  return jac_c2;
}
