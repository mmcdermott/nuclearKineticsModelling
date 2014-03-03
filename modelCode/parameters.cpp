#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <random>

typedef float float_T;
typedef float_T vec_T[2];
const float_T pi = 3.1415926535897;

//Type Parameters
const bool springsOn = true;
const bool translation = true;
const bool ONLY_COMMA = true;

//Parameters
const float_T Duration = 10;         //duration in minutes
const float_T Tau      = 1.0/4000.0; //time step in minutes

// MT Parameters
const unsigned MT_numb = 100; //# of MTs from one centrosome. 
vec_T MT_Pos_M[MT_numb];
vec_T MT_Pos_D[MT_numb];
bool MT_Growing_M[MT_numb];
bool MT_Growing_D[MT_numb];
float_T MT_GrowthVel_M[MT_numb];
float_T MT_GrowthVel_D[MT_numb];
//  Contact Parameters:
float_T contact_length = 400*Tau; //TODO: Find Better Estimate of Me. 
float_T MT_Contact_M[MT_numb];
float_T MT_Contact_D[MT_numb];
//  Growth Parameters
const float_T Vg   = 4*10; //growth velocity in mum/min
//const float_T Vg_c = 0;    //growth velocity after contact 
const float_T Vs_c = 12*10; //shortening velocity after contact
const float_T Vs   = 12*10; //shortening velocity in mum/min 
const float_T kc   = .1*60;   //catastrophe frequency in min-1 
const float_T kr   = .4*60;   //rescue frequency in min-1

float_T Pr_catastrophe = 1-exp(-kc*Tau);
float_T Pr_rescue      = 1-exp(-kr*Tau);
//   Force Parameters
const float_T F_MT = 1; //Force per MT in pN.
vec_T force_M;
vec_T force_D;
vec_T force;
float_T torque_M;
float_T torque_D;
float_T torque;
//   Spring Parameters
const float_T kM = 8;
const float_T kD = 16;

// Pronucleus Parameters
const float_T R1_max = 25;       //Embryo width in mum
const float_T R2_max = 15;       //Embryo length in mum
const float_T Prad   = 4;        //Pronucleous radius (mum)
const float_T Eta    = 1;        //Translational drag coeff (pN/(mum min))
const float_T Mu     = 5 * Eta;  //Rotational drag coeff ((pN mum)/min)
const float_T Eta2   = Eta/10;   //Trans. drag coeff of centrosome (pN/(mum min))
const float_T kbT    = .00414;   //pN mum
const float_T D      = kbT/Eta2; //Diffusion Coefficient for centrosome motion. (mum^2/min)

//  Starting Coordinates:
const float_T startPsi = pi/2.0;
const float_T startX   = 7.5;
const float_T startY   = 0;

//  General Coordinate Initializations: 
float_T psi;
vec_T basePosM;
vec_T basePosD;
vec_T proNucPos;

// Band parameters: 
//  regionAngles:
//   The regionAngles vector defines the endpoints of the angular span of
//   various regions along the cortex. It MUST start with 0 and it MUST end with
//   2*pi. Each region is presumed to be homogenous and all encompassing: In
//   particular, a single band of a protein ont he cortex does not totally
//   define a region unless it is the only protein defining connectivity
//   probability in that region. If two or more bands intersect, you must make
//   your regions various homogenous intersections such that the probability is
//   constant throughout any given region. Note that technically speaking, in
//   implementation a region is defined by two endpoints \theta_1, \theta_2 \in
//   [0,2\pi] such that a point on the cortex at angle \alpha is in that region
//   if \theta_1 <= \alpha < \theta_2
//  regionProbabilities: 
//    This defines the probabilities associated with the regions defined by the
//    regionAngles variable. It has length one less than the regionAngles
//    vector, as it is broken up into regions, not enpoints of regions. 
const int numRegions = 1;
float_T alpha = 3*pi/8.0;
//const float_T regionAngles[numRegions + 1] = {0, alpha, 2*pi - alpha, 2*pi};
//const float_T regionProbabilities[numRegions] = {0.5, 1, 0.5};
const float_T regionAngles[numRegions + 1] = {0, 2*pi};
const float_T regionProbabilities[numRegions] = {1};

// File Parameters
std::ofstream file;
std::string fileName  = "";
std::string fileDir   = "../data/";
std::string fileOrder = "t,proNucPos,psi,MT_Pos_M,MT_Pos_D,force_M,force_D,\
                         force,torque_M,torque_D,torque,basePosM,basePosD";

// Random Number Generation Parameters: 
std::normal_distribution<float_T> stdNormalDist;
std::default_random_engine generator;
