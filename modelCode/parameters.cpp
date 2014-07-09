#include <stdio.h>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <random>
#include "Vector.hpp"

//TODO: This typedef is unnecessary. The code in main.cpp should be updated to
//use Vector straight. 
typedef Vector vec_T;
const float_T pi = 3.1415926535897;

//This will define our centrosome options. 
enum MTOC {M_CENTROSOME, D_CENTROSOME};

//Debug parameters. 
bool spitValues = false;

//Type Parameters
const bool motherSpringOn = true;
const bool daughterSpringOn = false;
const bool translation = true;
const bool ONLY_COMMA = true;

//Parameters
const float_T Duration = 25;         //duration in minutes
const float_T Tau      = 1.0/4000.0; //time step in minutes

// MT Parameters
const unsigned MT_numb = 100; //# of MTs from one centrosome.  (maybe 1000)
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
const float_T F_MT   = 1; //Force per MT in pN.
const float_T Fratio = 1.0;// Env. killing: 0.32299363312512985;
vec_T force_M;
vec_T force_D;
vec_T force;
float_T torque_M;
float_T torque_D;
float_T torque;
//   Envelope Parameters
const float_T envWidthM = pi/2;
const float_T envelopeM[2] = {pi/2.0 - envWidthM/2.0, pi/2.0 + envWidthM/2.0}; // The envelope in which MTs from M can grow. 
const float_T envWidthD = pi/2;
const float_T envelopeD[2] = {pi/2.0 - envWidthD/2.0, pi/2.0 + envWidthD/2.0}; // The envelope in which MTs from M can grow. 
//   Spring Parameters
const float_T kM = 8;
const float_T kD = 16;

// Pronucleus Parameters
const float_T R1_max = 25;       //Embryo width in mum
const float_T R2_max = 15;       //Embryo length in mum
const float_T Prad   = 4;        //Pronucleous radius (mum)
const float_T Eta    = 2;        //Translational drag coeff (pN/(mum min))
const float_T Mu     = 5 * Eta;  //Rotational drag coeff ((pN mum)/min)
const float_T Eta2   = Eta/10;   //Trans. drag coeff of centrosome (pN/(mum min))
const float_T kbT    = .00414;   //pN mum
const float_T D      = kbT/Eta2; //Diffusion Coefficient for centrosome motion. (mum^2/min)

//  Starting Coordinates:
//   Centered Coordinates: 
//const float_T startPsi = pi/2.0;
//const float_T startX   = 0;
//const float_T startY   = 0;

//   Off-center Coordinates
const float_T startPsi = pi/2.0 + pi/8.0;
const float_T startX   = R1_max/5.0;
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
//  regionForceMultipliers: 
//    This defines the force multipliers associated with the regions defined by
//    the regionAngles variable. It has length one less than the regionAngles
//    vector, as it is broken up into regions, not enpoints of regions. 

//No Bands:
//const int numRegions = 1;
//const float_T regionAngles[numRegions + 1] = {0, 2*pi};
//const float_T regionProbabilities[numRegions] = {1};
//const float_T regionForceMultipliers[numRegions] = {1};

//Standard Bands:
const int numRegions = 5;
const float_T width = pi/4;
const float_T centerPos = 1.24287;
const float_T start = centerPos - width/2.0;
const float_T end = centerPos + width/2.0;
const float_T regionAngles[numRegions+1] = {0, start, end, 2*pi - end, 2*pi - start, 2*pi};
const float_T regionProbabilities[numRegions] = {1,1,1,1,1};
const float_T regionForceMultipliers[numRegions] = {1,-1,1,-1,1};

// MT density limitations
//  Only one contact per window, windows of length ~1
const size_t numberContactWindows = 128;
const float_T contactWindowAngles[numberContactWindows+1] = 
  {0., 0.039985, 0.079885, 0.119585, 0.159085, 0.198485, 0.237785, 0.277085,
    0.316485, 0.356085, 0.396085, 0.436485, 0.477585, 0.519485, 0.562285,
    0.606085, 0.651085, 0.697385, 0.745085, 0.794285, 0.845085, 0.897585,
    0.951785, 1.00768, 1.06538, 1.12468, 1.18558, 1.24798, 1.31168, 1.37648,
    1.44218, 1.50848, 1.57508, 1.64168, 1.70788, 1.77348, 1.83818, 1.90178,
    1.96398, 2.02468, 2.08378, 2.14119, 2.19689, 2.25089, 2.30309, 2.35369,
    2.40269, 2.45019, 2.49629, 2.54109, 2.58479, 2.62739, 2.66909, 2.71009,
    2.75049, 2.79039, 2.82999, 2.86939, 2.90869, 2.94799, 2.98739, 3.02699,
    3.06669, 3.10659, 3.14649, 3.18639, 3.22619, 3.26589, 3.30539, 3.34479,
    3.38409, 3.42339, 3.46279, 3.50239, 3.54239, 3.58289, 3.62409, 3.66609,
    3.70899, 3.75289, 3.79799, 3.84439, 3.89219, 3.94159, 3.99259, 4.04529,
    4.09969, 4.15579, 4.21359, 4.27309, 4.33419, 4.39679, 4.46069, 4.52569,
    4.59149, 4.65779, 4.72439, 4.79089, 4.85709, 4.92259, 4.98719, 5.05059,
    5.11269, 5.17329, 5.23219, 5.28939, 5.34489, 5.39869, 5.45069, 5.50109,
    5.54989, 5.59719, 5.64309, 5.68779, 5.73128, 5.77378, 5.81548, 5.85638,
    5.89668, 5.93658, 5.97608, 6.01538, 6.05468, 6.09398, 6.13338, 6.17298,
    6.21268, 6.25258, 6.28319};

bool contacts[numberContactWindows];

// File Parameters
std::ofstream file;
std::string fileName  = "";
std::string fileDir   = "../data/";
std::string fileOrder = "t,proNucPos,psi,MT_Pos_M,MT_Pos_D,force_M,force_D,\
                         force,torque_M,torque_D,torque,basePosM,basePosD";

// Random Number Generation Parameters: 
std::normal_distribution<float_T> stdNormalDist;
std::default_random_engine generator;
