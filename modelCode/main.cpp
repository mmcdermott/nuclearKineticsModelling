#include "initialization.cpp"
#include <iostream>
#include <chrono>
#include <random>
#include <fstream>
#include <assert.h>
#include <cstring>
#include <limits>

typedef std::numeric_limits< float_T >flt;

using namespace std;

float_T testStat() {
  /* testStat(): A function to generate a test statistic.
   * Input: None
   * Output: A uniform random float_T between 0 and 1
   */
  return ((float_T) rand())/((float_T) RAND_MAX);
}

ostream& operator<<(ostream& out, const vec_T& rhs) {
  /* operator<<(ostream, vec_T): A wrapper to write a vector to an output
   *   stream.
   * Inputs: 
   *   ostream& out: The stream to which rhs should be written. 
   *   const vec_T& rhs: The Vector that should be written.
   * Output: The ostream out. 
   */
  if (ONLY_COMMA) {
    //There are two modes of writing. This is by far the most common. 
    out << rhs[0] << "," << rhs[1];
  } else {
    out << rhs[0] << "|" << rhs[1];
  }
  return out;
}

ostream& operator<<(ostream& out, const vec_T (&rhs)[MT_numb]) {
  /* operator<<(ostream, vec_T[MT_numb]): A wrapper to write an array of vectors
   *   to an output stream.
   * Inputs: 
   *   ostream& out: The stream to which rhs should be written. 
   *   const vec_T (&rhs)[MT_numb]: The array of vectors that should be written
   * Output: The ostream out. 
   */
  out << rhs[0];
  for (size_t i = 1; i < MT_numb; ++i) {
    if (ONLY_COMMA) {
      //There are two modes of writing. This is by far the most common. 
      out << "," << rhs[i];
    } else {
      out << ";" << rhs[i];
    }
  }
  return out;
}

ostream& operator<<(ostream& out, const float_T (&rhs)[MT_numb]) {
  /* operator<<(ostream, float_T[MT_numb]): A wrapper to write an array of floats
   *   to an output stream.
   * Inputs: 
   *   ostream& out: The stream to which rhs should be written. 
   *   const vec_T (&rhs)[MT_numb]: The array of floats that should be written
   * Output: The ostream out. 
   */
  out << rhs[0];
  for (size_t i = 1; i < MT_numb; ++i) {
    if (ONLY_COMMA) {
      //There are two modes of writing. This is by far the most common. 
      out << "," << rhs[i];
    } else {
      out << ";" << rhs[i];
    }
  }
  return out;
}

void writeData(const float_T t) {
  /* writeData: A function that writes all relevant data to the data file
   * Inputs:
   *   const float_T t: The simulation time at which the points are being
   *     written.
   * Output: (none)
   */
  file << t <<","<< proNucPos <<","<< psi <<","<<  MT_Pos_M <<",";
  file << MT_Pos_D <<","<< force_M <<","<< force_D <<","<< force <<",";
  file << torque_M <<","<< torque_D <<","<< torque <<","<< basePosM <<",";
  file << basePosD << endl;
}

void writeFinalData() {
  /* writeFinalData: A function that writes the end position data to the data
   *   file.
   * Inputs: (none)
   * Output: (none)
   */
  file << proNucPos << "," << psi << endl;
}

void writePartialData(const float_T t) {
  /* writePartialData: A function that writes some relevant positional data to 
   *   the data file
   * Inputs:
   *   const float_T t: The simulation time at which the points are being
   *     written.
   * Output: (none)
   */
  file<<t<<","<<proNucPos <<","<< psi <<","<< basePosM <<","<< basePosD << endl;
}

void mtForceCalc(const vec_T (&mtEndPos)[MT_numb], const vec_T &basePos, const
    float_T (&mtContact)[MT_numb], vec_T &force, const float_T forceMag = F_MT) 
{
  /* mtForceCalc: A function which computes the force on an MTOC implied by the
   *   cortical pushing and pulling forces due to the MTs. 
   * Inputs: 
   *   const vec_T (&mtEndPos)[MT_numb]: A reference to an array of end mt
   *     positions. Basically, this tells the function where to look to find the
   *     MT end points, from which it can compute force if the MTs are making
   *     contact. Note that as this is a reference, the entire array is never
   *     passed to the function, just an address, making this function not quite
   *     as inefficient as it might seem. 
   *   const vec_T &basePos: A reference to the position of the MTOC. Without
   *     springs, this is easily calculatable, but with springs must be tracked,
   *     so we just keep it in all the time. 
   *   const float_T (&mtContact)[MT_numb]: A reference to an array flagging MT
   *     contacts, encoded via time remaining. If the time remaining is zero,
   *     the MT is not making contact currently. As we are passing a reference,
   *     it is again not so inefficient as it might seem. 
   *   vec_T &force: A reference to the particular force vector that stores the
   *     force on this MTOC (either force_M or force_D). This gets overwritten
   *     to store the result, and, as it is a reference, this change persists
   *     outside of this particular function. 
   *   const float_T forceMag: This is the magnitude of the force each MT should
   *     impart when being pulled by cortical dynein. Pushing is handled by
   *     using a multiplier, so this typically also is the magnitude of any
   *     cortical pushing forces as well. 
   * Output: (none)
   * Issues: (TODO) This function could probably be simplified so that not so
   *   many things need to be passed by reference here... maybe store the mt
   *   parameters in an array or something? It seems excessive, is all. 
   */
  //First we'll zero out the forces, because we're about to compute them from
  //scratch. 
  force.zero();
  for (size_t i=0; i < MT_numb; ++i) {
    //We only need to compute the force if the MT's making contact.
    if (mtContact[i] == 0) continue;

    //We need to search to see what multiplier we should assign, which requires
    //knowing the true cartesian angle of the vector. 
    //3D WARNING: This section would have to update when changing to 3d.
    float_T angle = atan2(mtEndPos[i][1],mtEndPos[i][0]);
    if (angle < 0) angle += 2*pi;

    //Now we find the multiplier: 
    float_T multiplier = 1;
    for (size_t i = 1; i <= numRegions; ++i) 
      if (angle < regionAngles[i] && angle >= regionAngles[i-1]) 
        multiplier = regionForceMultipliers[i-1];

    //Grabbing the unit vector associated with the mt...
    vec_T mtUnitVec = (mtEndPos[i] - basePos).normalize();
    //and adding it to the force vector. 
    force += multiplier*mtUnitVec;
  }
  //Finally, multiply the force vector by the force magnitude. 
  force *= forceMag;
}

void netMTForce(const char centrosome) {
  /* netMTForce: This is a wrapper for computing the net force on an MTOC. In
   *   reality, one only needs to specify for which MTOC they want to compute
   *   the forces in order to perform the computations in the function above, so
   *   this function takes that parameter, then fills in the remaining
   *   parameters automatically. 
   * Inputs: 
   *   const char centrosome: This parameter specifies which MTOC we care about
   *     here.
   * Outputs: (none, because we update the stored global force values instead)
   */
  //assert(centrosome == 'M' || centrosome == 'D')
  switch (centrosome) {
    case 'M':
      //If we care about the mother, then we need to use the M parameters, and
      //multiply the force magnitude by the Fratio parameter. 
      mtForceCalc(MT_Pos_M, basePosM, MT_Contact_M, force_M, Fratio*F_MT);
      break;
    case 'D':
      //If we care about the daughter, then we need to use the D parameters. 
      mtForceCalc(MT_Pos_D, basePosD, MT_Contact_D, force_D, F_MT);
      break;
  }
}

void updatePNPos() {
  /* updatePNPos: This functions updates the pronuclear position.
   * Inputs: (none)
   * Outputs: (none, though it changes various global parameters)
   * Issues: (TODO) This is rather inelegant, especially with all the repeated
   *   code. Maybe find a way to abstract some of it away? 
   */

  //First, we'll compute the radius vector of the Pronucleus, which arbitrarily
  //I've decided goes from the origin of the pronucleus to the mother spring
  //MTOC anchor point. 
  vec_T proNucRad({Prad*cos(psi), Prad*sin(psi)});
 // float_T cosinePrt  = Prad*cos(psi);
 // float_T sinePrt    = Prad*sin(psi);

  //Now, we zero out the force vector. 
  force.zero();

  //First to compute the net forces on the MTOCs (recall, these calls will
  //update the force_M and force_D values.
  netMTForce('M');
  netMTForce('D');
  //Now we need to care about which springs are on. 
  if (motherSpringOn) {
    //Now computing the position of the springAnchor...
    vec_T springAnchorM = proNucPos + proNucRad;
    //and the spring forces (which affect the centrosomes and
    //pronucleus evenly and oppositely)
    vec_T springForceM = kM*(springAnchorM - basePosM);
    //Adding the springForce into the MTOC force and pronuclear force...
    force_M += springForceM;
    force -= springForceM;
    
    //Updating Basepos: 
    // Base-pos obeys equation \eta_2 dx/dt = F + \xi,
    // where \xi is a random noise parameter on the order of \sqrt{2D \tau}. To
    // generate \xi, take a random number P from norm(0,1) and take 
    // P*\sqrt{2D\tau}. \eta_2 is a parameter to be set, but D = kB*T/\eta_2.
    // Right now, \eta_2 is one tenth the pronucleus drag coefficient. 
    vec_T randomVec({stdNormalDist(generator), stdNormalDist(generator)});
    vec_T xiM = sqrt(2*D*Tau)*randomVec;
    //Now we need to check and see if we must introduce a reflection force. 
    if (distBetween(proNucPos, basePosM) < Prad + 0.1) {
      //Reflection Force:
      vec_T radius = springAnchorM - proNucPos;
      //The reflection force is the opposite of the projectedForce of the MT:
      float_T forceProjCoeff = eucInnerProd(radius,force_M)/eucInnerProd(radius,radius);
      if (forceProjCoeff < 0) {
        vec_T reflectionForce = -(force_M.projectOn(radius));
        //The reflection force is impacted on the MTOC by the pronucleus.
        force_M += reflectionForce;
        force -= reflectionForce;
      }
      float_T randDispProjCoeff = eucInnerProd(radius,xiM)/eucInnerProd(radius,radius);
      if (randDispProjCoeff < 0) {
        vec_T reflectionDisp = -(xiM.projectOn(radius));
        //The reflection force is impacted on the MTOC by the pronucleus.
        xiM += reflectionDisp;
      }
    }
    basePosM += (1.0/Eta2)*Tau*force_M + xiM;

    //Calculating Torque:
    // To calculate the torque, note that we can see that the angle between
    // positive torque and the true x axis is \psi + \pi/2. Further, the angle
    // between positive torque and the true y is \psi. The daughter entails a
    // \pi switch, which merely negates things. Note that this could be
    // encapsulated via a cross product, but as we may move to 3d, that felt a
    // little silly. 
    //3D WARNING: This section would have to update when changing to 3d.
    torque_M = springForceM[0]*proNucRad[1] - springForceM[1]*proNucRad[0];
  } else {
    //With no springs, the MTOC's position is just a simple vector sum. 
    basePosM = proNucPos + proNucRad;
    if (spitValues) {
      cout << "Y-force on M: " << force_M[1] << endl;
    }
    force += force_M;

    torque_M = -force_M[0]*proNucRad[1] + force_M[1]*proNucRad[1];
  }
  //The next if clause is a near direct repeat of the mother clause above.
  //Hence, comments will be omitted. 
  if (daughterSpringOn) {
    vec_T springAnchorD = proNucPos - proNucRad;
    vec_T springForceD  = kD*(springAnchorD - basePosD);
    force_D += springForceD;
    force -= springForceD;

    vec_T randomVec({stdNormalDist(generator), stdNormalDist(generator)});
    vec_T xiD = sqrt(2*D*Tau)*randomVec;

    if (distBetween(proNucPos, basePosD) < Prad + 0.1) {
      vec_T radius = springAnchorD - proNucPos;
      float_T projectionCoefficient = eucInnerProd(radius,force_D)/eucInnerProd(radius,radius);
      if (projectionCoefficient < 0) {
        vec_T reflectionForce = -(force_D.projectOn(radius));
        force_D += reflectionForce;
        force -= reflectionForce;
      }
      float_T randDispProjCoeff = eucInnerProd(radius,xiD)/eucInnerProd(radius,radius);
      if (randDispProjCoeff < 0) {
        vec_T reflectionDisp = -(xiD.projectOn(radius));
        xiD += reflectionDisp;
      }
    }
    basePosD += (1.0/Eta2)*Tau*force_D + xiD;

    torque_D = -springForceD[0]*proNucRad[1] + springForceD[1]*proNucRad[0];
  } else {
    basePosD = proNucPos - proNucRad;
    if (spitValues) {
      cout << "Y-force on D: " << force_D[1] << endl;
    }
    force += force_D;
    torque_D = force_D[0]*proNucRad[1] - force_D[1]*proNucRad[0];
  }

  //Calculating Displacements:
  if (translation) 
    proNucPos += (1.0/Eta)*Tau*force;

  torque = torque_M + torque_D;
  psi += (1/Mu)*torque*Tau;
}

void advanceMT(const float_T vel, vec_T& vec, const float_T mag) {
  /* advanceMT: This grows or shrinks an MT according to velocity
   * Inputs: 
   *   const float_T vel: This is the growth velocity (can be negative)
   *   vec_T& vec: This is the vector of the MT, which must be updated via
   *     growing or shrinking.
   *   const float_T mag: This is a magnitude parameter. 
   * Output: (none, but vec is updated)
   */
  vec *= (1 + vel*Tau/mag);
}

float_T probContact(const float_T ang) {
  /* probContact: This computes the probability of successfully making contact
   *   at angle ang.
   * Inputs: 
   *   const float_T ang: This is the angle at which an MT is testing to make
   *     contact. 
   * Output: The probability of making contact at ang.
   */
  //First we check contact windows:
  for (size_t i = 1; i <= numberContactWindows; i++) {
    if (ang < contactWindowAngles[i]) {
      if (contacts[i-1]) return 0;
      else break;
    }
  }
  //Now we check the regions:
  for (size_t i = 1; i <= numRegions; i++) {
    if (ang < regionAngles[i]) {
      return regionProbabilities[i-1];
    }
  }
  //This is just a failsafe---the program should *never* get here. 
  cout << "This is a problem! I've worked through all regions and am still ";
  cout << "trying to determine the probability of contact at angle " << fixed;
  cout << ang << "! I'm going to just return 0!" << endl;
  return 0;
}

void addContact(const float_T angle) {
  /* addContact: This adds a contact at the specified angle.
   * Inputs:
   *   const float_T angle: This is the contact angle.
   * Output: (none, but globals are updated).
   */
  for (size_t i = 1; i <= numberContactWindows; i++) {
    if (angle < contactWindowAngles[i]){
      //We better not have already made contact! If the assert fails, the code
      //commented out below could help debugging.
      assert(!contacts[i-1]);
      //if (contacts[i-1]) {
      //  cout.precision(flt::digits10);
      //  cout << endl;
      //  cout << "Hey! I tried to add to an already contacted window?";
      //  cout << " Why is that!";
      //  cout << endl << "Angle of attempted contact: " << fixed << angle;
      //  cout << endl;
      //  cout << "Window in question: [" << fixed << contactWindowAngles[i-1];
      //  cout << ", " << contactWindowAngles[i] << "]." << endl;
      //  cout << "Setting it to true and moving on." << endl;
      //}
      //Making Contact!
      contacts[i-1] = true;
      break;
    }
  }
}

void removeContact(const float_T angle) {
  /* removeContact: This removes a contact at the specified angle.
   * Inputs:
   *   const float_T angle: This is the contact angle.
   * Output: (none, but globals are updated).
   */
  for (size_t i = 1; i <= numberContactWindows; i++) {
    if (angle < contactWindowAngles[i]){
      //We better have made contact here already! If the assert fails, the code
      //commented out below could help debugging.
      assert(contacts[i-1]);
      //if (!contacts[i-1]) {
      //  cout.precision(flt::digits10);
      //  cout << endl;
      //  cout << "Hey! I tried to remove a non-existent contact? Why is that!";
      //  cout << endl << "Angle of attempted contact: " << fixed << angle;
      //  cout << endl;
      //  cout << "Window in question: [" << fixed << contactWindowAngles[i-1];
      //  cout << ", " << contactWindowAngles[i] << "]." << endl;
      //  cout << "Setting it to false and moving on." << endl;
      //}
      //Remove the contact
      contacts[i-1] = false;
      break;
    }
  }
}


void mtContactTest(const char centrosome, const unsigned i) {
  /* mtContactTest: This performs a total test for MT contact for the ith MT off
   *   of the centrosome MTOC. Contact occurs if the test is successful, in this
   *   function. 
   * Inputs: 
   *   const char centrosome: This specifies from which MTOC the tested MT
   *     grows.
   *   const unsigned i: This specifies which MT is testing for contact.
   * Outputs: (none, but global values are changed). 
   */
  //First, declaring some angles we will specify. 
  float_T angleM, angleD;
  switch (centrosome) {
    //We need to determine which centrosome we care about here. 
    case 'M':
      angleM = atan2(MT_Pos_M[i][1],MT_Pos_M[i][0]);
      if (angleM < 0) angleM += 2*pi; //atan2 returns negative vals.
      if (testStat() < probContact(angleM)) {
        //cout << "Making Contact at Angle " << angleM << endl;
        MT_GrowthVel_M[i] = 0;
        MT_Contact_M[i]   = contact_length;
        MT_Growing_M[i]   = false;
        addContact(angleM);
      } else {
        MT_GrowthVel_M[i] = -Vs_c;
        MT_Growing_M[i]   = false;
      }
      break;
    case 'D':
      angleD = atan2(MT_Pos_D[i][1],MT_Pos_D[i][0]);
      if (angleD < 0) angleD += 2*pi; //atan2 returns negative vals.
      if (testStat() < probContact(angleD)) {
        //cout << "Making Contact at Angle " << angleD << endl;
        MT_GrowthVel_D[i] = 0;
        MT_Contact_D[i]   = contact_length;
        MT_Growing_D[i]   = false;
        addContact(angleD);
      } else {
        MT_GrowthVel_D[i] = -Vs_c;
        MT_Growing_D[i]   = false;
      }
      break;
  }
}

void respawnMTB(vec_T& vec, const float_T ang, const float_T envelope[2]) {
  float_T randR = testStat();  
  float_T randT = testStat();
  float_T r     = sqrt(randR);
  float_T t     = (envelope[1]-envelope[0])*randT + ang + envelope[0] - pi/2.0;
  vec[0]        = r*cos(t);
  vec[1]        = r*sin(t);
}

void respawnMT(const char centrosome, vec_T& vec, const unsigned i, const float_T envelope[2]) {
  //assert(centrosome == 'M' || centrosome == 'D')
  switch (centrosome) {
    case 'M':
      respawnMTB(vec, psi, envelope);
      MT_Growing_M[i]   = true;
      MT_GrowthVel_M[i] = Vg;
      MT_Contact_M[i]   = 0;
      break;
    case 'D':
      respawnMTB(vec, psi + pi, envelope);
      MT_Growing_D[i]   = true;
      MT_GrowthVel_D[i] = Vg;
      MT_Contact_D[i]   = 0;
      break;
  }
}

bool checkBoundary() {
    float_T dist = sqrt(pow(proNucPos[0]/(R1_max - Prad - .5),2) +
                        pow(proNucPos[1]/(R2_max-Prad - .5),2));
    return dist >= 1;
}

void runModel(bool writeAllData, bool writeTempData) {
  float_T nextWrite = 0;
  float_T writeInterval = 0.5;
  if (writeAllData) {
    writeData(0);
  } else if (writeTempData) {
    writePartialData(0);
    nextWrite += writeInterval;
  }
  for (float_T t=0; t <= Duration; t += Tau) {
    if (t > 400*Tau) {
      spitValues = false;
    }
    if (checkBoundary())
      break;
    while (psi >= 2*pi) {
      psi -= 2*pi;
    } 
    while (psi < 0) {
      psi += 2*pi;
    }
    //As MT_numb \to \infty, we'll have to calculate force/torque more and more
    //often. As such, we'll just go ahead and calculate it every time here. 
    updatePNPos();
    
    //Growing and Shrinking of MTs/State Change
    for (size_t i = 0; i < MT_numb; ++i) {
      //These are the relative positions of the MTs
      vec_T vecM, vecD;
      //Setting them properly. 
      vecM[0] = MT_Pos_M[i][0] - basePosM[0];
      vecM[1] = MT_Pos_M[i][1] - basePosM[1];
      vecD[0] = MT_Pos_D[i][0] - basePosD[0];
      vecD[1] = MT_Pos_D[i][1] - basePosD[1];

      //Grabbing their relative magnitudes. 
      float_T mag_M = sqrt(pow(vecM[0], 2) + pow(vecM[1], 2));
      float_T mag_D = sqrt(pow(vecD[0], 2) + pow(vecD[1], 2));

      //Respawning if necessary. First, we need to compute the angular position
      //of the MT relative to the pronucleus. 
      float_T angleM = atan2(vecM[1],vecM[0]);
      if (angleM < 0) angleM += 2*pi;
      float_T angleD = atan2(vecD[1],vecD[0]);
      if (angleD < 0) angleD += 2*pi;
      //angleM & angleD give us the cartesian angle of the MT relative to the
      //Pronucleus, but we really want to use a rotated axis system that aligns
      //with the pronucleus. The Mt base points are offset by $\pi/2$, and we
      //have $\psi$, so some simple geometry shows
      float_T thetaM = angleM - (psi - pi/2.0);
      float_T thetaD = angleD - (psi + pi/2.0);

      if (mag_M < 0.1) {
        //TODO if ((thetaM < envelopeM[0]) || (thetaM > envelopeM[1])))
        if (MT_Contact_M[i] > 0)
          removeContact(angleM);
        respawnMT('M', vecM, i, envelopeM);
      }
      if (mag_D < 0.1) {
        //if ((thetaD < envelopeD[0]) || (thetaD > envelopeD[1])))
        if (MT_Contact_D[i] > 0) 
          removeContact(angleD);
        respawnMT('D', vecD, i, envelopeD);
      }

      //Growing or Shrinking the MT. 
      advanceMT(MT_GrowthVel_M[i], vecM, mag_M);
      advanceMT(MT_GrowthVel_D[i], vecD, mag_D);

      //Updating the endpoint positions.
      MT_Pos_M[i][0] = basePosM[0] + vecM[0];
      MT_Pos_M[i][1] = basePosM[1] + vecM[1];
      MT_Pos_D[i][0] = basePosD[0] + vecD[0];
      MT_Pos_D[i][1] = basePosD[1] + vecD[1];

      //Now, for state updating, I'll need a random number. 
      float_T test = testStat();
      if (MT_Growing_M[i]) {
        if (test < Pr_catastrophe && mag_M > 0.4) { //TODO: Why is the mag_M thing here?
          MT_Growing_M[i]   = false;
          MT_GrowthVel_M[i] = -Vs;
        }
      } else {
        if (test < Pr_rescue && !(MT_Contact_M[i] > 0)) {
          MT_Growing_M[i]   = true;
          MT_GrowthVel_M[i] = Vg;
        }
      }
      test = testStat();
      if (MT_Growing_D[i]) {
        if (test < Pr_catastrophe && mag_D > 0.4) {
          MT_Growing_D[i]   = false;
          MT_GrowthVel_D[i] = -Vs;
        }
      } else {
        if (test < Pr_rescue && !(MT_Contact_D[i] > 0)) {
          MT_Growing_D[i]   = true;
          MT_GrowthVel_D[i] = Vg;
        }
      }

      //Computing their scaled real magnitudes (scaled via the ellipse). 
      float_T magScaledM=sqrt(pow(MT_Pos_M[i][0]/R1_max,2)+pow(MT_Pos_M[i][1]/R2_max,2));
      float_T magScaledD=sqrt(pow(MT_Pos_D[i][0]/R1_max,2)+pow(MT_Pos_D[i][1]/R2_max,2));

      //Updating Contact Indicator:
      if (MT_Contact_M[i] > 0) {
        MT_Contact_M[i] -= Tau;
        if (MT_Contact_M[i] <= 0) {
          MT_Contact_M[i] = 0;
          MT_GrowthVel_M[i] = -Vs_c;
          float_T angleM = atan2(MT_Pos_M[i][1],MT_Pos_M[i][0]);
          if (angleM < 0) angleM += 2*pi;
          removeContact(angleM);
        }
      } else  if (magScaledM >= 1) {
        mtContactTest('M',i);
      }
      if (MT_Contact_D[i] > 0) {
        MT_Contact_D[i] -= Tau;
        if (MT_Contact_D[i] <= 0) {
          MT_Contact_D[i] = 0;
          MT_GrowthVel_D[i] = -Vs_c;
          float_T angleD = atan2(MT_Pos_D[i][1],MT_Pos_D[i][0]);
          if (angleD < 0) angleD += 2*pi;
          removeContact(angleD);
        }
      } else if (magScaledD >= 1) {
        mtContactTest('D',i);
      }
    }
    if (writeAllData) {
      writeData(t+Tau);
    } else if (writeTempData && (t > nextWrite)) {
      writePartialData(t+Tau);
      nextWrite += writeInterval;
    }
  }
  if (!(writeAllData || writeTempData)) {
    writeFinalData();
  }
}

void usage() {
  cout << "Usage: " << endl;
  cout << "./main [number of runs = 1] [fileName = {MT#}-MT.csv]" << endl;
}

int numContacts() {
  unsigned int count = 0;
  for (size_t i = 0; i < numberContactWindows; i++)
    count += contacts[i];
  return count;
}

void test() {
  cout << "Currently, there are " << numContacts() << " contacts." << endl;
  setup();
  cout << "After setup(), there are " << numContacts() << " contacts." << endl;
  vector<float_T> anglesToAdd = {2,4,6,5,3,1,0.6};
  vector<float_T> anglesToRemove = {0.6,1,3,4,6,5,2};
  for (size_t i = 0; i < anglesToAdd.size(); i++) {
    float_T angle = anglesToAdd[i];
    cout << "Angle: " << angle << endl;
    cout << "Probability of Contact from M: " << probContact(angle) << endl;
    cout << "Probability of Contact from D: " << probContact(angle) << endl;
    addContact(angle);
    cout << "Number of Contacts after Making Contact: " << numContacts()<<endl;
  }
}

int main(int argc, const char* argv[]) {
  int numRuns;
  bool writeAllData = false;
  bool writeTempData = false;
  string writeAllDataKey = "all";
  string writeTempDataKey = "temp";
  if (argc == 2) {
    numRuns = atoi(argv[1]);
  } else if (argc == 3) {
    numRuns = atoi(argv[1]);
    fileName = argv[2];
  } else if (argc == 4) {
    numRuns = atoi(argv[1]);
    fileName = argv[2];
    if (argv[3] == writeAllDataKey) {
      writeAllData = true;
    } else if (argv[3] == writeTempDataKey) {
      cout << "Writing Temporary Data" << endl;
      writeTempData = true;
    }
  } else {
    usage();
    return 1;
  }
  setup();

  for (int run = 0; run < numRuns; ++run) {
    runModel(writeAllData,writeTempData);
    setToBasePos();
  }
  //test();
  return 0;
}
