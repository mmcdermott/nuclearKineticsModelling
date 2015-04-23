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

ostream& operator<<(ostream& out, const vec_T (&rhs)[MT_numb_M]) {
  /* operator<<(ostream, vec_T[MT_numb_M]): A wrapper to write an array of vectors
   *   to an output stream.
   * Inputs: 
   *   ostream& out: The stream to which rhs should be written. 
   *   const vec_T (&rhs)[MT_numb_M]: The array of vectors that should be written
   * Output: The ostream out. 
   * TODO: It is really dumb that we have two versions of this function. The MTs
   * should be encapsulated as objects, and then this would be MUCH cleaner and
   * easier to maintain.
   */
  out << rhs[0];
  for (size_t i = 1; i < MT_numb_M; ++i) {
    if (ONLY_COMMA) {
      //There are two modes of writing. This is by far the most common. 
      out << "," << rhs[i];
    } else {
      out << ";" << rhs[i];
    }
  }
  return out;
}
// TODO THIS IS REALLY BAD!!! THIS CODE NEEDS TO BE IMPROVED. ... but it should
// work. In the meantime, until it gets really fixed, if you need to change
// this, be careful you know what you're doing with pre-processor macros.
#if MT_numb_M != MT_numb_D
ostream& operator<<(ostream& out, const vec_T (&rhs)[MT_numb_D]) {
  /* operator<<(ostream, vec_T[MT_numb_D]): A wrapper to write an array of vectors
   *   to an output stream.
   * Inputs: 
   *   ostream& out: The stream to which rhs should be written. 
   *   const vec_T (&rhs)[MT_numb_D]: The array of vectors that should be written
   * Output: The ostream out. 
   * TODO: It is really dumb that we have two versions of this function. The MTs
   * should be encapsulated as objects, and then this would be MUCH cleaner and
   * easier to maintain.
   */
  out << rhs[0];
  for (size_t i = 1; i < MT_numb_D; ++i) {
    if (ONLY_COMMA) {
      //There are two modes of writing. This is by far the most common. 
      out << "," << rhs[i];
    } else {
      out << ";" << rhs[i];
    }
  }
  return out;
}
#endif

ostream& operator<<(ostream& out, const float_T (&rhs)[MT_numb_M]) {
  /* operator<<(ostream, float_T[MT_numb_M]): A wrapper to write an array of floats
   *   to an output stream.
   * Inputs: 
   *   ostream& out: The stream to which rhs should be written. 
   *   const vec_T (&rhs)[MT_numb_M]: The array of floats that should be written
   * Output: The ostream out. 
   * TODO: It is really dumb that we have two versions of this function. The MTs
   * should be encapsulated as objects, and then this would be MUCH cleaner and
   * easier to maintain.
   */
  out << rhs[0];
  for (size_t i = 1; i < MT_numb_M; ++i) {
    if (ONLY_COMMA) {
      //There are two modes of writing. This is by far the most common. 
      out << "," << rhs[i];
    } else {
      out << ";" << rhs[i];
    }
  }
  return out;
}
// TODO THIS IS REALLY BAD!!! THIS CODE NEEDS TO BE IMPROVED. ... but it should
// work. In the meantime, until it gets really fixed, if you need to change
// this, be careful you know what you're doing with pre-processor macros.
#if MT_numb_M != MT_numb_D
ostream& operator<<(ostream& out, const float_T (&rhs)[MT_numb_D]) {
  /* operator<<(ostream, float_T[MT_numb_D]): A wrapper to write an array of floats
   *   to an output stream.
   * Inputs: 
   *   ostream& out: The stream to which rhs should be written. 
   *   const vec_T (&rhs)[MT_numb_D]: The array of floats that should be written
   * Output: The ostream out. 
   * TODO: It is really dumb that we have two versions of this function. The MTs
   * should be encapsulated as objects, and then this would be MUCH cleaner and
   * easier to maintain.
   */
  out << rhs[0];
  for (size_t i = 1; i < MT_numb_D; ++i) {
    if (ONLY_COMMA) {
      //There are two modes of writing. This is by far the most common. 
      out << "," << rhs[i];
    } else {
      out << ";" << rhs[i];
    }
  }
  return out;
}
#endif

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

void mtForceCalcM()
{
  /* mtForceCalcM: A function which computes the force on the MTOC M implied by the
   *   cortical pushing and pulling forces due to the MTs. 
   * Inputs: 
   *   const float_T forceMag: This is the magnitude of the force each MT should
   *     impart when being pulled by cortical dynein. Pushing is handled by
   *     using a multiplier, so this typically also is the magnitude of any
   *     cortical pushing forces as well. 
   * Output: (none)
   * TODO: It is really dumb that we have two versions of this function. The MTs
   * should be encapsulated as objects, and then this would be MUCH cleaner and
   * easier to maintain.
   */
  //First we'll zero out the forces, because we're about to compute them from
  //scratch. 
  force_M.zero();
  for (size_t i=0; i < MT_numb_M; ++i) {
    //We only need to compute the force if the MT's making contact.
      
      //Continue if neither a pulling or pushing contact is not there
      if (MT_Contact_M[i] == 0 && MT_Contact_M_push[i] == 0) continue;
      
      float_T multiplier = 0; //Set this to null just in case
      
    //We need to search to see what multiplier we should assign, which requires
    //knowing the true cartesian angle of the vector. 
    //3D WARNING: This section would have to update when changing to 3d.
    float_T angle = atan2(MT_Pos_M[i][1],MT_Pos_M[i][0]);
      
    if (angle < 0) angle += 2*pi;

      //Now we find the multiplier:
      if (MT_Contact_M[i]>0)  multiplier = F_MT_pull; //Pull force multiplier
      if (MT_Contact_M_push[i]>0)  multiplier = -1*F_MT_push; //Push force multiplier
      
    //for (size_t i = 1; i <= numRegions; ++i)
     // if (angle < regionAngles[i] && angle >= regionAngles[i-1])
      //   multiplier = regionForceMultipliers[i-1];

    //Grabbing the unit vector associated with the mt...
    vec_T mtUnitVec = (MT_Pos_M[i] - basePosM).normalize();
    //and adding it to the force vector. 
    force_M += multiplier*mtUnitVec;
  }
  //Finally, multiply the force vector by the force magnitude. 
  //force_M *= forceMag; not needed any more
}
void mtForceCalcD()
{
  /* mtForceCalcM: A function which computes the force on the MTOC D implied by the
   *   cortical pushing and pulling forces due to the MTs. 
   * Inputs: 
   *   const float_T forceMag: This is the magnitude of the force each MT should
   *     impart when being pulled by cortical dynein. Pushing is handled by
   *     using a multiplier, so this typically also is the magnitude of any
   *     cortical pushing forces as well. 
   * Output: (none)
   * TODO: It is really dumb that we have two versions of this function. The MTs
   * should be encapsulated as objects, and then this would be MUCH cleaner and
   * easier to maintain.
   */
  //First we'll zero out the forces, because we're about to compute them from
  //scratch. 
  force_D.zero();
  for (size_t i=0; i < MT_numb_D; ++i) {
    //We only need to compute the force if the MT's making contact.
      
      //Continue if a pulling contact is not there
    if (MT_Contact_D[i] == 0 && MT_Contact_D_push[i] == 0) continue;
     
      float_T multiplier = 0;
    //We need to search to see what multiplier we should assign, which requires
    //knowing the true cartesian angle of the vector. 
    //3D WARNING: This section would have to update when changing to 3d.
    float_T angle = atan2(MT_Pos_D[i][1],MT_Pos_D[i][0]);
    if (angle < 0) angle += 2*pi;

    //Now we find the multiplier:
      if (MT_Contact_D[i]>0)  multiplier = F_MT_pull; //Pull force multiplier
      if (MT_Contact_D_push[i]>0)  multiplier = -1*F_MT_push; //Push force multiplier
      
    //Now check to see if we are in a push band region
    //for (size_t i = 1; i <= numRegions; ++i)
    //  if (angle < regionAngles[i] && angle >= regionAngles[i-1])
      //  multiplier = regionForceMultipliers[i-1];

    //Grabbing the unit vector associated with the mt...
    vec_T mtUnitVec = (MT_Pos_D[i] - basePosD).normalize();
    //and adding it to the force vector. 
    force_D += multiplier*mtUnitVec;
  }
  //Finally, multiply the force vector by the force magnitude.
    
  //force_D *= forceMag; not needed any more!
}

void netMTForce(const MTOC centrosome) {
  /* netMTForce: This is a wrapper for computing the net force on an MTOC. In
   *   reality, one only needs to specify for which MTOC they want to compute
   *   the forces in order to perform the computations in the function above, so
   *   this function takes that parameter, then fills in the remaining
   *   parameters automatically. 
   * Inputs: 
   *   const MTOC centrosome: This parameter specifies which MTOC we care about
   *     here.
   * Outputs: (none, because we update the stored global force values instead)
   */
  //assert(centrosome == 'M' || centrosome == 'D')
  switch (centrosome) {
    case M_CENTROSOME:
      //If we care about the mother, then we need to use the M parameters, and
      //multiply the force magnitude by the Fratio parameter. 
      mtForceCalcM();
      break;
    case D_CENTROSOME:
      //If we care about the daughter, then we need to use the D parameters. 
      mtForceCalcD();
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

  //First, we'll compute the radius vector of the Pronucleus, which I've 
  //arbitrarily decided goes from the origin of the pronucleus to the mother
  //spring MTOC anchor point. 
  vec_T proNucRad({Prad*cos(psi), Prad*sin(psi)});

  //Now, we zero out the force vector. 
  force.zero();

  //First to compute the net forces on the MTOCs (recall, these calls will
  //update the force_M and force_D values.
  netMTForce(M_CENTROSOME);
  netMTForce(D_CENTROSOME);
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

    torque_M = -force_M[0]*proNucRad[1] + force_M[1]*proNucRad[0];
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
  //Now we check the push band regions:
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
      // cout << "Hey! I tried to remove a non-existent contact? Why is that!";
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


void mtContactTest(const MTOC centrosome, const unsigned i) {
  /* mtContactTest: This performs a total test for MT contact for the ith MT off
   *   of the centrosome MTOC. Contact occurs if the test is successful, in this
   *   function. 
   * Inputs: 
   *   const MTOC centrosome: This specifies from which MTOC the tested MT
   *     grows.
   *   const unsigned i: This specifies which MT is testing for contact.
   * Outputs: (none, but global values are changed). 
   */

  //First, let's declare some angles we'll use. 
  float_T angleM, angleD;
  switch (centrosome) {
    //We need to determine which centrosome we care about here. 
    case M_CENTROSOME:
      //3D WARNING: The code below is used to determine the point of contact on
      //the cortex. It only works in 2D right now. 
      angleM = atan2(MT_Pos_M[i][1],MT_Pos_M[i][0]);
      if (angleM < 0) angleM += 2*pi; //atan2 returns negative vals.
      if (testStat() < probContact(angleM)) {
          
        MT_GrowthVel_M[i] = 0;
        MT_Contact_M[i] = contact_length_dynein;
        MT_Growing_M[i]   = false;
        addContact(angleM); // This counts the number of pull contacts? It should reject multiple pull requests in one window
          
      } else {
        
          MT_GrowthVel_M[i] = 0; //Stop growing but still hang on since MT is pushing
          MT_Contact_M_push[i] = contact_length_push;
          MT_Growing_M[i]   = false;
          
      }
      break;
    case D_CENTROSOME:
      //3D WARNING: The code below is used to determine the point of contact on
      //the cortex. It only works in 2D right now. 
      angleD = atan2(MT_Pos_D[i][1],MT_Pos_D[i][0]);
      if (angleD < 0) angleD += 2*pi; //atan2 returns negative vals.
      if (testStat() < probContact(angleD)) {
          
        MT_GrowthVel_D[i] = 0;
        MT_Contact_D[i] = contact_length_dynein;
        MT_Growing_D[i]   = false;
        addContact(angleD);
          
      } else {
        
          MT_GrowthVel_D[i] = 0; //Stop growing but still hang on
          MT_Contact_D_push[i] = contact_length_push;
          MT_Growing_D[i]   = false;
         
          
      }
      break;
  }
}

void respawnMTB(vec_T& vec, const float_T ang, const float_T envelope[2]) {
  /* respawnMTB: This function is the base respawn function for an MT. We have
   *   wrappers to call this for an explicit MT of MTOC M or D below. 
   * Inputs: 
   *   vec_T& vec: This is the MT vector that is to be respawned. 
   *   const float_T ang: This is the angular coordinate of the appropriate
   *     MTOC.
   *   const float_T envelope[2]: This stores the angular envelope for MT
   *     growth. 
   * Outputs: (none, but the reference vec is updated). 
   */

  float_T randR = testStat();  
  float_T randT = testStat();
  float_T r     = sqrt(randR);
  float_T t     = (envelope[1]-envelope[0])*randT + ang + envelope[0] - pi/2.0;
  vec           = Vector({r*cos(t),r*sin(t)});
}

void respawnMT(const MTOC centrosome, vec_T& vec, const unsigned int i,
               const float_T envelope[2]) {
  /* respawnMT: This function is the MTOC explicit wrapper function for
   *   respawnMTB above.  
   * Inputs: 
   *   const MTOC centrosome: This specifies the centrosome of the MT in
   *     question. 
   *   const unsigned int i: This specifies which MT of the many associated with
   *     MTOC centrosome is respawning. 
   *   const float_T envelope[2]: This specifies the allowed angular envelope of
   *     MT growth. 
   * Outputs: (none, but globals are updated). 
   */

  //The centrosome better be either M or D. 
  switch (centrosome) {
    case M_CENTROSOME:
      respawnMTB(vec, psi, envelope);
      MT_Growing_M[i]   = true;
      MT_GrowthVel_M[i] = Vg;
      if (MT_Contact_M[i]>0) MT_Contact_M[i]   = 0;
      if (MT_Contact_M_push[i]>0) MT_Contact_M_push[i] = 0;
      break;
    case D_CENTROSOME:
      respawnMTB(vec, psi + pi, envelope);
      MT_Growing_D[i]   = true;
      MT_GrowthVel_D[i] = Vg;
      if (MT_Contact_D[i]>0) MT_Contact_D[i]   = 0;
      if (MT_Contact_D_push[i]>0) MT_Contact_D_push[i] = 0;
      break;
  }
}

inline bool checkBoundary() {
  /* checkBoundary: This functions estimates if the pronucleus has impacted on
   *   the cortex. It estimates because that typically does a good enough job,
   *   and doing the exact calculation would be very difficult. 3D WARNING: This
   *   function is explicitly a 2D implementation. It's inline because its very
   *   short and called in the loop, so making it inline will make the program
   *   slightly more efficient. 
   * Inputs: (none)
   * Output: If the pronucleus has impacted the boundary.
   */
  float_T dist = sqrt(pow(proNucPos[0]/(R1_max - Prad - .5),2) +
                      pow(proNucPos[1]/(R2_max-Prad - .5),2));
  return dist >= 1;
}

void runModel(bool writeAllData, bool writeTempData) {
  /* runModel: This function actually runs the model and records the data as
   *   specified by the user. 
   * Inputs: 
   *   bool writeAllData: This specifies whether or not to write all data.
   *   bool writeTempData: This specifies whether or not to write data every 30
   *     seconds or not. 
   * Output: (none, but globals are altered)
   * Issues: The handling of data writing parameters is inefficient, but not
   *   enough to cause any issues, really. 
   */

  //First, we set up the data writing parameters: 
  float_T nextWrite = 0;
  float_T writeInterval = 0.5;
  if (writeAllData) {
    //We only need to write the full data if writeAllData is true. 
    writeData(0);
  } else if (writeTempData) {
    //If we're writing temporary data, then we can just write the partial data
    //and then note when we'll record data next. 
    writePartialData(0);
    nextWrite += writeInterval;
  }

  //Now we run the model. 
  for (float_T t=0; t <= Duration; t += Tau) {
    //Even if debugging info is desired, it will rarely be so after 400
    //iterations...
    if (t > 400*Tau) {
      spitValues = false;
    }

    //We need to ensure we're still in the cell; if not, we should exit. TODO:
    //improve these exit conditions, and/or induce a reflection force. 
    if (checkBoundary())
      break;

    //Now we normalize the angle. 
    while (psi >= 2*pi)
      psi -= 2*pi;
    while (psi < 0)
      psi += 2*pi;

    //As MT_numb \to \infty, we'll have to calculate force/torque more and more
    //often. As such, we'll just go ahead and calculate it every time here. 
    updatePNPos();
    
    //Growing and Shrinking of MTs/State Change
      
      //Loop through the mother cell MTs
    for (size_t i = 0; i < MT_numb_M; ++i) {
      //These are the relative positions of the MTs
      vec_T vecM;
      //Setting them properly. 
      vecM = MT_Pos_M[i] - basePosM;

      //Grabbing their relative magnitudes. 
      float_T mag_M = vecM.norm();

      //Respawning if necessary. First, we need to compute the angular position
      //of the MT relative to the pronucleus. 3D WARNING: This is explicitly 2D,
      //by relying on and computing angles as opposed to general cortex
      //positions. 
      float_T angleM = atan2(vecM[1],vecM[0]);
      if (angleM < 0) angleM += 2*pi;
      //angleM gives us the cartesian angle of the MT relative to the
      //Pronucleus, but we really want to use a rotated axis system that aligns
      //with the pronucleus. The Mt base points are offset by $\pi/2$, and we
      //have $\psi$, so some simple geometry shows
      //float_T thetaM = angleM - (psi - pi/2.0);

      //MT-ENV: Different MT growth envelope handling would go here, though it
      //currently isn't used beyond spawning new MTS.
      if (mag_M < 0.1) {
        if (MT_Contact_M[i] > 0)
          removeContact(angleM);
          
        respawnMT(M_CENTROSOME, vecM, i, envelopeM);
          
         // if (MT_Contact_M_push[i]>0)
        //respawnMT(M_CENTROSOME, vecM, i, envelopeM);
      }

      //Growing or Shrinking the MT. 
      advanceMT(MT_GrowthVel_M[i], vecM, mag_M);

      //Updating the endpoint positions.
      MT_Pos_M[i][0] = basePosM[0] + vecM[0];
      MT_Pos_M[i][1] = basePosM[1] + vecM[1];

      //Now, for state updating, I'll need a random number. 
      //For modularity, I'll declare the stat first. 
      float_T test;

      test = testStat();
      if (MT_Growing_M[i]) {
        if (test < Pr_catastrophe && mag_M > 0.4) { 
          //TODO: Is the length condition needed here? Why can't they
          //catastrophe if they are too short?
          MT_Growing_M[i]   = false;
          MT_GrowthVel_M[i] = -Vs;
        }
      } else {
        if (test < Pr_rescue && !(MT_Contact_M[i] > 0)) {
          MT_Growing_M[i]   = true;
          MT_GrowthVel_M[i] = Vg;
        }
      }

      //Computing their scaled real magnitudes (scaled via the ellipse), so as
      //to determine if a contact has occurred.  
      float_T magScaledM = sqrt(pow(MT_Pos_M[i][0]/R1_max,2) + 
                                pow(MT_Pos_M[i][1]/R2_max,2));

      //Updating Contact Indicator:
      if (MT_Contact_M[i] > 0) {
        //If we've already made contact, then we need to reduce the contact time
        //by the length of the time step. 
        MT_Contact_M[i] -= Tau;
          
        if (MT_Contact_M[i] <= 0) {
          //But, occasionally, this will actually drop our proposed contact time
          //down into the negatives. We need to correct for that. Moreover, we
          //can combine that with the case where we hit zero right on the head,
          //and just handle the overall contact-ended case here. 
          MT_Contact_M[i] = 0;
          MT_GrowthVel_M[i] = -Vs_c;
          float_T angleM = atan2(MT_Pos_M[i][1],MT_Pos_M[i][0]);
          if (angleM < 0) angleM += 2*pi;
          removeContact(angleM);
        }
      } else if (MT_Contact_M_push[i] > 0) {
          //If we've already made contact, then we need to reduce the contact time
          //by the length of the time step.
          MT_Contact_M_push[i] -= Tau;
          if (MT_Contact_M_push[i] <= 0) {
              //But, occasionally, this will actually drop our proposed contact time
              //down into the negatives. We need to correct for that. Moreover, we
              //can combine that with the case where we hit zero right on the head,
              //and just handle the overall contact-ended case here.
              MT_Contact_M_push[i] = 0;
              MT_GrowthVel_M[i] = -Vs_c;
              float_T angleM = atan2(MT_Pos_M[i][1],MT_Pos_M[i][0]);
              if (angleM < 0) angleM += 2*pi;
          }
      } else  if (magScaledM >= 1) {
        //If we haven't made contact, but are touching the cortex, we will test
        //for contact. 
        mtContactTest(M_CENTROSOME, i);
      }
    
        
        //Updating Contact Indicator:
      //  if (MT_Contact_M_push[i] > 0) {
            //If we've already made contact, then we need to reduce the contact time
            //by the length of the time step.
        //    MT_Contact_M_push[i] -= Tau;
          //  if (MT_Contact_M_push[i] <= 0) {
                //But, occasionally, this will actually drop our proposed contact time
                //down into the negatives. We need to correct for that. Moreover, we
                //can combine that with the case where we hit zero right on the head,
                //and just handle the overall contact-ended case here.
            //    MT_Contact_M_push[i] = 0;
              //  MT_GrowthVel_M[i] = -Vs_c;
              //  float_T angleM = atan2(MT_Pos_M[i][1],MT_Pos_M[i][0]);
               // if (angleM < 0) angleM += 2*pi;
            //}
            //} else  if (magScaledM >= 1) {
                //If we haven't made contact, but are touching the cortex, we will test
                //for contact. 
           //     mtContactTest(M_CENTROSOME, i);
            //}
    }
    
      
    //Loop through the daughter cell MTs
    for (size_t i = 0; i < MT_numb_D; ++i) {
      //These are the relative positions of the MTs
      vec_T vecD;
      //Setting them properly. 
      vecD = MT_Pos_D[i] - basePosD;

      //Grabbing their relative magnitudes. 
      float_T mag_D = vecD.norm();

      //Respawning if necessary. First, we need to compute the angular position
      //of the MT relative to the pronucleus. 3D WARNING: This is explicitly 2D,
      //by relying on and computing angles as opposed to general cortex
      //positions. 
      float_T angleD = atan2(vecD[1],vecD[0]);
      if (angleD < 0) angleD += 2*pi;
      //angleD gives us the cartesian angle of the MT relative to the
      //Pronucleus, but we really want to use a rotated axis system that aligns
      //with the pronucleus. The Mt base points are offset by $\pi/2$, and we
      //have $\psi$, so some simple geometry shows
      //float_T thetaD = angleD - (psi + pi/2.0);

      //MT-ENV: Different MT growth envelope handling would go here, though it
      //currently isn't used beyond spawning new MTS.
      if (mag_D < 0.1) {
          
        if (MT_Contact_D[i] > 0) 
          removeContact(angleD);
          respawnMT(D_CENTROSOME, vecD, i, envelopeD);
          
         // if (MT_Contact_D_push[i]>0)
         // respawnMT(D_CENTROSOME, vecD, i, envelopeD);
          
      }

      //Growing or Shrinking the MT. 
      advanceMT(MT_GrowthVel_D[i], vecD, mag_D);

      //Updating the endpoint positions.
      MT_Pos_D[i][0] = basePosD[0] + vecD[0];
      MT_Pos_D[i][1] = basePosD[1] + vecD[1];

      //Now, for state updating, I'll need a random number. 
      //For modularity, I'll declare the stat first. 
      float_T test;

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

      //Computing their scaled real magnitudes (scaled via the ellipse), so as
      //to determine if a contact has occurred.  
      float_T magScaledD = sqrt(pow(MT_Pos_D[i][0]/R1_max,2) + 
                                pow(MT_Pos_D[i][1]/R2_max,2));

      //Updating Contact Indicator:
      //The if clause below is identical to the one above, just for the D
      //centrosome over the M centrosome. 
        //Updating Contact Indicator:
        if (MT_Contact_D[i] > 0) {
            //If we've already made contact, then we need to reduce the contact time
            //by the length of the time step.
            MT_Contact_D[i] -= Tau;
            if (MT_Contact_D[i] <= 0) {
                //But, occasionally, this will actually drop our proposed contact time
                //down into the negatives. We need to correct for that. Moreover, we
                //can combine that with the case where we hit zero right on the head,
                //and just handle the overall contact-ended case here.
                MT_Contact_D[i] = 0;
                MT_GrowthVel_D[i] = -Vs_c;
                float_T angleD = atan2(MT_Pos_D[i][1],MT_Pos_D[i][0]);
                if (angleD < 0) angleD += 2*pi;
                removeContact(angleD);
            }
        } else if (MT_Contact_D_push[i] > 0){
            //If we've already made contact, then we need to reduce the contact time
            //by the length of the time step.
            MT_Contact_D_push[i] -= Tau;
            if (MT_Contact_D_push[i] <= 0) {
                //But, occasionally, this will actually drop our proposed contact time
                //down into the negatives. We need to correct for that. Moreover, we
                //can combine that with the case where we hit zero right on the head,
                //and just handle the overall contact-ended case here.
                MT_Contact_D_push[i] = 0;
                MT_GrowthVel_D[i] = -Vs_c;
                float_T angleD = atan2(MT_Pos_D[i][1],MT_Pos_D[i][0]);
                if (angleD < 0) angleD += 2*pi;
                
            }
        } else  if (magScaledD >= 1) {
                //If we haven't made contact, but are touching the cortex, we will test
                //for contact. 
                mtContactTest(D_CENTROSOME, i);
            }
        
        
        //Updating Contact Indicator:
       // if (MT_Contact_D_push[i] > 0) {
            //If we've already made contact, then we need to reduce the contact time
            //by the length of the time step.
         //   MT_Contact_D_push[i] -= Tau;
           // if (MT_Contact_D_push[i] <= 0) {
                //But, occasionally, this will actually drop our proposed contact time
                //down into the negatives. We need to correct for that. Moreover, we
                //can combine that with the case where we hit zero right on the head,
                //and just handle the overall contact-ended case here.
             //   MT_Contact_D_push[i] = 0;
               // MT_GrowthVel_D[i] = -Vs_c;
              //  float_T angleD = atan2(MT_Pos_D[i][1],MT_Pos_D[i][0]);
                //if (angleD < 0) angleD += 2*pi;
               
              //}
            //} else  if (magScaledD >= 1) {
                //If we haven't made contact, but are touching the cortex, we will test
                //for contact.
              //  mtContactTest(D_CENTROSOME, i);
            //}
    }
      
    //Now we just need to write our data for this run (if we're supposed to).
    
    if (writeAllData) {
      writeData(t + Tau);
    } else if (writeTempData && (t > nextWrite)) {
      writePartialData(t + Tau);
      nextWrite += writeInterval;
    }
   }
  //We've finished running the model; time to (potentially) write the final
  //data. 
  if (!(writeAllData || writeTempData)) {
    writeFinalData();
  }
}

void usage() {
  /* usage: This function prints a usage string. 
   * Inputs: (none)
   * Output: (none)
   */
  cout << "Usage: " << endl;
  //TODO: Update this usage string.
  cout << "./mtKineticModel <Number of Runs> <Output File> [temp/all]" << endl;
}

int numContacts() {
  /* numContacts: This function computes the number of total contacts made
   *   across the cortex. 
   * Inputs: (none)
   * Output: The number of contacts. 
   */
  unsigned int count = 0;
  for (size_t i = 0; i < numberContactWindows; i++)
    count += contacts[i];
  return count;
}

void test() {
  /* test: This is a dummy function useful for testing changes to the code. 
   * Inputs: (none, currently)
   * Output: (none, currently)
   */
}

int main(int argc, const char* argv[]) {
  /* main: This is the main function, run when you call the program as an
   *   executable from the command line. 
   * Inputs: 
   *   int argc: The number of arguments passed to the command line call to this
   *     function. 
   *   const char* argv[]: The arguments, encoded as an array of c-strings. 
   * Output: The return code for this program as a command line executable. 
   */
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
    //We want to run the model a set number of times, but after each run, we
    //must reset the parameters to the starting configuration. 
    runModel(writeAllData,writeTempData);
    setToBasePos();
  }
  //test();
  return 0;
}
