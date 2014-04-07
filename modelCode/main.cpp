#include "initialization.cpp"
#include <iostream>
#include <chrono>
#include <random>
#include <fstream>
#include <assert.h>
#include <cstring>

using namespace std;

float_T testStat() {
  return ((float_T) rand())/((float_T) RAND_MAX);
}

ostream& operator<<(ostream& out, const vec_T& rhs) {
  if (ONLY_COMMA) {
    out << rhs[0] << "," << rhs[1];
  } else {
    out << rhs[0] << "|" << rhs[1];
  }
  return out;
}

ostream& operator<<(ostream& out, const vec_T (&rhs)[MT_numb]) {
  out << rhs[0];
  for (size_t i = 1; i < MT_numb; ++i) {
    if (ONLY_COMMA) {
      out << "," << rhs[i];
    } else {
      out << ";" << rhs[i];
    }
  }
  return out;
}

ostream& operator<<(ostream& out, const float_T (&rhs)[MT_numb]) {
  out << rhs[0];
  for (size_t i = 1; i < MT_numb; ++i) {
    if (ONLY_COMMA) {
      out << "," << rhs[i];
    } else {
      out << ";" << rhs[i];
    }
  }
  return out;
}

void writeData(const float_T t) {
  file << t <<","<< proNucPos <<","<< psi <<","<<  MT_Pos_M <<",";
  file << MT_Pos_D <<","<< force_M <<","<< force_D <<","<< force <<",";
  file << torque_M <<","<< torque_D <<","<< torque <<","<< basePosM <<",";
  file << basePosD << endl;
}

void writeFinalData() {
  file << proNucPos << "," << psi << endl;
}

void writePartialData(const float_T t) {
  file<<t<<","<<proNucPos <<","<< psi <<","<< basePosM <<","<< basePosD << endl;
}

void updateRegionProbabilitiesD(const float_T angle) {
  //cout << "Calling updateRegionProbabilitiesD! effRegionAngles.size = ";
  //cout << effRegionAngles.size() << endl;
  //An MT Detached, now we need to remove some regions. 
  float_T low = angle - contactWindow;
  float_T high = angle + contactWindow;
  bool metHigh = false;
  if (low < 0) {
    low += 2*pi;
    high += 2*pi;
  }
  if (high > 2*pi) {
    //DO something different TODO;
  }
  vector<float_T>::iterator itA = effRegionAngles.begin()++;
  vector<float_T>::iterator itP = effRegionProbabilities.begin();
  for (; itA != effRegionAngles.end(); itA++, itP++) {
    if (low == *itA) {
      //low >= *(it-1), as otherwise you would have stopped there
      itA = effRegionAngles.erase(itA);
      break;
    }
  }
  for (; itA != effRegionAngles.end(); itA++,itP++) {
    if (high == *itA) {
      itP = effRegionProbabilities.erase(itP);
      itA = effRegionAngles.erase(itA);
      metHigh = true;
      break;
    }
    (*itP) =(1/probabilityFactor) * (*itP);
  }
  if (!metHigh) {
    itA = effRegionAngles.begin()++;
    itP = effRegionProbabilities.begin();
    for (; itA != effRegionAngles.end(); itA++,itP++) {
      if (high == *itA) {
        itP = effRegionProbabilities.erase(itP);
        itA = effRegionAngles.erase(itA);
        metHigh = true;
        break;
      }
      (*itP) =(1/probabilityFactor) * (*itP);
    }
  }
}

void insertAndScale(float_T low, float_T high) {
  //cout << "Calling InsertAndScale! effRegionAngles.size = ";
  //cout << effRegionAngles.size() << endl;
  if (low < 0) {
    low += 2*pi; 
    high += 2*pi;
  } 
  bool metHigh = false;
  vector<float_T>::iterator itA = effRegionAngles.begin()++;
  vector<float_T>::iterator itP = effRegionProbabilities.begin();
  for (; itA != effRegionAngles.end(); itA++, itP++) {
    if (low < *itA) {
      //low >= *(it-1), as otherwise you would have stopped there
      itA = effRegionAngles.insert(itA, low);
      break;
    }
  }
  itA++;
  for (; itA != effRegionAngles.end(); itA++,itP++) {
    if (high < *itA) {
      float_T oldProb = *itP;
      float_T newProb = probabilityFactor * (*itP);

      itP = effRegionProbabilities.insert(itP, newProb);
      itA = effRegionAngles.insert(itA,high);
      itP++;
      itP = effRegionProbabilities.insert(itP, oldProb);

      metHigh = true;
      break;
    }
    (*itP) = probabilityFactor * (*itP);
  }
  if (!metHigh) {
    itA = effRegionAngles.begin()++;
    itP = effRegionProbabilities.begin();
    for (; itA != effRegionAngles.end(); itA++,itP++) {
      if (high < *itA) {
        float_T oldProb = *itP;
        float_T newProb = probabilityFactor * (*itP);

        itP = effRegionProbabilities.insert(itP, newProb);
        itA = effRegionAngles.insert(itA,high);
        itP++;
        itP = effRegionProbabilities.insert(itP, oldProb);

        metHigh = true;
        break;
      }
      (*itP) = probabilityFactor * (*itP);
    }
  }
}

void updateRegionProbabilities(const float_T angle) {
  //An MT Attached, now we need to remove some regions. 
  insertAndScale(angle - contactWindow, angle + contactWindow);
  //float_T low = angle - contactWindow;
  //float_T high = angle + contactWindow;
  //if (low < 0) {
  //  low += 2*pi; 
  //  high += 2*pi;
  //} 
  //if (high > 2*pi) {
  //  //find and insert low
  //  for (; itA != effRegionAngles.end(); itA++, itP++) {
  //    if (low < *itA) {
  //      //low >= *(it-1), as otherwise you would have stopped there
  //      effRegionAngles.insert(itA, low);
  //      break;
  //    }
  //  }
  //  itA++;
  //  for (; itA != effRegionAngles.end(); itA++,itP++) {
  //    (*itP) = probabilityFactor * (*itP);
  //  } 
  //  //scale up through high
  //  for (; itA != effRegionAngles.end(); itA++,itP++) {
  //    if (high < *itA) {
  //      float_T oldProb = *itP;
  //      float_T newProb = probabilityFactor * (*itP);
  //  
  //      effRegionProbabilities.insert(itP, newProb);
  //      effRegionAngles.insert(itA,high);
  //      itP++;
  //      effRegionProbabilities.insert(itP, oldProb);
  //      break;
  //    }
  //    (*itP) = probabilityFactor * (*itP);
  //  }
  //}
  //vector<float_T>::iterator itA = effRegionAngles.begin()++;
  //vector<float_T>::iterator itP = effRegionProbabilities.begin();
  //for (; itA != effRegionAngles.end(); itA++, itP++) {
  //  if (low < *itA) {
  //    //low >= *(it-1), as otherwise you would have stopped there
  //    effRegionAngles.insert(itA, low);
  //    break;
  //  }
  //}
  //itA++;
  //for (; itA != effRegionAngles.end(); itA++,itP++) {
  //  if (high < *itA) {
  //    float_T oldProb = *itP;
  //    float_T newProb = probabilityFactor * (*itP);

  //    effRegionProbabilities.insert(itP, newProb);
  //    effRegionAngles.insert(itA,high);
  //    itP++;
  //    effRegionProbabilities.insert(itP, oldProb);
  //    break;
  //  }
  //  (*itP) = probabilityFactor * (*itP);
  //}
}

void mtForceCalc(const vec_T (&mtEndPos)[MT_numb], const vec_T &basePos, const
    float_T (&mtContact)[MT_numb], vec_T &force, float_T forceMag = F_MT) 
{
  force[0] = 0;
  force[1] = 0;
  for (size_t i=0; i < MT_numb; ++i) {
    if (mtContact[i] == 0) continue;
    float_T angle = atan2(mtEndPos[i][1],mtEndPos[i][0]);
    if (angle < 0) angle += 2*pi;
    float_T multiplier = 1;
    for (size_t i = 1; i <= numRegions; ++i) {
      if (angle < regionAngles[i] && angle >= regionAngles[i-1]) {
        //cout << "Region " << i << " with bounds [" << regionAngles[i-1];
        //cout << ", " << regionAngles[i] << ") yields multiplier ";
        multiplier = regionForceMultipliers[i-1];
        //cout << multiplier << endl;
      }
    }

    float_T mag = sqrt(pow(mtEndPos[i][0] - basePos[0],2) 
                       + pow(mtEndPos[i][1] - basePos[1],2));
    force[0] += multiplier*(mtEndPos[i][0] - basePos[0])/mag;
    force[1] += multiplier*(mtEndPos[i][1] - basePos[1])/mag;
  }
  force[0] *= forceMag;
  force[1] *= forceMag;
  //TODO: Include spring force F(x) = -kx
}

void netMTForce(char centrosome) {
  //assert(centrosome == 'M' || centrosome == 'D')
  switch (centrosome) {
    case 'M':
      mtForceCalc(MT_Pos_M, basePosM, MT_Contact_M, force_M, F_MT);
      break;
    case 'D':
      mtForceCalc(MT_Pos_D, basePosD, MT_Contact_D, force_D, F_MT);
      break;
  }
}

void updatePNPos() {
  //First, computing the springAnchor Points. These get used a lot, so we'll
  //start with them. 
  float_T cosinePrt  = Prad*cos(psi);
  float_T sinePrt    = Prad*sin(psi);

  if (springsOn) {
    //These compute the force on the centrosomes via the MT, which affects the
    //basePos
    netMTForce('M');
    netMTForce('D');

    vec_T springAnchorM, springAnchorD, springForceM, springForceD;
    springAnchorM[0] = proNucPos[0] + cosinePrt;
    springAnchorM[1] = proNucPos[1] + sinePrt;
    springAnchorD[0] = proNucPos[0] - cosinePrt;
    springAnchorD[1] = proNucPos[1] - sinePrt;
    //This computes the Spring Forces (which evenly affect the centrosomes and
    //pronucleus
    springForceM[0] = kM*(springAnchorM[0]-basePosM[0]);
    springForceM[1] = kM*(springAnchorM[1]-basePosM[1]);
    springForceD[0] = kD*(springAnchorD[0]-basePosD[0]);
    springForceD[1] = kD*(springAnchorD[1]-basePosD[1]);
    // Adding the springForce.
    force_M[0] += springForceM[0];
    force_M[1] += springForceM[1];
    force_D[0] += springForceD[0];
    force_D[1] += springForceD[1];
    //Now to compute the force on the ProNucleus. 
    force[0] = -(springForceM[0] + springForceD[0]);
    force[1] = -(springForceM[1] + springForceD[1]);
    //Note that we can see that the angle between positive torque and the true x
    //axis is \psi + \pi/2. Further, the angle between positive torque and the 
    //true y is \psi. The daughter entails a \pi switch, which merely negates
    //things. Thus, 
    torque_M          =  springForceM[0]*sinePrt - springForceM[1]*cosinePrt;
    torque_D          = -springForceD[0]*sinePrt + springForceD[1]*cosinePrt;
    torque            = torque_M + torque_D;

    // Updating Basepos: 
    // Base-pos obeys equation \eta_2 dx/dt = F + \xi,
    // where \xi is a random noise parameter on the order of \sqrt{2D \tau}. To
    // generate \xi, take a random number P from norm(0,1) and take 
    // P*\sqrt{2D\tau}. We don't know what \eta_2 is, that is a parameter to be
    // set, but D = kB*T/\eta_2. As an estimate, try setting \nu_2 to one 10th
    // the pronucleus drag coefficient. 
    float_T randNumXM = stdNormalDist(generator);
    float_T randNumYM = stdNormalDist(generator);
    float_T xiXM      = randNumXM*sqrt(2*D*Tau);
    float_T xiYM      = randNumYM*sqrt(2*D*Tau);

    float_T randNumXD = stdNormalDist(generator);
    float_T randNumYD = stdNormalDist(generator);
    float_T xiXD      = randNumXD*sqrt(2*D*Tau);
    float_T xiYD      = randNumYD*sqrt(2*D*Tau);

    basePosM[0]  += force_M[0]*(1.0/Eta2)*Tau + xiXM;
    basePosM[1]  += force_M[1]*(1.0/Eta2)*Tau + xiYM;
    basePosD[0]  += force_D[0]*(1.0/Eta2)*Tau + xiXD;
    basePosD[1]  += force_D[1]*(1.0/Eta2)*Tau + xiYD;
  } else {
    basePosM[0] = proNucPos[0] + cosinePrt;
    basePosM[1] = proNucPos[1] + sinePrt;
    basePosD[0] = proNucPos[0] - cosinePrt;
    basePosD[1] = proNucPos[1] - sinePrt;

    //These compute the force on the centrosomes via the MTs
    netMTForce('M');
    netMTForce('D');
    if (spitValues) {
      cout << "Y-force on M: " << force_M[1] << endl;
      cout << "Y-force on D: " << force_D[1] << endl;
    }

    force[0] = force_M[0] + force_D[0];
    force[1] = force_M[1] + force_D[1];

    torque_M          = -force_M[0]*sinePrt + force_M[1]*cosinePrt;
    torque_D          = force_D[0]*sinePrt - force_D[1]*cosinePrt;
    torque            = torque_M + torque_D;
  }

  //Calculating Displacements:
  if (translation) {
    proNucPos[0] += force[0]*(1.0/Eta)*Tau;
    proNucPos[1] += force[1]*(1.0/Eta)*Tau;
  }
  psi += (1/Mu)*torque*Tau;
}

void advanceMT(const float_T vel, vec_T& vec, const float_T mag) {
  vec[0] *= (1 + vel*Tau/mag);
  vec[1] *= (1 + vel*Tau/mag);
}

float_T probContact(float_T ang, const char centrosome) {
  //for (size_t i = 1; i <= numRegions; ++i)
  for (size_t i = 1; i < effRegionAngles.size(); i++) 
    if (ang < effRegionAngles[i] && ang >= effRegionAngles[i-1])
      return regionProbabilities[i-1];
  //    if (regionProbabilities[i-1] != 1)
  //      return regionProbabilities[i-1];
  return 1;
}

void mtContactTest(const char centrosome, const unsigned i) {
  float_T angleM, angleD;
  switch (centrosome) {
    case 'M':
      angleM = atan2(MT_Pos_M[i][1],MT_Pos_M[i][0]);
      if (angleM < 0) angleM += 2*pi; //atan2 returns negative vals.
      if (testStat() < probContact(angleM,'M')) {
        MT_GrowthVel_M[i] = 0;
        MT_Contact_M[i]   = contact_length;
        MT_Growing_M[i]   = false;
        updateRegionProbabilities(angleM);
      } else {
        MT_GrowthVel_M[i] = -Vs_c;
        MT_Growing_M[i]   = false;
      }
      break;
    case 'D':
      angleD = atan2(MT_Pos_D[i][1],MT_Pos_D[i][0]);
      if (angleD < 0) angleD += 2*pi; //atan2 returns negative vals.
      if (testStat() < probContact(angleD,'D')) {
        MT_GrowthVel_D[i] = 0;
        MT_Contact_D[i]   = contact_length;
        MT_Growing_D[i]   = false;
        updateRegionProbabilities(angleD);
      } else {
        MT_GrowthVel_D[i] = -Vs_c;
        MT_Growing_D[i]   = false;
      }
      break;
  }
}

void respawnMTB(vec_T& vec, const float_T ang, const vec_T& envelope) {
  float_T randR = testStat();  
  float_T randT = testStat();
  float_T r     = sqrt(randR);
  float_T t     = (envelope[1]-envelope[0])*randT + ang + envelope[0] - pi/2.0;
  vec[0]        = r*cos(t);
  vec[1]        = r*sin(t);
}

void respawnMT(const char centrosome, vec_T& vec, const unsigned i, const vec_T& envelope) {
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
    float_t dist = sqrt(pow(proNucPos[0]/(R1_max - Prad - .5),2) +
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
  for (float_t t=0; t <= Duration; t += Tau) {
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

      if (mag_M < 0.2 || ((thetaM < envelopeM[0]) || (thetaM > envelopeM[1])))
        respawnMT('M', vecM, i, envelopeM);
      if (mag_D < 0.2 || ((thetaD < envelopeD[0]) || (thetaD > envelopeD[1])))
        respawnMT('D', vecD, i, envelopeD);

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
        if (test < Pr_catastrophe && mag_M > 0.4) {
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
      float_T magScaledM = sqrt(pow(MT_Pos_M[i][0]/R1_max, 2) + pow(MT_Pos_M[i][1]/R2_max, 2));
      float_T magScaledD = sqrt(pow(MT_Pos_D[i][0]/R1_max, 2) + pow(MT_Pos_D[i][1]/R2_max, 2));

      //Updating Contact Indicator:
      if (MT_Contact_M[i] > 0) {
        MT_Contact_M[i] -= Tau;
        if (MT_Contact_M[i] <= 0) {
          MT_Contact_M[i] = 0;
          MT_GrowthVel_M[i] = -Vs_c;
          float_T angleM = atan2(MT_Pos_M[i][1],MT_Pos_M[i][0]);
          updateRegionProbabilitiesD(angleM);
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
          updateRegionProbabilitiesD(angleD);
        }
      } else  if (magScaledD >= 1) {
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
  return 0;
}
