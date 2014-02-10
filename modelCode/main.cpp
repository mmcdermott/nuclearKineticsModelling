#include "initialization.cpp"
#include <iostream>
#include <chrono>
#include <random>
#include <fstream>
#include <assert.h>

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
  file << torque_M <<","<< torque_D <<","<< torque <<","<< MT_Contact_M;
  file <<","<< MT_Contact_D << "," << basePosM << "," << basePosD << endl;
}

void writeFinalData() {
  file << proNucPos << "," << psi << endl;
}

void mtForceCalc(const vec_T (&mtEndPos)[MT_numb], const vec_T &basePos, const
    float_T (&mtContact)[MT_numb], vec_T &force, float_T forceMag = F_MT) 
{
  force[0] = 0;
  force[1] = 0;
  for (size_t i=0; i < MT_numb; ++i) {
    if (mtContact[i] == 0) continue;
    float_T mag = sqrt(pow(mtEndPos[i][0] - basePos[0],2) 
                       + pow(mtEndPos[i][1] - basePos[1],2));
    force[0] += (mtEndPos[i][0] - basePos[0])/mag;
    force[1] += (mtEndPos[i][1] - basePos[1])/mag;
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

  vec_T springAnchorM, springAnchorD, springForceM, springForceD;
  springAnchorM[0] = proNucPos[0] + cosinePrt;
  springAnchorM[1] = proNucPos[1] + sinePrt;
  springAnchorD[0] = proNucPos[0] - cosinePrt;
  springAnchorD[1] = proNucPos[1] - sinePrt;
  //cout << "springAnchorM = [" << springAnchorM << "]" << endl;
  //This computes the Spring Forces (which evenly affect the centrosomes and
  //pronucleus
  springForceM[0] = kM*(springAnchorM[0]-basePosM[0]);
  springForceM[1] = kM*(springAnchorM[1]-basePosM[1]);
  springForceD[0] = kD*(springAnchorD[0]-basePosD[0]);
  springForceD[1] = kD*(springAnchorD[1]-basePosD[1]);
  //cout << "springForceM = [" << springAnchorM << "]" << endl;
  //These compute the force on the centrosomes via the MT, which affects the
  //basePos
  netMTForce('M');
  netMTForce('D');
  // Adding the springForce.
  force_M[0] += springForceM[0];
  force_M[1] += springForceM[1];
  force_D[0] += springForceD[0];
  force_D[1] += springForceD[1];
  //cout << "force_M = [" << force_M << "]" << endl;
  // Updating Basepos: 
  // Base-pos obeys equation \eta_2 dx/dt = F + \xi,
  // where \xi is a random noise parameter on the order of \sqrt{2D \tau}. To
  // generate \xi, take a random number P from norm(0,1) and take 
  // \sqrt{2DP\tau}. We don't know what \eta_2 or D is, those are
  // parameters to be set. As an estimate, try setting \nu_2 to one 10th the
  // pronucleus drag coefficient. 
  float_T randNumXM = stdNormalDist(generator);
  float_T randNumYM = stdNormalDist(generator);
  //cout << "randNumXM = " << randNumXM << endl;
  //cout << "randNumYM = " << randNumYM << endl;
  //TODO: Check if this is right. 
  float_T xiXM      = randNumXM*sqrt(2*D*Tau);
  float_T xiYM      = randNumYM*sqrt(2*D*Tau);
  //cout << "xiXM = " << xiXM << endl;
  //cout << "xiYM = " << xiYM << endl;

  float_T randNumXD = stdNormalDist(generator);
  float_T randNumYD = stdNormalDist(generator);
  float_T xiXD      = randNumXD*sqrt(2*D*Tau);
  float_T xiYD      = randNumYD*sqrt(2*D*Tau);

  force_M[0] += xiXM;
  force_M[1] += xiYM;
  force_D[0] += xiXD;
  force_D[1] += xiYD;
  //cout << "force_M = [" << force_M << "]" << endl;
  //exit(1);
  
  //Now to compute the force on the ProNucleus. 
  // This is old stuff, when we didn't have springs. 
  //  force[0] = force_M[0] + force_D[0];
  //  force[1] = force_M[1] + force_D[1];
  force[0] = -(springForceM[0] + springForceD[0]);
  force[1] = -(springForceM[1] + springForceD[1]);
  //Note that we can see that the angle between positive torque and the true x
  //axis is \psi + \pi/2. Further, the angle between positive torque and the 
  //true y is \psi. The daughter entails a \pi switch, which merely negates
  //things. Thus, 
  torque_M          = +springForceM[0]*sinePrt - springForceM[1]*cosinePrt;
  torque_D          = -springForceD[0]*sinePrt + springForceD[1]*cosinePrt;
  torque            = torque_M + torque_D;
  //Calculating Displacements:
  if (translation) {
    proNucPos[0] += force[0]*(1.0/Eta)*Tau;
    proNucPos[1] += force[1]*(1.0/Eta)*Tau;
    basePosM[0]  += force_M[0]*(1.0/Eta2)*Tau;
    basePosM[1]  += force_M[1]*(1.0/Eta2)*Tau;
    basePosD[0]  += force_D[0]*(1.0/Eta2)*Tau;
    basePosD[1]  += force_D[1]*(1.0/Eta2)*Tau;
  }
  psi += (1/Mu)*torque*Tau;
}

void advanceMT(const float_T vel, vec_T& vec, const float_T mag) {
  vec[0] *= (1 + vel*Tau/mag);
  vec[1] *= (1 + vel*Tau/mag);
}

void respawnMTB(vec_T& vec, const float_T ang) {
  float_T randR = testStat();  
  float_T randT = testStat();
  float_T r     = sqrt(randR);
  float_T t     = pi*randT + ang - pi/2.0;
  vec[0]        = r*cos(t);
  vec[1]        = r*sin(t);
}

float_T probContact(float_T ang, const char centrosome) {
  for (size_t i = 1; i <= numRegions; ++i)
    if (ang < regionAngles[i] && ang >= regionAngles[i-1])
      if (regionProbabilities[i-1] != 1)
        return regionProbabilities[i-1];
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
      } else {
        MT_GrowthVel_D[i] = -Vs_c;
        MT_Growing_D[i]   = false;
      }
      break;
  }
}

void respawnMT(const char centrosome, vec_T& vec, const unsigned i) {
  //assert(centrosome == 'M' || centrosome == 'D')
  switch (centrosome) {
    case 'M':
      respawnMTB(vec, psi);
      MT_Growing_M[i]   = true;
      MT_GrowthVel_M[i] = Vg;
      MT_Contact_M[i]   = 0;
      break;
    case 'D':
      respawnMTB(vec, psi + pi);
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

void runModel(bool writeAllData) {
  if (writeAllData) {
    writeData(0);
  }
  for (float_t t=0; t <= Duration; t += Tau) {
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

      //Respawning if necessary. 
      if (mag_M < 0.2) 
        respawnMT('M', vecM, i);
      if (mag_D < 0.2) 
        respawnMT('D', vecD, i);

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
        }
      } else  if (magScaledM >= 1) {
        mtContactTest('M',i);
      }
      if (MT_Contact_D[i] > 0) {
        MT_Contact_D[i] -= Tau;
        if (MT_Contact_D[i] <= 0) {
          MT_Contact_D[i] = 0;
          MT_GrowthVel_D[i] = -Vs_c;
        }
      } else  if (magScaledD >= 1) {
        mtContactTest('D',i);
      }
    }
    if (writeAllData) {
      writeData(t+Tau);
    }
  }
  if (!writeAllData) {
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
  if (argc == 1) {
    numRuns = 1;
    writeAllData = true;
  } else if (argc == 2) {
    numRuns = atoi(argv[1]);
  } else if (argc == 3) {
    numRuns = atoi(argv[1]);
    fileName = argv[2];
  } else if (argc == 4) {
    numRuns = atoi(argv[1]);
    fileName = argv[2];
    writeAllData = true;
  } else {
    usage();
    return 1;
  }
  setup();

  for (int run = 0; run < numRuns; ++run) {
    runModel(writeAllData);
    setToBasePos();
  }
  return 0;
}
