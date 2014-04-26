#include "parameters.cpp"
#include <sstream>
#include <chrono>
#include <random>

void writeParams() {
  //file << fileOrder << std::endl;
  //  Non-plotting params: 
  file << F_MT << "," << contact_length << "," << Eta << "," << Mu << ",";
  file << Vg << "," << Vs_c << "," << Vs << "," << kc << "," << kr << ",";
  file << translation << "," << startX << "," << startY << "," << startPsi;
  file << std::endl;

  //  Plotting Params: 
  file << MT_numb << "," << R1_max << "," << R2_max << "," << Prad << ",";
  file << Duration << "," << Tau << std::endl;
  file << numRegions << std::endl;
  file << 0;
  for (int i = 1; i <= numRegions; ++i) {
    file << "," << regionAngles[i];
  }
  file << std::endl;
  file << regionProbabilities[0];
  for (int i = 1; i < numRegions; ++i) {
    file << "," << regionProbabilities[i];
  }
  file << std::endl;
  file << regionForceMultipliers[0];
  for (int i = 1; i < numRegions; ++i) {
    file << "," << regionForceMultipliers[i];
  }
  file << std::endl;
}

void setToBasePos() {
  //Initialize:
  force_M[0] = 0;
  force_M[1] = 0;
  force_D[0] = 0;
  force_D[1] = 0;
  force[0]   = 0;
  force[1]   = 0;
  torque_M   = 0;
  torque_D   = 0;
  torque     = 0;

  //Base System: 
  psi = startPsi;
  proNucPos[0] = startX;
  proNucPos[1] = startY;

  float_T cosinePrt  = Prad * cos(psi);
  float_T sinePrt    = Prad * sin(psi);

  basePosM[0] = proNucPos[0] + cosinePrt;
  basePosM[1] = proNucPos[1] + sinePrt;
  basePosD[0] = proNucPos[0] - cosinePrt;
  basePosD[1] = proNucPos[1] - sinePrt;

  for (unsigned i = 0; i < MT_numb; ++i) {
    //MTs:
    //Constructing the 4 random numbers necessary per MT.
    // Suffix *R means radius, *T theta. 
    // Suffix *M* means mother, *D* means daughter.
    //
    // TODO: Implement Distribution (Boost ?)
    float_T randMR = 100*((float_T) rand())/((float_T) RAND_MAX);  //Random Values
    float_T randMT = ((float_T) rand())/((float_T) RAND_MAX);  //Random Values
    float_T randDR = 100*((float_T) rand())/((float_T) RAND_MAX);  //Random Values
    float_T randDT = ((float_T) rand())/((float_T) RAND_MAX);  //Random Values
    //std::cout << "randMR: " << randMR << std::endl;
    //std::cout << "randMT: " << randMT << std::endl;
    //std::cout << "randDR: " << randDR << std::endl;
    //std::cout << "randDT: " << randDT << std::endl;
    //Using these RVs to construct the position variables of the MTs.
    float_T r_M    = sqrt(randMR); //Random Radius of MT;
    //std::cout << r_M << std::endl;
    float_T t_M    = (startPsi - pi/2) + (envelopeM[1]-envelopeM[0])*randMT + envelopeM[0];    //Random Theta of MT;
    float_T r_D    = sqrt(randDR); //Random Radius of MT;
    float_T t_D    = (startPsi + pi/2) + (envelopeD[1]-envelopeD[0])*(randDT) + envelopeD[0];//Random Theta of MT;
    //Assigning the positions. 
    MT_Pos_M[i][0] = basePosM[0]+ r_M*cos(t_M);
    MT_Pos_M[i][1] = basePosM[1]+ r_M*sin(t_M);
    MT_Pos_D[i][0] = basePosD[0]+ r_D*cos(t_D);
    MT_Pos_D[i][1] = basePosD[1]+ r_D*sin(t_D);

    //Growing or Shrinking?
    MT_Growing_M[i] = true;
    MT_Growing_D[i] = true;
    MT_GrowthVel_M[i] = Vg;
    MT_GrowthVel_D[i] = Vg;

    //Made Contact?
    MT_Contact_M[i] = 0;
    MT_Contact_D[i] = 0;
  }
}

void setup() {
  srand(time(NULL));
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  generator = std::default_random_engine(seed);
  //TODO: Use one random system. 

  //Set up file: 
  
  std::stringstream fileNameSS;
  if (fileName == "") {
    fileNameSS << fileDir << "NO-FILENAME-SET-MT-" << seed << ".csv";
  } else { 
    fileNameSS << fileDir << fileName << ".csv";
  } 

  fileName = fileNameSS.str();
  file.open(fileName);
  writeParams();
  setToBasePos();
  for (size_t i = 0; i < numberContactWindows; i++)
    contacts[i] = false;
}
