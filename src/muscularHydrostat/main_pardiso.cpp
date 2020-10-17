//##################################################################################
//
// rigid body-elastic body interaction
//
// Copyright (c) 2016- Mechanical and Bioengineering Systems Lab.,
//                     Department of Mechanical Science and Bioengineering,
//                     Graduate School of Engineering Science,
//                     Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   main_pardiso.cpp
 * @author T. Otani
 * @detail Fem solid analysis with displacement control.
 */

#include "muscularHydrostat.h"

#if _GLOG
#include "glog/logging.h"
#endif

using namespace std;

int main(int argc,char *argv[])
{
  #if _GLOG
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
  #endif

  muscularHydrostat::Muscle muscle;

  if(argc!=2){
    cout << "Invalid input" << endl;
    return -1;
  }

  //read tp file
  std::string input_file = argv[1];
  int ierror;
  if ((ierror = muscle.tp.read(input_file)) != TP_NO_ERROR) {
    printf("\tError at reading '%s' file\n", input_file.c_str());
    return 1;
  }

  cout << "---------preprocess start----------" << endl;
  muscle.preprocess();
  cout << "---------preprocess completed----------" << endl << endl;
  cout << "---------main loop start----------" << endl;
  muscle.femSolidAnalysis();

  return 0;
}
