#ifndef _INFO_H_
#define _INFO_H_

/**
 * @file info.h
 * @brief  solverInfo Class Header
 * @brief  outputInfo Class Headerer.h"

 * @author Tomohiro Otani
 */

#include <string>
#include "TextPars"
class outputInfo{
 public:
  outputInfo(){};
  ~outputInfo(){};
  std::string outputDir;
  std::string format;
  std::string filepath;
  void getOutputInfo(TextParser &tp);
};

class solverInfo{
 public:
  solverInfo(){};
  ~solverInfo(){};

  double dt;
  int maxIteration;
  int ompNumThreads;
  std::string LinearSolverlibrary;
  std::string solver;
  int precondition;
  int LSMaxIteration;
  double tolerance;
  std::string initx;
  std::string scale;
  void getSolverInfo(TextParser &tp);

 private:
};

#endif //_INFO_H_