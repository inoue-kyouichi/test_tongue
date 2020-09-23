/**
 * @file minimumComplianceProblem.cpp
 * @brief minimumComplianceProblem class
 * @author T. Otani
 */

#include "minimumComplianceProblem.h"

using namespace std;

// #################################################################
/**
 * @brief rigid body-elastic body interaction problem
 */
void MinimumComplianceProblem::minimumComplianceProblem::mainLoop()
{
  int output_iter=1;
  string output;
  FILE *fp;

  for (int loop = 1; loop <= 1000; loop++)
  {
    mainProblem.mainLoop();

    //self-adjoint problem
    #pragma omp parallel for
    for(int i=0;i<mainProblem.elasticBody.numOfNode;i++){
      for (int j = 0; j < 3; j++) adjointV(i, j) = mainProblem.elasticBody.U(i, j);
    }

    // f0 = costFunction();
    // calcShapeGradient();
    // H1GradientMethod(PARDISO_for_H1GradientMethod, dt);

    // updateSurface();

    // fp = fopen("history.dat", "a+");
    // fprintf(fp, "%d %e\n", loop, f0);
    // fclose(fp);

    if (loop % 10 == 0)
    {
        string fileName = mainProblem.outputDir + "/test_shapeGradient_" + to_string(loop / 10) + ".vtu";
        // exportVTU_boundaryShapeGradient(fileName);
        // exportVTU_ShapeGradient(fileName);
    }
  }

}

// #################################################################
/**
 * @brief initialize and preprocess
 */
void MinimumComplianceProblem::minimumComplianceProblem::preprocess()
{
  mainProblem.preprocess();

  adjointV.allocate(mainProblem.elasticBody.numOfNode,3);
}

