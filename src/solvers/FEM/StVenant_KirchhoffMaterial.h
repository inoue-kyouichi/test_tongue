#ifndef _StVENANT_KIRCHHOFF_
#define _StVENANT_KIRCHHOFF_

#include "fem.h"


class StVenantKirchhoffMaterial : public Fem{
 public:
    StVenantKirchhoffMaterial(){}
    StVenantKirchhoffMaterial(const double young,const double poisson){
        Young = young;
        Poisson = poisson;
        lambda = young * poisson / ((1e0+poisson) * (1e0-2e0*poisson));
        mu = 5e-1 * young / (1e0+poisson);
    }
    ~StVenantKirchhoffMaterial(){}
    double Young,Poisson;
    void calcStressTensor();
    void inputMaterialParameters(TextParser &tp);
 private:
    double lambda,mu;
    void calcStressTensor_SantVenant_element_spatialForm(const int ic,DOUBLEARRAY2D &U_tmp,const bool option);
    double SantVenant_inGaussIntegral(DOUBLEARRAY2D &dNdr,DOUBLEARRAY2D &x_current,DOUBLEARRAY2D &x_ref,DOUBLEARRAY2D &dNdx,const int numOfNodeInElm,
    const double weight,const int ic,const bool option);
};

#endif