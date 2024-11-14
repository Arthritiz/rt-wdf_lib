/*
 ==============================================================================

 This file is part of the RT-WDF library.
 Copyright (c) 2015,2016 - Maximilian Rest, Ross Dunkel, Kurt Werner.

 Permission is granted to use this software under the terms of either:
 a) the GPL v2 (or any later version)
 b) the Affero GPL v3

 Details of these licenses can be found at: www.gnu.org/licenses

 RT-WDF is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 -----------------------------------------------------------------------------
 To release a closed-source product which uses RT-WDF, commercial licenses are
 available: write to rt-wdf@e-rm.de for more information.

 ==============================================================================

 rt-wdf_nlModels.cpp
 Created: 2 Dec 2015 4:10:47pm
 Author:  mrest

 ==============================================================================
*/

#include "rt-wdf_nlModels.h"


//==============================================================================
// Parent class for nlModels
//==============================================================================
nlModel::nlModel( int numPorts ) : numPorts (numPorts) {

}

nlModel::~nlModel( ) {

}

//----------------------------------------------------------------------
int nlModel::getNumPorts( ) {
    return numPorts;
}


//==============================================================================
// Diode Models according to Kurt Werner et al
// ("An Improved and Generalized Diode Clipper Model for Wave Digital Filters")
//==============================================================================
#define Is_DIODE    2.52e-9
#define VT_DIODE    0.02585

diodeModel::diodeModel() : nlModel( 1 ) {

}

//----------------------------------------------------------------------
void diodeModel::calculate( Wvec* fNL,
                            Wmat* JNL,
                            Wvec* x,
                            Wvec* lastX,
                            int* currentPort ) {

    const double vd = (*x)(*currentPort);
    const double arg1 = vd/VT_DIODE;

    (*fNL)(*currentPort) = Is_DIODE*(exp(arg1)-1);
    (*JNL)(*currentPort,*currentPort) = (Is_DIODE/VT_DIODE)*exp(arg1);

    (*currentPort) = (*currentPort)+getNumPorts();
}

//==============================================================================
diodeApModel::diodeApModel( ) : nlModel( 1 ) {

}

//----------------------------------------------------------------------
void diodeApModel::calculate( Wvec* fNL,
                              Wmat* JNL,
                              Wvec* x,
                              Wvec* lastX,
                              int* currentPort) {

    const double vd = (*x)(*currentPort);
    const double arg1 = vd/VT_DIODE;

    (*fNL)(*currentPort) = Is_DIODE*(exp(arg1)-1)-Is_DIODE*(exp(-arg1)-1);
    (*JNL)(*currentPort,*currentPort) = (Is_DIODE/VT_DIODE)*(exp(arg1)+exp(-arg1));

    (*currentPort) = (*currentPort)+getNumPorts();
}


//==============================================================================
// Transistor Models using Ebers-Moll equations
// ("Large-signal behavior of junction transistors")
//==============================================================================
#define VT_BJT      0.02585

npnEmModel::npnEmModel() : nlModel( 2 ) {

}

//----------------------------------------------------------------------
void npnEmModel::calculate( Wvec* fNL,
                            Wmat* JNL,
                            Wvec* x,
                            Wvec* lastX,
                            int* currentPort) {

    const double Is_BJT = 5.911e-15;
    const double BETAF  = 1.434e3;
    const double BETAR  = 1.262;
    const double ALPHAF = (BETAF/(1.0+BETAF));
    const double ALPHAR = (BETAR/(1.0+BETAR));

    const double vBC = (*x)(*currentPort);
    const double vBE = (*x)((*currentPort)+1);

    const double vBC_o_VT_BJT = vBC/VT_BJT;
    const double vBE_o_VT_BJT = vBE/VT_BJT;
    const double Is_BJT_o_VT_BJT = Is_BJT/VT_BJT;
    const double Is_BJT_o_ALPHAR = Is_BJT/ALPHAR;
    const double Is_BJT_o_ALPHAF = Is_BJT/ALPHAF;

    // i_bc
    (*fNL)(*currentPort) = -Is_BJT*(exp(vBE_o_VT_BJT )-1)+(Is_BJT_o_ALPHAR)*(exp(vBC_o_VT_BJT)-1);

    // dI_bc/dVbc
    (*JNL)((*currentPort),(*currentPort)) = (Is_BJT_o_ALPHAR/VT_BJT)*exp(vBC_o_VT_BJT);
    // dI_bc/dVbe
    (*JNL)((*currentPort),((*currentPort)+1)) = (-Is_BJT_o_VT_BJT)*exp(vBE_o_VT_BJT );

    // i_be
    (*fNL)((*currentPort)+1) = (Is_BJT_o_ALPHAF)*(exp(vBE_o_VT_BJT )-1)-Is_BJT*(exp(vBC_o_VT_BJT)-1);
    // dI_be/dVbc
    (*JNL)(((*currentPort)+1),(*currentPort)) = (-Is_BJT_o_VT_BJT)*exp(vBC_o_VT_BJT);
    // dI_be/dVbe
    (*JNL)(((*currentPort)+1),((*currentPort)+1)) = (Is_BJT_o_ALPHAF/VT_BJT)*exp(vBE_o_VT_BJT );

    (*currentPort) = (*currentPort)+getNumPorts();
}

const std::map<std::string, pnpEmModel::ModelSpec> pnpEmModel::modelSpecs = {
    {"default", ModelSpec(5.911e-15, 1.434e3, 1.262)},
    {"AC128-A", ModelSpec(85.8e-9, 85, 20)},
    {"AC128-B", ModelSpec(120.8e-9, 120, 20)},
};

double pnpEmModel::limitStep(double vnew, double vold) {
    double arg;

    if ((vnew > vcrit) && (abs(vnew - vold) > (VT_BJT + VT_BJT))) {
        if (vold > 0) {
            arg = 1 + (vnew - vold)/VT_BJT;
            if (arg > 0) {
                vnew = vold + VT_BJT * std::log(arg);
            } else {
                vnew = vcrit;
            }
        } else {
            vnew = VT_BJT * std::log(vnew/VT_BJT);
        }
    }

    return vnew;
}

pnpEmModel::pnpEmModel(std::string modelName): nlModel(2)
{
    auto iter = pnpEmModel::modelSpecs.find(modelName);
    if (iter == pnpEmModel::modelSpecs.end())
    {
        throw;
    }

    this->modelSpec = iter->second;

    vcrit = VT_BJT * std::log(VT_BJT/(std::sqrt(2)*this->modelSpec.Is_BJT));
}

void pnpEmModel::calculate( Wvec* fNL,
                            Wmat* JNL,
                            Wvec* x,
                            Wvec* lastX,
                            int* currentPort)
{
    for (int i = *currentPort; i < *currentPort + getNumPorts(); i++)
    {
        (*x)(i) = limitStep((*x)(i), (*lastX)(i));
    }

    const double vEB = (*x)(*currentPort);
    const double vCB = (*x)((*currentPort)+1);

    const double vEB_o_VT_BJT = vEB/VT_BJT;
    const double vCB_o_VT_BJT = vCB/VT_BJT;
    const double Is_BJT_o_VT_BJT = modelSpec.Is_BJT/VT_BJT;
    const double Is_BJT_o_ALPHAR = modelSpec.Is_BJT/modelSpec.ALPHAR;
    const double Is_BJT_o_ALPHAF = modelSpec.Is_BJT/modelSpec.ALPHAF;

    // i_eb
    (*fNL)(*currentPort) = (Is_BJT_o_ALPHAF)*(exp(vEB_o_VT_BJT)-1) - modelSpec.Is_BJT*(exp(vCB_o_VT_BJT)-1);
    // dI_eb/dVeb
    (*JNL)((*currentPort),(*currentPort)) = (Is_BJT_o_ALPHAF/VT_BJT)*exp(vEB_o_VT_BJT);
    // dI_eb/dVcb
    (*JNL)((*currentPort),((*currentPort)+1)) = (-Is_BJT_o_VT_BJT)*exp(vCB_o_VT_BJT);

    // i_cb (-i_c)
    (*fNL)((*currentPort)+1) = -modelSpec.Is_BJT*(exp(vEB_o_VT_BJT)-1) + (Is_BJT_o_ALPHAR)*(exp(vCB_o_VT_BJT)-1);
    // dI_cb/dVeb
    (*JNL)(((*currentPort)+1),(*currentPort)) = (-Is_BJT_o_VT_BJT)*exp(vEB_o_VT_BJT);
    // dI_cb/dVcb
    (*JNL)(((*currentPort)+1),((*currentPort)+1)) = (Is_BJT_o_ALPHAR/VT_BJT)*exp(vCB_o_VT_BJT);

    (*currentPort) = (*currentPort)+getNumPorts();
}

//==============================================================================
// Triode model according to Dempwolf et al
// ("A physically-motivated triode model for circuit simulations")
//==============================================================================
triDwModel::triDwModel() : nlModel( 2 ) {


}

//----------------------------------------------------------------------
void triDwModel::calculate( Wvec* fNL,
                            Wmat* JNL,
                            Wvec* x,
                            Wvec* lastX,
                            int* currentPort) {

    const double G = 2.242E-3;
    const double C = 3.40;
    const double mu = 103.2;
    const double y = 1.26;

    const double Gg = 6.177E-4;
    const double Cg = 9.901;
    const double E = 1.314;
    const double Ig0 = 8.025E-8;

    const double vAC_mu = (*x)(*currentPort) / mu;
    const double vGC = (*x)((*currentPort)+1);


    const double exp_Cg_vGC = exp( Cg * vGC );
    const double log_1_exp_Cg_vGC_Cg = (log( 1 + exp_Cg_vGC ) / Cg);


    // Ig
    (*fNL)((*currentPort)+1) = Gg * pow( log_1_exp_Cg_vGC_Cg, E ) + Ig0;

    // dIg / dvAC
    (*JNL)(((*currentPort)+1),(*currentPort)) = 0;
    // dIg / dvGC
    (*JNL)(((*currentPort)+1),((*currentPort)+1)) = ( Gg * E * exp_Cg_vGC *
                                                      pow( log_1_exp_Cg_vGC_Cg, (E-1)) ) /
                                                    (1 + exp_Cg_vGC);


    const double exp_C_vAC_mu_vGC = exp( C * ( vAC_mu + vGC ));
    const double pow_log_1_exp_C_vAC_mu_vGC_C_y_1 = pow( (log(1 + exp_C_vAC_mu_vGC) / C), (y-1));

    // Ik
    (*fNL)(*currentPort) = G * pow( log( 1 + exp_C_vAC_mu_vGC ) / C , y ) - (*fNL)((*currentPort)+1);

    // dIk / dvAC
    (*JNL)((*currentPort),(*currentPort)) = ( G * y * exp_C_vAC_mu_vGC *
                                              pow_log_1_exp_C_vAC_mu_vGC_C_y_1 ) /
                                            (mu * (1 + exp_C_vAC_mu_vGC));
    // dIk / dvGC
    (*JNL)((*currentPort),((*currentPort)+1)) = ( G * y * exp_C_vAC_mu_vGC *
                                                  pow_log_1_exp_C_vAC_mu_vGC_C_y_1 ) /
                                                (1 + exp_C_vAC_mu_vGC) - (*JNL)(((*currentPort)+1),((*currentPort)+1));


    (*currentPort) = (*currentPort)+getNumPorts();

}
