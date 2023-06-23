/**
 *  @file PengRobinsonGasPhase.cpp
 *   ThermoPhase object for the Peng-Robinson gas equation of
 * state - workhorse for %Cantera (see \ref thermoprops
 * and class \link Cantera::PengRobinsonGasPhase PengRobinsonGasPhase\endlink).
 */

#include "cantera/thermo/PengRobinsonGasPhase.h"
#include "cantera/base/utilities.h"

#include <iostream>

using namespace std;

namespace Cantera
{
PengRobinsonGasPhase::PengRobinsonGasPhase() : m_p0(-1.0) {}

PengRobinsonGasPhase::PengRobinsonGasPhase(const std::string& inputFile,
                                           const std::string& id_)
    : m_p0(-1.0)
{
    initThermoFile(inputFile, id_);
}

PengRobinsonGasPhase::PengRobinsonGasPhase(XML_Node& phaseRef,
                                           const std::string& id_)
    : m_p0(-1.0)
{
    initThermoXML(phaseRef, id_);
}

// Molar Thermodynamic Properties of the Solution ------------------

doublereal PengRobinsonGasPhase::enthalpy_mole() const
{
    _updateThermoRealFluid();
    double h0 = GasConstant * temperature() * mean_X(enthalpy_RT_ref());
    double departure = -GasConstant * temperature() +
                       K1 * (Am - temperature() * dAmdT) +
                       pressure() * molarVolume();
    // printf("h0 = %f, departure = %f\n", h0, departure);
    // printf("rt = %f, v = %f, pv - rt = %f\n", GasConstant * temperature(),
    // molarVolume(), pressure() * molarVolume() - GasConstant * temperature());
    return h0 + departure;
}

doublereal PengRobinsonGasPhase::entropy_mole() const
{
    // this is still the ideal gas value
    return GasConstant * (mean_X(entropy_R_ref()) - sum_xlogx() -
                          std::log(pressure() / m_spthermo->refPressure()));
}

doublereal PengRobinsonGasPhase::cp_mole() const
{
    _updateThermoRealFluid();
    doublereal cp0 = GasConstant * mean_X(cp_R_ref()); // ideal gas part
    doublereal departure =
        -GasConstant - K1 * temperature() * d2AmdT2 -
        temperature() * pow(dPdT, 2) / dPdV; // departure function
    //    printf("cp0 = %f, departure = %f\n", cp0, departure);
    //    printf("dPdT = %f, dPdV = %f, d2AmdT2 = %f, K1 = %f\n", dPdT, dPdV,
    //    d2AmdT2, K1);
    // return cp0;
    return cp0 + departure;
}

doublereal PengRobinsonGasPhase::cv_mole() const
{
    _updateThermoRealFluid();
    doublereal cv0 =
        GasConstant * mean_X(cp_R_ref()) - GasConstant; // ideal gas part
    doublereal departure =
        (-temperature() * d2AmdT2 * K1); // departure function
    return cv0 + departure;
}

// TODO: this is currently the same as in IdealGasPhase
doublereal PengRobinsonGasPhase::standardConcentration(size_t k) const
{
    return pressure() / (GasConstant * temperature());
}

// TODO: this is currently the same as in IdealGasPhase
void PengRobinsonGasPhase::getActivityCoefficients(doublereal* ac) const
{
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

// TODO: this is currently the same as in IdealGasPhase
void PengRobinsonGasPhase::getStandardChemPotentials(doublereal* muStar) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), muStar, RT());
    double tmp = log(pressure() / m_spthermo->refPressure()) * RT();
    for (size_t k = 0; k < m_kk; k++) {
        muStar[k] += tmp; // add RT*ln(P/P_0)
    }
}

//  Partial Molar Properties of the Solution --------------

// TODO: this is currently the same as in IdealGasPhase
void PengRobinsonGasPhase::getChemPotentials(doublereal* mu) const
{
    getStandardChemPotentials(mu);
    doublereal rt = temperature() * GasConstant;
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] += rt * (log(xx));
    }
}

void PengRobinsonGasPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    //    const vector_fp& _h = enthalpy_RT_ref();
    //    doublereal rt = GasConstant * temperature();
    //    scale(_h.begin(), _h.end(), hbar, rt);

    _updateThermoRealFluid();
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal rt = GasConstant * temperature();
    scale(_h.begin(), _h.end(), &hbar[0], rt);
    doublereal temp = Am - temperature() * dAmdT;
    for (size_t k = 0; k < m_kk; k++) {
        if (IsCrit[k] == 0) continue;
        hbar[k] += -rt + dK1dN[k] * temp +
                   K1 * (dAmdN[k] - temperature() * d2AmdTdN[k]) +
                   pressure() * dVdN[k];
        // hbar[k] = hbar0[k];
        // printf("T = %.5g, k = %d, X = %.5g, hbar = %.5g, hbar0 = %.5g,
        // departure
        // = %.5g\n", temperature(), k, moleFraction(k), hbar[k], hbar0[k],
        // hbar[k]-hbar0[k]);
        // printf("rt = %f, dvdn = %f, pv - rt = %f\n", rt, dVdN[k], pressure()
        // *
        // dVdN[k] - rt);
    }
}

// TODO: this is currently the same as in IdealGasPhase
void PengRobinsonGasPhase::getPartialMolarEntropies(doublereal* sbar) const
{
    const vector_fp& _s = entropy_R_ref();
    scale(_s.begin(), _s.end(), sbar, GasConstant);
    doublereal logp = log(pressure() / m_spthermo->refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        doublereal xx = std::max(SmallNumber, moleFraction(k));
        sbar[k] += GasConstant * (-logp - log(xx));
    }
}

// TODO: this is currently the same as in IdealGasPhase
void PengRobinsonGasPhase::getPartialMolarIntEnergies(doublereal* ubar) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal rt = GasConstant * temperature();
    for (size_t k = 0; k < m_kk; k++) {
        ubar[k] = rt * (_h[k] - 1.0);
    }
}

// TODO: this is currently the same as in IdealGasPhase
void PengRobinsonGasPhase::getPartialMolarCp(doublereal* cpbar) const
{
    const vector_fp& _cp = cp_R_ref();
    scale(_cp.begin(), _cp.end(), cpbar, GasConstant);
}

// TODO: this is currently the same as in IdealGasPhase
void PengRobinsonGasPhase::getPartialMolarVolumes(doublereal* vbar) const
{
    double vol = 1.0 / molarDensity();
    for (size_t k = 0; k < m_kk; k++) {
        vbar[k] = vol;
    }
}

// Properties of the Standard State of the Species in the Solution --
// TODO: these are all currently the same as in IdealGasPhase

void PengRobinsonGasPhase::getEnthalpy_RT(doublereal* hrt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    copy(_h.begin(), _h.end(), hrt);
}

void PengRobinsonGasPhase::getEntropy_R(doublereal* sr) const
{
    const vector_fp& _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), sr);
    double tmp = log(pressure() / m_spthermo->refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        sr[k] -= tmp;
    }
}

void PengRobinsonGasPhase::getGibbs_RT(doublereal* grt) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt);
    double tmp = log(pressure() / m_spthermo->refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] += tmp;
    }
}

void PengRobinsonGasPhase::getPureGibbs(doublereal* gpure) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), gpure, RT());
    double tmp = log(pressure() / m_spthermo->refPressure());
    tmp *= RT();
    for (size_t k = 0; k < m_kk; k++) {
        gpure[k] += tmp;
    }
}

void PengRobinsonGasPhase::getIntEnergy_RT(doublereal* urt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - 1.0;
    }
}

void PengRobinsonGasPhase::getCp_R(doublereal* cpr) const
{
    const vector_fp& _cpr = cp_R_ref();
    copy(_cpr.begin(), _cpr.end(), cpr);
}

void PengRobinsonGasPhase::getStandardVolumes(doublereal* vol) const
{
    double tmp = 1.0 / molarDensity();
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = tmp;
    }
}

// Thermodynamic Values for the Species Reference States ---------

void PengRobinsonGasPhase::getEnthalpy_RT_ref(doublereal* hrt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    copy(_h.begin(), _h.end(), hrt);
}

void PengRobinsonGasPhase::getGibbs_RT_ref(doublereal* grt) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt);
}

void PengRobinsonGasPhase::getGibbs_ref(doublereal* g) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), g, RT());
}

void PengRobinsonGasPhase::getEntropy_R_ref(doublereal* er) const
{
    const vector_fp& _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), er);
}

void PengRobinsonGasPhase::getIntEnergy_RT_ref(doublereal* urt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - 1.0;
    }
}

void PengRobinsonGasPhase::getCp_R_ref(doublereal* cprt) const
{
    const vector_fp& _cpr = cp_R_ref();
    copy(_cpr.begin(), _cpr.end(), cprt);
}

void PengRobinsonGasPhase::getStandardVolumes_ref(doublereal* vol) const
{
    doublereal tmp = RT() / m_p0;
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = tmp;
    }
}

// Real fluid stuff ---------

doublereal PengRobinsonGasPhase::pressure() const
{
    _updateThermoRealFluid();
    return (GasConstant * temperature()) / (molarVolume() - Bm) -
           Am / (pow(molarVolume(), 2) + 2.0 * molarVolume() * Bm - pow(Bm, 2));
}

void PengRobinsonGasPhase::setPressure(doublereal p)
{
    // SetRealFluidConstants(); // only depends on mass fractions
    // SetRealFluidThermodynamics(); // also depends on temperature except for
    // Bm
    _updateThermoRealFluid();
    setDensity(meanMolecularWeight() /
               GetVolumeFromPressureTemperature(p, temperature()));
    //SetRealFluidThermodynamics();
    //printf("density = %g\n", density());
}

void PengRobinsonGasPhase::ReadCriticalProperties() const
{
    //! \brief Reads critical properties based on JP's papaer
    for (size_t k = 0; k < m_kk; k++) {
        if (speciesName(k) == "H2") {
            IsCrit[k] = 1;
            Tcrit[k] = 33.145;    // K
            Pcrit[k] = 1.2964e+6; // Pa
            Vcrit[k] = 64.4834e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = -0.219;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "O2") {
            IsCrit[k] = 1;
            Tcrit[k] = 154.5800; // K
            Pcrit[k] = 5.0430e+6; // Pa
            Vcrit[k] = 73.3682e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0222;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "H2O") {
            IsCrit[k] = 1;
            Tcrit[k] = 647.10;    // K
            Pcrit[k] = 22.064e+6; // Pa
            Vcrit[k] = 55.95e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.3443;
            sigma[k] = 0.0;
            dipole[k] = 1.855;
        } else if (speciesName(k) == "O") {
            IsCrit[k] = 0;
            Tcrit[k] = 105.28;   // K
            Pcrit[k] = 7.088e+6; // Pa
            Vcrit[k] = 41.21e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "H") {
            IsCrit[k] = 0;
            Tcrit[k] = 190.82; // K
            Pcrit[k] = 31.013e+6; // Pa
            Vcrit[k] = 17.07e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "OH") {
            IsCrit[k] = 0;
            Tcrit[k] = 105.28;   // K
            Pcrit[k] = 7.088e+6; // Pa
            Vcrit[k] = 41.21e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "H2O2") {
            IsCrit[k] = 0;
            Tcrit[k] = 141.34;   // K
            Pcrit[k] = 4.786e+6; // Pa
            Vcrit[k] = 81.93e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "HO2") {
            IsCrit[k] = 0;
            Tcrit[k] = 141.34;   // K
            Pcrit[k] = 4.786e+6; // Pa
            Vcrit[k] = 81.93e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "N2") {
            IsCrit[k] = 1;
            Tcrit[k] = 126.19;    // K
            Pcrit[k] = 3.3958e+6; // Pa
            Vcrit[k] = 89.41e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0372;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "NC12H26") {
            IsCrit[k] = 1;
            Tcrit[k] = 658.10;    // K
            Pcrit[k] = 1.817e+6; // Pa
            Vcrit[k] = 751.8797e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.57639;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "CH4") {
            IsCrit[k] = 1;
            Tcrit[k] = 190.564;    // K
            Pcrit[k] = 4.5992e+6; // Pa
            Vcrit[k] = 98.6291e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.01142;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "CO") {
            IsCrit[k] = 1;
            Tcrit[k] = 132.86;    // K
            Pcrit[k] = 3.4935e+6; // Pa
            Vcrit[k] = 92.1659e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.05;
            sigma[k] = 0.0;
            dipole[k] = 0.1;
        } else if (speciesName(k) == "CO2") {
            IsCrit[k] = 1;
            Tcrit[k] = 304.1282;    // K
            Pcrit[k] = 7.3773e+6; // Pa
            Vcrit[k] = 94.1185e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.22394;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "C2H6") {
            IsCrit[k] = 1;
            Tcrit[k] = 305.33;    // K
            Pcrit[k] = 	4.8718e+6; // Pa
            Vcrit[k] = 145.5604e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0993;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "C2H4") {
            IsCrit[k] = 1;
            Tcrit[k] = 282.35;    // K
            Pcrit[k] = 	5.0418e+6; // Pa
            Vcrit[k] = 130.9414e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0866;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "CH3OH") { // from here
            IsCrit[k] = 1;
            Tcrit[k] = 513;    // K
            Pcrit[k] =  8.01e+06; // Pa
            Vcrit[k] = 0.1179; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.556;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "C2H2") {
            IsCrit[k] = 1;
            Tcrit[k] = 309.2;    // K
            Pcrit[k] =  6.25e+06; // Pa
            Vcrit[k] = 0.1129; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.19;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "AC3H4") {
            IsCrit[k] = 1;
            Tcrit[k] = 385.4;    // K
            Pcrit[k] =  5.32e+06; // Pa
            Vcrit[k] = 0.162; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.313;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "C2H5CHO") {
            IsCrit[k] = 1;
            Tcrit[k] = 515.3;    // K
            Pcrit[k] =  6.33e+06; // Pa
            Vcrit[k] = 0.176; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.313;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "NC4H8") {
            IsCrit[k] = 1;
            Tcrit[k] = 419.6;    // K
            Pcrit[k] =  4.02e+06; // Pa
            Vcrit[k] = 0.240; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.191;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "C4H6") {
            IsCrit[k] = 1;
            Tcrit[k] = 443.7;    // K
            Pcrit[k] =  4.49e+06; // Pa
            Vcrit[k] = 0.219; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.255;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "C3H6") {
            IsCrit[k] = 1;
            Tcrit[k] = 365.0;    // K
            Pcrit[k] =  4.61e+06; // Pa
            Vcrit[k] = 0.1813; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.144;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "NC5H10") {
            IsCrit[k] = 1;
            Tcrit[k] = 464.7;    // K
            Pcrit[k] =  3.55e+06; // Pa
            Vcrit[k] = 0.3026; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.233;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "NC7H14") {
            IsCrit[k] = 1;
            Tcrit[k] = 537.2;    // K
            Pcrit[k] =  2.85e+06; // Pa
            Vcrit[k] = 0.4389; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.358;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "NC10H20") {
            IsCrit[k] = 1;
            Tcrit[k] = 615.0;    // K
            Pcrit[k] =  2.2e+06; // Pa
            Vcrit[k] = 0.650; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.491;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "C3H7CHO") {
            IsCrit[k] = 1;
            Tcrit[k] = 545.4;    // K
            Pcrit[k] =  5.38e+06; // Pa
            Vcrit[k] = 0.2176; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.352;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "NO") {
            IsCrit[k] = 1;
            Tcrit[k] = 180.0;    // K
            Pcrit[k] =  6.48e+06; // Pa
            Vcrit[k] = 0.058; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.582;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "N2O") {
            IsCrit[k] = 1;
            Tcrit[k] = 309.6;    // K
            Pcrit[k] =  7.255e+06; // Pa
            Vcrit[k] = 0.097; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.165;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "NH3") {
            IsCrit[k] = 1;
            Tcrit[k] = 405.5;    // K
            Pcrit[k] =  11.35e+06; // Pa
            Vcrit[k] = 0.0725; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.250;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else if (speciesName(k) == "NO2") {
            IsCrit[k] = 1;
            Tcrit[k] = 431.4;    // K
            Pcrit[k] =  1.013e+06; // Pa
            Vcrit[k] = 0.169; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
        } else {
            IsCrit[k] = 0;
            Tcrit[k] = 126.19;    // K
            Pcrit[k] = 3.3958e+6; // Pa
            Vcrit[k] = 89.41e-3; // m3/kmol
            rhocrit[k] = molecularWeight(k) / Vcrit[k]; // kg/m3
            Zcrit[k] = (Pcrit[k] * Vcrit[k]) / (GasConstant * Tcrit[k]);
            omega[k] = 0.0372;
            sigma[k] = 0.0;
            dipole[k] = 0.0;
            cout << " Unknown or non-major species : " << speciesName(k)
                 << ". All critical properties were set to N2 values with IsCrit 0." << endl;
        }
    }
}

void PengRobinsonGasPhase::SetRealFluidConstants() const
{
    //! \brief Assign arrays of critical point data
    for (size_t k = 0; k < m_kk; k++) {
        if (IsCrit[k] == 0) continue;
        for (size_t l = 0; l < m_kk; l++) {
            if (IsCrit[l] == 0) continue;

            int apos = k * m_kk + l;

            // binary interation parameter
            double tmp_k = 0.0;
            //if (k == l) tmp_k = 0.0;

            Tcrit_IJ[apos] = sqrt(Tcrit[l] * Tcrit[k]) * (1.0 - tmp_k);
            Vcrit_IJ[apos] = pow(pow(Vcrit[l], 1.0 / 3.0) + pow(Vcrit[k], 1.0 / 3.0), 3.0) / 8.0;
            Zcrit_IJ[apos] = 0.5 * (Zcrit[l] + Zcrit[k]);
            Pcrit_IJ[apos] = Zcrit_IJ[apos] * GasConstant * Tcrit_IJ[apos] / Vcrit_IJ[apos];
            omega_IJ[apos] = 0.5 * (omega[l] + omega[k]);
        }
    }

    //! \brief Assign constants used for computing real fluid effects
    for (size_t k = 0; k < m_kk; k++) {
        if (IsCrit[k] == 0) continue;

        cst_b[k] = 0.077796 * GasConstant * Tcrit[k] / Pcrit[k];
        for (size_t l = 0; l < m_kk; l++) {
            if (IsCrit[l] == 0) continue;

            int apos = k * m_kk + l;
            cst_a[apos] = 0.457236 * pow(GasConstant * Tcrit_IJ[apos], 2.0) / Pcrit_IJ[apos];
            cst_c[apos] = 0.37464 + 1.54226 * omega_IJ[apos] - 0.26992 * pow(omega_IJ[apos], 2);

            cst_aa[apos] = cst_a[apos] * pow(cst_c[apos], 2) / Tcrit_IJ[apos];
            cst_bb[apos] = cst_a[apos] * 2.0 * (cst_c[apos] + pow(cst_c[apos], 2)) / sqrt(Tcrit_IJ[apos]);
            cst_cc[apos] = cst_a[apos] * (1.0 + 2.0 * cst_c[apos] + pow(cst_c[apos], 2));
        }
    }

    // compute Bm
    Bm = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        if (IsCrit[k] == 0) continue;

        Bm += moleFraction(k) * cst_b[k];
    }
}

void PengRobinsonGasPhase::SetRealFluidThermodynamics() const
{
    //! \brief Computes the thermodynamic properties of the mixture given a molar volume and temperature
    Am = 0.0;
    dAmdT = 0.0;
    d2AmdT2 = 0.0;
    double temp = molarVolume() * molarVolume() + 2.0 * Bm * molarVolume() - Bm * Bm;
    for (size_t k = 0; k < m_kk; k++) {
        if (IsCrit[k] == 0) continue;

        dAmdN[k] = 0.0;
        d2AmdTdN[k] = 0.0;
        for (size_t l = 0; l < m_kk; l++) {
            if (IsCrit[l] == 0) continue;
            int apos = k * m_kk + l;

            double X_X = moleFraction(l) * moleFraction(k);
            double A_IJ = cst_a[apos] * pow(1.0 + cst_c[apos] * (1.0 - sqrt(temperature() / Tcrit_IJ[apos])), 2);
            double G = cst_c[apos] * sqrt(temperature() / Tcrit_IJ[apos]) / (1.0 + cst_c[apos] * (1.0 - sqrt(temperature() / Tcrit_IJ[apos])));
            double D = cst_c[apos] * (1.0 + cst_c[apos]) * Tcrit_IJ[apos] / Pcrit_IJ[apos] * sqrt(Tcrit_IJ[apos] / temperature());

            Am += X_X * A_IJ;
            dAmdT -= X_X * A_IJ * G;
            d2AmdT2 += X_X * D;

            dAmdN[k] += moleFraction(l) * A_IJ;
            d2AmdTdN[k] += moleFraction(l) * A_IJ * G;
        }
        dAmdN[k] *= 2.0;
        d2AmdTdN[k] *= -2.0 / temperature();
    }
    dAmdT /= temperature();
    d2AmdT2 *= 0.457236 * pow(GasConstant, 2) / (2.0 * temperature());
    dPdT = GasConstant / (molarVolume() - Bm) - dAmdT / (pow(molarVolume(), 2) + 2.0 * molarVolume() * Bm - pow(Bm, 2));
    double arg = GasConstant * temperature() * (molarVolume() + Bm) * pow(((molarVolume()) / (molarVolume() - Bm) + Bm / (molarVolume() + Bm)), 2);
    dPdV = -GasConstant * temperature() / pow((molarVolume() - Bm), 2) * (1.0 - 2.0 * Am / arg);

    //printf("rho = %g, V = %g, R = %g, T = %g, dPdV = %g\n", density(), molarVolume(), GasConstant, temperature(), dPdV);

    //expansivity = -(dPdT) / (molarVolume() * dPdV); //ideal gas: equal to 3.34E-3 (1/K)
    K1 = 1.0 / (sqrt(8.0) * Bm) * log((molarVolume() + (1 - sqrt(2.0)) * Bm) / (molarVolume() + (1 + sqrt(2.0)) * Bm));

    for (size_t k = 0; k < m_kk; k++) {
        if (IsCrit[k] == 0) continue;

        dPdN[k] = GasConstant * temperature() / (molarVolume() - Bm) + GasConstant * temperature() * cst_b[k] / pow((molarVolume() - Bm), 2)
            - dAmdN[k] / temp + 2.0 * Am * cst_b[k] * (molarVolume() - Bm) / pow(temp, 2);
        dVdN[k] = -dPdN[k] / dPdV;
        dK1dN[k] = 1.0 / temp * dVdN[k] - cst_b[k] / Bm * (K1 + molarVolume() / temp);
    }
}

doublereal PengRobinsonGasPhase::GetTemperatureFromPressureDensity_anal(
    doublereal p_in, doublereal rho_in)
{
    double aa = 0.0;
    double bb = 0.0;
    double cc = 0.0;

    Bm = 0.0;
    for (size_t k = 0; k < m_kk; k++)
        Bm += moleFraction(k) * cst_b[k];

    // set rho (density should already been set)
    //setDensity(rho_in);

    for (size_t k = 0; k < m_kk; k++) {
        for (size_t l = 0; l < m_kk; l++) {
            int apos = k * m_kk + l;
            double X_X = moleFraction(l) * moleFraction(k);
            aa -= X_X * cst_aa[apos];
            bb += X_X * cst_bb[apos];
            cc -= X_X * cst_cc[apos];
            //printf("k = %d, l = %d, X_X = %g, cst_aa = %g, cst_bb = %g, cst_cc = %g\n", k, l, X_X, cst_aa[apos], cst_bb[apos], cst_cc[apos]);
        }
    }
    double tmp = 1.0 / (pow(molarVolume(), 2) + 2.0 * Bm * molarVolume() - pow(Bm, 2));
    //printf("tmp = %g\n", tmp);
    aa *= tmp;
    bb *= tmp;
    cc *= tmp;
    aa += GasConstant / (molarVolume() - Bm);
    cc -= p_in;

    //printf("p_in = %g, rho_in = %g, v = %g, R = %g, Bm = %g\n", p_in, rho_in, molarVolume(), GasConstant, Bm);
    //printf("aa = %g, bb = %g, cc = %g\n", aa, bb, cc);

    return pow((-bb + sqrt(pow(bb, 2) - 4.0 * aa * cc)) / 2.0 / aa, 2);
}

doublereal PengRobinsonGasPhase::GetVolumeFromPressureTemperature(
    doublereal p_in, doublereal T_in) const
{
    // Computation for Z directly
    doublereal Amat = Am * p_in / pow(GasConstant * T_in, 2);
    doublereal Bmat = Bm * p_in / (GasConstant * T_in);

    doublereal coef_v0 = (pow(Bmat, 3) + pow(Bmat, 2) - Amat * Bmat);
    doublereal coef_v1 = (-3.0 * pow(Bmat, 2) - 2.0 * Bmat + Amat);
    doublereal coef_v2 = (Bmat - 1.0);
    return GasConstant * T_in * (GetCubicRoots(coef_v0, coef_v1, coef_v2)) / p_in;
}

doublereal PengRobinsonGasPhase::GetCubicRoots(doublereal a0, doublereal a1,
                                               doublereal a2) const
{
    //! \brief Solves the cubic root equation given three parameters
    const double eps = 1.0e-12;
    double p, q, u, v, Det;
    double arg, signArg;
    double phi;

    double Z  = 0.0;
    double Z1 = 0.0;
    double Z2 = 0.0;
    double Z3 = 0.0;

    p = (3.0 * a1 - pow(a2, 2)) / 3.0;
    q = a0 + 2 * pow(a2, 3) / 27.0 - a2 * a1 / 3.0;

    Det = pow(p / 3.0, 3) + pow(q / 2.0, 2);

    if (Det > 0) { // only one real root
        arg     = -q/2.0 + sqrt(Det); 
        signArg = (arg>0)?1.0:-1.0; 
        u       = signArg*pow(fabs(arg), 1.0/3.0);
        arg     = -q/2.0 - sqrt(Det); 
        signArg = (arg>0)?1.0:-1.0; 
        v       = signArg*pow(fabs(arg), 1.0/3.0);
        Z       = -a2/3.0 + u + v;
    } else if (fabs(Det) <= eps) {
        arg     = -q/2.0;
        signArg = (arg>0)?1.0:-1.0;
        u       = signArg*pow(fabs(arg), 1.0/3.0);
        arg     = -q/2.0;
        signArg = (arg>0)?1.0:-1.0;
        v       = signArg*pow(fabs(arg), 1.0/3.0);
        Z       = -a2/3.0 + u + v;
    } else {
        arg = -q / (2.0 * sqrt(pow(fabs(p) / 3.0, 3.0)));
        phi = acos(arg);

        Z1 = -a2 / 3.0 + 2.0 * sqrt(fabs(p) / 3.0) * cos((phi) / 3.0);
        Z2 = -a2 / 3.0 - 2.0 * sqrt(fabs(p) / 3.0) * cos((phi - Pi) / 3.0);
        Z3 = -a2 / 3.0 - 2.0 * sqrt(fabs(p) / 3.0) * cos((phi + Pi) / 3.0);

        //jph For fully consistent results we should impose a fugacity constraint
        Z = min(min(Z1, Z2), Z3);
        if (Z < 0) {
            Z = max(max(Z1, Z2), Z3);
        }
    }
    // double sanity = pow(Z,3)+pow(Z,2)*a2+Z*a1+a0;
    //  printf("Sanity = %e\n ", sanity);
    return Z;
}

doublereal PengRobinsonGasPhase::critTemperature() const
{
    return mean_X(Tcrit);
}

doublereal PengRobinsonGasPhase::critPressure() const
{
    return mean_X(Pcrit);
}

doublereal PengRobinsonGasPhase::critVolume() const
{
    return mean_X(Vcrit);
}

void PengRobinsonGasPhase::getCritTemperature(doublereal* const tc) const
{
    copy(Tcrit.begin(), Tcrit.end(), tc);
}

void PengRobinsonGasPhase::getCritPressure(doublereal* const pc) const
{
    copy(Pcrit.begin(), Pcrit.end(), pc);
}

void PengRobinsonGasPhase::getCritVolume(doublereal* const vc) const
{
    copy(Vcrit.begin(), Vcrit.end(), vc);
}

void PengRobinsonGasPhase::getCritCompressibility(doublereal* const zc) const
{
    copy(Zcrit.begin(), Zcrit.end(), zc);
}

void PengRobinsonGasPhase::getAcentricFactor(doublereal* const om) const
{
    copy(omega.begin(), omega.end(), om);
}

void PengRobinsonGasPhase::getDipoleMoment(doublereal* const di) const
{
    copy(dipole.begin(), dipole.end(), di);
}

doublereal PengRobinsonGasPhase::critCompressibility() const
{
    return mean_X(Zcrit);
}

void PengRobinsonGasPhase::initThermo()
{
    ThermoPhase::initThermo();

    // read in critical properties (currently hard-coded)
    ReadCriticalProperties();

    Tcrit_IJ.resize(m_kk * m_kk);
    Pcrit_IJ.resize(m_kk * m_kk);
    Vcrit_IJ.resize(m_kk * m_kk);
    Zcrit_IJ.resize(m_kk * m_kk);
    omega_IJ.resize(m_kk * m_kk);
    fill(Tcrit_IJ.begin(), Tcrit_IJ.end(), 0);
    fill(Pcrit_IJ.begin(), Pcrit_IJ.end(), 0);
    fill(Vcrit_IJ.begin(), Vcrit_IJ.end(), 0);
    fill(Zcrit_IJ.begin(), Zcrit_IJ.end(), 0);
    fill(omega_IJ.begin(), omega_IJ.end(), 0);

    cst_a.resize(m_kk * m_kk);
    cst_b.resize(m_kk);
    cst_c.resize(m_kk * m_kk);
    fill(cst_a.begin(), cst_a.end(), 0);
    fill(cst_b.begin(), cst_b.end(), 0);
    fill(cst_c.begin(), cst_c.end(), 0);

    cst_aa.resize(m_kk * m_kk);
    cst_bb.resize(m_kk * m_kk);
    cst_cc.resize(m_kk * m_kk);
    fill(cst_aa.begin(), cst_aa.end(), 0);
    fill(cst_bb.begin(), cst_bb.end(), 0);
    fill(cst_cc.begin(), cst_cc.end(), 0);

    dPdN.resize(m_kk);
    dVdN.resize(m_kk);
    dAmdN.resize(m_kk);
    d2AmdTdN.resize(m_kk);
    dK1dN.resize(m_kk);
}

bool PengRobinsonGasPhase::addSpecies(shared_ptr<Species> spec)
{
    bool added = ThermoPhase::addSpecies(spec);
    if (added) {
        if (m_kk == 1) {
            m_p0 = refPressure();
        }
        m_h0_RT.push_back(0.0);
        m_g0_RT.push_back(0.0);
        m_expg0_RT.push_back(0.0);
        m_cp0_R.push_back(0.0);
        m_s0_R.push_back(0.0);
        m_pp.push_back(0.0);

        // real fluid stuff
        IsCrit.push_back(0);
        Tcrit.push_back(0.0);
        Pcrit.push_back(0.0);
        rhocrit.push_back(0.0);
        Vcrit.push_back(0.0);
        Zcrit.push_back(0.0);
        omega.push_back(0.0);
        sigma.push_back(0.0);
        dipole.push_back(0.0);
    }
    return added;
}

void PengRobinsonGasPhase::setToEquilState(const doublereal* mu_RT)
{
    const vector_fp& grt = gibbs_RT_ref();

    /*
     * Within the method, we protect against inf results if the
     * exponent is too high.
     *
     * If it is too low, we set
     * the partial pressure to zero. This capability is needed
     * by the elemental potential method.
     */
    doublereal pres = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        double tmp = -grt[k] + mu_RT[k];
        if (tmp < -600.) {
            m_pp[k] = 0.0;
        } else if (tmp > 300.0) {
            double tmp2 = tmp / 300.;
            tmp2 *= tmp2;
            m_pp[k] = m_p0 * exp(300.) * tmp2;
        } else {
            m_pp[k] = m_p0 * exp(tmp);
        }
        pres += m_pp[k];
    }
    // set state
    setState_PX(pres, &m_pp[0]);
}

void PengRobinsonGasPhase::_updateThermo() const
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    doublereal tnow = temperature();

    // If the temperature has changed since the last time these
    // properties were computed, recompute them.
    if (cached.state1 != tnow) {
        m_spthermo->update(tnow, &m_cp0_R[0], &m_h0_RT[0], &m_s0_R[0]);
        cached.state1 = tnow;
        // update the species Gibbs functions
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
    }
}

void PengRobinsonGasPhase::_updateThermoRealFluid() const
{
    static const int cacheId = m_cache.getId();
    CachedArray cached = m_cache.getArray(cacheId);

    doublereal tnow = temperature();
    doublereal rhonow = density();
    vector_fp ynow(m_kk);
    getMassFractions(&ynow[0]);

    if (cached.value != ynow) {
        // update everything since mass fractions have changed
        SetRealFluidConstants();
        SetRealFluidThermodynamics();
        cached.value = ynow;
        cached.state1 = tnow;
        cached.state2 = rhonow;
    } else if (cached.state1 != tnow || cached.state2 != rhonow) {
        // only update thermodynamics since only temperature/density has changed
        SetRealFluidThermodynamics();
        cached.state1 = tnow;
        cached.state2 = rhonow;
    }
}

}
