/**
 *  @file ConstDensityThermo.cpp
 * Declarations for a Thermo manager for incompressible ThermoPhases
 * (see \ref thermoprops and \link Cantera::ConstDensityThermo ConstDensityThermo
\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/mix_defs.h"
#include "cantera/thermo/ConstDensityThermo.h"
#include "cantera/base/ctml.h"

namespace Cantera
{

ConstDensityThermo::ConstDensityThermo(const ConstDensityThermo& right)
{
    *this = right;
}

ConstDensityThermo& ConstDensityThermo::operator=(const ConstDensityThermo& right)
{
    if (&right == this) {
        return *this;
    }

    m_h0_RT = right.m_h0_RT;
    m_cp0_R = right.m_cp0_R;
    m_g0_RT = right.m_g0_RT;
    m_s0_R = right.m_s0_R;

    return *this;

}

ThermoPhase* ConstDensityThermo::duplMyselfAsThermoPhase() const
{
    return new ConstDensityThermo(*this);
}

int ConstDensityThermo::eosType() const
{
    warn_deprecated("ConstDensityThermo::eosType",
                    "To be removed after Cantera 2.3.");
    return cIncompressible;
}

doublereal ConstDensityThermo::enthalpy_mole() const
{
    doublereal p0 = m_spthermo->refPressure();
    return RT() * mean_X(enthalpy_RT()) + (pressure() - p0)/molarDensity();
}

doublereal ConstDensityThermo::entropy_mole() const
{
    return GasConstant * (mean_X(entropy_R()) - sum_xlogx());
}

doublereal ConstDensityThermo::cp_mole() const
{
    return GasConstant * mean_X(cp_R());
}

doublereal ConstDensityThermo::cv_mole() const
{
    return cp_mole();
}

doublereal ConstDensityThermo::pressure() const
{
    return m_press;
}

void ConstDensityThermo::setPressure(doublereal p)
{
    m_press = p;
}

void ConstDensityThermo::getActivityConcentrations(doublereal* c) const
{
    getConcentrations(c);
}

void ConstDensityThermo::getActivityCoefficients(doublereal* ac) const
{
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

doublereal ConstDensityThermo::standardConcentration(size_t k) const
{
    return molarDensity();
}

void ConstDensityThermo::getChemPotentials(doublereal* mu) const
{
    doublereal vdp = (pressure() - m_spthermo->refPressure())/
                     molarDensity();
    const vector_fp& g_RT = gibbs_RT();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] = RT()*(g_RT[k] + log(xx)) + vdp;
    }
}


void ConstDensityThermo::getStandardChemPotentials(doublereal* mu0) const
{
    getPureGibbs(mu0);
}

bool ConstDensityThermo::addSpecies(shared_ptr<Species> spec)
{
    bool added = ThermoPhase::addSpecies(spec);
    if (added) {
        m_h0_RT.push_back(0.0);
        m_g0_RT.push_back(0.0);
        m_cp0_R.push_back(0.0);
        m_s0_R.push_back(0.0);
    }
    return added;
}

void ConstDensityThermo::_updateThermo() const
{
    doublereal tnow = temperature();
    if (m_tlast != tnow) {
        m_spthermo->update(tnow, &m_cp0_R[0], &m_h0_RT[0],
                           &m_s0_R[0]);
        m_tlast = tnow;
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
        m_tlast = tnow;
    }
}

void ConstDensityThermo::setParametersFromXML(const XML_Node& eosdata)
{
    eosdata._require("model","Incompressible");
    doublereal rho = getFloat(eosdata, "density", "toSI");
    setDensity(rho);
}

}
