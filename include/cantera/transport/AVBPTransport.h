/**
 *
 *  @file AVBPTransport.h
 *   Header file defining class AVBPTransport which implements
 *   the simplified transport model used in the solver AVBP
 */

/* $Author: B. Franzelli (v. 1.7) $
 * $Revision: A. Felden (v 2.1-2.3) $
 * $Date: 01/2018 $
 */


#ifndef CT_AVBPTRAN_H
#define CT_AVBPTRAN_H

#include "GasTransport.h"
#include "cantera/numerics/DenseMatrix.h"

// STL includes
#include <vector>
#include <string>
#include <map>
#include <numeric>
#include <algorithm>

namespace Cantera
{

class GasTransportParams;

/**
 * Class AVBPTransport
 * Constant Sch for each species, Pr, and simplified viscosity
 */
class AVBPTransport : public GasTransport
{

public:
    //! Default constructor.
    AVBPTransport();

    //! Return the model id for transport
    /*!
     * @return cAVBPAverage
     */
    virtual int model() const {
        warn_deprecated("AVBPTransport::model",
                        "To be removed after Cantera 2.3.");
        //return cAVBPTransport;
        return 0;
    }

    virtual std::string transportType() const {
        return "AVBP";
    }

    //! Return the thermal diffusion coefficients
    //virtual void getThermalDiffCoeffs(doublereal* const dt);

    //! Returns the mixture thermal conductivity
    virtual doublereal thermalConductivity();

    //! Get the Electrical mobilities (m^2/V/s).
    //virtual void getMobilities(doublereal* const mobil);

    //! Update the internal parameters whenever the temperature has changed
    virtual void update_T();

    //! Update the internal parameters whenever the concentrations have changed
    virtual void update_C();

    //virtual void getSpeciesFluxes(size_t ndim,
    //                              const doublereal* const grad_T,
    //                              size_t ldx,
    //                              const doublereal* const grad_X,
    //                              size_t ldf, doublereal* const fluxes);

    //! Initialize the transport object
    //virtual bool initGas(GasTransportParams& tr);
    virtual void init(thermo_t* thermo, int mode=0, int log_level=0);

    //! Viscosity of the mixture
    virtual doublereal viscosity();

    virtual void getSpeciesViscosities(doublereal* const visc) {
        update_T();
        updateViscosity_T();
        std::copy(m_visc.begin(), m_visc.end(), visc);
    }

    //! Mixture diffusion coefficients [m^2/s].
    virtual void getMixDiffCoeffs(doublereal* const d);

    virtual void read_mixture(std::string s);

    size_t avbp_ipea;
    vector_fp avbp_pea_coeffs;

private:

    //! Calculate the pressure from the ideal gas law
    doublereal pressure_ig() const {
        return (m_thermo->molarDensity() * GasConstant *
                m_thermo->temperature());
    }

    doublereal m_lambda;
    //bool m_debug;

    // AVBP variables
    vector_fp avbp_Sch;
    vector_fp avbp_Le;
    doublereal avbp_Prandtl;
    doublereal avbp_mu0;
    doublereal avbp_T0;
    doublereal avbp_beta;
    std::string avbp_fuel;

};
}
#endif
