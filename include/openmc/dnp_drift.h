#ifndef OPENMC_DNP_DRIFT_H
#define OPENMC_DNP_DRIFT_H

#include "openmc/particle.h"
#include "openmc/particle_data.h"

namespace openmc {

//! Initialize the connection to the external transport library for the explicit
//! transport of Delayed Neutron Precursors.
void initialize_dnp_drift();

//! Explicit transport of Delayed Neutron Precursor (DNP).
//!
//! \param[in] dnp_decay_time Decay time of the DNP
//! \param[in,out] site Fission site corresponding to the emitted delayed
//! neutron
//! \param[in, out] p Particle at the origin of the fission event passed for
//! the random seed
//! \return true if the DNP is still in the modeled system, false otherwise.
bool transport_dnp(double dnp_decay_time, SourceSite* site, Particle& p);

//! Free memory associated with the external transport library.
void finalize_dnp_drift();

} // namespace openmc

#endif // OPENMC_DNP_DRIFT_H
