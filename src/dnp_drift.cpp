#include "openmc/dnp_drift.h"

#include <dlfcn.h>
#include <string>

#include "openmc/simulation.h"

namespace openmc {

void initialize_dnp_drift()
{
  typedef void (*dnp_init_handle)(std::string, std::string, std::string, double,
    double, std::map<std::string, std::vector<int>>);

  simulation::dnp_library =
    dlopen(settings::dnp_drift_library_path.c_str(), RTLD_LAZY);

  if (!simulation::dnp_library) {
    fatal_error("Error loading external library for precursor drift!");
  }

  auto dnp_init = reinterpret_cast<dnp_init_handle>(
    dlsym(simulation::dnp_library, "initialize"));
  auto dlsym_error = dlerror();
  if (dlsym_error) {
    dlclose(simulation::dnp_library);
    fatal_error(dlsym_error);
  }
  dnp_init(settings::nekrs_re2_path, settings::nekrs_fld_path,
    settings::dnp_drift_method, settings::dnp_drift_dt,
    settings::dnp_drift_external_time, settings::dnp_drift_bcs);

  simulation::dnp_transport =
    reinterpret_cast<simulation::dnp_transport_handle>(
      dlsym(simulation::dnp_library, "transport_site"));
  dlsym_error = dlerror();
  if (dlsym_error) {
    dlclose(simulation::dnp_library);
    fatal_error(dlsym_error);
  }
}

bool transport_dnp(double dnp_decay_time, SourceSite* site, Particle& p)
{
  bool available = simulation::dnp_transport(
    site->r.x, site->r.y, site->r.z, dnp_decay_time, *p.current_seed());

  if (!available) {
    return false;
  }
  return true;
}

void finalize_dnp_drift()
{
  if (simulation::dnp_library) {
    dlclose(simulation::dnp_library);
  }
}

} // namespace openmc
