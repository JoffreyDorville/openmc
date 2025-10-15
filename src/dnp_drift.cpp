#include "openmc/dnp_drift.h"

#include <dlfcn.h>
#include <string>

#include "openmc/simulation.h"

namespace openmc {

void initialize_dnp_drift()
{
  if (settings::dnp_drift_method == "streamline") {
    typedef void (*dnp_init_handle)(std::string, std::string, std::string, double,
      double, std::map<std::string, std::vector<int>>, bool);

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
    dnp_init(settings::dnp_drift_mesh_path, settings::dnp_drift_field_path,
      settings::dnp_drift_integration_method, settings::dnp_drift_dt,
      settings::dnp_drift_external_time, settings::dnp_drift_bcs,
      settings::dnp_drift_recycling);

    simulation::dnp_transport =
      reinterpret_cast<simulation::dnp_transport_handle>(
        dlsym(simulation::dnp_library, "transport_site"));
    dlsym_error = dlerror();
    if (dlsym_error) {
      dlclose(simulation::dnp_library);
      fatal_error(dlsym_error);
    }
  }
}

bool transport_dnp(double dnp_decay_time, SourceSite* site, Particle& p)
{
  if (settings::dnp_drift_model == "msre") {
    return transport_dnp_msre(dnp_decay_time, site, p);
  } else {
    double time;
    return simulation::dnp_transport(
      site->r.x, site->r.y, site->r.z, time, dnp_decay_time, *p.current_seed());
  }
}

bool transport_dnp_msre(double dnp_decay_time, SourceSite* site, Particle& p)
{
  // MSRE characteristics
  double h_channel = 170.311;
  double mean_velocity_channel = settings::dnp_drift_dt; // Use dritf_dt settings as the mean velocity channel
  double h_upper_head = 17.13;
  double mean_velocity_upper_head = h_upper_head / 3.9;

  // Method selection
  //std::string channel_method = "residence_time";
  std::string channel_method = "explicit";
  std::string upper_head_method = "random";
  //std::string upper_head_method = "linear_z_random_x_y";

  // Initialization
  double x = site->r.x;
  double y = site->r.y;
  double z = site->r.z;
  bool available = true;
  double remaining_time = dnp_decay_time;
  double dist;
  double time;

  // Define location from coordinates
  std::string location;
  if (z <= h_channel) {
    location = "channel";
  } else if (z <= (h_channel + h_upper_head)) {
    location = "upper_head";
  } else {
    throw std::runtime_error {"Unknown location!"};
  }

  // Transport
  while (remaining_time > 0.){

    // Travelling in the channel
    if (location == "channel") {

      // First approach: residence time
      if (channel_method == "residence_time") {

        dist = h_channel - z;
        time = dist / mean_velocity_channel;

        // Decay in the channel for this iteration
        if (time > remaining_time) {
          z += remaining_time * mean_velocity_channel;
          remaining_time = 0.;
          break;

        // Continue to next location
        } else {
          location = "upper_head";
          z = h_channel;
          remaining_time -= time;
        }

      // Second approach: explicit transport using the external transport library
      } else if (channel_method == "explicit") {

        if (!transport_dnp(remaining_time, time, site, p)) {

          x = site->r.x;
          y = site->r.y;
          z = site->r.z;
          remaining_time = time;

          // If there is still time, we need to continue to the upper head
          if (time > 0) {
            location = "upper_head";
          // The particle stopped right at the outlet
          } else {
            break;
          }
        } else {
          x = site->r.x;
          y = site->r.y;
          z = site->r.z;
          remaining_time = time;
          break;
        }
      }

    // Travelling in the upper head
    } else if (location == "upper_head") {

      dist = (h_channel + h_upper_head) - z;
      time = dist / mean_velocity_upper_head;

      // Decay in the upper head for this iteration
      if (time > remaining_time) {

        std::string upper_head_method = "random";

        // First approach: random
        if (upper_head_method == "random") {

          remaining_time = 0.;

          // Sample x and y for mixing
          double a = 2.54 * sqrt(2) / 2.;

          // Sample uniformly in an axis-aligned representation of the upper head
          double new_x = -a + prn(p.current_seed()) * 2 * a;
          double new_y = -a + prn(p.current_seed()) * 2 * a;
          double new_z = h_channel + prn(p.current_seed()) * h_upper_head;

          // Rotate 45 degree
          double coef = sqrt(2)/2.;
          double new_x_r = new_x * coef - new_y * coef;
          double new_y_r = new_x * coef + new_y * coef;

          // Store new coordinates
          double x = new_x_r;
          double y = new_y_r;
          double z = new_z;
          break;

        // Second approach: adjust z linearily, randomly sample x and y
        } else if (upper_head_method == "linear_z_random_x_y") {

          z += remaining_time * mean_velocity_upper_head;
          remaining_time = 0.;

          // Sample x and y for mixing
          double a = 2.54 * sqrt(2) / 2.;

          // Sample uniformly in an axis-aligned representation of the upper head
          double new_x = -a + prn(p.current_seed()) * 2 * a;
          double new_y = -a + prn(p.current_seed()) * 2 * a;

          // Rotate 45 degree
          double coef = sqrt(2)/2.;
          double new_x_r = new_x * coef - new_y * coef;
          double new_y_r = new_x * coef + new_y * coef;

          // Store new coordinates
          double x = new_x_r;
          double y = new_y_r;
          break;
        
        // Approach not implemented
        } else {
          throw std::runtime_error {"Not implemented!!!"};
        }

      // Continue to next location
      } else {
        location = "outside";
        z = h_channel + h_upper_head;
        remaining_time -= time;
      }

    // Travelling outside
    } else if (location == "outside") {

      // Decay outside during this iteration
      if (remaining_time < settings::dnp_drift_external_time) {
        available = false;
        remaining_time = 0.;
        break;

      // Continue to next location
      } else {
        location = "channel";
        remaining_time -= settings::dnp_drift_external_time;
        z = 0.;
      }

    // Error
    } else{ 
      throw std::runtime_error {"Unknown location!!!"};
    }
  }

  site->r.x = x;
  site->r.y = y;
  site->r.z = z;

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
