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
  double h_channel = settings::dnp_drift_msre_h_channel;
  double h_upper_head = settings::dnp_drift_msre_h_upper_head;
  double mean_velocity_upper_head = settings::dnp_drift_msre_v_upper_head;

  // Method selection
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
  if (settings::dnp_drift_msre_representation == "channel") {
    if (z <= h_channel && z >= 0.) {
      location = "channel";
    } else {
      throw std::runtime_error {"Unknown location!"};
    }
  } else if (settings::dnp_drift_msre_representation == "channel_upper_head") {
    if (z <= h_channel) {
      location = "channel";
    } else if (z <= (h_channel + h_upper_head)) {
      location = "upper_head";
    } else {
      throw std::runtime_error {"Unknown location!"};
    }
  }

  // Transport
  while (remaining_time > 0.){

    // Travelling in the channel
    if (location == "channel") {

      // First approach: residence time
      if (settings::dnp_drift_method == "residence-time") {

        dist = h_channel - z;
        time = dist / settings::dnp_drift_msre_v_channel;

        // Decay in the channel for this iteration
        if (time > remaining_time) {
          z += remaining_time * settings::dnp_drift_msre_v_channel;
          remaining_time = 0.;
          break;

        // Continue to next location
        } else {
          if (settings::dnp_drift_msre_representation == "channel") {
            location = "outside";
          } else if (settings::dnp_drift_msre_representation == "channel_upper_head") {
            location = "upper_head";
          } else {
            fatal_error("Not implemented yet.");
          }
          z = h_channel;
          remaining_time -= time;
        }

      // Second approach: explicit transport using the external transport library
      } else if (settings::dnp_drift_method == "streamline") {

        // Transport
        if (!simulation::dnp_transport(
          x, y, z, time, remaining_time, *p.current_seed())) {

          remaining_time = time;

          // If there is still time, we need to continue to the next part
          if (remaining_time > 0) {
            if (settings::dnp_drift_msre_representation == "channel") {
              location = "outside";
            } else if (settings::dnp_drift_msre_representation == "channel_upper_head") {
              location = "upper_head";
            } else {
              fatal_error("Not implemented yet.");
            }
          // The particle stopped right at the outlet
          } else {
            break;
          }
        } else {
          remaining_time = time;
          break;
        }
      }

    // Travelling in the upper head
    } else if (location == "upper_head") {

      // Upper head is residence-time only by default.
      // Because it has the largest cross-section on the XY plane, we do not verify
      // that the point is inside the XY cross-section.

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
          x = new_x_r;
          y = new_y_r;
          z = new_z;
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

      // Continue to the channel
      } else {
        location = "channel";
        remaining_time -= settings::dnp_drift_external_time;
        z = 0.;

        // If x and y are not in the channel cross-section, we resample x and y
        if (!is_inside_msre_channel_2d(x, y)) {
          resample_msre_channel_2d(x, y, p);
        }
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

bool is_inside_msre_channel_2d(double x, double y)
{
  double IN_TO_CM = 2.54;

  // The 2D cross-section of the MSRE channel is composed of a rectangle and two
  // disc halves.

  // Verify if point is in rectangle
  double rect_x_min = -0.2 * IN_TO_CM;
  double rect_x_max = 0.2 * IN_TO_CM;
  double rect_y_min = -0.4 * IN_TO_CM;
  double rect_y_max = 0.4 * IN_TO_CM;

  if ((x >= rect_x_min) && (x <= rect_x_max) && (y >= rect_y_min) &&
      (y <= rect_y_max)) {
    return true;
  }

  // Verify if point is in the first circle
  double c1_x = 0.0;
  double c1_y = -0.4 * IN_TO_CM;
  double c1_r = 0.2 * IN_TO_CM;

  if (sqrt(pow(x - c1_x, 2) + pow(y - c1_y, 2)) <= c1_r) {
    return true;
  }

  // Verify if point is in the second circle
  double c2_x = 0.0;
  double c2_y = 0.4 * IN_TO_CM;
  double c2_r = 0.2 * IN_TO_CM;

  if (sqrt(pow(x - c2_x, 2) + pow(y - c2_y, 2)) <= c2_r) {
    return true;
  }

  // Return false otherwise
  return false;
}

void resample_msre_channel_2d(double& x, double& y, Particle& p)
{
  double IN_TO_CM = 2.54;

  // Declare axis-aligned bounding box of the channel cross-section
  double aabb_min_x = -0.2 * IN_TO_CM;
  double aabb_max_x = +0.2 * IN_TO_CM;
  double aabb_min_y = -0.6 * IN_TO_CM;
  double aabb_max_y = +0.6 * IN_TO_CM;

  // Initialization
  bool found = false;
  double sampled_x;
  double sampled_y;
  int iter = 0;

  while (!found) {

    if (iter > 1000) {
      fatal_error("Could not sample a point inside the MSRE channel 2D cross-section!");
    }

    // Sample positions
    sampled_x = aabb_min_x + (aabb_max_x - aabb_min_x) * prn(p.current_seed());
    sampled_y = aabb_min_y + (aabb_max_y - aabb_min_y) * prn(p.current_seed());

    // Check if the sampled positions are inside the channel cross-section
    if (is_inside_msre_channel_2d(sampled_x, sampled_y)) {
      found = true;
      x = sampled_x;
      y = sampled_y;
      break;
    }

    iter++;
  }
}

void finalize_dnp_drift()
{
  if (simulation::dnp_library) {
    dlclose(simulation::dnp_library);
  }
}

} // namespace openmc
