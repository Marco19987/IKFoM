#include <Eigen/Eigen>
#include "../include/test_ikfom/use-ikfom.hpp"

#define MAXIMUM_ITER 10  // check the meaning
int main(int argc, char** argv)
{
  // define variables
  double dt = 0.01;                                                 // time step
  state_ikfom init_state;                                           // state vector (struct)
  input_ikfom u;                                                    // input vector (struct)
  measurement_ikfom z;                                              // measurement vector (struct)
  MTK::get_cov<process_noise_ikfom>::type Q = process_noise_cov();  // process noise covariance: Q, an Eigen matrix
  MTK::get_cov<measurement_noise_ikfom>::type R =
      measurement_noise_cov();  // measurement noise covariance: R, an Eigen matrix
  esekfom::esekf<state_ikfom, process_noise_dof, input_ikfom, measurement_ikfom, measurement_noise_dof>::cov init_P;

  // initalize variables
  init_state.bpo << 0.0, 0.0, 0.0;                // initial position
  init_state.bQo = MTK::SO3<double>::Identity();  // initial orientation
  u.ovo << 0.0, 0.0, 0.0;                         // input velocity in body frame
  u.oomegao << 0.0, 0.0, 0.0;                     // input angular velocity in body frame
  z.bpo_meas << 0.0, 0.0, 0.1;                    // measurement position in base frame
  z.bQo_meas = MTK::SO3<double>::Identity();      // measurement orientation in base frame

  // create the filter object
  esekfom::esekf<state_ikfom, process_noise_dof, input_ikfom, measurement_ikfom, measurement_noise_dof> kf(init_state,
                                                                                                           init_P);

  double epsi[state_dof] = { 0.001 };
  std::fill(epsi, epsi + state_dof,
            0.001);  // if the absolute of innovation of ekf update is smaller than epso, the update iteration is
                     // converged
  kf.init(get_f, df_dx, df_dw, get_h, dh_dx, dh_dv, MAXIMUM_ITER, epsi);

  // predict
  kf.predict(dt, Q, u);

  // update
  kf.update_iterated(z, R);

  // print the state
  state_ikfom cur_state_aft = kf.get_x();

  std::cout << "Estimated state: \n" << cur_state_aft << std::endl;

  return 0;
}
