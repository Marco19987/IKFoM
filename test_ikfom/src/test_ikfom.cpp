#include <Eigen/Eigen>
#include "../include/test_ikfom/use-ikfom.hpp"

#define MAXIMUM_ITER 10 // check the meaning
int main(int argc, char **argv)
{
  // define variables
  double dt = 0.01;                                                // time step
  state_ikfom ekf_state;                                           // state vector (struct)
  input_ikfom u;                                                   // input vector (struct)
  measurement_ikfom z;                                             // measurement vector (struct)
  MTK::get_cov<process_noise_ikfom>::type Q = process_noise_cov(); // process noise covariance: Q, an Eigen matrix
  MTK::get_cov<measurement_noise_ikfom>::type R =
      measurement_noise_cov(); // measurement noise covariance: R, an Eigen matrix
  esekfom::esekf<state_ikfom, process_noise_dof, input_ikfom, measurement_ikfom, measurement_noise_dof>::cov init_P;

  // initalize variables
  ekf_state.bpo << 0.0, 0.0, 0.0;               // initial position
  ekf_state.bQo = MTK::SO3<double>::Identity(); // initial orientation
  u.ovo << 0.1, 0.0, 0.0;                       // input velocity in body frame
  u.oomegao << 0.1, 0.0, 0.0;                   // input angular velocity in body frame
  // z.bpo_meas << 0.0, 0.0, 0.1;                  // measurement position in base frame
  // z.bQo_meas = MTK::SO3<double>::Identity();    // measurement orientation in base frame
  bool valid = false;
  z = get_h(ekf_state, valid);

  // create the filter object
  esekfom::esekf<state_ikfom, process_noise_dof, input_ikfom, measurement_ikfom, measurement_noise_dof> kf(ekf_state,
                                                                                                           init_P);

  double epsi[state_dof] = {0.001};
  std::fill(epsi, epsi + state_dof,
            0.001); // if the absolute of innovation of ekf update is smaller than epso, the update iteration is
                    // converged
  kf.init(get_f, df_dx, df_dw, get_h, dh_dx, dh_dv, MAXIMUM_ITER, epsi);

  // setup simulator real process
  state_ikfom true_state;
  true_state.bpo << 0.3, 0.2, 0.1;                                  // initial position
  true_state.bQo = Eigen::Quaternion<double>(0.7071, 0.7071, 0, 0); // initial orientation
  double tf = 10;
  double t = 0.0; // simulation time

  while (t < tf)
  {
    // predict
    kf.predict(dt, Q, u);

    // update
    kf.update_iterated(z, R);

    // // print the state
    // state_ikfom cur_state_aft = kf.get_x();

    // simulate real system
    true_state.boxplus(get_f(true_state, u), dt);
    z = get_h(true_state, valid);

    std::cout << "Real state: \n"
              << true_state << std::endl;

    std::cout << "Estimated state: \n"
              << kf.get_x() << std::endl;

    t = t + dt;
  }

  return 0;
}
