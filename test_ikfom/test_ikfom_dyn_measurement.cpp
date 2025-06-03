#include <Eigen/Eigen>
// #include "use-ikfom.hpp"
// #include  <boost/preprocessor/repetition/repeat_from_to.hpp>
#include "robot_spring_object_system_ikfom.hpp"

#define MAXIMUM_ITER 10  // check the meaning
// MTK_BUILD_MANIFOLD(measurement_dyn, ((vect3, bpo_meas))((SO3, bQo_meas)));

// #define DECL_SE_VEC_MANIFOLD(z, n, text) ((vect3, po_meas ## n))((SO3, Qo_meas ## n))

// #define BUILD_SE_VEC_MANIFOLD(z1, n1, text1) \
//   MTK_BUILD_MANIFOLD(SE_vec ## n1, ((SO3, bQo2_meas)) BOOST_PP_REPEAT( n1, DECL_SE_VEC_MANIFOLD, 0)); \


// BOOST_PP_REPEAT_FROM_TO(1, 12, BUILD_SE_VEC_MANIFOLD, 0)

// MTK_BUILD_MANIFOLD(measurement_dyn2, ((SE_vec1, sevec1))((SO3, bQo2_meas)));

// measurement_dyn6<3> a__;

int main(int argc, char** argv)
{
  //  define variables
  double dt = 0.01;       // time step
  state_ikfom ekf_state;  // state vector (struct)
  input_ikfom u;          // input vector (struct)
  MTK::get_cov<process_noise_ikfom>::type Q =
      Eigen::Matrix<double, state_dof, state_dof>::Identity();  // process noise covariance: Q, an Eigen matrix
  esekfom::esekf<state_ikfom, process_noise_dof, input_ikfom>::cov init_P;

  // initalize variables
  // ekf_state.bpo << 0.0, 0.0, 0.0;                // initial position
  // ekf_state.bQo = MTK::SO3<double>::Identity();  // initial orientation
  // u.ovo << 0.1, 0.0, 0.0;                        // input velocity in body frame
  // u.oomegao << 0.1, 0.0, 0.0;                    // input angular velocity in body frame
  bool valid = false;

  // create the filter object
  esekfom::esekf<state_ikfom, process_noise_dof, input_ikfom> kf(ekf_state, init_P);

  double epsi[state_dof] = { 0.001 };
  std::fill(epsi, epsi + state_dof,
            0.001);  // if the absolute of innovation of ekf update is smaller than epso, the update iteration is
                     // converged

  kf.init_dyn_runtime_share(get_f, df_dx, df_dw, MAXIMUM_ITER, epsi);

  // setup simulator real process
  // state_ikfom true_state;
  // true_state.bpo << 0.3, 0.2, 0.1;                                   // initial position
  // true_state.bQo = Eigen::Quaternion<double>(0.7071, 0.7071, 0, 0);  // initial orientation
  // double tf = 10;
  // double t = 0.0;  // simulation time

  // while (t < tf)
  // {
  //   // predict
  kf.predict(dt, Q, u);

  //   // update

  // define manifold measure
  MTK::MyManifold<double, 3> z();
  kf.update_iterated_dyn_runtime_share(z, h_dyn_runtime_share<3>);

  MTK::MyManifold<double, 7> z1();
  kf.update_iterated_dyn_runtime_share(z1, h_dyn_runtime_share<7>);
                                       

  //   // // print the state
  //   // state_ikfom cur_state_aft = kf.get_x();

  //   // simulate real system
  //   true_state.boxplus(get_f(true_state, u), dt);
  //   z = get_h(true_state, valid);

  //   std::cout << "Real state: \n" << true_state << std::endl;

    std::cout << "Estimated state: \n" << kf.get_x() << std::endl;

  //   t = t + dt;
  // }

  return 0;
}
