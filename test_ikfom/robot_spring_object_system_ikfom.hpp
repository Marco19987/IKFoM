
#ifndef USE_IKFOM_H
#define USE_IKFOM_H

#include <IKFoM_toolkit/esekfom/esekfom.hpp>
#include "my_manifold.hpp"

typedef MTK::vect<3, double> vect3;
typedef MTK::SO3<double> SO3;

MTK_BUILD_MANIFOLD(state_ikfom, ((vect3, bpo))((SO3, bQo))((vect3, ovo))((vect3, oomegao))((vect3, b1pe1))(
									(SO3, b1Qe1))((vect3, b2pe2))((SO3, b2Qe2))((vect3, b2pb1))((SO3, b2Qb1))(
									(vect3, opg1))((SO3, oQg1))((vect3, opg2))((SO3, oQg2)));

MTK_BUILD_MANIFOLD(input_ikfom, ((vect3, e1ve1))((vect3, e1omegae1))((vect3, e2ve2))((vect3, e2omegae2)));

MTK_BUILD_MANIFOLD(process_noise_ikfom,
				   ((vect3, w_ovo))((vect3, w_oomegao))((vect3, w_linear_acc))((vect3, w_angular_acc))((
					   vect3, w_e1ve1))((vect3, w_e1omegae1))((vect3, w_e2ve2))((vect3, w_e2omegae2))((vect3, w_b2vb1))(
					   (vect3, w_b2omegab1))((vect3, w_ovg1))((vect3, w_oomegag1))((vect3, w_ovg2))((vect3,
																									 w_oomegag2)));

struct ProcessParams
{
  vect3 bg;							 // gravity vector in base frame
  vect3 bpb1;						 // position of robot1 base frame in base frame
  SO3 bQb1;							 // quaternion of robot1 base frame in base frame
  Eigen::Matrix<double, 6, 6> M;	 // inertia matrix
  Eigen::Matrix<double, 6, 6> K1;	 // robot1-object elastic stiffness
  Eigen::Matrix<double, 6, 6> B1;	 // robot1-object damping springs
  Eigen::Matrix<double, 6, 6> K2;	 // robot2-object elastic stiffness
  Eigen::Matrix<double, 6, 6> B2;	 // robot2-object damping stiffness
  Eigen::Matrix<double, 6, 6> Bair;	 // air damping matrix
};

inline ProcessParams& get_process_params()
{
  static ProcessParams params;
  return params;
}

#define state_dof 42
#define process_noise_dof 42
#define measurement_dof(num_pose_measure) (12 + 12 + 6 * (num_pose_measure))	//! zk = h(xk,vk) the dof is not the dim of zk

Eigen::Matrix<double, state_dof, 1> get_f(state_ikfom& s, const input_ikfom& in)
{
  Eigen::Matrix<double, state_dof, 1> res = Eigen::Matrix<double, state_dof, 1>::Identity();
  return res;
}

Eigen::Matrix<double, state_dof, state_dof> df_dx(state_ikfom& s, const input_ikfom& in)
{
  Eigen::Matrix<double, state_dof, state_dof> cov = Eigen::Matrix<double, state_dof, state_dof>::Identity();
  return cov;
}

Eigen::Matrix<double, state_dof, process_noise_dof> df_dw(state_ikfom& s, const input_ikfom& in)
{
  Eigen::Matrix<double, state_dof, process_noise_dof> cov =
	  Eigen::Matrix<double, state_dof, process_noise_dof>::Identity();
  return cov;
}

template <int num_aruco>
MTK::MyManifold<double, num_aruco> get_h(state_ikfom& s, bool& valid)
{
  MTK::MyManifold<double, num_aruco> h_;
  return h_;
}

template <int num_aruco>
Eigen::Matrix<double, measurement_dof(num_aruco), state_dof> dh_dx(state_ikfom& s, bool& valid)
{
  Eigen::Matrix<double, measurement_dof(num_aruco), state_dof> cov =
	  Eigen::Matrix<double, measurement_dof(num_aruco), state_dof>::Identity();

  return cov;
}

template <int num_aruco>
Eigen::Matrix<double, measurement_dof(num_aruco), measurement_dof(num_aruco)> dh_dv(state_ikfom& s, bool& valid)
{
  Eigen::Matrix<double, measurement_dof(num_aruco), measurement_dof(num_aruco)> cov =
	  Eigen::Matrix<double, measurement_dof(num_aruco), measurement_dof(num_aruco)>::Identity();
  return cov;
}

template <int num_aruco>
MTK::MyManifold<double, num_aruco>
h_dyn_runtime_share(state_ikfom& s, esekfom::dyn_runtime_share_datastruct<double>& dyn_runtime_share_data)
{
  if (dyn_runtime_share_data.converge)
  {
  }	 // this value is true means iteration is converged

  dyn_runtime_share_data.h_x =
	  dh_dx<num_aruco>(s, dyn_runtime_share_data.valid);  // H_x is the result matrix of the first differentiation
  dyn_runtime_share_data.h_v =
	  dh_dv<num_aruco>(s, dyn_runtime_share_data.valid);  // H_v is the result matrix of the second differentiation
  dyn_runtime_share_data.R =
	  Eigen::Matrix<double, measurement_dof(num_aruco), measurement_dof(num_aruco)>::Identity();  // R is the
																								  // measurement noise
																								  // covariance

  MTK::MyManifold<double, num_aruco> h_();

  return h_;
}

#endif