
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
  vect3 bg;                          // gravity vector in base frame
  vect3 bpb1;                        // position of robot1 base frame in base frame
  SO3 bQb1;                          // quaternion of robot1 base frame in base frame
  Eigen::Matrix<double, 6, 6> M;     // inertia matrix
  Eigen::Matrix<double, 6, 6> K1;    // robot1-object elastic stiffness
  Eigen::Matrix<double, 6, 6> B1;    // robot1-object damping springs
  Eigen::Matrix<double, 6, 6> K2;    // robot2-object elastic stiffness
  Eigen::Matrix<double, 6, 6> B2;    // robot2-object damping stiffness
  Eigen::Matrix<double, 6, 6> Bair;  // air damping matrix
  int num_pose_measure_robot1;       //!   (num_pose_measure_robot1 + num_pose_measure_robot2) = num_pose_measure
  int num_pose_measure_robot2;  // num measure aruco from robot-i range from 0 to N (with N number arucos on the object)
};

#define state_dof 42
#define process_noise_dof 42
#define measurement_dof(num_pose_measure)                                                                              \
  (12 + 12 + 6 * (num_pose_measure))  //! zk = h(xk,vk) the dof is not the dim of zk

// Helper functions
inline void printProcessParams(const ProcessParams& process_params)
{
  std::cout << "Process Params:" << std::endl;
  std::cout << "bg: " << process_params.bg.transpose() << std::endl;
  std::cout << "bpb1: " << process_params.bpb1.transpose() << std::endl;
  std::cout << "bQb1: " << process_params.bQb1.coeffs().transpose() << std::endl;
  std::cout << "M: " << std::endl << process_params.M << std::endl;
  std::cout << "K1: " << std::endl << process_params.K1 << std::endl;
  std::cout << "B1: " << std::endl << process_params.B1 << std::endl;
  std::cout << "K2: " << std::endl << process_params.K2 << std::endl;
  std::cout << "B2: " << std::endl << process_params.B2 << std::endl;
  std::cout << "Bair: " << std::endl << process_params.Bair << std::endl;
}

inline void skew_matrix(const Eigen::Matrix<double, 3, 1>& v, Eigen::Matrix<double, 3, 3>& skew_mat)
{
  skew_mat << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
}

inline void get_grasp_matrix(const state_ikfom& s, Eigen::Matrix<double, 6, 12>& W_)
{
  Eigen::Matrix<double, 6, 6> Wg1, Wg2;
  Eigen::Matrix3d skew_opg1, skew_opg2;

  skew_matrix(s.opg1, skew_opg1);
  skew_matrix(s.opg2, skew_opg2);

  Wg1 << Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero(), skew_opg1, Eigen::Matrix3d::Identity();
  Wg2 << Eigen::Matrix3d::Identity(), Eigen::Matrix3d::Zero(), skew_opg2, Eigen::Matrix3d::Identity();

  W_.block<6, 6>(0, 0) = Wg1;
  W_.block<6, 6>(0, 6) = Wg2;
}

inline void get_Rbar(const state_ikfom& s, Eigen::Matrix<double, 12, 12>& Rbar_)
{
  Eigen::Matrix3d oRg1 = s.oQg1.toRotationMatrix();
  Eigen::Matrix3d oRg2 = s.oQg2.toRotationMatrix();
  Rbar_.setZero();
  Rbar_.block<3, 3>(0, 0) = oRg1;
  Rbar_.block<3, 3>(3, 3) = oRg1;
  Rbar_.block<3, 3>(6, 6) = oRg2;
  Rbar_.block<3, 3>(9, 9) = oRg2;
}

inline void pose_to_matrix(const vect3& position, const SO3& quaternion, Eigen::Matrix<double, 4, 4>& T)
{
  T = Eigen::Matrix<double, 4, 4>::Identity();  // Ensure T is initialized as a 4x4 matrix
  T.template block<3, 3>(0, 0) = quaternion.toRotationMatrix();
  T.template block<3, 1>(0, 3) = position;
}

Eigen::Matrix<double, 3, 1> angle_axis_error(Eigen::Matrix<double, 3, 3> R)
{
  // error = r*sin(theta); r = 0.5 * (1/sin(theta)) * (axis_versor)
  Eigen::Matrix<double, 3, 1> axis_versor;
  axis_versor << R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1);
  axis_versor = axis_versor / 2.0;
  return axis_versor;
}

inline void spring_model(state_ikfom& s, ProcessParams& process_params, Eigen::Matrix<double, 12, 1>& out,
                         bool compute_viscous, bool rotate_in_grasp_frame, const input_ikfom& in = input_ikfom())
{
  Eigen::Matrix<double, 4, 4> bTo = Eigen::Matrix<double, 4, 4>::Identity();
  pose_to_matrix(s.bpo, s.bQo, bTo);
  Eigen::Matrix<double, 4, 4> b1Te1 = Eigen::Matrix<double, 4, 4>::Identity();
  pose_to_matrix(s.b1pe1, s.b1Qe1, b1Te1);
  Eigen::Matrix<double, 4, 4> b2Te2 = Eigen::Matrix<double, 4, 4>::Identity();
  pose_to_matrix(s.b2pe2, s.b2Qe2, b2Te2);
  Eigen::Matrix<double, 4, 4> b2Tb1 = Eigen::Matrix<double, 4, 4>::Identity();
  pose_to_matrix(s.b2pb1, s.b2Qb1, b2Tb1);
  Eigen::Matrix<double, 4, 4> oTg1 = Eigen::Matrix<double, 4, 4>::Identity();
  pose_to_matrix(s.opg1, s.oQg1, oTg1);
  Eigen::Matrix<double, 4, 4> oTg2 = Eigen::Matrix<double, 4, 4>::Identity();
  pose_to_matrix(s.opg2, s.oQg2, oTg2);
  Eigen::Matrix<double, 4, 4> bTb1 = Eigen::Matrix<double, 4, 4>::Identity();
  pose_to_matrix(process_params.bpb1, process_params.bQb1, bTb1);

  Eigen::Matrix<double, 4, 4> g1Te1 = oTg1.inverse() * bTo.inverse() * bTb1 * b1Te1;
  Eigen::Matrix<double, 4, 4> g2Te2 = oTg2.inverse() * bTo.inverse() * bTb1 * b2Tb1.inverse() * b2Te2;

  Eigen::Matrix<double, 3, 3> g1Re1 = g1Te1.template block<3, 3>(0, 0);
  Eigen::Matrix<double, 3, 3> g2Re2 = g2Te2.template block<3, 3>(0, 0);
  Eigen::Matrix<double, 3, 3> bRo = bTo.template block<3, 3>(0, 0);

  // Compute wrench spring 1
  Eigen::Matrix<double, 3, 3> oRg1 = oTg1.template block<3, 3>(0, 0);
  Eigen::Matrix<double, 3, 3> g1Rb = oRg1.transpose() * bRo.transpose();
  Eigen::Matrix<double, 3, 1> bpe1_b = bTb1.template block<3, 1>(0, 3) + bTb1.template block<3, 3>(0, 0) * s.b1pe1;
  Eigen::Matrix<double, 3, 1> bpg1_b = s.bpo + bRo * oTg1.template block<3, 1>(0, 3);
  Eigen::Matrix<double, 3, 1> g1fg1_e = process_params.K1.template block<3, 3>(0, 0) * g1Rb * (bpe1_b - bpg1_b);

  Eigen::Matrix<double, 3, 3> e1Rb =
      s.b1Qe1.toRotationMatrix().transpose() * bTb1.template block<3, 3>(0, 0).transpose();
  Eigen::Matrix<double, 3, 1> e1fe1_e = -process_params.K1.template block<3, 3>(0, 0) * e1Rb * (bpe1_b - bpg1_b);

  Eigen::Matrix<double, 3, 1> bpe1_dot = bTb1.template block<3, 3>(0, 0) * s.b1Qe1.toRotationMatrix() * in.e1ve1;
  Eigen::Matrix<double, 3, 3> skew_b_omega_o;
  skew_matrix(bRo * s.oomegao, skew_b_omega_o);
  Eigen::Matrix<double, 3, 1> bpg1_dot = bRo * s.ovo + skew_b_omega_o * bRo * oTg1.template block<3, 1>(0, 3);
  Eigen::Matrix<double, 3, 1> e1fe1_beta = -process_params.B1.template block<3, 3>(0, 0) * e1Rb * (bpe1_dot - bpg1_dot);

  Eigen::Matrix<double, 3, 1> e1_tau_e1_e =
      process_params.K1.template block<3, 3>(3, 3) * angle_axis_error(g1Re1.transpose());

  Eigen::Matrix<double, 3, 1> e1_tau_e1_beta =
      -process_params.B1.template block<3, 3>(3, 3) * e1Rb *
      (bTb1.template block<3, 3>(0, 0) * s.b1Qe1.toRotationMatrix() * in.e1omegae1 - bRo * s.oomegao);

  Eigen::Matrix<double, 6, 1> e1he1;

  if (compute_viscous)
  {
    e1he1 << e1fe1_e + e1fe1_beta, e1_tau_e1_e + e1_tau_e1_beta;
  }
  else
  {
    e1he1 << e1fe1_e, e1_tau_e1_e;
  }

  // Compute wrench spring 2
  Eigen::Matrix<double, 3, 3> oRg2 = oTg2.template block<3, 3>(0, 0);
  Eigen::Matrix<double, 3, 3> g2Rb = oRg2.transpose() * bRo.transpose();
  Eigen::Matrix<double, 4, 4> b1Tb2 = b2Tb1.inverse();
  Eigen::Matrix<double, 3, 1> bpe2_b =
      bTb1.template block<3, 1>(0, 3) +
      bTb1.template block<3, 3>(0, 0) * (b1Tb2.template block<3, 1>(0, 3) + b1Tb2.template block<3, 3>(0, 0) * s.b2pe2);
  Eigen::Matrix<double, 3, 1> bpg2_b = s.bpo + bRo * oTg2.template block<3, 1>(0, 3);
  Eigen::Matrix<double, 3, 3> e2Rb = s.b2Qe2.toRotationMatrix().transpose() *
                                     b1Tb2.template block<3, 3>(0, 0).transpose() *
                                     bTb1.template block<3, 3>(0, 0).transpose();
  Eigen::Matrix<double, 3, 1> e2fe2_e = -process_params.K2.template block<3, 3>(0, 0) * e2Rb * (bpe2_b - bpg2_b);

  Eigen::Matrix<double, 3, 1> bpe2_dot =
      bTb1.template block<3, 3>(0, 0) * b1Tb2.template block<3, 3>(0, 0) * s.b2Qe2.toRotationMatrix() * in.e2ve2;
  Eigen::Matrix<double, 3, 1> bpg2_dot = bRo * s.ovo + skew_b_omega_o * bRo * oTg2.template block<3, 1>(0, 3);
  Eigen::Matrix<double, 3, 1> g2fg2_beta = process_params.B2.template block<3, 3>(0, 0) * g2Rb * (bpe2_dot - bpg2_dot);
  Eigen::Matrix<double, 3, 1> e2fe2_beta = -process_params.B2.template block<3, 3>(0, 0) * e2Rb * (bpe2_dot - bpg2_dot);

  Eigen::Matrix<double, 3, 1> e2_tau_e2_e =
      process_params.K2.template block<3, 3>(3, 3) * angle_axis_error(g2Re2.transpose());

  Eigen::Matrix<double, 3, 1> e2_tau_e2_beta =
      -process_params.B2.template block<3, 3>(3, 3) * e2Rb *
      (bTb1.template block<3, 3>(0, 0) * b1Tb2.template block<3, 3>(0, 0) * s.b2Qe2.toRotationMatrix() * in.e2omegae2 -
       bRo * s.oomegao);

  Eigen::Matrix<double, 6, 1> e2he2;

  if (compute_viscous)
  {
    e2he2 << e2fe2_e + e2fe2_beta, e2_tau_e2_e + e2_tau_e2_beta;
  }
  else
  {
    e2he2 << e2fe2_e, e2_tau_e2_e;
  }

  out << e1he1, e2he2;

  if (rotate_in_grasp_frame)
  {
    Eigen::Matrix<double, 12, 12> Rbar_grasp = Eigen::Matrix<double, 12, 12>::Zero();
    Rbar_grasp.template block<3, 3>(0, 0) = g1Re1;
    Rbar_grasp.template block<3, 3>(3, 3) = g1Re1;
    Rbar_grasp.template block<3, 3>(6, 6) = g2Re2;
    Rbar_grasp.template block<3, 3>(9, 9) = g2Re2;

    Eigen::Vector<double, 3> g1pe1 = g1Te1.template block<3, 1>(0, 3);
    Eigen::Vector<double, 3> g2pe2 = g2Te2.template block<3, 1>(0, 3);
    Eigen::Matrix<double, 6, 6> Wg1_e1;
    Eigen::Matrix<double, 6, 6> Wg2_e2;

    Eigen::Matrix<double, 3, 3> skew_g1pe1, skew_g2pe2;
    skew_matrix(g1pe1, skew_g1pe1);
    skew_matrix(g2pe2, skew_g2pe2);

    Wg1_e1 << Eigen::Matrix<double, 3, 3>::Identity(), Eigen::Matrix<double, 3, 3>::Zero(), skew_g1pe1,
        Eigen::Matrix<double, 3, 3>::Identity();
    Wg2_e2 << Eigen::Matrix<double, 3, 3>::Identity(), Eigen::Matrix<double, 3, 3>::Zero(), skew_g2pe2,
        Eigen::Matrix<double, 3, 3>::Identity();

    Eigen::Matrix<double, 12, 12> W_grasp = Eigen::Matrix<double, 12, 12>::Zero();
    W_grasp.template block<6, 6>(0, 0) = Wg1_e1;
    W_grasp.template block<6, 6>(6, 6) = Wg2_e2;

    out = -W_grasp * Rbar_grasp * out;
  }
}

Eigen::Matrix<double, state_dof, 1> get_f(state_ikfom& s, const input_ikfom& in, ProcessParams& process_params)
{
  Eigen::Matrix<double, state_dof, 1> res = Eigen::Matrix<double, state_dof, 1>::Zero();

  Eigen::Matrix<double, 6, 6> bRo_bar = Eigen::Matrix<double, 6, 6>::Zero();
  bRo_bar.template block<3, 3>(0, 0) = s.bQo.toRotationMatrix();
  bRo_bar.template block<3, 3>(3, 3) = s.bQo.toRotationMatrix();

  Eigen::Matrix<double, 12, 1> ghg = Eigen::Matrix<double, 12, 1>::Zero();  // spring model output in the grasp frame
  spring_model(s, process_params, ghg, true, true, in);                     // spring model

  Eigen::Matrix<double, 6, 12> oWg = Eigen::Matrix<double, 6, 12>::Zero();
  get_grasp_matrix(s, oWg);
  Eigen::Matrix<double, 12, 12> oRg_bar = Eigen::Matrix<double, 12, 12>::Identity();
  get_Rbar(s, oRg_bar);

  Eigen::Matrix<double, 6, 1> oh = oWg * oRg_bar * ghg;  // resulting wrench in the object frame

  Eigen::Matrix<double, 6, 1> o_twist_o;
  o_twist_o << s.ovo, s.oomegao;  // object twist in the object frame

  Eigen::Matrix<double, 3, 3> skew_oomegao;
  Eigen::Matrix<double, 3, 3> skew_ovo;
  skew_matrix(s.oomegao, skew_oomegao);
  skew_matrix(s.ovo, skew_ovo);
  Eigen::Matrix<double, 6, 6> double_skew = Eigen::Matrix<double, 6, 6>::Zero();
  double_skew.template block<3, 3>(0, 0) = skew_oomegao;
  double_skew.template block<3, 3>(3, 0) = skew_ovo;
  double_skew.template block<3, 3>(3, 3) = skew_oomegao;

  Eigen::Matrix<double, 6, 1> bg_ext = Eigen::Matrix<double, 6, 1>::Zero();
  bg_ext.template block<3, 1>(0, 0) = process_params.bg;  // gravity vector in base frame

  Eigen::Matrix<double, 6, 1> o_twist_dot_o =
      process_params.M.inverse() * (oh - process_params.Bair * o_twist_o - double_skew * process_params.M * o_twist_o) +
      bRo_bar.transpose() * bg_ext;

  res.template block<3, 1>(0, 0) = s.bQo.toRotationMatrix() * s.ovo;
  res.template block<3, 1>(3, 0) = s.oomegao;
  res.template block<6, 1>(6, 0) = o_twist_dot_o;
  res.template block<3, 1>(12, 0) = s.b1Qe1.toRotationMatrix() * in.e1ve1;  // b1pe1_dot
  res.template block<3, 1>(15, 0) = in.e1omegae1;
  res.template block<3, 1>(18, 0) = s.b2Qe2.toRotationMatrix() * in.e2ve2;  // b2pe2_dot
  res.template block<3, 1>(21, 0) = in.e2omegae2;
  std::cout << "res: " << res.transpose() << std::endl;

  return res;
}

Eigen::Matrix<double, state_dof, state_dof> df_dx(state_ikfom& s, const input_ikfom& in, ProcessParams& process_params)
{
  Eigen::Matrix<double, state_dof, state_dof> cov = Eigen::Matrix<double, state_dof, state_dof>::Identity();
  return cov;
}

Eigen::Matrix<double, state_dof, process_noise_dof> df_dw(state_ikfom& s, const input_ikfom& in,
                                                          ProcessParams& process_params)
{
  Eigen::Matrix<double, state_dof, process_noise_dof> cov =
      Eigen::Matrix<double, state_dof, process_noise_dof>::Identity();
  return cov;
}

template <int num_aruco>
MTK::MyManifold<double, num_aruco> get_h(state_ikfom& s, ProcessParams& process_params)
{
  MTK::MyManifold<double, num_aruco> h_(wrench_measure(), wrench_measure(), pose_measure(), pose_measure(),
                                        std::vector<pose_measure>(num_aruco, pose_measure()));

  Eigen::Matrix<double, 4, 4> bTo = Eigen::Matrix<double, 4, 4>::Identity();
  pose_to_matrix(s.bpo, s.bQo, bTo);

  Eigen::Matrix<double, 4, 4> b2Tb1 = Eigen::Matrix<double, 4, 4>::Identity();
  pose_to_matrix(s.b2pb1, s.b2Qb1, b2Tb1);

  Eigen::Matrix<double, 4, 4> bTb1 = Eigen::Matrix<double, 4, 4>::Identity();
  pose_to_matrix(process_params.bpb1, process_params.bQb1, bTb1);

  Eigen::Matrix<double, 3, 1> b1po;
  Eigen::Matrix<double, 3, 1> b2po;

  // measured object position in the robots base frame b1 and b2
  b1po = bTb1.template block<3, 3>(0, 0).transpose() * (s.bpo - process_params.bpb1);
  b2po = b2Tb1.template block<3, 1>(0, 3) + b2Tb1.template block<3, 3>(0, 0) * b1po;

  // measured object orientation in the robots base frame b1 and b2
  Eigen::Quaternion<double> b1Qo = process_params.bQb1.inverse() * s.bQo;
  Eigen::Quaternion<double> b2Qo = s.b2Qb1 * b1Qo;

  // out is [b1po;b1Qo;b1po;b1Qo; ... ;b2po;b2Qo; .. ; b2po;b2Qo]
  for (int i = 0; i < process_params.num_pose_measure_robot1; i++)
  {
    h_.object_pose[i].position = b1po;
    h_.object_pose[i].quaternion = b1Qo;
  }
  for (int i = process_params.num_pose_measure_robot1; i < process_params.num_pose_measure_robot2; i++)
  {
    h_.object_pose[i].position = b2po;
    h_.object_pose[i].quaternion = b2Qo;
  }
  // output of the wrenches exerted by the robots on the object
  Eigen::Matrix<double, 12, 1> ehe;  // wrench in the end effector frame applied by the object on the robots
  spring_model(s, process_params, ehe, false, false);
  h_.e1he1_e1.force = ehe.template block<3, 1>(0, 0);
  h_.e1he1_e1.torque = ehe.template block<3, 1>(3, 0);
  h_.e2he2_e2.force = ehe.template block<3, 1>(6, 0);
  h_.e2he2_e2.torque = ehe.template block<3, 1>(9, 0);

  // output of the fkine of the robots
  h_.b1_pose_e1.position = s.b1pe1;    // position of robot1 end effector in the base frame
  h_.b1_pose_e1.quaternion = s.b1Qe1;  // orientation of robot1 end effector in the base frame
  h_.b2_pose_e2.position = s.b2pe2;    // position of robot2 end effector in the base frame
  h_.b2_pose_e2.quaternion = s.b2Qe2;  // orientation of robot2 end effector in the base frame

  std::cout << "h_.b1_pose_e1.position: " << h_.b1_pose_e1.position.transpose() << std::endl;
  std::cout << "h_.b1_pose_e1.orientation: " << h_.b1_pose_e1.quaternion.coeffs().transpose() << std::endl;
  std::cout << "h_.b2_pose_e2.position: " << h_.b2_pose_e2.position.transpose() << std::endl;
  std::cout << "h_.b2_pose_e2.orientation: " << h_.b2_pose_e2.quaternion.coeffs().transpose() << std::endl;
  std::cout << "h_.e1he1_e1 force: " << h_.e1he1_e1.force.transpose() << std::endl;
  std::cout << "h_.e1he1_e1 torque: " << h_.e1he1_e1.torque.transpose() << std::endl;
  std::cout << "h_.e2he2_e2 force: " << h_.e2he2_e2.force.transpose() << std::endl;
  std::cout << "h_.e2he2_e2 torque: " << h_.e2he2_e2.torque.transpose() << std::endl;

  std::cout << "h_.object_pose.size(): " << h_.object_pose.size() << std::endl;
  for (const auto& pose : h_.object_pose)
  {
    std::cout << "pose.position: " << pose.position.transpose() << std::endl;
    std::cout << "pose.orientation: " << pose.quaternion.coeffs().transpose() << std::endl;
  }

  return h_;
}

template <int num_aruco>
Eigen::Matrix<double, measurement_dof(num_aruco), state_dof> dh_dx(state_ikfom& s, bool& valid,
                                                                   ProcessParams& process_params)
{
  Eigen::Matrix<double, measurement_dof(num_aruco), state_dof> cov =
      Eigen::Matrix<double, measurement_dof(num_aruco), state_dof>::Identity();

  return cov;
}

template <int num_aruco>
Eigen::Matrix<double, measurement_dof(num_aruco), measurement_dof(num_aruco)> dh_dv(state_ikfom& s, bool& valid,
                                                                                    ProcessParams& process_params)
{
  Eigen::Matrix<double, measurement_dof(num_aruco), measurement_dof(num_aruco)> cov =
      Eigen::Matrix<double, measurement_dof(num_aruco), measurement_dof(num_aruco)>::Identity();
  return cov;
}

template <int num_aruco>
MTK::MyManifold<double, num_aruco>
h_dyn_runtime_share(state_ikfom& s, esekfom::dyn_runtime_share_datastruct<double>& dyn_runtime_share_data,
                    ProcessParams& process_params)
{
  if (process_params.num_pose_measure_robot1 + process_params.num_pose_measure_robot2 != num_aruco)
  {
    throw std::runtime_error("Error: num_pose_measure_robot1 + num_pose_measure_robot2 != num_aruco");
  }

  if (dyn_runtime_share_data.converge)
  {
  }  // this value is true means iteration is converged

  dyn_runtime_share_data.h_x = dh_dx<num_aruco>(
      s, dyn_runtime_share_data.valid, process_params);  // H_x is the result matrix of the first differentiation
  dyn_runtime_share_data.h_v = dh_dv<num_aruco>(
      s, dyn_runtime_share_data.valid, process_params);  // H_v is the result matrix of the second differentiation
  dyn_runtime_share_data.R =
      Eigen::Matrix<double, measurement_dof(num_aruco), measurement_dof(num_aruco)>::Identity();  // R is the
                                                                                                  // measurement noise
                                                                                                  // covariance

  MTK::MyManifold<double, num_aruco> h_(wrench_measure(), wrench_measure(), pose_measure(), pose_measure(),
                                        std::vector<pose_measure>(num_aruco, pose_measure()));

  return get_h<num_aruco>(s, process_params);
}

#endif