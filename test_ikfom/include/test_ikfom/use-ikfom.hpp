
#ifndef USE_IKFOM_H
#define USE_IKFOM_H

#include "../IKFoM_toolkit/esekfom/esekfom.hpp"


typedef MTK::vect<3, double> vect3;
typedef MTK::SO3<double> SO3;


MTK_BUILD_MANIFOLD(state_ikfom,
((vect3, bpo)) 
((SO3, bQo))
);

MTK_BUILD_MANIFOLD(input_ikfom,
((vect3, ovo))
((vect3, oomegao))
);

MTK_BUILD_MANIFOLD(process_noise_ikfom,
((vect3, n_ovo))
((vect3, n_oomegao))
);

MTK_BUILD_MANIFOLD(measurement_ikfom,
((vect3, bpo_meas)) 
((SO3, bQo_meas))
);

MTK_BUILD_MANIFOLD(measurement_noise_ikfom,
((vect3, v_bpo))
((vect3, v_bQo))
);

MTK::get_cov<process_noise_ikfom>::type process_noise_cov()
{
	MTK::get_cov<process_noise_ikfom>::type cov = MTK::get_cov<process_noise_ikfom>::type::Zero();
	MTK::setDiagonal<process_noise_ikfom, vect3, 0>(cov, &process_noise_ikfom::n_ovo, 0.0001);// 0.03
	MTK::setDiagonal<process_noise_ikfom, vect3, 3>(cov, &process_noise_ikfom::n_oomegao, 0.0001); // *dt 0.01 0.01 * dt * dt 0.05
	
	return cov;
}

MTK::get_cov<measurement_noise_ikfom>::type measurement_noise_cov()
{
	MTK::get_cov<measurement_noise_ikfom>::type cov = MTK::get_cov<measurement_noise_ikfom>::type::Zero();
	MTK::setDiagonal<measurement_noise_ikfom, vect3, 0>(cov, &measurement_noise_ikfom::v_bpo, 0.0001);// 0.03
	MTK::setDiagonal<measurement_noise_ikfom, vect3, 3>(cov, &measurement_noise_ikfom::v_bQo, 0.0001); // *dt 0.01 0.01 * dt * dt 0.05
	
	return cov;
}

#define state_dof 6
#define process_noise_dof 6
#define measurement_dof 6 //! zk = h(xk,vk) the dof is not the dim of zk
#define measurement_noise_dof 6


Eigen::Matrix<double, state_dof, 1> get_f(state_ikfom &s, const input_ikfom &in)
{
	Eigen::Matrix<double, state_dof, 1> res = Eigen::Matrix<double, state_dof, 1>::Zero();
	res.template block<3,1>(0,0) = (s.bQo * Eigen::Quaternion<double>(0, in.ovo[0], in.ovo[1], in.ovo[2]) * s.bQo.conjugate()).vec();	
	res.template block<3,1>(3,0) = in.oomegao;
	return res;
}

Eigen::Matrix<double, state_dof, state_dof> df_dx(state_ikfom &s, const input_ikfom &in)
{
	Eigen::Matrix<double, state_dof, state_dof> cov = Eigen::Matrix<double, state_dof, state_dof>::Identity();
	// cov.template block<3, 3>(0, 6) = Eigen::Matrix3d::Identity();
	// vect3 acc_;
	// in.acc.boxminus(acc_, s.ba);
	// vect3 omega;
	// in.gyro.boxminus(omega, s.bg);
	// cov.template block<3, 3>(6, 3) = -s.rot.toRotationMatrix()*MTK::hat(acc_);
	// cov.template block<3, 3>(6, 12) = -s.rot.toRotationMatrix();
	// Eigen::Matrix<state_ikfom::scalar, 2, 1> vec = Eigen::Matrix<state_ikfom::scalar, 2, 1>::Zero();
	// Eigen::Matrix<state_ikfom::scalar, 3, 2> grav_matrix;
	// s.S2_Mx(grav_matrix, vec, 15);
	// cov.template block<3, 2>(6, 15) =  grav_matrix; 
	// cov.template block<3, 3>(3, 9) = -Eigen::Matrix3d::Identity(); 
	return cov;
}


Eigen::Matrix<double, state_dof, process_noise_dof> df_dw(state_ikfom &s, const input_ikfom &in)
{
	Eigen::Matrix<double, state_dof, process_noise_dof> cov = Eigen::Matrix<double, state_dof, process_noise_dof>::Identity();
	// cov.template block<3, 3>(6, 3) = -s.rot.toRotationMatrix();
	// cov.template block<3, 3>(3, 0) = -Eigen::Matrix3d::Identity();
	// cov.template block<3, 3>(9, 6) = Eigen::Matrix3d::Identity();
	// cov.template block<3, 3>(12, 9) = Eigen::Matrix3d::Identity();
	return cov;
}

measurement_ikfom get_h(state_ikfom &s, bool &valid)
{
	measurement_ikfom h_;
	h_.bpo_meas = s.bpo;	
	h_.bQo_meas = s.bQo;
	return h_;
}

Eigen::Matrix<double, measurement_dof, state_dof> dh_dx(state_ikfom &s, bool &valid)
{
	Eigen::Matrix<double, measurement_dof, state_dof> cov = Eigen::Matrix<double, measurement_dof, state_dof>::Identity();
	
	return cov;
}

Eigen::Matrix<double, measurement_dof, measurement_noise_dof> dh_dv(state_ikfom &s, bool &valid)
{
	Eigen::Matrix<double, measurement_dof, measurement_noise_dof> cov = Eigen::Matrix<double, measurement_dof, measurement_noise_dof>::Identity();
	return cov;
}



#endif