
#ifndef USE_IKFOM_H
#define USE_IKFOM_H

#include "../include/IKFoM_toolkit/esekfom/esekfom.hpp"

typedef MTK::vect<3, double> vect3;
typedef MTK::SO3<double> SO3;

MTK_BUILD_MANIFOLD(state_ikfom,
				   ((vect3, bpo))((SO3, bQo)));

MTK_BUILD_MANIFOLD(input_ikfom,
				   ((vect3, ovo))((vect3, oomegao)));

MTK_BUILD_MANIFOLD(process_noise_ikfom,
				   ((vect3, n_ovo))((vect3, n_oomegao)));

MTK_BUILD_MANIFOLD(measurement_ikfom,
				   ((vect3, bpo_meas))((SO3, bQo_meas)));

MTK_BUILD_MANIFOLD(measurement_noise_ikfom,
				   ((vect3, v_bpo))((vect3, v_bQo)));

MTK::get_cov<process_noise_ikfom>::type process_noise_cov()
{
	MTK::get_cov<process_noise_ikfom>::type cov = MTK::get_cov<process_noise_ikfom>::type::Zero();
	MTK::setDiagonal<process_noise_ikfom, vect3, 0>(cov, &process_noise_ikfom::n_ovo, 0.0001);	   // 0.03
	MTK::setDiagonal<process_noise_ikfom, vect3, 3>(cov, &process_noise_ikfom::n_oomegao, 0.0001); // *dt 0.01 0.01 * dt * dt 0.05

	return cov;
}

MTK::get_cov<measurement_noise_ikfom>::type measurement_noise_cov()
{
	MTK::get_cov<measurement_noise_ikfom>::type cov = MTK::get_cov<measurement_noise_ikfom>::type::Zero();
	MTK::setDiagonal<measurement_noise_ikfom, vect3, 0>(cov, &measurement_noise_ikfom::v_bpo, 0.0001); // 0.03
	MTK::setDiagonal<measurement_noise_ikfom, vect3, 3>(cov, &measurement_noise_ikfom::v_bQo, 0.0001); // *dt 0.01 0.01 * dt * dt 0.05

	return cov;
}

#define state_dof 6
#define process_noise_dof 6
#define measurement_dof 6 //! zk = h(xk,vk) the dof is not the dim of zk
#define measurement_noise_dof 6

// Suggestion: Use static/global variables or a singleton class to hold external parameters.
// Example using a static struct for parameters:

struct ExternalParams {
	Eigen::Vector3d bias{0.0, 0.5, 0.2}; // Example parameter
	double scale = 1.0;
};

inline ExternalParams& get_external_params() {
	static ExternalParams params;
	return params;
}

// std::bind(f,_,_,std::ref(param)) , const &params
Eigen::Matrix<double, state_dof, 1> get_f(state_ikfom &s, const input_ikfom &in)
{
	const ExternalParams& params = get_external_params(); // a way to pass extra arguments to the function
	Eigen::Matrix<double, state_dof, 1> res = Eigen::Matrix<double, state_dof, 1>::Zero();
	// Example usage of external parameters:
	Eigen::Vector3d adjusted_ovo = in.ovo;
	res.template block<3, 1>(0, 0) = (s.bQo * Eigen::Quaternion<double>(0, adjusted_ovo[0], adjusted_ovo[1], adjusted_ovo[2]) * s.bQo.conjugate()).vec();
	res.template block<3, 1>(3, 0) = in.oomegao;
	return res;
}

Eigen::Matrix<double, state_dof, state_dof> df_dx(state_ikfom &s, const input_ikfom &in)
{
	Eigen::Matrix<double, state_dof, state_dof> cov = Eigen::Matrix<double, state_dof, state_dof>::Zero();
	cov.template block<3, 3>(0, 3) = -s.bQo.toRotationMatrix() * MTK::hat(in.ovo);
	return cov;
}

Eigen::Matrix<double, state_dof, process_noise_dof> df_dw(state_ikfom &s, const input_ikfom &in)
{
	Eigen::Matrix<double, state_dof, process_noise_dof> cov = Eigen::Matrix<double, state_dof, process_noise_dof>::Zero();
	cov.template block<3, 3>(0, 0) = s.bQo.toRotationMatrix();
	cov.template block<3, 3>(3, 3) = Eigen::Matrix3d::Identity();
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



// template<typename measurement_runtime, typename measurementModel_dyn_runtime_share>
// measurement h_dyn_runtime_share(state &s, esekfom::dyn_runtime_share_datastruct<double> &dyn_runtime_share_data) 
// {
// 	if(dyn_runtime_share_data.converge) {} // this value is true means iteration is converged 
// 	if(condition) dyn_runtime_share_data.valid = false; // the iteration stops before convergence when this value is false, if conditions other than convergence is satisfied
// 	dyn_runtime_share_data.h_x = H_x; // H_x is the result matrix of the first differentiation 
// 	dyn_runtime_share_data.h_v = H_v; // H_v is the result matrix of the second differentiation
// 	dyn_runtime_share_data.R = R; // R is the measurement noise covariance
	
// 	measurement h_;
// 	h_.pos = s.pos;
// 	return h_;
// }

#endif