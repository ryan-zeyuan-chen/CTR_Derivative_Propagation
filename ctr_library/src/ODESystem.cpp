// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "ODESystem.hpp"

// default constructor
ODESystem::ODESystem()
{
	m_u_ast_x = 0.00;
	m_u_ast_y = 0.00;
	m_EI = 0.00;
	m_GJ = 0.00;
	m_e3 = 0.00;
	m_f = 0.00;
	m_V = 0.00;
	m_E = 0.00;
	m_Vs = 0.00;
	m_Es = 0.00;
	m_e3 = {0.00, 0.00, 1.00};
}

// copy constructor
ODESystem::ODESystem(const ODESystem &rhs) : m_u_ast_x(rhs.m_u_ast_x), m_u_ast_y(rhs.m_u_ast_y), m_EI(rhs.m_EI), m_GJ(rhs.m_GJ), m_e3(rhs.m_e3), m_f(rhs.m_f), m_V(rhs.m_V), m_E(rhs.m_E), m_Vs(rhs.m_Vs), m_Es(rhs.m_Es) {}

// move constructor
ODESystem::ODESystem(ODESystem &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_u_ast_x = std::move(rhs.m_u_ast_x);
		this->m_u_ast_y = std::move(rhs.m_u_ast_y);
		this->m_EI = std::move(rhs.m_EI);
		this->m_GJ = std::move(rhs.m_GJ);
		this->m_e3 = std::move(rhs.m_e3);
		this->m_f = std::move(rhs.m_f);
		this->m_V = std::move(rhs.m_V);
		this->m_E = std::move(rhs.m_E);
		this->m_Vs = std::move(rhs.m_Vs);
		this->m_Es = std::move(rhs.m_Es);
	}
}

// ODESystem destructor
ODESystem::~ODESystem()
{
	// nothing to be done
}

// Copy assignment operator
ODESystem &ODESystem::operator=(const ODESystem &rhs)
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_u_ast_x = rhs.m_u_ast_x;
		this->m_u_ast_y = rhs.m_u_ast_y;
		this->m_EI = rhs.m_EI;
		this->m_GJ = rhs.m_GJ;
		this->m_e3 = rhs.m_e3;
		this->m_f = rhs.m_f;
		this->m_V = rhs.m_V;
		this->m_E = rhs.m_E;
		this->m_Vs = rhs.m_Vs;
		this->m_Es = rhs.m_Es;
	}
	return *this;
}

// move assignment operator
ODESystem &ODESystem::operator=(ODESystem &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_u_ast_x = std::move(rhs.m_u_ast_x);
		this->m_u_ast_y = std::move(rhs.m_u_ast_y);
		this->m_EI = std::move(rhs.m_EI);
		this->m_GJ = std::move(rhs.m_GJ);
		this->m_e3 = std::move(rhs.m_e3);
		this->m_f = std::move(rhs.m_f);
		this->m_V = std::move(rhs.m_V);
		this->m_E = std::move(rhs.m_E);
		this->m_Vs = std::move(rhs.m_Vs);
		this->m_Es = std::move(rhs.m_Es);
	}
	return *this;
}

// functor that implements the system of ODEs governing a three-tube CTR
void ODESystem::operator()(const state_type &y, state_type &dyds, const double s)
{
	// 1st element of y computes the curvature of the first (innermost) tube along the x direction
	// 2nd element of y computes the curvature of the first (innermost) tube along the y direction
	// next 3 elements of y are the torsional curvatures for the three tubes, e.g., y = [u1_z  u2_z  u3_z]
	// next 3 elements of y are twist angles, theta_i = [theta_1  theta_2  theta_3]
	// last 7 elements are r(position) and h(quaternion-orientations) of the local frame, respectively at each arc-length s

	double dtheta_2 = y[3UL] - y[2UL];
	double dtheta_3 = y[4UL] - y[2UL];

	blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> R1(mathOp::rotz(y[5UL])), R2(mathOp::rotz(y[6UL])), R3(mathOp::rotz(y[7UL]));

	// implementing curvature equation u_i = transpose(R_z(theta_i))*u_1 + \dot{theta_i}*e3
	blaze::StaticVector<double, 3UL> u1, u2, u3;

	u1 = blaze::subvector<0UL, 3UL>(y);
	u2 = blaze::trans(R2) * u1 + (dtheta_2 * m_e3);
	u3 = blaze::trans(R3) * u1 + (dtheta_3 * m_e3);

	// computing the twist curvatures (uz_i) and twist angles (theta_i)
	auto computeTwists = [&](size_t idx, const blaze::StaticVector<double, 3UL> &u)
	{
		if (m_GJ[idx] != 0.00)
		{
			// uz_i = ( (E_i * I_i) / (G_i * J_i) ) * (ux_i * uy_ast - uy_i * ux_ast)
			dyds[2 + idx] = (m_EI[idx] / m_GJ[idx]) * (u[0UL] * m_u_ast_y[idx] - u[1UL] * m_u_ast_x[idx]);
			// dtheta_i = uz_i - uz_1
			dyds[5 + idx] = u[2UL] - u1[2UL];
		}
		else
			dyds[2 + idx] = dyds[5 + idx] = 0.00;
	};

	computeTwists(0UL, u1);
	computeTwists(1UL, u2);
	computeTwists(2UL, u3);

	// computing curvature of the first tube along the x and y directions
	blaze::StaticMatrix<double, 3UL, 3UL> dR2_dtheta, dR3_dtheta; // dR1_dtheta, dR2_dtheta, dR3_dtheta;
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> K_inv, K1, K2, K3;
	K_inv(0UL, 0UL) = K_inv(1UL, 1UL) = 1.00 / blaze::sum(m_EI);
	K_inv(2UL, 2UL) = 1.00 / blaze::sum(m_GJ);

	// dR1_dtheta = mathOp::dRz_dTheta(y[5UL]);
	dR2_dtheta = mathOp::dRz_dTheta(y[6UL]);
	dR3_dtheta = mathOp::dRz_dTheta(y[7UL]);

	// du_i = R(theta_i) * (K_i * d_Theta_i * d_R(theta_i)' * u_1 + u_i^ * K_i * (u_i - u_i_ast)
	blaze::StaticVector<double, 3UL> du_1, du_2, du_3, u1_ast, u2_ast, u3_ast, Du;
	// partial component due to tube 1
	blaze::diagonal(K1) = {m_EI[0UL], m_EI[0UL], m_GJ[0UL]};
	u1_ast = {m_u_ast_x[0UL], m_u_ast_y[0UL], 0.0};
	du_1 = mathOp::rotz(y[5UL]) * mathOp::hatPreMultiply(u1, K1) * (u1 - u1_ast);

	// partial component due to tube 2
	blaze::diagonal(K2) = {m_EI[1UL], m_EI[1UL], m_GJ[1UL]};
	u2_ast = {m_u_ast_x[1UL], m_u_ast_y[1UL], 0.0};
	du_2 = mathOp::rotz(y[6UL]) * (K2 * dyds[6UL] * blaze::trans(dR2_dtheta) * u1 + mathOp::hatPreMultiply(u2, K2) * (u2 - u2_ast));

	// partial component due to tube 3
	blaze::diagonal(K3) = {m_EI[2UL], m_EI[2UL], m_GJ[2UL]};
	u3_ast = {m_u_ast_x[2UL], m_u_ast_y[2UL], 0.0};
	du_3 = mathOp::rotz(y[7UL]) * (K3 * dyds[7UL] * blaze::trans(dR3_dtheta) * u1 + mathOp::hatPreMultiply(u3, K3) * (u3 - u3_ast));

	// R (orientation) of the local frame at arc-length s
	mathOp::getSO3(blaze::subvector<11UL, 4UL>(y), R1);

	// Equilibrium bending curvature along x, y
	Du = -K_inv * ((du_1 + du_2 + du_3) + mathOp::hatPreMultiply(m_e3, blaze::trans(R1)) * m_f);

	// curvature of tube 1 along the x and y directions
	blaze::subvector<0UL, 2UL>(dyds) = blaze::subvector<0UL, 2UL>(Du);

	// spatial derivative of the quaternion representation h_dot
	blaze::subvector<11UL, 4UL>(dyds) = mathOp::quaternionDiff(u1, blaze::subvector<11UL, 4UL>(y));

	// calculating r_dot
	blaze::subvector<8UL, 3UL>(dyds) = blaze::column<2UL>(R1); // r_dot = R * e3

	/*
	==================================================================================================================================
	********************************* computing the matrices for the derivative propagation approach *********************************
	==================================================================================================================================
	*/

	// Unpacking the matrix V
	blaze::row<0UL>(this->m_V) = blaze::trans(blaze::subvector<15UL, 17UL>(y));
	blaze::row<1UL>(this->m_V) = blaze::trans(blaze::subvector<32UL, 17UL>(y));
	blaze::row<2UL>(this->m_V) = blaze::trans(blaze::subvector<49UL, 17UL>(y));
	blaze::row<3UL>(this->m_V) = blaze::trans(blaze::subvector<66UL, 17UL>(y));
	blaze::row<4UL>(this->m_V) = blaze::trans(blaze::subvector<83UL, 17UL>(y));
	blaze::row<5UL>(this->m_V) = blaze::trans(blaze::subvector<100UL, 17UL>(y));
	blaze::row<6UL>(this->m_V) = blaze::trans(blaze::subvector<117UL, 17UL>(y));

	// Unpacking the matrix E
	blaze::row<0UL>(this->m_E) = blaze::trans(blaze::subvector<134UL, 17UL>(y));
	blaze::row<1UL>(this->m_E) = blaze::trans(blaze::subvector<151UL, 17UL>(y));
	blaze::row<2UL>(this->m_E) = blaze::trans(blaze::subvector<168UL, 17UL>(y));
	blaze::row<3UL>(this->m_E) = blaze::trans(blaze::subvector<185UL, 17UL>(y));
	blaze::row<4UL>(this->m_E) = blaze::trans(blaze::subvector<202UL, 17UL>(y));
	blaze::row<5UL>(this->m_E) = blaze::trans(blaze::subvector<219UL, 17UL>(y));

	// ==>> MATRIX V_s <<==
	blaze::StaticMatrix<double, 3UL, 3UL> dR1_dxk;
	blaze::StaticMatrix<double, 3UL, 6UL> O_I;
	O_I(0UL, 3UL) = O_I(1UL, 4UL) = O_I(2UL, 5UL) = 1.00;

	blaze::StaticVector<double, 3UL> du1_dxk, du2_dxk, du3_dxk, df_dfi;
	double dTheta2Dot_dxk, dTheta3Dot_dxk, du1z_dxk, du2z_dxk, du3z_dxk;

	auto setDyDs = [&](size_t col, const blaze::StaticVector<double, 3UL> &df_dfi)
	{
		dR1_dxk = mathOp::hatPreMultiply(O_I * blaze::column(m_E, col), R1);														// d_R1/d_xk
		dTheta2Dot_dxk = m_V(5UL, col) - m_V(4UL, col);																				// d_theta2Dot/d_xk
		dTheta3Dot_dxk = m_V(6UL, col) - m_V(4UL, col);																				// d_theta3Dot/d_xk
		du1_dxk = {m_V(2UL, col), m_V(3UL, col), m_V(4UL, col)};																	// d_u1/d_xk
		du2_dxk = m_V(0UL, col) * blaze::trans(dR2_dtheta) * u1 + blaze::trans(R2) * du1_dxk + dTheta2Dot_dxk * m_e3;				// d_u2/d_xk
		du3_dxk = m_V(1UL, col) * blaze::trans(dR3_dtheta) * u1 + blaze::trans(R3) * du1_dxk + dTheta3Dot_dxk * m_e3;				// d_u3/d_xk
		du1z_dxk = (m_GJ[0UL] > 0.00) ? (m_EI[0UL] / m_GJ[0UL]) * (du1_dxk[0UL] * u1_ast[1UL] - du1_dxk[1UL] * u1_ast[0UL]) : 0.00; // du1z/d_xk
		du2z_dxk = (m_GJ[1UL] > 0.00) ? (m_EI[1UL] / m_GJ[1UL]) * (du2_dxk[0UL] * u2_ast[1UL] - du2_dxk[1UL] * u2_ast[0UL]) : 0.00; // du2z/d_xk
		du3z_dxk = (m_GJ[2UL] > 0.00) ? (m_EI[2UL] / m_GJ[2UL]) * (du3_dxk[0UL] * u3_ast[1UL] - du3_dxk[1UL] * u3_ast[0UL]) : 0.00; // du3z/d_xk

		du_1 = dR1_dxk * mathOp::hatPreMultiply(u1, K1) * (u1 - u1_ast) + R1 * (mathOp::hatPreMultiply(du1_dxk, K1) * (u1 - u1_ast) + mathOp::hatPreMultiply(u1, K1) * du1_dxk);
		du_2 = dTheta2Dot_dxk * dR2_dtheta * (K2 * dyds[6UL] * blaze::trans(dR2_dtheta) * u1 + mathOp::hatPreMultiply(u2, K2) * (u2 - u2_ast)) + R2 * (mathOp::hatPreMultiply(du2_dxk, K2) * (u2 - u2_ast) + mathOp::hatPreMultiply(u2, K2) * du2_dxk);
		du_3 = dTheta3Dot_dxk * dR3_dtheta * (K3 * dyds[7UL] * blaze::trans(dR3_dtheta) * u1 + mathOp::hatPreMultiply(u3, K3) * (u3 - u3_ast)) + R3 * (mathOp::hatPreMultiply(du3_dxk, K3) * (u3 - u3_ast) + mathOp::hatPreMultiply(u3, K3) * du3_dxk);

		Du = -K_inv * ((du_1 + du_2 + du_3) + mathOp::hatPreMultiply(m_e3, blaze::trans(dR1_dxk)) * m_f + mathOp::hatPreMultiply(m_e3, blaze::trans(R1)) * df_dfi); // du1Dot_xy/d_xk

		blaze::column(m_Vs, col) = {dTheta2Dot_dxk, dTheta3Dot_dxk, Du[0UL], Du[1UL], du1z_dxk, du2z_dxk, du3z_dxk};
	};

	//>> ============================ dy(s) / d_beta_1  ============================
	setDyDs(0UL, df_dfi);
	//>> ============================ dy(s) / d_beta_2  ============================
	setDyDs(1UL, df_dfi);
	//>> ============================ dy(s) / d_beta_3  ============================
	setDyDs(2UL, df_dfi);
	//>> ============================ dy(s) / d_alpha_1  ============================
	setDyDs(3UL, df_dfi);
	//>> ============================ dy(s) / d_alpha_2  ============================
	setDyDs(4UL, df_dfi);
	//>> ============================ dy(s) / d_alpha_3  ============================
	setDyDs(5UL, df_dfi);
	//>> ============================ dy(s) / d_wf_1  ============================
	df_dfi = {1.0, 0.0, 0.0};
	setDyDs(6UL, df_dfi);
	//>> ============================ dy(s) / d_wf_2  ============================
	df_dfi = {0.0, 1.0, 0.0};
	setDyDs(7UL, df_dfi);
	//>> ============================ dy(s) / d_wf_3  ============================
	df_dfi = {0.0, 0.0, 1.0};
	setDyDs(8UL, df_dfi);

	// must set df_dfi to zero after force derivatives have ben considered
	df_dfi = 0.0;

	//>> ============================ dy(s) / d_wm_1  ============================
	setDyDs(9UL, df_dfi);
	//>> ============================ dy(s) / d_wm_2  ============================
	setDyDs(10UL, df_dfi);
	//>> ============================ dy(s) / d_wm_3  ============================
	setDyDs(11UL, df_dfi);
	//>> =========================== dy(s) / d_u1x(0)  ===========================
	setDyDs(12UL, df_dfi);
	//>> =========================== dy(s) / d_u1y(0)  ===========================
	setDyDs(13UL, df_dfi);
	//>> =========================== dy(s) / d_u1z(0)  ===========================
	setDyDs(14UL, df_dfi);
	//>> =========================== dy(s) / d_u2z(0)  ===========================
	setDyDs(15UL, df_dfi);
	//>> =========================== dy(s) / d_u3z(0)  ===========================
	setDyDs(16UL, df_dfi);

	// ==>> MATRIX E_s <<==
	blaze::StaticMatrix<double, 6UL, 6UL> Ad_g;
	blaze::StaticMatrix<double, 6UL, 17UL> dXi_dyTimesV;

	// building the adjoint matrix
	blaze::submatrix<0UL, 0UL, 3UL, 3UL>(Ad_g) = R1;
	blaze::submatrix<3UL, 3UL, 3UL, 3UL>(Ad_g) = R1;
	blaze::submatrix<0UL, 3UL, 3UL, 3UL>(Ad_g) = mathOp::hatPreMultiply(blaze::subvector<8UL, 3UL>(y), R1);
	// building the matrix dXi_dyTimesV
	blaze::submatrix<3UL, 0UL, 3UL, 17UL>(dXi_dyTimesV) = blaze::submatrix<2UL, 0UL, 3UL, 17UL>(this->m_V);
	// computing Es = Ad_g * dXi_dyTimesV
	this->m_Es = Ad_g * dXi_dyTimesV;

	// Packing the matrix Vs
	blaze::subvector<15UL, 17UL>(dyds) = blaze::trans(blaze::row<0UL>(this->m_Vs));
	blaze::subvector<32UL, 17UL>(dyds) = blaze::trans(blaze::row<1UL>(this->m_Vs));
	blaze::subvector<49UL, 17UL>(dyds) = blaze::trans(blaze::row<2UL>(this->m_Vs));
	blaze::subvector<66UL, 17UL>(dyds) = blaze::trans(blaze::row<3UL>(this->m_Vs));
	blaze::subvector<83UL, 17UL>(dyds) = blaze::trans(blaze::row<4UL>(this->m_Vs));
	blaze::subvector<100UL, 17UL>(dyds) = blaze::trans(blaze::row<5UL>(this->m_Vs));
	blaze::subvector<117UL, 17UL>(dyds) = blaze::trans(blaze::row<6UL>(this->m_Vs));

	// Packing the matrix Es
	blaze::subvector<134UL, 51UL>(dyds) = 0.00;
	blaze::subvector<185UL, 17UL>(dyds) = blaze::trans(blaze::row<3UL>(this->m_Es));
	blaze::subvector<202UL, 17UL>(dyds) = blaze::trans(blaze::row<4UL>(this->m_Es));
	blaze::subvector<219UL, 17UL>(dyds) = blaze::trans(blaze::row<5UL>(this->m_Es));
}

// setter method for updating the parameters forward kinematics computation
void ODESystem::setEquationParameters(const blaze::StaticVector<double, 3UL> &u_ast_x,
									  const blaze::StaticVector<double, 3UL> &u_ast_y,
									  const blaze::StaticVector<double, 3UL> &EI,
									  const blaze::StaticVector<double, 3UL> &GJ,
									  const blaze::StaticVector<double, 3UL> &force)
{
	this->m_u_ast_x = u_ast_x;
	this->m_u_ast_y = u_ast_y;
	this->m_EI = EI;
	this->m_GJ = GJ;
	this->m_f = force;
}