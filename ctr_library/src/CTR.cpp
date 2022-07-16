// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include <boost/config.hpp>
#ifdef BOOST_MSVC
#pragma warning(disable : 4996)
#endif

#include "CTR.hpp"
#include <execution>

// overloaded class constructor
CTR::CTR(const std::array<std::shared_ptr<Tube>, 3UL> &Tb, blaze::StaticVector<double, 6UL> &q, double Tol, mathOp::rootFindingMethod method) : m_accuracy(Tol), m_method(method), m_Tubes(Tb), m_beta(blaze::subvector<0UL, 3UL>(q)), m_q(q)
{
	this->m_theta_0 = {0, q[4UL] - q[3UL], q[5UL] - q[4UL]};
	this->m_e3 = {0, 0, 1};
	this->m_segment = std::make_unique<Segment>(this->m_Tubes, this->m_beta);
	this->m_r_0 = 0.0;
	this->m_h_0 = {1.0, 0.0, 0.0, 0.0};

	this->m_stateEquations = std::make_unique<ODESystem>();
	this->m_stateObserver = std::make_unique<Observer>(this->m_y, this->m_s);
}

// copy constructor
CTR::CTR(const CTR &rhs) : m_accuracy(rhs.m_accuracy), m_method(rhs.m_method), m_Tubes(rhs.m_Tubes),
						   m_beta(rhs.m_beta), m_q(rhs.m_q), m_theta_0(rhs.m_theta_0), m_e3(rhs.m_e3),
						   m_r_0(rhs.m_r_0), m_h_0(rhs.m_h_0), m_y(rhs.m_y), m_s(rhs.m_s)
{
	this->m_segment = std::make_unique<Segment>(this->m_Tubes, this->m_beta);
	this->m_stateEquations = std::make_unique<ODESystem>();
	this->m_stateObserver = std::make_unique<Observer>(this->m_y, this->m_s);
}

// move constructor
CTR::CTR(CTR &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_accuracy = rhs.m_accuracy;
		this->m_method = rhs.m_method;
		this->m_Tubes = std::move(rhs.m_Tubes);
		this->m_beta = std::move(rhs.m_beta);
		this->m_q = std::move(rhs.m_q);
		this->m_theta_0 = std::move(rhs.m_theta_0);
		this->m_e3 = std::move(rhs.m_e3);
		this->m_segment = std::move(rhs.m_segment);
		this->m_r_0 = std::move(rhs.m_r_0);
		this->m_h_0 = std::move(rhs.m_h_0);
		this->m_y = std::move(rhs.m_y);
		this->m_s = std::move(rhs.m_s);
		this->m_stateEquations = std::move(rhs.m_stateEquations);
		this->m_stateObserver = std::move(rhs.m_stateObserver);
	}
}

// CTR destructor
CTR::~CTR()
{
	// nothing to be done as smart pointers are being used
}

// copy assignment operator
CTR &CTR::operator=(const CTR &rhs)
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_accuracy = rhs.m_accuracy;
		this->m_method = rhs.m_method;
		this->m_Tubes = rhs.m_Tubes;
		this->m_beta = rhs.m_beta;
		this->m_q = rhs.m_q;
		this->m_theta_0 = rhs.m_theta_0;
		this->m_e3 = rhs.m_e3;
		this->m_segment = std::make_unique<Segment>(this->m_Tubes, this->m_beta);
		this->m_r_0 = rhs.m_r_0;
		this->m_h_0 = rhs.m_h_0;
		this->m_y = rhs.m_y;
		this->m_s = rhs.m_s;
		this->m_stateEquations = std::make_unique<ODESystem>();
		this->m_stateObserver = std::make_unique<Observer>(this->m_y, this->m_s);
	}

	return *this;
}

// move assignment operator
CTR &CTR::operator=(CTR &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_accuracy = rhs.m_accuracy;
		this->m_method = rhs.m_method;
		this->m_Tubes = std::move(rhs.m_Tubes);
		this->m_beta = std::move(rhs.m_beta);
		this->m_q = std::move(rhs.m_q);
		this->m_theta_0 = std::move(rhs.m_theta_0);
		this->m_e3 = std::move(rhs.m_e3);
		this->m_segment = std::move(rhs.m_segment);
		this->m_r_0 = std::move(rhs.m_r_0);
		this->m_h_0 = std::move(rhs.m_h_0);
		this->m_y = std::move(rhs.m_y);
		this->m_s = std::move(rhs.m_s);
		this->m_stateEquations = std::move(rhs.m_stateEquations);
		this->m_stateObserver = std::move(rhs.m_stateObserver);
	}

	return *this;
}

// function that resets the initial parameters for the ODESolver
void CTR::reset(const blaze::StaticVector<double, 5UL> &initGuess)
{
	vec3d uz_0 = {initGuess[2UL], initGuess[3UL], initGuess[4UL]};
	// alpha1_0 =  alpha_1 - beta_1 * uz_1(0)
	double alpha1_0 = this->m_q[3UL] - this->m_beta[0UL] * uz_0[0UL];

	// clearing the observer's containers
	if (!this->m_s.empty())
	{
		this->m_y.clear();
		this->m_y.reserve(300UL);
		this->m_s.clear();
		this->m_s.reserve(300UL);
	}

	// theta_i(0) = alpha_1 - alpha_i - (beta_i * uz_i(0) - beta_1 * uz_1(0))
	this->m_theta_0 = {0.0, this->m_q[4UL] - this->m_beta[1UL] * uz_0[1UL] - alpha1_0, this->m_q[5UL] - this->m_beta[2UL] * uz_0[2UL] - alpha1_0};

	// origin of the local frame at s = 0
	this->m_r_0 = 0.0;
	// transforming proximal orientation to quaternion representation
	mathOp::euler2Quaternion(0.0, alpha1_0, 0.0, this->m_h_0);

	// ==== adjusting the initial condition for the derivative propagation matrices ====

	//>> MATRIX V(0)

	// m_V0 = 0.0;
	blaze::row<0UL>(this->m_V0) = {-m_beta[1UL] * this->m_V0(5UL, 0UL) + uz_0[0UL] + m_beta[0UL] * this->m_V0(4UL, 0UL), // d_theta2/d_beta1
								   -uz_0[1UL] - m_beta[1UL] * this->m_V0(5UL, 1UL) + m_beta[0UL] * this->m_V0(4UL, 1UL), // d_theta2/d_beta2
								   -m_beta[1UL] * this->m_V0(5UL, 2UL) + m_beta[0UL] * this->m_V0(4UL, 2UL),	   		 // d_theta2/d_beta3
								   -1.0 - m_beta[1UL] * this->m_V0(5UL, 3UL) + m_beta[0UL] * this->m_V0(4UL, 3UL), 		 // d_theta2/d_alpha1
								   1.0 - m_beta[1UL] * this->m_V0(5UL, 4UL) + m_beta[0UL] * this->m_V0(4UL, 4UL), 		 // d_theta2/d_alpha2
								   -m_beta[1UL] * this->m_V0(5UL, 5UL) + m_beta[0UL] * this->m_V0(4UL, 5UL),		 	 // d_theta2/d_alpha3
								   -m_beta[1UL] * this->m_V0(5UL, 6UL) + m_beta[0UL] * this->m_V0(4UL, 6UL),	   		 // d_theta2/d_wf1
								   -m_beta[1UL] * this->m_V0(5UL, 7UL) + m_beta[0UL] * this->m_V0(4UL, 7UL),	   		 // d_theta2/d_wf2
								   -m_beta[1UL] * this->m_V0(5UL, 8UL) + m_beta[0UL] * this->m_V0(4UL, 8UL),	   		 // d_theta2/d_wf3
								   -m_beta[1UL] * this->m_V0(5UL, 9UL) + m_beta[0UL] * this->m_V0(4UL, 9UL),	 		 // d_theta2/d_wm1
								   -m_beta[1UL] * this->m_V0(5UL, 10UL) + m_beta[0UL] * this->m_V0(4UL, 10UL),			 // d_theta2/d_wm2
								   -m_beta[1UL] * this->m_V0(5UL, 11UL) + m_beta[0UL] * this->m_V0(4UL, 11UL),			 // d_theta2/d_wm3
								   -m_beta[1UL] * this->m_V0(5UL, 12UL) + m_beta[0UL] * this->m_V0(4UL, 12UL),			 // d_theta2/d_u1x(0)
								   -m_beta[1UL] * this->m_V0(5UL, 13UL) + m_beta[0UL] * this->m_V0(4UL, 13UL),			 // d_theta2/d_u1y(0)
								   -m_beta[1UL] * this->m_V0(5UL, 14UL) + m_beta[0UL],			 						 // d_theta2/d_u1z(0)
								   -m_beta[1UL] + m_beta[0UL] * this->m_V0(4UL, 15UL),						   			 // d_theta2/d_u2z(0)
								   -m_beta[1UL] * this->m_V0(5UL, 16UL) + m_beta[0UL] * this->m_V0(4UL, 16UL)}; 		 // d_theta2/d_u3z(0)

	blaze::row<1UL>(this->m_V0) = {-m_beta[2UL] * this->m_V0(6UL, 0UL) + uz_0[0UL] + m_beta[0UL] * this->m_V0(4UL, 0UL), // d_theta3/d_beta1
								   -m_beta[2UL] * this->m_V0(6UL, 1UL) + m_beta[0UL] * this->m_V0(4UL, 1UL),			 // d_theta3/d_beta2
								   -m_beta[2UL] * this->m_V0(6UL, 2UL) + m_beta[0UL] * this->m_V0(4UL, 2UL),	   		 // d_theta3/d_beta3
								   -1.0 - m_beta[2UL] * this->m_V0(6UL, 3UL) + m_beta[0UL] * this->m_V0(4UL, 3UL), 		 // d_theta3/d_alpha1
								   -m_beta[2UL] * this->m_V0(6UL, 4UL) + m_beta[0UL] * this->m_V0(4UL, 4UL),		 	 // d_theta3/d_alpha2
								   1.0 - m_beta[2UL] * this->m_V0(6UL, 5UL) + m_beta[0UL] * this->m_V0(4UL, 5UL), 		 // d_theta3/d_alpha3
								   -m_beta[2UL] * this->m_V0(6UL, 6UL) + m_beta[0UL] * this->m_V0(4UL, 6UL),	   		 // d_theta3/d_wf1
								   -m_beta[2UL] * this->m_V0(6UL, 7UL) + m_beta[0UL] * this->m_V0(4UL, 7UL),	   		 // d_theta3/d_wf2
								   -m_beta[2UL] * this->m_V0(6UL, 8UL) + m_beta[0UL] * this->m_V0(4UL, 8UL),	   		 // d_theta3/d_wf3
								   - m_beta[2UL] * this->m_V0(6UL, 9UL) + m_beta[0UL] * this->m_V0(4UL, 9UL), 			 // d_theta3/d_wm1
								   - m_beta[2UL] * this->m_V0(6UL, 10UL) + m_beta[0UL] * this->m_V0(4UL, 10UL), 		 // d_theta3/d_wm2
								   - m_beta[2UL] * this->m_V0(6UL, 11UL) + m_beta[0UL] * this->m_V0(4UL, 11UL), 		 // d_theta3/d_wm3
								   - m_beta[2UL] * this->m_V0(6UL, 12UL) + m_beta[0UL] * this->m_V0(4UL, 12UL), 		 // d_theta3/d_u1x(0)
								   - m_beta[2UL] * this->m_V0(6UL, 13UL) + m_beta[0UL] * this->m_V0(4UL, 13UL), 		 // d_theta3/d_u1y(0)
								   -m_beta[2UL] * this->m_V0(6UL, 14UL) + m_beta[0UL],			 						 // d_theta3/d_u1z(0)
								   -m_beta[2UL] * this->m_V0(6UL, 15UL) + m_beta[0UL] * this->m_V0(4UL, 15UL), 			 // d_theta3/d_u2z(0)
								   -m_beta[2UL] + m_beta[0UL] * this->m_V0(4UL, 16UL)};						   			 // d_theta3/d_u3z(0)

	blaze::row<2UL>(this->m_V0) = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0}; // d_u1x(0)/d_xk

	blaze::row<3UL>(this->m_V0) = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}; // d_u1y(0)/d_xk

	blaze::row<4UL>(this->m_V0) = {(m_V0(1UL,0UL)+m_V0(0UL,0UL)+m_beta[1UL]*m_V0(5UL,0UL)+m_beta[2UL]*m_V0(6UL,0UL))/(2*m_beta[0UL])
									-(m_theta_0[2UL]+m_theta_0[1UL]+2*m_q[3UL]-m_q[4UL]-m_q[5UL]+m_beta[1UL]*uz_0[1UL]+m_beta[2UL]*uz_0[2UL])/(2*m_beta[0UL]*m_beta[0UL]),  // d_u1z/d_beta1
									(m_V0(1UL,1UL)+m_V0(0UL,1UL)+uz_0[1UL]+m_beta[1UL]*m_V0(5UL,1UL)+m_beta[2UL]*m_V0(6UL,1UL))/(2*m_beta[0UL]), 							// d_u1z/d_beta2
									(m_V0(1UL,2UL)+m_V0(0UL,2UL)+m_beta[1UL]*m_V0(5UL,2UL)+uz_0[2UL]+m_beta[2UL]*m_V0(6UL,2UL))/(2*m_beta[0UL]), 							// d_u1z/d_beta3
									(m_V0(1UL,3UL)+m_V0(0UL,3UL)+2.0+m_beta[1UL]*m_V0(5UL,3UL)+m_beta[2UL]*m_V0(6UL,3UL))/(2*m_beta[0UL]), 									// d_u1z/d_alpha1
									(m_V0(1UL,4UL)+m_V0(0UL,4UL)-1.0+m_beta[1UL]*m_V0(5UL,4UL)+m_beta[2UL]*m_V0(6UL,4UL))/(2*m_beta[0UL]), 									// d_u1z/d_alpha2
									(m_V0(1UL,5UL)+m_V0(0UL,5UL)-1.0+m_beta[1UL]*m_V0(5UL,5UL)+m_beta[2UL]*m_V0(6UL,5UL))/(2*m_beta[0UL]), 									// d_u1z/d_alpha3
									(m_V0(1UL,6UL)+m_V0(0UL,6UL)+m_beta[1UL]*m_V0(5UL,6UL)+m_beta[2UL]*m_V0(6UL,6UL))/(2*m_beta[0UL]), 										// d_u1z/d_wf1
									(m_V0(1UL,7UL)+m_V0(0UL,7UL)+m_beta[1UL]*m_V0(5UL,7UL)+m_beta[2UL]*m_V0(6UL,7UL))/(2*m_beta[0UL]), 										// d_u1z/d_wf2
									(m_V0(1UL,8UL)+m_V0(0UL,8UL)+m_beta[1UL]*m_V0(5UL,8UL)+m_beta[2UL]*m_V0(6UL,8UL))/(2*m_beta[0UL]), 										// d_u1z/d_wf3
									(m_V0(1UL,9UL)+m_V0(0UL,9UL)+m_beta[1UL]*m_V0(5UL,9UL)+m_beta[2UL]*m_V0(6UL,9UL))/(2*m_beta[0UL]), 										// d_u1z/d_wm1
									(m_V0(1UL,10UL)+m_V0(0UL,10UL)+m_beta[1UL]*m_V0(5UL,10UL)+m_beta[2UL]*m_V0(6UL,10UL))/(2*m_beta[0UL]), 									// d_u1z/d_wm2
									(m_V0(1UL,11UL)+m_V0(0UL,11UL)+m_beta[1UL]*m_V0(5UL,11UL)+m_beta[2UL]*m_V0(6UL,11UL))/(2*m_beta[0UL]), 									// d_u1z/d_wm3
									0.0,																																	// d_u1z/d_u1x(0)
									0.0,																																	// d_u1z/d_u1y(0)
									1.0,																																	// d_u1z/d_u1z(0)
									(m_V0(1UL,15UL)+m_V0(0UL,15UL)+m_beta[1UL]+m_beta[2UL]*m_V0(6UL,15UL))/(2*m_beta[0UL]), 												// d_u1z/d_u2z(0)
									(m_V0(1UL,16UL)+m_V0(0UL,16UL)+m_beta[1UL]*m_V0(5UL,16UL)+m_beta[2UL])/(2*m_beta[0UL])}; 												// d_u1z/d_u3z(0)

	blaze::row<5UL>(this->m_V0) = {(-m_V0(0UL,0UL)+uz_0[0UL]+m_beta[0UL]*m_V0(4UL,0UL))/m_beta[1UL],  																		// d_u2z/d_beta1
									(m_beta[0UL]*m_V0(4UL,1UL)-m_V0(0UL,1UL))/m_beta[1UL] 
									-(m_q[4UL]-m_q[3UL]-m_theta_0[1UL]+m_beta[0UL]*uz_0[0UL])/(m_beta[1UL]*m_beta[1UL]),													// d_u2z/d_beta2
									(m_beta[0UL]*m_V0(4UL,2UL)-m_V0(0UL,2UL))/m_beta[1UL], 																					// d_u2z/d_beta3
									(-1.0+m_beta[0UL]*m_V0(4UL,3UL)-m_V0(0UL,3UL))/m_beta[1UL], 																			// d_u2z/d_alpha1
									(1.0+m_beta[0UL]*m_V0(4UL,4UL)-m_V0(0UL,4UL))/m_beta[1UL], 																				// d_u2z/d_alpha2
									(m_beta[0UL]*m_V0(4UL,5UL)-m_V0(0UL,5UL))/m_beta[1UL], 																					// d_u2z/d_alpha3
									(m_beta[0UL]*m_V0(4UL,6UL)-m_V0(0UL,6UL))/m_beta[1UL], 																					// d_u2z/d_wf1
									(m_beta[0UL]*m_V0(4UL,7UL)-m_V0(0UL,7UL))/m_beta[1UL], 																					// d_u2z/d_wf2
									(m_beta[0UL]*m_V0(4UL,8UL)-m_V0(0UL,8UL))/m_beta[1UL], 																					// d_u2z/d_wf3
									(m_beta[0UL]*m_V0(4UL,9UL)-m_V0(0UL,9UL))/m_beta[1UL], 																					// d_u2z/d_wm1
									(m_beta[0UL]*m_V0(4UL,10UL)-m_V0(0UL,10UL))/m_beta[1UL], 																				// d_u2z/d_wm2
									(m_beta[0UL]*m_V0(4UL,11UL)-m_V0(0UL,11UL))/m_beta[1UL], 																				// d_u2z/d_wm3
									0.0,																																	// d_u2z/d_u1x(0)
									0.0,																																	// d_u2z/d_u1y(0)
									(-m_V0(0UL,14UL)+m_beta[0UL]+m_beta[0UL]*m_V0(4UL,14UL))/m_beta[1UL],																	// d_u2z/d_u1z(0)
									1.0, 																																	// d_u2z/d_u2z(0)
									(-m_V0(0UL,16UL)+m_beta[0UL]*m_V0(4UL,16UL))/m_beta[1UL]}; 																				// d_u2z/d_u3z(0)

	blaze::row<6UL>(this->m_V0) = {(-m_V0(1UL,0UL)+uz_0[0UL]+m_beta[0UL]*m_V0(4UL,0UL))/m_beta[2UL],  																		// d_u3z/d_beta1
									(m_beta[0UL]*m_V0(4UL,1UL)-m_V0(1UL,1UL))/m_beta[2UL], 																					// d_u3z/d_beta2
									(m_beta[0UL]*m_V0(4UL,2UL)-m_V0(1UL,2UL))/m_beta[2UL] 
									-(m_q[5UL]-m_q[3UL]-m_theta_0[2UL]+m_beta[0UL]*uz_0[0UL])/(m_beta[2UL]*m_beta[2UL]),													// d_u3z/d_beta3
									(-1.0+m_beta[0UL]*m_V0(4UL,3UL)-m_V0(1UL,3UL))/m_beta[2UL], 																			// d_u3z/d_alpha1
									(m_beta[0UL]*m_V0(4UL,4UL)-m_V0(1UL,4UL))/m_beta[2UL], 																					// d_u3z/d_alpha2
									(1.0+m_beta[0UL]*m_V0(4UL,5UL)-m_V0(1UL,5UL))/m_beta[2UL],																				// d_u3z/d_alpha3
									(m_beta[0UL]*m_V0(4UL,6UL)-m_V0(1UL,6UL))/m_beta[2UL], 																					// d_u3z/d_wf1
									(m_beta[0UL]*m_V0(4UL,7UL)-m_V0(1UL,7UL))/m_beta[2UL], 																					// d_u3z/d_wf2
									(m_beta[0UL]*m_V0(4UL,8UL)-m_V0(1UL,8UL))/m_beta[2UL], 																					// d_u3z/d_wf3
									(m_beta[0UL]*m_V0(4UL,9UL)-m_V0(1UL,9UL))/m_beta[2UL], 																					// d_u3z/d_wm1
									(m_beta[0UL]*m_V0(4UL,10UL)-m_V0(1UL,10UL))/m_beta[2UL], 																				// d_u3z/d_wm2
									(m_beta[0UL]*m_V0(4UL,11UL)-m_V0(1UL,11UL))/m_beta[2UL], 																				// d_u3z/d_wm3
									0.0,																																	// d_u3z/d_u1x(0)
									0.0,																																	// d_u3z/d_u1y(0)
									(-m_V0(1UL,14UL)+m_beta[0UL]+m_beta[0UL]*m_V0(4UL,14UL))/m_beta[2UL],																	// d_u3z/d_u1z(0)
									(-m_V0(1UL,15UL)+m_beta[0UL]*m_V0(4UL,15UL))/m_beta[2UL], 																				// d_u3z/d_u2z(0)
									1.0}; 																																	// d_u3z/d_u3z(0)
		
	std::cout << "m_V0 =\n" << m_V0 << std::endl;

	//>> MATRIX E(0)
	blaze::StaticMatrix<double, 4UL, 4UL> g_inv, dH_dxk;
	blaze::StaticMatrix<double, 3UL, 3UL> R0;
	mathOp::getSO3(this->m_h_0, R0);
	blaze::submatrix<0UL, 0UL, 3UL, 3UL>(g_inv) = blaze::trans(R0);
	g_inv(3UL, 3UL) = 1.0;

	double c(cos(alpha1_0)), s(sin(alpha1_0)), scale;

	// dH/d_beta_1
	scale = uz_0[0UL] + m_beta[0UL] * this->m_V0(4UL, 0UL);
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<0UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);
	// dH/d_beta_2
	scale = m_beta[0UL] * this->m_V0(4UL, 1UL);
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<1UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);
	// dH/d_beta_3
	scale = m_beta[0UL] * this->m_V0(4UL, 2UL);
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<2UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);
	// dH/d_alpha_1
	scale = 1 - m_beta[0UL] * this->m_V0(4UL, 3UL);
	dH_dxk = {{-scale * s, -scale * c, 0.0, 0.0}, {scale * c, -scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<3UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);
	// dH/d_alpha_2
	scale = -m_beta[0UL] * this->m_V0(4UL, 4UL);
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<4UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);
	// dH/d_alpha_3
	scale = -m_beta[0UL] * this->m_V0(4UL, 5UL);
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<5UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);

	// dH/d_wf_1
	scale = m_beta[0UL] * this->m_V0(4UL, 6UL);
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<6UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);
	// dH/d_wf_2
	scale = m_beta[0UL] * this->m_V0(4UL, 7UL);
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<7UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);
	// dH/d_wf_3
	scale = m_beta[0UL] * this->m_V0(4UL, 8UL);
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<8UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);
	// dH/d_wm_1
	scale = m_beta[0UL] * this->m_V0(4UL, 9UL);
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<9UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);
	// dH/d_wm_2
	scale = m_beta[0UL] * this->m_V0(4UL, 10UL);
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<10UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);
	// dH/d_wm_3
	scale = m_beta[0UL] * this->m_V0(4UL, 11UL);
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<11UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);

	// dH/d_u1x(0)
	scale = m_beta[0UL] * this->m_V0(4UL, 12UL);
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<12UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);
	// dH/d_u1y(0)
	scale = m_beta[0UL] * this->m_V0(4UL, 13UL);
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<13UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);
	// dH/d_u1z(0)
	scale = m_beta[0UL];
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<14UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);
	// dH/d_u2z(0)
	scale = m_beta[0UL] * this->m_V0(4UL, 15UL);
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<15UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);
	// dH/d_u3z(0)
	scale = m_beta[0UL] * this->m_V0(4UL, 16UL);
	dH_dxk = {{scale * s, scale * c, 0.0, 0.0}, {-scale * c, scale * s, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}, {0.0, 0.0, 0.0, 0.0}};
	blaze::column<16UL>(this->m_E0) = mathOp::wedgeOperator(dH_dxk * g_inv);

	std::cout << "m_E0 =\n" << m_E0 << std::endl;
}

// function that solves (integrates) the CTR ode (state) equations
blaze::StaticVector<double, 5UL> CTR::ODESolver(const blaze::StaticVector<double, 5UL> &initGuess)
{
	// initGuess = [u1_x(0) u1_y(0) u1_z(0) u2_z(0) u3_z(0)]  --> vector of initial guesses for solving the BVP
	this->reset(initGuess); // resets CTR parameters and variables for a new iteration of the ode-solver

	// retrieving the bending & torsional stiffness and precurvatures in all segments of the CTR in the current configuration
	auto tuple = this->m_segment->returnParameters();

	blaze::HybridMatrix<double, 3UL, 18UL, blaze::columnMajor> EI(std::get<0UL>(tuple)), GJ(std::get<1UL>(tuple)), U_x(std::get<2UL>(tuple)), U_y(std::get<3UL>(tuple));
	std::vector<double> S(std::get<4UL>(tuple));

	// ##################################################### NUMERICAL METHODS FOR ODE INTEGRATION #####################################################

	// ********************************  8-th ORDER ADAPTIVE ADAMS-BASHFORTH-MOULTON STEPPER ********************************
	boost::numeric::odeint::adaptive_adams_bashforth_moulton<8UL, state_type, double, state_type, double,
															 boost::numeric::odeint::vector_space_algebra, boost::numeric::odeint::default_operations, boost::numeric::odeint::initially_resizer>
		abm8_stepper;

	// ********************************  4-th ORDER CLASSIC RUNGE-KUTTA STEPPER ********************************
	// typedef boost::numeric::odeint::runge_kutta4_classic<state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra,
	// 	boost::numeric::odeint::default_operations, boost::numeric::odeint::initially_resizer> rk4_stepper;

	// ********************************  5-th ORDER CASH-KARP STEPPER ********************************
	// typedef boost::numeric::odeint::runge_kutta_cash_karp54<state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra,
	// 	boost::numeric::odeint::default_operations, boost::numeric::odeint::initially_resizer> rkk54_stepper;

	// ********************************  5-th ORDER DORMAND-PRINCE RUNGE-KUTTA ********************************
	// typedef boost::numeric::odeint::runge_kutta_dopri5<state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra> rkd5_stepper;

	// ********************************  BULIRSCH-STOER STEPPER ********************************
	// typedef boost::numeric::odeint::bulirsch_stoer<state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra,
	// 	boost::numeric::odeint::default_operations, boost::numeric::odeint::initially_resizer> blstr_stepper;

	// ******************************** RUUNGE-KUTTA-FEHLBERG (RKF78) STEPPER ********************************
	// typedef boost::numeric::odeint::runge_kutta_fehlberg78<state_type, double, state_type, double, boost::numeric::odeint::vector_space_algebra,
	// 	boost::numeric::odeint::default_operations, boost::numeric::odeint::initially_resizer> rk78_stepper;

	// #################################################################################################################################################

	// initial guess for the BC (twist curvatures)
	vec2d u1_xy_BC = {initGuess[0UL], initGuess[1UL]};
	vec3d uz_BC, theta_BC(this->m_theta_0);
	uz_BC = {initGuess[2UL], initGuess[3UL], initGuess[4UL]};
	blaze::StaticMatrix<double, 7UL, 17UL> V_BC(this->m_V0);
	blaze::StaticMatrix<double, 6UL, 17UL> E_BC(this->m_E0);

	// start and end points, in terms of arc-length s, of each CTR segment and initial step-size for integration (ds)
	double s_start, s_end, ds;

	// instantiating variables needed for enforcing the continuity of internal moments and Jacobians across CTR segments
	vec3d u1, u2, u3, u1_next, u1_ast, u2_ast, u3_ast, u1_ast_next, u2_ast_next, u3_ast_next;
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> K1, K2, K3, K1_next, K2_next, K3_next;
	blaze::StaticMatrix<double, 3UL, 3UL> R, R_theta_2, R_theta_3;
	double d_theta2, d_theta3;

	blaze::StaticMatrix<double, 7UL, 17UL> V_next;
	blaze::StaticMatrix<double, 6UL, 17UL> E_next;
	blaze::StaticVector<double, 7UL> ydot, ydot_next;
	blaze::StaticMatrix<double, 6UL, 6UL> Ad_g;
	blaze::StaticVector<double, 6UL> xi, xi_next, Ad_gXi, Ad_gXi_next;
	blaze::StaticVector<double, 3UL> du_1, du_2, du_3;

	// instantiating the vector of initial conditions for solving the state equations (236 x 1)
	state_type y_0;

	// iterating through the tube segments comprising the CTR
	size_t len_seg = S.size() - 1;
	for (size_t seg = 0; seg < len_seg; ++seg)
	{
		// 1st element of y computes the curvature of the first (innermost) tube along the x direction
		// 2nd element of y computes the curvature of the first (innermost) tube along the y direction
		// next 3 elements of y are the torsional curvatures for the three tubes, e.g., ui_z = [u1_z  u2_z  u3_z]
		// next 3 elements of y are the twist angles, theta_i = [theta_1  theta_2  theta_3]
		// last 12 elements are r(position) and R(orientations) of the local frames, respectively at each arc-length s

		// initializing the initial conditions vector (236 x 1)
		blaze::subvector<0UL, 15UL>(y_0) = {u1_xy_BC[0UL], u1_xy_BC[1UL],
											uz_BC[0UL], uz_BC[1UL], uz_BC[2UL],
											theta_BC[0UL], theta_BC[1UL], theta_BC[2UL],
											m_r_0[0UL], m_r_0[1UL], m_r_0[2UL],
											m_h_0[0UL], m_h_0[1UL], m_h_0[2UL], m_h_0[3UL]};
		// packs the V_BC and E_BC matrices into the state vector
		this->packVMatrix(V_BC, y_0);
		this->packEMatrix(E_BC, y_0);

		// specifying the interval of integration (in terms of tube segment arc-lengths)
		s_start = S[seg];
		s_end = S[seg + 1];
		ds = (s_end - s_start) / 20; // 20 points per segment

		// passing the tube parameters in the segment to the state equation method
		this->m_stateEquations->setEquationParameters(blaze::column(U_x, seg), blaze::column(U_y, seg), blaze::column(EI, seg), blaze::column(GJ, seg), this->m_wf);

		// ##################################################### NUMERICAL INTEGRATION #####################################################
		// Employs the selected stepper (Numerical method) and integrates the system of ODEs along the segment considered
		boost::numeric::odeint::integrate_adaptive(abm8_stepper, *this->m_stateEquations, y_0, s_start, s_end, ds, *this->m_stateObserver);

		// rate of twist of tube 2 at the end of the i - th segment (y0 is changed to the approximative solution f(x,t') at the end of integration)
		d_theta2 = y_0[3UL] - y_0[2UL];
		// rate of twist of tube 3 at the end of the i - th segment
		d_theta3 = y_0[4UL] - y_0[2UL];

		// setting the boundary(initial conditions) for the next segment
		uz_BC = blaze::subvector<2UL, 3UL>(y_0);
		theta_BC = blaze::subvector<5UL, 3UL>(y_0);
		m_r_0 = blaze::subvector<8UL, 3UL>(y_0);
		m_h_0 = blaze::subvector<11UL, 4UL>(y_0);
		u1 = blaze::subvector<0UL, 3UL>(y_0);

		if (seg < len_seg - 1)
		{
			// stiffness matrices in the current segment
			blaze::diagonal(K1) = {EI(0UL, seg), EI(0UL, seg), GJ(0UL, seg)};
			blaze::diagonal(K2) = {EI(1UL, seg), EI(1UL, seg), GJ(1UL, seg)};
			blaze::diagonal(K3) = {EI(2UL, seg), EI(2UL, seg), GJ(2UL, seg)};
			// precurvature vectors in the current segment
			u1_ast = {U_x(0, seg), U_y(0, seg), 0.0};
			u2_ast = {U_x(1, seg), U_y(1, seg), 0.0};
			u3_ast = {U_x(2, seg), U_y(2, seg), 0.0};
			// stiffness matrices in the next segment
			blaze::diagonal(K1_next) = {EI(0UL, seg + 1), EI(0UL, seg + 1), GJ(0UL, seg + 1)};
			blaze::diagonal(K2_next) = {EI(1UL, seg + 1), EI(1UL, seg + 1), GJ(1UL, seg + 1)};
			blaze::diagonal(K3_next) = {EI(2UL, seg + 1), EI(2UL, seg + 1), GJ(2UL, seg + 1)};
			// precurvature vectors in the next segment
			u1_ast_next = {U_x(0, seg + 1), U_y(0, seg + 1), 0.0};
			u2_ast_next = {U_x(1, seg + 1), U_y(1, seg + 1), 0.0};
			u3_ast_next = {U_x(2, seg + 1), U_y(2, seg + 1), 0.0};

			// orientation of tube 2's local frame at the end of the i-th segment
			R_theta_2 = mathOp::rotz(theta_BC[1UL]);
			// orientation of tube 3's local frame at the end of the i-th segment
			R_theta_3 = mathOp::rotz(theta_BC[2UL]);

			// curvatures of tubes 2 and 3 at the end of the current segment
			u2 = mathOp::transposePreMultiply(R_theta_2, u1) + d_theta2 * m_e3;
			u3 = mathOp::transposePreMultiply(R_theta_3, u1) + d_theta3 * m_e3;

			// estimating curvature of next segment through continuity of internal bending moments
			// * * * # --- # ====>> For reference, check: Computing Jacobians and Compliance Matrices(Rucker, 2011) pg. 5 <<==== # --- # * * *
			u1_next = blaze::inv(K1_next + K2_next + K3_next) * (K1 * (u1 - u1_ast) + R_theta_2 * K2 * (u2 - u2_ast) + R_theta_3 * K3 * (u3 - u3_ast) + K1_next * u1_ast_next + R_theta_2 * K2_next * (u2_ast_next - d_theta2 * m_e3) + R_theta_3 * K3_next * (u3_ast_next - d_theta3 * m_e3));

			// effectively setting the initial condition for u1_dot along the x, y directions
			u1_xy_BC = blaze::subvector<0UL, 2UL>(u1_next);

			/*
			==================================================================================================================================
			********************************** Enforcing Continuity of the Derivative Propagation Matrices **********************************
			==================================================================================================================================
			*/

			// ydot = [theta2_dot, theta3_dot, u1x_dot, u1y_dot, u1z_dot, u2z_dot, u3z_dot]

			// ====>> ydot
			// theta2_dot
			ydot[0UL] = (K2(2UL, 2UL) > 0.0) ? d_theta2 : 0.0; // u2[2UL] - u1[2UL];

			// theta3_dot
			ydot[1UL] = (K3(2UL, 2UL) > 0.0) ? d_theta3 : 0.0; // u3[2UL] - u1[2UL];

			// u1xy_dot
			du_1 = mathOp::rotz(y_0[5UL]) * mathOp::hatPreMultiply(u1, K1) * (u1 - u1_ast);
			du_2 = mathOp::rotz(y_0[6UL]) * (K2 * d_theta2 * mathOp::rotz_dot_transpose(theta_BC[1UL]) * u1 + mathOp::hatPreMultiply(u2, K2) * (u2 - u2_ast));
			du_3 = mathOp::rotz(y_0[7UL]) * (K3 * d_theta3 * mathOp::rotz_dot_transpose(theta_BC[2UL]) * u1 + mathOp::hatPreMultiply(u3, K3) * (u3 - u3_ast));
			// R (orientation) of the local frame at arc-length s
			mathOp::getSO3(this->m_h_0, R);
			blaze::subvector<2UL, 3UL>(ydot) = -blaze::inv(K1 + K2 + K3) * ((du_1 + du_2 + du_3) + mathOp::hatPreMultiply(m_e3, blaze::trans(R)) * this->m_wf);

			// u1z_dot
			ydot[4UL] = (K1(2UL, 2UL) > 0.0) ? (K1(0UL, 0UL) / K1(2UL, 2UL)) * (u1[0UL] * u1_ast[1UL] - u1[1UL] * u1_ast[0UL]) : 0.0;

			// u2z_dot
			ydot[5UL] = (K2(2UL, 2UL) > 0.0) ? (K2(0UL, 0UL) / K2(2UL, 2UL)) * (u2[0UL] * u2_ast[1UL] - u2[1UL] * u2_ast[0UL]) : 0.0;

			// u3z_dot
			ydot[6UL] = (K3(2UL, 2UL) > 0.0) ? (K3(0UL, 0UL) / K3(2UL, 2UL)) * (u3[0UL] * u3_ast[1UL] - u3[1UL] * u3_ast[0UL]) : 0.0;

			// ------------------------------------------------------------------------------------------------------------------------------------------------
			// ====>> ydot_next
			// theta2_dot
			ydot_next[0UL] = (K2_next(2UL, 2UL) > 0.0) ? d_theta2 : 0.0; // u2[2UL] - u1[2UL];

			// theta3_dot
			ydot_next[1UL] = (K3_next(2UL, 2UL) > 0.0) ? d_theta3 : 0.0; // u3[2UL] - u1[2UL];

			// u1xy_dot
			u2 = blaze::trans(R_theta_2) * u1_next + d_theta2 * m_e3; // u2_next
			u3 = blaze::trans(R_theta_3) * u1_next + d_theta3 * m_e3; // u3_next
			du_1 = mathOp::rotz(y_0[5UL]) * mathOp::hatPreMultiply(u1_next, K1_next) * (u1_next - u1_ast_next);
			du_2 = mathOp::rotz(y_0[6UL]) * (K2_next * d_theta2 * mathOp::rotz_dot_transpose(theta_BC[1UL]) * u1_next + mathOp::hatPreMultiply(u2, K2_next) * (u2 - u2_ast_next));
			du_3 = mathOp::rotz(y_0[7UL]) * (K3_next * d_theta3 * mathOp::rotz_dot_transpose(theta_BC[2UL]) * u1_next + mathOp::hatPreMultiply(u3, K3_next) * (u3 - u3_ast_next));
			// R (orientation) of the local frame at arc-length s
			mathOp::getSO3(this->m_h_0, R);
			blaze::subvector<2UL, 3UL>(ydot_next) = -blaze::inv(K1_next + K2_next + K3_next) * ((du_1 + du_2 + du_3) + mathOp::hatPreMultiply(m_e3, blaze::trans(R)) * this->m_wf);

			// u1z_dot
			ydot_next[4UL] = (K1_next(2UL, 2UL) > 0.0) ? (K1_next(0UL, 0UL) / K1_next(2UL, 2UL)) * (u1_next[0UL] * u1_ast_next[1UL] - u1_next[1UL] * u1_ast_next[0UL]) : 0.0;

			// u2z_dot
			ydot_next[5UL] = (K2_next(2UL, 2UL) > 0.0) ? (K2_next(0UL, 0UL) / K2_next(2UL, 2UL)) * (u2[0UL] * u2_ast_next[1UL] - u2[1UL] * u2_ast_next[0UL]) : 0.0;

			// u3z_dot
			ydot_next[6UL] = (K3_next(2UL, 2UL) > 0.0) ? (K3_next(0UL, 0UL) / K3_next(2UL, 2UL)) * (u3[0UL] * u3_ast_next[1UL] - u3[1UL] * u3_ast_next[0UL]) : 0.0;

			xi = {0.0, 0.0, 1.0, u1[0UL], u1[1UL], u1[2UL]};
			xi_next = {0.0, 0.0, 1.0, u1_next[0UL], u1_next[1UL], u1_next[2UL]};

			// building the adjoint matrix
			blaze::submatrix<0UL, 0UL, 3UL, 3UL>(Ad_g) = R;
			blaze::submatrix<3UL, 3UL, 3UL, 3UL>(Ad_g) = R;
			blaze::submatrix<0UL, 3UL, 3UL, 3UL>(Ad_g) = mathOp::hatPreMultiply(m_r_0, R);

			Ad_gXi = Ad_g * xi;
			Ad_gXi_next = Ad_g * xi_next;

			// unpacking the V matrix for enforcing continuity of the matrix differential equations across boundaries
			this->m_V = this->unpackVMatrix(y_0);			

			// unpacking the E matrix for enforcing continuity of the matrix differential equations across boundaries
			this->m_E = this->unpackEMatrix(y_0);

			for (size_t col = 0UL; col < 17UL; ++col)
			{
				blaze::column(V_BC, col) = blaze::column(this->m_V, col) + ydot - ydot_next;
				blaze::column(E_BC, col) = blaze::column(this->m_E, col) + Ad_gXi - Ad_gXi_next;
			}

			// std::cout << "V_BC =\n" << V_BC << std::endl;
			// std::cout << "E_BC =\n" << E_BC << std::endl;
		}
		else
		{
			// stiffness matrices in the current segment
			blaze::diagonal(K1) = {EI(0UL, seg), EI(0UL, seg), GJ(0UL, seg)};
			u1_ast = {U_x(0, seg), U_y(0, seg), 0.0};
			mathOp::getSO3(this->m_h_0, R);
			this->m_V = this->unpackVMatrix(y_0);
			this->m_E = this->unpackEMatrix(y_0);
		}
	}

	//
	//	****************  #####  -------------- ___ DISTAL BOUNDARY CONDITIONS ___ --------------  #####  ****************
	//			1) internal moment at the tip of tube 1 must equal the external moment applied
	//			   Namely, R1K1(u1(L1) - u1 * (L1)) = Wm = > u1(l1) - u1 * (l1)-K1\R1'Wm = 0
	//
	//			2) at the distal ends of the remaining tubes, the axial component of the internal moments must equal zero
	//			   Namely, ui, z(Li) - ui, z* (Li) = 0

	// external moment applied at the distal end of the CTR
	blaze::StaticVector<double, 3UL> axialMoment = blaze::inv(K1) * blaze::trans(R) * this->m_wm;

	// Residue vector due to infringment of the distal boundary conditions
	blaze::StaticVector<double, 5UL> b;

	// computing the violation of the distal boundary conditions
	b = {u1[0UL] - u1_ast[0UL] - axialMoment[0UL],
		 u1[1UL] - u1_ast[1UL] - axialMoment[1UL],
		 u1[2UL] - axialMoment[2UL],
		 0.0,
		 0.0};	

	// lambda function that finds the u_z curvatures at the distal ends of tubes 2 and 3
	auto computeResidue = [&](double distalEnd, size_t index)->size_t
	{
		// must use some tolerance when comparing floating points
		auto itt = std::find_if(this->m_s.begin(), this->m_s.end(), [&](double x)
								{ return (std::abs(x - distalEnd) <= 1e-7) ? true : false; }); // finds where tube ends (with a 0.0001mm tolerance)

		size_t id = std::distance(this->m_s.begin(), itt);
		// ui_z at the distal end of the i-th tube
		b[2UL + index] = this->m_y[id][2UL + index];

		// return the index of discretized arc-length at which the i-th tube ends
		return id;
	};

	// Computing the Residues associated to the twist curvatures (rate of twist) at the distal ends of tubes 2 and 3
	vec3d distEnd(this->m_segment->getDistEnd()); // arc-lengths at the distal ends of each one of the tubes

	size_t distalTb2 = computeResidue(distEnd[1UL], 1UL);
	size_t distalTb3 = computeResidue(distEnd[2UL], 2UL);

	// computes the Jacobian of residues
	// ==>> MATRIX B <<==
	blaze::StaticMatrix<double, 3UL, 3UL> dR1_dxk;
	blaze::StaticMatrix<double, 3UL, 6UL> O_I;
	O_I(0UL, 3UL) = O_I(1UL, 4UL) = O_I(2UL, 5UL) = 1.0;

	blaze::StaticVector<double, 3UL> du1_dxk, dm_dmi, res;

	auto db_dxk = [&](const size_t col, const size_t distalTb2, const size_t distalTb3, const blaze::StaticVector<double, 3UL> &dm_dmi)
	{
		dR1_dxk = mathOp::hatPreMultiply(O_I * blaze::column(m_E, col), R); // d_R1/d_xk

		du1_dxk = {m_V(2UL, col), m_V(3UL, col), m_V(4UL, col)};
		res = du1_dxk - blaze::inv(K1) * (blaze::trans(dR1_dxk) * this->m_wm + blaze::trans(R) * dm_dmi);

		blaze::column(this->m_B, col) = {res[0UL], res[1UL], res[2UL], this->m_y[distalTb2][99UL + col], this->m_y[distalTb3][116UL + col]};
	};

	//>> ============================ db / d_beta_1  ============================
	db_dxk(0UL, distalTb2, distalTb3, dm_dmi);
	//>> ============================ db / d_beta_2  ============================
	db_dxk(1UL, distalTb2, distalTb3, dm_dmi);
	//>> ============================ db / d_beta_3  ============================
	db_dxk(2UL, distalTb2, distalTb3, dm_dmi);
	//>> ============================ db / d_alpha_1  ============================
	db_dxk(3UL, distalTb2, distalTb3, dm_dmi);
	//>> ============================ db / d_alpha_2  ============================
	db_dxk(4UL, distalTb2, distalTb3, dm_dmi);
	//>> ============================ db / d_alpha_3  ============================
	db_dxk(5UL, distalTb2, distalTb3, dm_dmi);
	//>> ============================ db / d_wf_1  ============================
	db_dxk(6UL, distalTb2, distalTb3, dm_dmi);
	//>> ============================ db / d_wf_2  ============================
	db_dxk(7UL, distalTb2, distalTb3, dm_dmi);
	//>> ============================ db / d_wf_3  ============================
	db_dxk(8UL, distalTb2, distalTb3, dm_dmi);
	//>> ============================ db / d_wm_1  ============================
	dm_dmi = {1.0, 0.0, 0.0};
	db_dxk(9UL, distalTb2, distalTb3, dm_dmi);
	//>> ============================ db / d_wm_2  ============================
	dm_dmi = {0.0, 1.0, 0.0};
	db_dxk(10UL, distalTb2, distalTb3, dm_dmi);
	//>> ============================ db / d_wm_3  ============================
	dm_dmi = {0.0, 0.0, 1.0};
	db_dxk(11UL, distalTb2, distalTb3, dm_dmi);

	dm_dmi = 0.0;
	//>> =========================== db / d_u1x(0)  ===========================
	db_dxk(12UL, distalTb2, distalTb3, dm_dmi);
	//>> =========================== db / d_u1y(0)  ===========================
	db_dxk(13UL, distalTb2, distalTb3, dm_dmi);
	//>> =========================== db / d_u1z(0)  ===========================
	db_dxk(14UL, distalTb2, distalTb3, dm_dmi);
	//>> =========================== db / d_u2z(0)  ===========================
	db_dxk(15UL, distalTb2, distalTb3, dm_dmi);
	//>> =========================== db / d_u3z(0)  ===========================
	db_dxk(16UL, distalTb2, distalTb3, dm_dmi);

	std::cout << "initGuess = " << blaze::trans(initGuess);
	std::cout << "ODESolver residue = " << blaze::trans(b) << std::endl;

	return b;
}

// function that computes the finite-differences Jacobian for solving the BVP
blaze::StaticMatrix<double, 5UL, 5UL> CTR::jac_BVP()
{
	blaze::StaticMatrix<double, 5UL, 5UL> Bu = blaze::submatrix<0UL, 12UL, 5UL, 5UL>(m_B);

	return Bu;
}

// function that computes the finite-differences Jacobian wrt actuation inputs
std::pair<blaze::StaticMatrix<double, 6UL, 6UL>, blaze::StaticMatrix<double, 6UL, 6UL>> CTR::jacobian()
{
	blaze::StaticMatrix<double, 6UL, 6UL> J, C;

	auto Eq = blaze::submatrix<0UL, 0UL, 6UL, 6UL>(m_E);
	auto Ew = blaze::submatrix<0UL, 6UL, 6UL, 6UL>(m_E);
	auto Eu = blaze::submatrix<0UL, 12UL, 6UL, 5UL>(m_E);

	auto Bq = blaze::submatrix<0UL, 0UL, 5UL, 6UL>(m_B);
	auto Bw = blaze::submatrix<0UL, 6UL, 5UL, 6UL>(m_B);
	auto Bu = blaze::submatrix<0UL, 12UL, 5UL, 5UL>(m_B);

	auto EuBu_inv = Eu * mathOp::pInv(Bu);
	
	J = Eq - EuBu_inv * Bq;
	C = Ew - EuBu_inv * Bw;

	return std::make_pair(J, C);
}

// function that implements Powell's Dog Leg Method (Nonlinear root-finding method for solving the BVP)
bool CTR::PowellDogLeg(blaze::StaticVector<double, 5UL> &initGuess)
{
	bool found;
	size_t k = 0UL;
	const size_t k_max = 300UL;
	double alpha, beta, delta, eps1, eps2, rho, c;
	blaze::StaticVector<double, 5UL> g, f, f_new, x_new, h_sd, h_gn, h_dl;
	blaze::StaticMatrix<double, 5UL, 5UL> J;

	// initializing parameters
	delta = 1;
	eps1 = eps2 = 1e-22;

	f = this->ODESolver(initGuess);
	J = this->jac_BVP();
	g = blaze::trans(J) * f;

	// checking if the initial guess satisfies the BVP without the need of any further refinement
	found = ((blaze::linfNorm(f) <= this->m_accuracy) || (blaze::linfNorm(g) <= eps1)) ? true : false;

	while (!found && (k < k_max))
	{
		k++;

		alpha = blaze::sqrNorm(g) / blaze::sqrNorm(J * g);
		h_sd = -alpha * g;			 // steepest descend (this is a direction, not a step!)
		h_gn = -mathOp::pInv(J) * f; // Gauss-Newton step (Least Square solution)

		// two candidates for the step to take from this point, a = alpha*h_sd & b = h_gn

		// computing the dog leg direction
		if (blaze::norm(h_gn) <= delta)
			h_dl = h_gn;
		else
		{
			if (blaze::norm(h_sd) >= delta)
				h_dl = delta * blaze::normalize(h_sd);
			else
			{
				c = blaze::trans(h_sd) * (h_gn - h_sd);

				if (c <= 0.0)
				{
					beta = (-c + sqrt(c * c + blaze::sqrNorm(h_gn - h_sd) * (delta * delta - blaze::sqrNorm(h_sd)))) / blaze::sqrNorm(h_gn - h_sd);
				}
				else
				{
					beta = (delta * delta - blaze::sqrNorm(h_sd)) / (c + sqrt(c * c + blaze::sqrNorm(h_gn - h_sd) * (delta * delta - blaze::sqrNorm(h_sd))));
				}

				h_dl = h_sd + beta * (h_gn - h_sd); // Dog Leg step
			}
		}

		if (blaze::norm(h_dl) <= eps2 * (blaze::norm(initGuess) + eps2))
			found = true;
		else
		{
			x_new = initGuess + h_dl;
			f_new = this->ODESolver(x_new);
			rho = (blaze::sqrNorm(f) - blaze::sqrNorm(f_new)) / (0.5 * blaze::trans(h_dl) * ((delta * h_dl) - g));

			if (rho > 0.0)
			{
				initGuess = std::move(x_new);
				f = std::move(f_new);
				J = this->jac_BVP();
				g = blaze::trans(J) * f;

				if ((blaze::linfNorm(f) <= this->m_accuracy) || (blaze::linfNorm(g) <= eps1))
					found = true;
			}

			if (rho > 0.75)
				delta = std::max(delta, 3 * blaze::norm(h_dl));
			else
			{
				if (rho < 0.25)
					delta *= 0.5;
			}

			if (delta < eps2 * (blaze::norm(initGuess) + eps2))
				found = true;
		}
	}

	return found;
}

// function that implements the Levenberg-Marquardt Method (Nonlinear root-finding method for solving the BVP)
bool CTR::Levenberg_Marquardt(blaze::StaticVector<double, 5UL> &initGuess)
{
	size_t k = 0UL;
	const size_t k_max = 300UL;
	blaze::StaticVector<double, 5UL> h, g, f, f_new;
	blaze::StaticMatrix<double, 5UL, 5UL> J, A;
	blaze::IdentityMatrix<double> I(5UL);
	double rho, nu = 2, mu, tau = 1e-3, e1 = 1e-18, e2 = 1e-25;
	bool found;

	// computing the residue and residue Jacobian associated to initGuess
	f = this->ODESolver(initGuess);
	J = this->jac_BVP();
	A = blaze::trans(J) * J;
	g = blaze::trans(J) * f;
	found = (blaze::linfNorm(g) <= e1) ? true : false;
	mu = tau * blaze::max(blaze::diagonal(A));

	// starting the iterative minimization loop
	while ((!found) && (k < k_max))
	{
		k++;
		blaze::solve(blaze::declsym(A + (mu * I)), h, -g);

		f_new = this->ODESolver(initGuess + h);
		rho = (blaze::sqrNorm(f) - blaze::sqrNorm(f_new)) / (0.5 * blaze::trans(h) * ((mu * h) - g));

		if (rho > 0.0)
		{
			// accept the decrease in the function
			initGuess += h;
			// computing the residue Jacobian at the new initial guess
			J = this->jac_BVP();
			A = blaze::trans(J) * J;
			f = std::move(f_new);
			g = blaze::trans(J) * f;
			found = (blaze::linfNorm(g) <= e1) ? true : false;
			mu = mu * std::max(0.33333333, 1 - blaze::pow(2 * rho - 1, 3));
			nu = 2;
		}
		else
		{
			mu = mu * nu;
			nu = 2 * nu;
		}

		// checking if the tolerance has been satisfied
		if (blaze::linfNorm(f) <= this->m_accuracy)
			found = true;
	}

	return found;
}

// function that implements the Broyden (Nonlinear root-finding method for solving the BVP)
bool CTR::Broyden(blaze::StaticVector<double, 5UL> &initGuess)
{
	// found: returns true (false) when the root-finding method converges (does not converge) within k_max iterations
	bool found;

	// initial Hessian matrix --> computed via finite differences
	blaze::StaticMatrix<double, 5UL, 5UL> JacInv, JacInvNew;

	// setting up and starting my handmadeBFGS method
	blaze::StaticVector<double, 5UL> F, Fold, X, Xold, deltaX, deltaF; // staticVectors are automatically initialized to 0

	// Residue yielded by the initial guess for the CTR BVP
	F = this->ODESolver(initGuess); // F(x_k)	: residue
	X = std::move(initGuess);		// x_k		: initial guess
	JacInvNew = JacInv = mathOp::pInv(this->jac_BVP());

	// checking if the initial guess already satisfies the BVP
	found = (blaze::linfNorm(F) <= this->m_accuracy) ? true : false;

	size_t k = 0UL;
	const size_t k_max = 300UL;
	while (!found && (k < k_max))
	{
		k++;

		deltaX = X - Xold; // dX := x_k - x_k-1
		deltaF = F - Fold; // dF := F(x_k) - F(x_k-1)

		JacInv = std::move(JacInvNew);
		if ((blaze::norm(deltaX) > 0.0) && (blaze::norm(deltaF) > 0.0))
			JacInvNew = JacInv + ((deltaX - JacInv * deltaF) / (blaze::trans(deltaX) * JacInv * deltaF)) * blaze::trans(deltaX) * JacInv;
		else
			JacInvNew = JacInv;

		Xold = std::move(X);
		Fold = std::move(F);

		// update the initial guess
		X = Xold - JacInv * F;
		F = this->ODESolver(X);

		while (blaze::isnan(F))
		{
			if (blaze::isnan(X))
				X = 0.0;
			else
				X /= blaze::max(blaze::abs(X));

			F = this->ODESolver(X);
			JacInv = JacInvNew = mathOp::pInv(this->jac_BVP());
			Xold = std::move(X);
			X = Xold - JacInv * F;
		}

		if (k % 10 == 0.0)
		{
			JacInv = JacInvNew = mathOp::pInv(this->jac_BVP());
			X = Xold - JacInv * F;
		}

		if (blaze::linfNorm(F) <= this->m_accuracy)
			found = true;
	}

	initGuess = std::move(X);
	return found;
}

// function that implements Broyden's Nonlinear root-finding method for solving the BVP
bool CTR::Broyden_II(blaze::StaticVector<double, 5UL> &initGuess)
{
	bool found;

	// initial Hessian matrix --> computed via finite differences
	blaze::StaticMatrix<double, 5UL, 5UL> Jac, JacNew;

	// setting up and starting my handmadeBFGS method
	blaze::StaticVector<double, 5UL> F, Fold, X, Xold, deltaX, deltaF; // staticVectors are automatically initialized to 0

	// Residue yielded by the initial guess for the CTR BVP
	F = this->ODESolver(initGuess); // F(x_k)	: residue
	X = std::move(initGuess);		// x_k		: initial guess
	JacNew = this->jac_BVP();

	// checking if the initial guess already satisfies the BVP
	found = (blaze::linfNorm(F) <= this->m_accuracy) ? true : false;

	size_t k = 0UL;
	const size_t k_max = 300UL;
	while (!found && (k < k_max))
	{
		k++;

		deltaX = X - Xold; // dX := x_k - x_k-1
		deltaF = F - Fold; // dF := F(x_k) - F(x_k-1)

		Jac = std::move(JacNew);
		if (blaze::sqrNorm(deltaX) > 0.0)
			JacNew = Jac + blaze::sqrNorm(deltaX) * (deltaF - (Jac * deltaX)) * blaze::trans(deltaX);
		else
			JacNew = Jac;

		Xold = std::move(X);
		Fold = std::move(F);

		// update the initial guess
		X = Xold - mathOp::pInv(Jac) * F;
		F = this->ODESolver(X);

		while (blaze::isnan(F))
		{
			if (blaze::isnan(X))
				X = 0.0;
			else
				X /= blaze::max(blaze::abs(X));

			F = this->ODESolver(X);
			JacNew = this->jac_BVP();
			Xold = std::move(X);
			X = Xold - mathOp::pInv(Jac) * F;
		}

		if (k % 10 == 0.0)
		{
			JacNew = this->jac_BVP();
			X = Xold - mathOp::pInv(Jac) * F;
		}

		if (blaze::linfNorm(F) <= this->m_accuracy)
			found = true;
	}

	initGuess = std::move(X);
	return found;
}

// function that implements the Newton-Raphson method (Nonlinear root-finding method for solving the BVP)
bool CTR::Newton_Raphson(blaze::StaticVector<double, 5UL> &initGuess)
{
	bool found;
	// setting up and starting my handmade Newton-Raphson method
	blaze::StaticVector<double, 5UL> Residue, Residue_new, d_Residue, dGuess; // staticVectors are automatically initialized to 0

	// Residue of the unperturbed initial guess for the CTR
	Residue = this->ODESolver(initGuess);

	found = (blaze::linfNorm(Residue) <= this->m_accuracy) ? true : false;

	//  Weighing matrices for adjusting the initial guess iteratively (Implementing a PD regulator)
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 5UL, 5UL, blaze::rowMajor>> Kp, Kd;
	blaze::StaticMatrix<double, 5UL, 5UL, blaze::columnMajor> jac_bvp;
	blaze::diagonal(Kp) = 0.3;	// 0.45 | 0.6  | 0.3
	blaze::diagonal(Kd) = 2e-3; // 3e-3 | 5e-3 | 2e-3

	size_t k = 0UL;
	const size_t k_max = 300UL;

	// starting iterations for adjusting the initial guess "u_guess ~ initGuess"
	while (!found && (k < k_max))
	{
		k++;
		jac_bvp = this->jac_BVP();
		// error equation(globally asymptotically stable)
		dGuess = mathOp::pInv(jac_bvp) * (Kp * Residue + Kd * d_Residue);
		// updating the initial guess(weighted negative gradient of the cost function)
		initGuess -= dGuess;
		// computing the new cost associated to the newly readjusted initial guess
		Residue_new = this->ODESolver(initGuess);

		// Checking if the Jacobian has large elements beyond machine precision
		while (blaze::isnan(Residue_new))
		{
			if (blaze::isnan(initGuess))
				initGuess = 0.0;
			else
				initGuess /= blaze::max(blaze::abs(initGuess));

			Residue_new = this->ODESolver(initGuess);
			d_Residue = Residue_new - Residue;
			Residue = std::move(Residue_new);
			continue;
		}

		// cost variation due to initial guess refinement
		d_Residue = Residue_new - Residue;
		// updating the cost
		Residue = std::move(Residue_new);

		if (blaze::linfNorm(Residue) <= this->m_accuracy)
			found = true;
	}

	return found;
}

// function that implements the Modified, globally convergent Newton-Raphson method (Nonlinear root-finding method for solving the BVP)
bool CTR::Modified_Newton_Raphson(blaze::StaticVector<double, 5UL> &initGuess)
{
	/*
		Algorithm extracted from page 309 of Introduction to Numerical Analysis 3rd edition by Josef Stoer & Roland Bulirsch
	*/

	bool found;
	// computes the residue associated to the initial guess
	blaze::StaticVector<double, 5UL> f(this->ODESolver(initGuess)), d;
	blaze::StaticVector<double, 5UL, blaze::rowVector> Dh;
	blaze::StaticMatrix<double, 5UL, 5UL> D, D_inv;
	double h, h_0, lambda, gamma, improvementFactor, d_norm, Dh_norm;
	size_t j = 0UL, k = 0UL;
	const size_t k_max = 1500UL;
	std::vector<double> h_k; // vector to store all h_k's
	h_k.reserve(k_max);
	std::vector<double>::iterator result; // iterator for determining the min element in the vector h_k

	found = (blaze::linfNorm(f) <= this->m_accuracy) ? true : false;

	while (!found && (k < k_max))
	{
		k++;
		// computing the residue Jacobian
		D = this->jac_BVP();
		D_inv = mathOp::pInv(D);

		// search direction (directional derivative)
		d = D_inv * f;
		gamma = 1 / (blaze::norm(D_inv) * blaze::norm(D)); // gamma := 1/cond(Df)
		h_0 = blaze::sqrNorm(f);						   // h := f'f
		// Dh := D(f'f) = 2f'Df
		Dh = 2 * blaze::trans(f) * D;
		d_norm = blaze::norm(d);
		Dh_norm = blaze::norm(Dh);

		while (true)
		{
			f = this->ODESolver(initGuess - blaze::pow(0.5, j) * d);
			// std::cout << "Modified_Newton_Raphson -- j = : " << j << " | residue = " << blaze::trans(f);
			while (blaze::isnan(f))
			{
				j++;
				f = this->ODESolver(initGuess - blaze::pow(0.5, j) * d);
			}
			h = blaze::sqrNorm(f);
			improvementFactor = blaze::pow(0.5, j) * 0.25 * gamma * d_norm * Dh_norm;
			// storig the value of h_k to determine step size posteriorly
			h_k.push_back(h);

			if (h <= (h_0 - improvementFactor))
				break;
			else
				j++;
		}

		// determining the value of the step-size lambda
		result = std::min_element(h_k.begin(), h_k.end());
		// retrieving the minimum h_k
		lambda = blaze::pow(0.5, std::distance(h_k.begin(), result));
		initGuess -= lambda * d;
		h_k.clear();

		// resets the exponent variable j
		j = 0UL;

		// compute the residue associated to the newly refined initGuess
		f = this->ODESolver(initGuess);

		// checking the terminating condition
		if (blaze::linfNorm(f) <= this->m_accuracy)
		{
			found = true;
		}
		else
		{
			if (blaze::isnan(f))
			{
				if (blaze::isnan(initGuess))
					initGuess = 0.0;
				else
					initGuess /= blaze::max(blaze::abs(initGuess));

				// recompute the residue of the readjusted initial guess
				f = this->ODESolver(initGuess);
			}
		}
	}

	if (!found)
	{
		found = this->PowellDogLeg(initGuess);

		if (!found)
			found = this->Newton_Raphson(initGuess);

		if (!found)
			found = this->Levenberg_Marquardt(initGuess);
	}

	return found;
}

// function that implements the CTR actuation for any inputs joint values q
bool CTR::actuate_CTR(blaze::StaticVector<double, 5UL> &initGuess, const blaze::StaticVector<double, 6UL> &q_input)
{
	// boolean flag for indicating convergence (1: zero found | 0: zero not found)
	bool found = false;

	// updating the CTR joints for desired input values
	this->setConfiguration(q_input);

	// recalculates the CTR transition points and segments
	m_segment->recalculateSegments(this->m_Tubes, this->m_beta);

	// initial guess for proximal boundary condition--[u_x(0) u_y(0) u1_z(0) u2_z(0) u3_z(0)]
	switch (this->m_method)
	{
	case mathOp::rootFindingMethod::NEWTON_RAPHSON:
		found = this->Newton_Raphson(initGuess);
		break;
	case mathOp::rootFindingMethod::LEVENBERG_MARQUARDT:
		found = this->Levenberg_Marquardt(initGuess);
		break;
	case mathOp::rootFindingMethod::POWELL_DOG_LEG:
		found = this->PowellDogLeg(initGuess);
		break;
	case mathOp::rootFindingMethod::MODIFIED_NEWTON_RAPHSON:
		found = this->Modified_Newton_Raphson(initGuess);
		break;
	case mathOp::rootFindingMethod::BROYDEN:
		found = this->Broyden(initGuess);
		break;
	case mathOp::rootFindingMethod::BROYDEN_II:
		found = this->Broyden_II(initGuess);
		break;
	}

	return found;
}

// sets the external point force acting on the distal end of the CTR
void CTR::setDistalForce(const blaze::StaticVector<double, 3UL> &f)
{
	this->m_wf = f;
}

// sets the external pure moment acting on the distal end of the CTR
void CTR::setDistalMoment(const blaze::StaticVector<double, 3UL> &m)
{
	this->m_wm = m;
}

// sets the external force and pure moment acting on the distal end of the CTR
void CTR::setDistalWrench(const blaze::StaticVector<double, 3UL> &f, const blaze::StaticVector<double, 3UL> &m)
{
	this->m_wf = f;
	this->m_wf = m;
}

// function that implements the position control ==> returns [u_0, Jac, q_min, timeout]
std::tuple<blaze::StaticMatrix<double, 3UL, 6UL>, blaze::StaticVector<double, 6UL>, bool> CTR::posCTRL(blaze::StaticVector<double, 5UL> &initGuess, const vec3d &target, const double posTol)
{
	double t = 0.05;											 		// time base
	double minError = 1e3;										 		// minimum distance to target
	bool status;												 		// status = TRUE (FALSE) indicates convergence (lack thereof)
	blaze::StaticMatrix<double, 6UL, 6UL, blaze::columnMajor> J, J_inv; // Jacobian matrix
	// proportional gain for control
	blaze::DiagonalMatrix<blaze::StaticMatrix<double, 6UL, 6UL>> Kp;
	blaze::diagonal(Kp) = 1.0;

	// Capturing the CTR's current joint configuration
	blaze::StaticVector<double, 6UL> dqdt, q_min(this->m_q), q(this->m_q);
	// Capturing the proximal BC for the minimum distance
	blaze::StaticVector<double, 5UL> initGuessMin(initGuess);
	// Calculate the CTR Jacobian in the present configuration and retrieves convergence status
	status = this->actuate_CTR(initGuess, q);
	if (!status)
	{
		// timeout = 1UL; // Newton-Raphson convergence failure
		return std::make_tuple(J, q, status);
	}

	// Capturing the position of the end-effector
	vec3d x_CTR;
	blaze::StaticVector<double, 6UL> poseError;
	auto tipError = blaze::subvector<3UL, 3UL>(poseError);

	x_CTR = this->getTipPos();
	// Current position error
	tipError = target - x_CTR;

	// Euclidean distance to target
	double dist2Tgt = blaze::norm(tipError);

	if (dist2Tgt < minError)
	{
		minError = dist2Tgt;
		q_min = q;

		if (dist2Tgt <= posTol)
			return std::make_tuple(J, q_min, status);
	}

	// function to implement actuators sigularity avoidance
	blaze::StaticVector<double, 6UL> f;
	// clearance between linear actuators
	double Clr = 5e-3, deltaBar = 0.0;
	// lengths of straight sections of the CTR tubes
	vec3d ls, L;
	L = {this->m_Tubes[0UL]->getTubeLength(), this->m_Tubes[1UL]->getTubeLength(), this->m_Tubes[2UL]->getTubeLength()};
	ls = {this->m_Tubes[0UL]->getStraightLen(), this->m_Tubes[1UL]->getStraightLen(), this->m_Tubes[2UL]->getStraightLen()};
	// lower and upper bounds on prismatic joint limits
	vec3d betaMax, betaMin;

	size_t N_itr = 0UL;			   // iterations counter
	const size_t maxIter = 1500UL; // maximum admissible number of iterations in the position control loop

	blaze::IdentityMatrix<double, blaze::rowMajor> I(6UL); // 6 x 6 Identity matrix

	// parameters for local optimization (joint limits avoidance)
	auto f1 = blaze::subvector<3UL, 3UL>(f);
	double ke = 15;

	// pointers to the revolute joints of the CTR
	auto revJoints = blaze::subvector<3UL, 3UL>(q);

	// position control loop
	while ((dist2Tgt > posTol) && (N_itr < maxIter))
	{
		// incrementing the number of iterations
		N_itr++;

		// compute the Jacobian in the present configuration
		auto pair = this->jacobian();
		J = pair.first;
		// C = pair.second;
		// Pseudo-inverse of Jacobian for resolving CTR joint motion rates
		J_inv = mathOp::pInv(J);

		// Nullspace control (collision and actuation limits)
		betaMin = {std::max(-ls[0UL] + deltaBar, L[1UL] + this->m_beta[1UL] - L[0UL] + deltaBar),
				   std::max({-ls[1UL] + deltaBar, L[2UL] + this->m_beta[2UL] - L[1UL] + deltaBar, this->m_beta[0UL] + Clr}),
				   std::max(-ls[2UL] + deltaBar, this->m_beta[1UL] + Clr)};

		betaMax = {m_beta[1UL] - Clr,
				   std::min(L[0UL] + this->m_beta[0UL] - deltaBar - L[1UL], this->m_beta[2UL] - Clr),
				   std::min(L[1UL] + this->m_beta[1UL] - deltaBar - L[2UL], -deltaBar)};

		// penalty function for local optimization (actuator collision avoidance)
		// Had to add an infinitesimal (1e-10) to the denominador (betaMax - betaMin + 1e-12) to avoid blow-up when the actuator interval degenerates to a point
		f1 = blaze::pow(blaze::abs((betaMax + betaMin - 2 * this->m_beta) / (betaMax - betaMin + 1e-10)), ke) * blaze::sign(this->m_beta - (betaMax + betaMin) * 0.5);

		// Inverse kinematics
		blaze::subvector<0UL, 3UL>(f) = 1e-4 * f1;
		// resolved rates -- Nullspacec local optimization (joint limit avoidance)
		dqdt = (J_inv * Kp * poseError + (I - blaze::trans(J_inv * J)) * (-f)) * t;
		// dqdt = J_inv * Kp * poseError * t;

		auto q_dot = [&] { // rescaling joint variables for limit avoidance
			for (size_t i = 0UL; i < 3UL; ++i)
			{
				if (this->m_beta[i] + dqdt[i] > betaMax[i])
					dqdt[i] = (betaMax[i] - this->m_beta[i]) * 0.5;

				if (this->m_beta[i] + dqdt[i] < betaMin[i])
					dqdt[i] = (betaMin[i] - this->m_beta[i]) * 0.5;
			}
			return dqdt;
		};

		dqdt = q_dot();

		// updating the CTR joints->q: [beta, theta]
		q += dqdt;

		// wrapping the actuation angles to the [-Pi,Pi) interval
		revJoints = blaze::map(revJoints, [](double theta)
							   {
				static constexpr double TWO_PI = 2 * M_PI;
				double wrappedAngle;

				// wrappedAngle = fmod(theta, TWO_PI);
				wrappedAngle = remainder(theta, TWO_PI);

				return wrappedAngle; });

		// actuate the CTR to new configuration and retrieve execution timeout status
		status = this->actuate_CTR(initGuess, q);

		// interrupts the loop execution if actuation fails
		if (!status)
		{
			initGuess = 0.0;
			status = this->actuate_CTR(initGuess, q);

			if (!status)
			{
				std::cerr << "############################################# posCTR() #############################################" << std::endl;
				std::cerr << "		posCTRL() ==> Nonlinear root-finder did not converge when solving the BVP problem!" << std::endl;
				std::cerr << "############################################# posCTR() #############################################" << std::endl
						  << std::endl;
			}
		}

		// CTR tip position after control adjustment
		x_CTR = this->getTipPos();

		tipError = target - x_CTR;
		dist2Tgt = blaze::norm(tipError);

		if (dist2Tgt < minError)
		{
			minError = dist2Tgt;
			q_min = q;
			initGuessMin = initGuess;
		}

		// stops the control loop when the position update becomes significantly small
		if (blaze::norm(dqdt) <= 1e-6)
		{
			std::cout << "Exited out of position control loop due small incremental threshold!" << std::endl;
			return std::make_tuple(J, q_min, status);
		}

		// std::cout << "Tip error after adustment: " << dist2Tgt << " Tip Position: " << blaze::trans(x_CTR);
	}

	// Actuating the CTR to the configuration which yields the minimum position error
	initGuess = std::move(initGuessMin);
	this->actuate_CTR(initGuess, q_min);

	// std::cout << "CTR Inverse Kinematics ended in " << N_itr << " iterations. Minimum error: " << minError << std::endl;

	return std::make_tuple(J, q_min, status);
}

// returns the arc-length of the backbone which is the closest point to the calyx / renal pyramid
double CTR::returnArcLength(const vec3d &calyx)
{
	// lambda comparator: compute the point on the backbone closest to the renal calyx
	auto dist = [&](state_type x, state_type y)
	{
		return blaze::norm(blaze::subvector<8UL, 3UL>(x) - calyx) < blaze::norm(blaze::subvector<8UL, 3UL>(y) - calyx);
	};

	auto it = std::min_element(this->m_y.begin(), this->m_y.end(), dist);
	// index of the nearest neighbor on the array
	size_t idx = std::distance(this->m_y.begin(), it);

	// returning the arc-length value
	return this->m_s[idx];
}

// getter method for retrieving the torsional curvatures at proximal end of the CTR (s = 0)
vec3d CTR::getProximalTwist()
{
	return blaze::subvector<2UL, 3UL>(this->m_y.front());
}

// getter method for retrieving the torsional curvatures at the last available CTR segment
vec3d CTR::getDistalTwist()
{
	return blaze::subvector<2UL, 3UL>(this->m_y.back());
}

// function that returns the Vector of tubes comprising the CTR
std::array<std::shared_ptr<Tube>, 3UL> CTR::getTubes()
{
	return this->m_Tubes;
}

// function that returns the current linear joint values of the CTR
vec3d CTR::getBeta()
{
	return this->m_beta;
}

// function that returns the current joint values of the CTR
blaze::StaticVector<double, 6UL> CTR::getConfiguration()
{
	return this->m_q;
}

// function that returns the position of the CTR tip
vec3d CTR::getTipPos()
{
	vec3d pos;
	if (!this->m_y.empty())
		pos = blaze::subvector<8UL, 3UL>(this->m_y.back());

	return pos;
}

// function that returns the arc-lenghts at each tube's distal end
vec3d CTR::getDistalEnds()
{
	return this->m_segment->getDistalEnds();
}

// function that returns the individual tube shapes
std::tuple<blaze::HybridMatrix<double, 3UL, 200UL, blaze::columnMajor>, blaze::HybridMatrix<double, 3UL, 200UL, blaze::columnMajor>, blaze::HybridMatrix<double, 3UL, 200UL, blaze::columnMajor>> CTR::getTubeShapes()
{
	blaze::HybridMatrix<double, 3UL, 200UL, blaze::columnMajor> Tb_1, Tb_2, Tb_3;
	// arc-lengths at the distal ends of each tube
	vec3d distal_idx, distalEnds(this->m_segment->getDistalEnds());

	// lambda retrieves the shape of the ith tube, transforms to patient coord frame and rescales to mm
	auto tubeShape = [&](size_t tube_index, blaze::HybridMatrix<double, 3UL, 200UL, blaze::columnMajor> &M)
	{
		double element;
		size_t numOfPoints;

		// find the index in the arc-length vector at which each tube ends
		element = distalEnds[tube_index];
		std::vector<double>::iterator it = std::find_if(this->m_s.begin(), this->m_s.end(), [&element](double x)
														{ return (std::abs(x - element) <= 1e-7) ? true : false; }); // finds where tube ends (0.0001mm tolerance)

		numOfPoints = std::distance(this->m_s.begin(), it);

		// resizes the columns of the hybrid matrix M to the exact number of points
		M.resize(3UL, numOfPoints);

		for (size_t col = 0UL; col < numOfPoints; ++col)
			blaze::column(M, col) = {this->m_y[col][8UL], this->m_y[col][9UL], this->m_y[col][10UL] };
	};

	tubeShape(0UL, Tb_1);
	tubeShape(1UL, Tb_2);
	tubeShape(2UL, Tb_3);

	// returns the tuple containing the shape of the tubes
	return std::make_tuple(Tb_1, Tb_2, Tb_3);
}

// function that returns a vector with the CTR shape
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> CTR::getShape()
{
	std::vector<double> r_x, r_y, r_z;
	r_x.reserve(300UL);
	r_y.reserve(300UL);
	r_y.reserve(300UL);

	if (this->m_y.size() > 0UL)
	{
		for (size_t i = 0UL; i < this->m_y.size(); ++i)
		{
			r_x.push_back(this->m_y[i][8UL]);
			r_y.push_back(this->m_y[i][9UL]);
			r_z.push_back(this->m_y[i][10UL]);
		}
	}

	/*printf("s = \n");
	for (double d : m_s)
		printf("%lf, ", d);
	printf("\n\n");*/

	return std::make_tuple(r_x, r_y, r_z);
}

// setter method for setting the actuation joint values (without actuating the CTR) <--> used for computing the Jacobian
void CTR::setConfiguration(const blaze::StaticVector<double, 6UL> &q)
{
	this->m_q = q;
	this->m_beta = blaze::subvector<0UL, 3UL>(this->m_q);
}

// function that sets which method to use for solving the BVP
void CTR::setBVPMethod(mathOp::rootFindingMethod mthd)
{
	this->m_method = mthd;
}

// function that returns the orientation of the local frame at the distal end
void CTR::getDistalSO3()
{
	blaze::StaticVector<double, 4UL> h_distal = blaze::subvector<11UL, 4UL>(this->m_y.back());
	blaze::StaticMatrix<double, 3UL, 3UL> R;
	mathOp::getSO3(h_distal, R);
	std::cout << "At the distal end, R =\n"
			  << R << std::endl;
}

// function that packs the derivative propagation matrix V into the state vector
void CTR::packVMatrix(const blaze::StaticMatrix<double, 7UL, 17UL> &V, state_type &s)
{
	blaze::subvector<15UL, 17UL>(s) = blaze::trans(blaze::row<0UL>(V));
	blaze::subvector<32UL, 17UL>(s) = blaze::trans(blaze::row<1UL>(V));
	blaze::subvector<49UL, 17UL>(s) = blaze::trans(blaze::row<2UL>(V));
	blaze::subvector<66UL, 17UL>(s) = blaze::trans(blaze::row<3UL>(V));
	blaze::subvector<83UL, 17UL>(s) = blaze::trans(blaze::row<4UL>(V));
	blaze::subvector<100UL, 17UL>(s) = blaze::trans(blaze::row<5UL>(V));
	blaze::subvector<117UL, 17UL>(s) = blaze::trans(blaze::row<6UL>(V));
}

// function that unpacks the derivative propagation matrix V from the state vector
blaze::StaticMatrix<double, 7UL, 17UL> CTR::unpackVMatrix(const state_type &s)
{
	blaze::StaticMatrix<double, 7UL, 17UL> V;

	blaze::row<0UL>(V) = blaze::trans(blaze::subvector<15UL, 17UL>(s));
	blaze::row<1UL>(V) = blaze::trans(blaze::subvector<32UL, 17UL>(s));
	blaze::row<2UL>(V) = blaze::trans(blaze::subvector<49UL, 17UL>(s));
	blaze::row<3UL>(V) = blaze::trans(blaze::subvector<66UL, 17UL>(s));
	blaze::row<4UL>(V) = blaze::trans(blaze::subvector<83UL, 17UL>(s));
	blaze::row<5UL>(V) = blaze::trans(blaze::subvector<100UL, 17UL>(s));
	blaze::row<6UL>(V) = blaze::trans(blaze::subvector<117UL, 17UL>(s));

	return V;
}

// function that packs the derivative propagation matrix E into the state vector
void CTR::packEMatrix(const blaze::StaticMatrix<double, 6UL, 17UL> &E, state_type &s)
{
	blaze::subvector<134UL, 51UL>(s) = 0.0;
	blaze::subvector<185UL, 17UL>(s) = blaze::trans(blaze::row<3UL>(E));
	blaze::subvector<202UL, 17UL>(s) = blaze::trans(blaze::row<4UL>(E));
	blaze::subvector<219UL, 17UL>(s) = blaze::trans(blaze::row<5UL>(E));
}

// function that unpacks the derivative propagation matrix E from the state vector
blaze::StaticMatrix<double, 6UL, 17UL> CTR::unpackEMatrix(const state_type &s)
{
	blaze::StaticMatrix<double, 6UL, 17UL> E;

	blaze::row<0UL>(E) = blaze::trans(blaze::subvector<134UL, 17UL>(s));
	blaze::row<1UL>(E) = blaze::trans(blaze::subvector<151UL, 17UL>(s));
	blaze::row<2UL>(E) = blaze::trans(blaze::subvector<168UL, 17UL>(s));
	blaze::row<3UL>(E) = blaze::trans(blaze::subvector<185UL, 17UL>(s));
	blaze::row<4UL>(E) = blaze::trans(blaze::subvector<202UL, 17UL>(s));
	blaze::row<5UL>(E) = blaze::trans(blaze::subvector<219UL, 17UL>(s));

	return E;
}