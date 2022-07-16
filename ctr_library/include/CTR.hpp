// ***************************************************************************************** //
// *  This file is part of NIH_CSTAR, a CTR kinematics library for Concentric Tube Robots  * //
// *																					   * //
// *  ----------- # Copyright (C) 2021 Filipe C. Pedrosa <fpedrosa@uwo.ca> # -----------   * //
// *																					   * //
// *  Project developed under the supervision of Prof Dr Rajni Patel <rvpatel@uwo.ca>	   * //
// *			  CSTAR (Canadian Surgical Technologies & Advanced Robotics)			   * //
// *				   Western University, London, ON, Canada							   * //
// ***************************************************************************************** //

#pragma once

#include "Tube.hpp"
#include "Segment.hpp"
#include "ODESystem.hpp"
#include "Observer.hpp"
#include "mathOperations.hpp"
#include <boost/numeric/odeint.hpp>
#include <memory>
#include <tuple>
#include <chrono>

// 2-dimensional vector for u_{x,y} BC
typedef blaze::StaticVector<double, 2UL> vec2d;

class CTR
{
public:
	// removes the default constructor
	CTR() = delete;

	// overloaded constructor
	CTR(const std::array<std::shared_ptr<Tube>, 3UL> &Tb, blaze::StaticVector<double, 6UL> &q, const double Tol, const mathOp::rootFindingMethod method);

	// copy constructor
	CTR(const CTR &rhs);

	// move constructor
	CTR(CTR &&rhs)
	noexcept;

	// CTR destructor
	~CTR();

	// copy assignment operator
	CTR &operator=(const CTR &rhs);

	// move assignment operator
	CTR &operator=(CTR &&rhs) noexcept;

	// function that resets the initial parameters for the ODESolver
	void reset(const blaze::StaticVector<double, 5UL> &initGuess);

	// function that solves (integrates) the CTR ode (state) equations
	blaze::StaticVector<double, 5UL> ODESolver(const blaze::StaticVector<double, 5UL> &initGuess);

	// function that computes the Jacobian of residues for solving the BVP
	blaze::StaticMatrix<double, 5UL, 5UL> jac_BVP();

	// function that computes the spatial velocities Jacobian wrt actuation inputs
	std::pair<blaze::StaticMatrix<double, 6UL, 6UL>, blaze::StaticMatrix<double, 6UL, 6UL>> jacobian();

	// function that implements Powell's Dog Leg Method (Nonlinear root-finding method for solving the BVP)
	bool PowellDogLeg(blaze::StaticVector<double, 5UL> &initGuess);

	// function that implements the Levenberg-Marquardt Method (Nonlinear root-finding method for solving the BVP)
	bool Levenberg_Marquardt(blaze::StaticVector<double, 5UL> &initGuess);

	// function that implements Broyden's Nonlinear root-finding method for solving the BVP (Jacobian inverse is estimated)
	bool Broyden(blaze::StaticVector<double, 5UL> &initGuess);

	// function that implements Broyden's Nonlinear root-finding method for solving the BVP
	bool Broyden_II(blaze::StaticVector<double, 5UL> &initGuess);

	// function that implements the Newton-Raphson method (Nonlinear root-finding method for solving the BVP)
	bool Newton_Raphson(blaze::StaticVector<double, 5UL> &initGuess);

	// function that implements the Modified, globally convergent Newton-Raphson method (Nonlinear root-finding method for solving the BVP)
	bool Modified_Newton_Raphson(blaze::StaticVector<double, 5UL> &initGuess);

	// function that implements the CTR actuation for any inputs joint values q
	bool actuate_CTR(blaze::StaticVector<double, 5UL> &initGuess, const blaze::StaticVector<double, 6UL> &q_input);

	// sets the external point force acting on the distal end of the CTR
	void setDistalForce(const blaze::StaticVector<double, 3UL> &f);

	// sets the external pure moment acting on the distal end of the CTR
	void setDistalMoment(const blaze::StaticVector<double, 3UL> &m); 

	// sets the external force and pure moment acting on the distal end of the CTR
	void setDistalWrench(const blaze::StaticVector<double, 3UL> &f, const blaze::StaticVector<double, 3UL> &m);

	// function that implements the position control ==> returns [u_0, Jac, q_min, timeout]
	std::tuple<blaze::StaticMatrix<double, 3UL, 6UL>, blaze::StaticVector<double, 6UL>, bool> posCTRL(blaze::StaticVector<double, 5UL> &initGuess, const vec3d &target, const double Tol);
	
	// returns the arc-length of the backbone which is the closest point to the renal pyramic
	double returnArcLength(const vec3d &calyx);

	// getter method for retrieving the torsional curvatures at proximal end of the CTR (s = 0)
	vec3d getProximalTwist();

	// getter method for retrieving the torsional curvatures at the last available CTR segment
	vec3d getDistalTwist();

	// function that returns the Vector of tubes comprising the CTR
	std::array<std::shared_ptr<Tube>, 3UL> getTubes();

	// function that returns the current linear joint values of the CTR
	vec3d getBeta();

	// function that returns the current joint values of the CTR
	blaze::StaticVector<double, 6UL> getConfiguration();

	// function that returns the position of the CTR tip
	vec3d getTipPos();

	// function that returns the arc-lenghts of each tube's distal end
	vec3d getDistalEnds(); 

	// function that returns the individual tube shapes
	std::tuple<blaze::HybridMatrix<double, 3UL, 200UL, blaze::columnMajor>, blaze::HybridMatrix<double, 3UL, 200UL, blaze::columnMajor>, blaze::HybridMatrix<double, 3UL, 200UL, blaze::columnMajor>> getTubeShapes();

	// function that returns a vector with the CTR shape
	std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> getShape();

	// setter method for setting the actuation joint values (without actuating the CTR) <--> used for computing the Jacobian
	void setConfiguration(const blaze::StaticVector<double, 6UL> &q);

	// function that sets which method to use for solving the BVP
	void setBVPMethod(mathOp::rootFindingMethod mthd);

	// function that returns the orientation of the local frame at the distal end
	void getDistalSO3();

	// function that packs the derivative propagation matrix V into the state vector
	void packVMatrix(const blaze::StaticMatrix<double, 7UL, 17UL>& V, state_type& s);

	// function that unpacks the derivative propagation matrix V from the state vector
	blaze::StaticMatrix<double, 7UL, 17UL> unpackVMatrix(const state_type& s);

	// function that packs the derivative propagation matrix E into the state vector
	void packEMatrix(const blaze::StaticMatrix<double, 6UL, 17UL>& E, state_type& s);

	// function that unpacks the derivative propagation matrix E from the state vector
	blaze::StaticMatrix<double, 6UL, 17UL> unpackEMatrix(const state_type& s);

private:
	double m_accuracy;									// accuracy to which the BVP must be solved
	mathOp::rootFindingMethod m_method;					// methods available: Newton-Raphson, Levenberg-Marquardt, Broyden
	std::array<std::shared_ptr<Tube>, 3UL> m_Tubes; 	// Vector of tubes comprising the CTR
	vec3d m_beta;										// linear actuation
	blaze::StaticVector<double, 6UL> m_q;				// joint actuation values
	vec3d m_theta_0;									// initial twist angle for all tubes at s = 0
	vec3d m_wf;											// external force at the CTR tip (force component of the external wrench)
	vec3d m_wm;											// external moment at the CTR tip (moment component of the external wrench
	vec3d m_e3;											// third canonical basis of R^3
	std::unique_ptr<Segment> m_segment;					// segments between transition points in the CTR
	vec3d m_r_0;										// origin of the local frame at s = 0 or at the end of the i-th segment (for BC)
	blaze::StaticVector<double, 4UL> m_h_0;				// initial orientation of local frame at s = 0 (or at the end of the i-th segment (for BC)) ==>> QUATERNION
	std::vector<state_type> m_y;						// stores the CTR state vector at each integration step
	std::vector<double> m_s;							// stores the arc-length points along the backbone
	std::unique_ptr<ODESystem> m_stateEquations;		// implements the state differential equations for a three-tube CTR
	std::unique_ptr<Observer> m_stateObserver;			// implements the state observer for Boost::odeInt
	// matrices for implementing the derivative propagation approach
	blaze::StaticMatrix<double, 6UL, 17UL> m_E, m_E0;
	blaze::StaticMatrix<double, 5UL, 17UL> m_B;
	blaze::StaticMatrix<double, 7UL, 17UL> m_V, m_V0;
};