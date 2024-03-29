// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "Tube.hpp"
#include <math.h> // pow()

// defining fractions of the constant pi
constexpr double Tube::pi_64 = M_PI / 64;
constexpr double Tube::pi_32 = M_PI / 32;

// default class constructor
Tube::Tube()
{
	m_OD = 0.0;
	m_ID = 0.0;
	m_E = 0.0;
	m_I = 0.0;
	m_G = 0.0;
	m_J = 0.0;
	m_K = 0.0;
	m_ls = 0.0;
	m_lc = 0.0;
	m_u_ast = 0.0;
}

// overloaded class constructor
Tube::Tube(double OD, double ID, double E, double G, double ls, double lc, const vec3d &u_ast) : m_OD(OD), m_ID(ID), m_E(E), m_G(G), m_ls(ls), m_lc(lc), m_u_ast(u_ast)
{
	m_I = pi_64 * (pow(OD, 4) - pow(ID, 4));
	m_J = pi_32 * (pow(OD, 4) - pow(ID, 4));
	m_K(0, 0) = m_K(1, 1) = m_E * m_I;
	m_K(2, 2) = m_G * m_J;
}

// copy constructor
Tube::Tube(const Tube &rhs) : m_OD(rhs.m_OD), m_ID(rhs.m_ID), m_E(rhs.m_E), m_I(rhs.m_I), m_G(rhs.m_G),
							  m_J(rhs.m_J), m_K(rhs.m_K), m_ls(rhs.m_ls), m_lc(rhs.m_lc), m_u_ast(rhs.m_u_ast){};

// Tube destructor
Tube::~Tube()
{
	// nothing to be done
}

// move constructor
Tube::Tube(Tube &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_OD = rhs.m_OD;
		this->m_ID = rhs.m_ID;
		this->m_E = rhs.m_E;
		this->m_I = rhs.m_I;
		this->m_G = rhs.m_G;
		this->m_J = rhs.m_J;
		this->m_K = std::move(rhs.m_K);
		this->m_ls = rhs.m_ls;
		this->m_lc = rhs.m_lc;
		this->m_u_ast = std::move(rhs.m_u_ast);
	}
}

// Copy assignment operator
Tube &Tube::operator=(const Tube &rhs)
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_OD = rhs.m_OD;
		this->m_ID = rhs.m_ID;
		this->m_E = rhs.m_E;
		this->m_I = rhs.m_I;
		this->m_G = rhs.m_G;
		this->m_J = rhs.m_J;
		this->m_K = rhs.m_K;
		this->m_ls = rhs.m_ls;
		this->m_lc = rhs.m_lc;
		this->m_u_ast = rhs.m_u_ast;
	}

	return *this;
}

// move assignment operator
Tube &Tube::operator=(Tube &&rhs) noexcept
{
	// handling self assignment
	if (this != &rhs)
	{
		this->m_OD = rhs.m_OD;
		this->m_ID = rhs.m_ID;
		this->m_E = rhs.m_E;
		this->m_I = rhs.m_I;
		this->m_G = rhs.m_G;
		this->m_J = rhs.m_J;
		this->m_K = std::move(rhs.m_K);
		this->m_ls = rhs.m_ls;
		this->m_lc = rhs.m_lc;
		this->m_u_ast = std::move(rhs.m_u_ast);
	}

	return *this;
}

// // get method for retrieving the tube parameters
std::tuple<double, double, double, double, double, vec3d> Tube::getTubeParameters()
{
	return std::make_tuple(this->m_OD, this->m_ID, this->m_E, this->m_ls, this->m_lc, this->m_u_ast);
}

// get method for retrieving the tube overall length
double Tube::getTubeLength()
{
	return this->m_ls + this->m_lc;
}

// get method for retrieving the tube precurvature vector
vec3d Tube::get_u_ast()
{
	return this->m_u_ast;
}

// get method for retrieving the precurvature along x or y directions
double Tube::get_u_ast(const size_t id)
{
	switch (id)
	{
	case 1:
		return this->m_u_ast[0];
		break;
	case 2:
		return this->m_u_ast[1];
		break;
	default:
		return 0.0;
	}
}

// set method for updating the tube precurvature vector
void Tube::set_u_ast(const vec3d &u_ast)
{
	this->m_u_ast = u_ast;
}

// set method for updadting the tube's precurvature along x or y directions
void Tube::set_u_ast(const size_t id, const double u)
{
	switch (id)
	{
	case 1:
		this->m_u_ast[0] = u;
		break;
	case 2:
		this->m_u_ast[1] = u;
		break;
	default:
		std::cerr << "Invalid entry. Cannot update tube's precurvature!";
	}
}

// get method for retrieving the length of the straight section
double Tube::getStraightLen()
{
	return this->m_ls;
}

// get method for retrieving the length of the tube curved section
double Tube::getCurvLen()
{
	return this->m_lc;
}

// set method for updating the length of the straight section (if negative length is provided, the length will be automatically set to zero)
void Tube::setStraightLen(double ls)
{
	this->m_ls = (ls > 0) ? ls : 0.0;
}

// set method for updating the legnth of the curved section (if negative length is provided, the length will be automatically set to zero)
void Tube::setCurvLen(double lc)
{
	this->m_lc = (lc > 0) ? lc : 0.0;
}

// get method for retrieving the stiffness matrix
blaze::DiagonalMatrix<blaze::StaticMatrix<double, 3UL, 3UL, blaze::rowMajor>> Tube::getK_Matrix()
{
	return this->m_K;
}

// get method for retrieving the ith entry of the main diagonal
double Tube::getK(int i)
{
	switch (i)
	{
	case 1:
		return this->m_K(0, 0);
		break;
	case 2:
		return this->m_K(1, 1);
		break;
	case 3:
		return this->m_K(2, 2);
		break;
	default:
		return 0.0;
	}
}