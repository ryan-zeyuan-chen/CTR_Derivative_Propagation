#pragma once

#include <blaze/Math.h>
#include <iostream>

namespace mathOp
{
	// enumeration for defining which root finding method to be used
	enum class rootFindingMethod
	{
		NEWTON_RAPHSON,
		LEVENBERG_MARQUARDT,
		POWELL_DOG_LEG,
		MODIFIED_NEWTON_RAPHSON,
		BROYDEN,
		BROYDEN_II
	};

	// Function for converting from degrees to radians
	inline double deg2Rad(const double degree)
	{
		static constexpr double pi_180 = M_PI / 180.00;
		return (degree * pi_180);
	}

	// Function that reduces any angle to its congruent in the [-Pi, Pi) interval
	inline double wrapTo2Pi(const double theta)
	{
		static constexpr double TWO_PI = 2.00 * M_PI;
		double wrappedAngle;

		// wrappedAngle = fmod(theta, TWO_PI);
		wrappedAngle = remainder(theta, TWO_PI);

		// if (wrappedAngle < 0.0)
		// 	wrappedAngle += TWO_PI;

		return wrappedAngle;
	}

	// method for computing the Penrose Pseudoinverse via SVD decomposition
	inline blaze::HybridMatrix<double, 6UL, 6UL> pInv(const blaze::HybridMatrix<double, 6UL, 6UL> &M)
	{
		// declaring the auxiliary matrices for pInv computation
		blaze::HybridMatrix<double, 6UL, 6UL> U; // The matrix for the left singular vectors
		blaze::HybridVector<double, 6UL> s;		 // The vector for the singular values
		blaze::HybridMatrix<double, 6UL, 6UL> V; // The matrix for the right singular vectors

		// SVD decomposition of the matrix M
		try
		{
			blaze::svd(M, U, s, V);
		}
		catch (std::exception &e)
		{
			std::cerr << "Blaze SVD has failed: " << e.what() << std::endl;
			// std::cout << "CTR joints: " << blaze::trans(this->m_q) << std::endl;
			std::cerr << M << std::endl;
		}

		// computing the pseudoinverse of M via SVD decomposition
		blaze::DiagonalMatrix<blaze::HybridMatrix<double, 6UL, 6UL>> S_inv(s.size(), s.size());

		// Creating a reference to the diagonal of matrix S
		auto diag = blaze::diagonal(S_inv);
		// applies a "damping factor" to the zero singular values of the matrix M
		diag = blaze::map(s, [](double d)
						  { return d / ((d * d) + 1.00E-25); }); // Damped least squares -- SVD pseudo inverse

		return blaze::trans(U * S_inv * V);
	}

	inline void SO3_To_Quaternion(blaze::StaticVector<double, 4UL> &h, const blaze::StaticMatrix<double, 3UL, 3UL> &R)
	// inline void toQuaternion(blaze::StaticVector<double, 4UL> &h, const blaze::StaticMatrix<double, 3UL, 3UL> &R)
	{
		double n4;					 // the norm of quaternion multiplied by 4
		double tr = blaze::trace(R); // trace of martix

		if (tr > 0.00)
		{
			h[0UL] = tr + 1.00;
			h[1UL] = R(1UL, 2UL) - R(2UL, 1UL);
			h[2UL] = R(2UL, 0UL) - R(0UL, 2UL);
			h[3UL] = R(0UL, 1UL) - R(1UL, 0UL);
			n4 = h[0UL];
		}
		else if (R(0UL, 0UL) > std::max(R(1UL, 1UL), R(2UL, 2UL)))
		{
			h[0UL] = R(1UL, 2UL) - R(2UL, 1UL);
			h[1UL] = 1.00 + R(0UL, 0UL) - R(1UL, 1UL) - R(2UL, 2UL);
			h[2UL] = R(1UL, 0UL) + R(0UL, 1UL);
			h[3UL] = R(2UL, 0UL) + R(0UL, 2UL);
			n4 = h[1UL];
		}
		else if (R(1UL, 1UL) > R(2UL, 2UL))
		{
			h[0UL] = R(2UL, 0UL) - R(0UL, 2UL);
			h[1UL] = R(1UL, 0UL) + R(0UL, 1UL);
			h[2UL] = 1.00 + R(1UL, 1UL) - R(0UL, 0UL) - R(2UL, 2UL);
			h[3UL] = R(2UL, 1UL) + R(1UL, 2UL);
			n4 = h[2UL];
		}
		else
		{
			h[0UL] = R(0UL, 1UL) - R(1UL, 0UL);
			h[1UL] = R(2UL, 0UL) + R(0UL, 2UL);
			h[2UL] = R(2UL, 1UL) + R(1UL, 2UL);
			h[3UL] = 1.00 + R(2UL, 2UL) - R(0UL, 0UL) - R(1UL, 1UL);
			n4 = h[3UL];
		}

		h *= 1.00 / (2.00 * sqrt(n4)); // R(UL, UL)
	}

	inline void euler2Quaternion(const double heading, const double attitude, const double bank, blaze::StaticVector<double, 4UL> &h)
	{
		// Gotta convert the angles to radians first.
		double theta(0.50 * heading), phi(0.50 * attitude), psi(0.50 * bank);
		double c1(cos(theta)), s1(sin(theta)), c2(cos(phi)), s2(sin(phi)), c3(cos(psi)), s3(sin(psi)), c1c2, s1s2;

		c1c2 = c1 * c2;
		s1s2 = s1 * s2;

		h[0UL] = c1c2 * c3 - s1s2 * s3;
		h[1UL] = c1c2 * s3 + s1s2 * c3;
		h[2UL] = s1 * c2 * c3 + c1 * s2 * s3;
		h[3UL] = c1 * s2 * c3 - s1 * c2 * s3;
	}

	// function that returns a rotation matrix in SO(3) from a set of non-unity quaternions
	inline void getSO3(const blaze::StaticVector<double, 4UL> &h, blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> &R)
	{
		blaze::IdentityMatrix<double> I(3UL);

		R = {{-h[2UL] * h[2UL] - h[3UL] * h[3UL], h[1UL] * h[2UL] - h[3UL] * h[0UL], h[1UL] * h[3UL] + h[2UL] * h[0UL]},
			 {h[1UL] * h[2UL] + h[3UL] * h[0UL], -h[1UL] * h[1UL] - h[3UL] * h[3UL], h[2UL] * h[3UL] - h[1UL] * h[0UL]},
			 {h[1UL] * h[3UL] - h[2UL] * h[0UL], h[2UL] * h[3UL] + h[1UL] * h[0UL], -h[1UL] * h[1UL] - h[2UL] * h[2UL]}};

		R *= 2.00 / blaze::sqrNorm(h);
		R += I;
	}

	// function that returns the rotation matrix Rz of any angle theta
	inline blaze::StaticMatrix<double, 3UL, 3UL> rotz(const double &theta)
	{
		blaze::StaticMatrix<double, 3UL, 3UL, blaze::columnMajor> R;
		double c(cos(theta)), s(sin(theta));

		R = {{c, -s, 0.00},
			 {s, c, 0.00},
			 {0.00, 0.00, 1.00}};

		return R;
	}

	// function that returns the derivative of the rotation matrix Rz
	inline blaze::StaticMatrix<double, 3UL, 3UL> rotz_dot_transpose(const double &theta)
	{
		blaze::StaticMatrix<double, 3UL, 3UL> dR;
		double c(cos(theta)), s(sin(theta));

		dR = {{-s, c, 0.00},
			  {-c, -s, 0.00},
			  {0.00, 0.00, 0.00}};

		return dR;
	}

	// function that returns the derivative of the rotation matrix Rz
	inline blaze::StaticMatrix<double, 3UL, 3UL> dRz_dTheta(const double &theta)
	{
		blaze::StaticMatrix<double, 3UL, 3UL> dR;
		double c(cos(theta)), s(sin(theta));

		dR = {{-s, -c, 0.00},
			  {c, -s, 0.00},
			  {0.00, 0.00, 0.00}};

		return dR;
	}

	// function that implements the wedge operator (reshapes "flattens out" a matrix into a vector)
	inline blaze::StaticVector<double, 6UL> wedgeOperator(blaze::StaticMatrix<double, 4UL, 4UL> M)
	{
		// M =  [\hat{w} v;    ==> \wedge{M} = [v; w]
		//        0 0] }

		blaze::StaticVector<double, 6UL> v = {M(0UL, 3UL), M(1UL, 3UL), M(2UL, 3UL), M(2UL, 1UL), M(0UL, 2UL), M(1UL, 0UL)};

		return v;
	}

	// function that implements the wedge operator (reshapes "flattens out" a matrix into a vector)
	inline blaze::StaticMatrix<double, 3UL, 3UL> hatOperator(const blaze::StaticVector<double, 3UL> &v)
	{
		blaze::StaticMatrix<double, 3UL, 3UL> hat;

		hat = {{0.00, -v[2UL], v[1UL]},
			   {v[2UL], 0.00, -v[0UL]},
			   {-v[1UL], v[0UL], 0.00}};

		return hat;
	}

	// function that computes the premultiplication of a matrix M by v^, i.e., v^ * M
	inline blaze::StaticMatrix<double, 3UL, 3UL> hatPreMultiply(const blaze::StaticVector<double, 3UL> &v, const blaze::StaticMatrix<double, 3UL, 3UL> &M)
	{
		blaze::StaticMatrix<double, 3UL, 3UL> Res = {{-M(1UL, 0UL) * v[2UL] + M(2UL, 0UL) * v[1UL], -M(1UL, 1UL) * v[2UL] + M(2UL, 1UL) * v[1UL], -M(1UL, 2UL) * v[2UL] + M(2UL, 2UL) * v[1UL]},
													 {M(0UL, 0UL) * v[2UL] - M(2UL, 0UL) * v[0UL], M(0UL, 1UL) * v[2UL] - M(2UL, 1UL) * v[0UL], M(0UL, 2UL) * v[2UL] - M(2UL, 2UL) * v[0UL]},
													 {-M(0UL, 0UL) * v[1UL] + M(1UL, 0UL) * v[0UL], -M(0UL, 1UL) * v[1UL] + M(1UL, 1UL) * v[0UL], -M(0UL, 2UL) * v[1UL] + M(1UL, 2UL) * v[0UL]}};

		return Res;
	}

	// function that computes the posmultiplication of a matrix M by v^, i.e., M * v^
	inline blaze::StaticMatrix<double, 3UL, 3UL> hatPostMultiply(const blaze::StaticMatrix<double, 3UL, 3UL> &M, const blaze::StaticVector<double, 3UL> &v)
	{
		blaze::StaticMatrix<double, 3UL, 3UL> Res = {{M(0UL, 1UL) * v[2UL] - M(0UL, 2UL) * v[1UL], -M(0UL, 0UL) * v[2UL] + M(0UL, 2UL) * v[0UL], M(0UL, 0UL) * v[1UL] - M(0UL, 1UL) * v[0UL]},
													 {M(1UL, 1UL) * v[2UL] - M(1UL, 2UL) * v[1UL], -M(1UL, 0UL) * v[2UL] + M(1UL, 2UL) * v[0UL], M(1UL, 0UL) * v[1UL] - M(1UL, 1UL) * v[0UL]},
													 {M(2UL, 1UL) * v[2UL] - M(2UL, 2UL) * v[1UL], -M(2UL, 0UL) * v[2UL] + M(2UL, 2UL) * v[0UL], M(2UL, 0UL) * v[1UL] - M(2UL, 1UL) * v[0UL]}};

		return Res;
	}

	// function that computes the differential quaternion evolution
	inline blaze::StaticVector<double, 4UL> quaternionDiff(const blaze::StaticVector<double, 3UL> &u, const blaze::StaticVector<double, 4UL> &h)
	{
		blaze::StaticVector<double, 4UL> hs;

		hs[0UL] = 0.50 * (-u[0UL] * h[1UL] - u[1UL] * h[2UL] - u[2UL] * h[3UL]);
		hs[1UL] = 0.50 * (u[0UL] * h[0UL] + u[2UL] * h[2UL] - u[1UL] * h[3UL]);
		hs[2UL] = 0.50 * (u[1UL] * h[0UL] - u[2UL] * h[1UL] + u[0UL] * h[3UL]);
		hs[3UL] = 0.50 * (u[2UL] * h[0UL] + u[1UL] * h[1UL] - u[0UL] * h[2UL]);

		return hs;
	}
}