

#ifndef ITERMETHOD_H
#define ITERMETHOD_H

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <vector>
#include <iostream>
#include <fstream>

//////////////////////////////////////////////////
//              Interface                       //
//////////////////////////////////////////////////

using namespace std;
using namespace Eigen;
typedef Eigen::SparseMatrix<double> SpMat;


//////////////////////////////////////////////////////////////////////
//                        GMRES related                             //
//////////////////////////////////////////////////////////////////////

/**
 * @brief obtain the cos and sin from coordinate (x, y)
 */
template <class Type>
void rotmat(const Type& x, const Type& y,
	Type* c, Type* s) {
	if (y == 0) {
		*c = 1;
		*s = 0;
	}
	else if (fabs(y) > fabs(x)) {
		double tmp = x / y;
		*s = 1.0 / sqrt(1.0 + tmp * tmp);
		*c = tmp * (*s);
	}
	else {
		double tmp = y / x;
		*c = 1.0 / sqrt(1.0 + tmp * tmp);
		*s = tmp * (*c);
	}
}

template <class Adotx>
std::tuple<VectorXd, std::vector<double>, int>
Gmres0(Adotx Ax, const VectorXd& b, const VectorXd& x0,
	const int restart,
	const int maxit, const double rtol) {

	/* initial setup */
	const int N = b.size();
	int M = restart;
	VectorXd x = x0;
	double bnrm2 = b.norm();
	std::vector<double> errVec;	/* error at each iteration */

	/* initialize workspace */
	MatrixXd V(MatrixXd::Zero(N, M + 1)); // orthogonal vectors in Arnoldi iteration
	MatrixXd H(MatrixXd::Zero(M + 1, M)); // Hessenberg matrix in Arnoldi iteration

	ArrayXd C(ArrayXd::Zero(M + 1));  // cos of Givens rotation
	ArrayXd S(ArrayXd::Zero(M + 1));  // sin of Givens rotation

	/* outer iteration */
	for (size_t iter = 0; iter < maxit; iter++) {
		/* obtain residule */
		VectorXd r = b - Ax(x);
		double rnorm = r.norm();
		double err = rnorm / bnrm2;
		if (err < rtol) return std::make_tuple(x, errVec, 0);

		V.col(0) = r / rnorm;	// obtain V_1

		VectorXd g(VectorXd::Zero(M + 1));// vector g = ||r|| * e1
		g(0) = rnorm;

		/* inner iteration : Arnoldi Iteration */
		for (size_t i = 0; i < M; i++) {
			// form the V_{i+1} and H(:, i)
			V.col(i + 1) = Ax(V.col(i));
			for (size_t j = 0; j <= i; j++) {
				H(j, i) = V.col(i + 1).dot(V.col(j));
				V.col(i + 1) -= H(j, i) * V.col(j);
			}
			H(i + 1, i) = V.col(i + 1).norm();
			if (H(i + 1, i) != 0) V.col(i + 1) /= H(i + 1, i);

			/* Givens Rotation
			 *  | c, s| * |x| = |r|
			 *  |-s, c|   |y|   |0|
			 */
			for (size_t j = 0; j < i; j++) {
				double tmp = C(j) * H(j, i) + S(j) * H(j + 1, i);
				H(j + 1, i) = -S(j) * H(j, i) + C(j) * H(j + 1, i);
				H(j, i) = tmp;
			}

			/* rotate the last element */
			rotmat(H(i, i), H(i + 1, i), &C(i), &S(i));
			H(i, i) = C(i) * H(i, i) + S(i) * H(i + 1, i);
			H(i + 1, i) = 0;

			/* rotate g(i) and g(i+1)
			 *  | c, s| * |g_i| = |c * g_i |
			 *  |-s, c|   | 0 |   |-s * g_i|
			 */
			g(i + 1) = -S(i) * g(i); // be careful about the order
			g(i) = C(i) * g(i);

			/* after all above operations, we obtain dimensions
			 *  H : [i+2, i+1]
			 *  V : [N,   i+2]
			 *  g : [i+2,   1]
			 *
			 * Here, it is better denote H as R since it is transfored
			 * into upper triangular form.
			 * residual = || \beta * e_1 - H * y|| = ||g - R*y||
			 * since the last row of R is zero, then residual =
			 * last element of g, namely g(i+1)
			 */
			double err = fabs(g(i + 1)) / bnrm2;
			errVec.push_back(err);
			if (err < rtol) {
				VectorXd y = H.topLeftCorner(i + 1, i + 1).lu().solve(g.head(i + 1));
				x += V.leftCols(i + 1) * y;
				return std::make_tuple(x, errVec, 0);
			}
		}

		// the inner loop finishes => has not converged
		// so need to update x, and go to outer loop
		VectorXd y = H.topLeftCorner(M, M).lu().solve(g.head(M));
		x += V.leftCols(M) * y;
	}

	// if the outer loop finished => not converged
	return std::make_tuple(x, errVec, 1);

}

/**
 * @brief a wrapper of Gmres0.
 *
 * This function uses an explicit form of matrix A.
 */
template <typename Mat>
std::tuple<VectorXd, std::vector<double>, int>
Gmres(const Mat& A, const VectorXd& b, const VectorXd& x0,
	const int restart,
	const int maxit, const double rtol) {

	return Gmres0([&A](VectorXd x) {return A * x; }, b, x0, restart, maxit, rtol);
}

#endif	/* ITERMETHOD_H */