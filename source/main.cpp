/*
Title: Matrix Mathematics
File Name: main.cpp
Copyright © 2016
Author: Andrew Litfin
Written under the supervision of David I. Schwartz, Ph.D., and
supported by a professional development seed grant from the B. Thomas
Golisano College of Computing & Information Sciences
(https://www.rit.edu/gccis) at the Rochester Institute of Technology.
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.
This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// The primary objects of study in linear algebra are matrices.
// This tutorial series will explore the applications of matrices to computer games and simulation,
//  especially in the realm of physical transformations.
// The exposition follows that of Eric Lengyel in "Foundations of Game Engine Development" (Volume 1).
// We have included the Vector structs from the previous series, and introduced Matrix structs that act similarly.
// These structs are based upon and largely follow code samples given in FGED.
//  As before, Matrix2D is heavily annotated, with other structs being annotated in places of difference.

// This tutorial delves deeper into projection.

#include "../header/Matrix4D.h"
#include "../header/tests.h"
#include "../header/helpers.h"

#include <iostream>
#include <ctime>

int main()
{
	// Required for some helper functions
	srand((unsigned)time(0));

	// Previously we examined the idea that vectors can be projected onto by using a matrix operation,
	//  where the matrix is constructed by taking the outer product of the vector we want to project onto with itself.
	// We proved this by showing that, if the vector we project onto is a unit vector, then Dot(b, a)*b == Outer(b, b)*a
	// Suppose we could extend this idea to more than one vector?

	// Recall that when a vector space is equipped with a basis, any element of the space can be uniquely written
	//  as a linear combination of the basis vectors.
	// This is also true for subspaces.

	// Suppose that u_1, u_2, ..., u_k is an orthonormal basis for a subspace U of Rn.
	// We then have by definition that every u in U can be written as a unique linear combination of u_1, u_2, ..., u_k.
	// To project onto this subspace we can take a similar approach to the 1D subspace case (projecting onto a single vector):
	//  each vector v can be projected onto the kD subspace as Sum_{i=1}^{k} Project(v, u_i) = Dot(u_i, v) * u_i = Outer(u_i, u_i) * v
	// Factoring out the v we arrive at the fact that to project onto any subspace spanned by an orthonormal basis,
	//  we just need to take the outer product of each basis vector and sum the results.

	// But we can get even faster than that!
	// It turns out that our friend the outer product which we initially defined only for vectors can be extended to matrices with similar meaning.

	// Let the matrix A have columns u_1, u_2, ..., u_k.
	// Then define P := A*Transpose(A).  P is then the projection matrix onto the subspace spanned by u_i.

	// To show that the two are equal, let's look at a low-dimensional example.
	// Suppose we are in R3 and want to project onto a 2D plane spanned by u and v, such that |u| = |v| = 1 and Dot(u, v) = 0.
	//  (That is, u and v form an orthonormal basis).
	// Define A as
	// A = [ u v ]
	// =>
	// [ u_1 v_1 ]
	// [ u_2 v_2 ]
	// [ u_3 v_3 ]
	// Then P = A * Transpose(A) = 
	// [ u_1^2+v_1^2  u_12+v_12  u_13+v_13 ]
	// [ u_12+v_12  u_2^2+v_2^2  u_23+v_23 ]
	// [ u_13+v_13  u_23+v_23  u_3^2+v_3^2 ]
	// = Outer(u, u) + Outer(v, v)

	// However, in the name of brevity at the expense of completeness, we don't have a mat2x3 or mat3x2 struct to work with.
	// (GLM, on the other hand, does, so if you're using GLM you needn't concern yourself with the following.)
	// So how can we make this calculation with nothing but square matrices?
	// The answer is we pad our matrix with zeroes.
	// So for example, in R3, we would instead have A = [ u v 0 ].

	// Consider the following:
	// Suppose we want to project onto the span of u and v:

	Vector3D u(1 / sqrtf(2), 1 / sqrtf(2), 0);
	Vector3D v(0, 0, 1);

	// Then

	Matrix3D A(u, v, Vector3D());

	// And

	Matrix3D P = A * Transpose(A);

	std::cout << "P =\n" << P;

	// !!REMEMBER!! The basis vectors are assumed to be unit length!
	// If they are not, then there will be distortion in the image of the matrix.
	// For example,

	Vector3D u_nonunit(1, 1, 0);
	Vector3D v_nonunit(0, 0, 2);
	Matrix3D A_nonunit(u_nonunit, v_nonunit, Vector3D());
	Matrix3D P_nonunit = A_nonunit * Transpose(A_nonunit);

	std::cout << "P_nonunit =\n" << P_nonunit;

	// To see the distortion, consider a = (0, 0, 1).
	// a is clearly already in the subspace spanned by u and v, so P is the identity on a.
	// (Indeed, P is the identity on its image (P^2 = P))

	std::cout << "P*(0, 0, 1) = " << P * Vector3D(0, 0, 1) << "\n";

	// However,

	std::cout << "P_nonunit*(0, 0, 1) = " << P_nonunit * Vector3D(0, 0, 1) << "\n";

	// A stretching by a factor of 4!
	// Clearly this is undesirable.

	// This can be fixed in one of two ways:
	// 1) instead of P = A*Transpose(A), define P' = A * Inverse(Transpose(A)*A) * Transpose(A)
	//    Clearly this is very computationally intensive and undesirable.
	//    IN ADDITION THIS WILL NOT WORK IF THE ZERO VECTOR IN INCLUDED AS IN THE "TRICK" ABOVE
	// 2) Normalize the vectors *before* constructing the projection matrix.

	std::cout << "\nPress Enter to exit . . . ";
	std::cin.get();
	return 0;
}
