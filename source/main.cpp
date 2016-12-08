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

	////////////////////////////
	// Projection as a matrix //
	////////////////////////////
	{
		/*
		In [math-vector-projection], we examined the act of projecting one vector onto another.
		We showed that proj is a function taking two vectors and returning another, defined as
		 proj_b(a) = (Dot(a, b) / Dot(b, b)) * b
		
		Suppose we want to project MANY vectors onto a single vector.
		In particular, could we represent it as a matrix?
		The answer is yes!
		
		First we need to introduce yes another vector operation called the Outer product.
		 (As opposed to the Inner product (dot product)).
		Let u be an m by 1 column vector and v be an n by 1 column vector.
		Then Outer(u, v) := u * Transpose(v), yielding an m by n matrix where the (i, j) element equals u_i * v_j
		
		In addition, recall that Dot(u, v) can be realized as Transpose(u) * v.
		
		Then,
		proj_b(a) = ((Dot(a, b) / Dot(b, b)) * b
		= 1/|b|^2 * Dot(a, b) * b (Dot(b, b) = 1/|b|^2)
		= 1/|b|^2 * b * Dot(a, b) (Scalar multiplication is commutative in vectors over a field)
		= 1/|b|^2 * b * Dot(b, a) (Dot(a, b) = Dot(b, a) when the underlying field is the reals)
		= 1/|b|^2 * b * (Transpose(b) * a) (Dot(a, b) = Transpose(a) * b)
		= 1/|b|^2 * (b * Transpose(b)) * a (Associativity of matrix multiplication)
		= 1/|b|^2 * Outer(b, b) * a (Definition of Outer product)
		
		Then we can realize proj_b as a matrix by itself that equals Outer(b, b) / Dot(b, b) !
		
		This has the benefit that we can project many vectors onto a single vector simply with matrix multiplication.
		*/

		Vector3D b(randFloat(0, 1), randFloat(0, 1), randFloat(0, 1));
		std::cout << "b = " << b << "\n";

		Matrix3D proj_b = MakeProjection(b);
		std::cout << "proj_b =\n" << proj_b << "\n";

		for (int i = 0; i < 5; i++)
		{
			Vector3D a(randFloat(0, 1), randFloat(0, 1), randFloat(0, 1));
			std::cout << "a = " << a << "\nproj_b(a) (matrix) = " << proj_b * a << "\nproj_b(a) (vector) = " << Project(a, b) << "\n";
			if (MagSquared((proj_b * a) - Project(a, b)) < FLT_EPSILON)
			{
				std::cout << "The two formulations are equal!\n\n";
			}
		}

		// Note, however, that when projecting just one vector onto another vector, the dot product version is actually faster.

		/*
		Projection as currying (optional):

		When considering proj as a function, the vector formulation takes in to vectors and returns another.
		Formally, we write this as proj : V x V -> V
		We said earlier that a matrix represents a function between vector spaces.
		Formally, if A is an m by n matrix over the reals, then A : Rn -> Rm.
		So the matrix formulation of proj, being that it returns a matrix, in essence returns a function.
		Formally, proj : V -> (V -> V).
		Note that the number of parameters has decreased by one but our function returns a function that takes one parameter.
		This is an example of a process called "currying," and if you ever want to study functional programming, you'll need to know it.
		*/
	}
	std::cout << "Press Enter to continue . . .\n";
	std::cin.get();

	////////////////////////////////
	// Projection onto a subspace //
	////////////////////////////////
	{
		/*
		Previously we examined the idea that vectors can be projected onto by using a matrix operation,
		Suppose we could extend this idea to more than one vector?

		Recall that when a vector space is equipped with a basis, any element of the space can be uniquely written
		 as a linear combination of the basis vectors.
		This is also true for subspaces.

		Suppose that u_1, u_2, ..., u_k is an orthonormal basis for a subspace U of Rn.
		We then have by definition that every u in U can be written as a unique linear combination of u_1, u_2, ..., u_k.
		To project onto this subspace we can take a similar approach to the 1D subspace case (projecting onto a single vector):
		 each vector v can be projected onto the kD subspace as Sum_{i=1}^{k} Project(v, u_i) = Dot(u_i, v) * u_i = Outer(u_i, u_i) * v
		Factoring out the v we arrive at the fact that to project onto any subspace spanned by an orthonormal basis,
		 we just need to take the outer product of each basis vector and sum the results.

		But we can get even faster than that!
		It turns out that our friend the outer product which we initially defined only for vectors can be extended to matrices with similar meaning.

		Let the matrix A have columns u_1, u_2, ..., u_k.
		Then define P := A*Transpose(A).  P is then the projection matrix onto the subspace spanned by u_1, u_2, ..., u_k.

		To show that the two are equal, let's look at a low-dimensional example.
		Suppose we are in R3 and want to project onto a 2D plane spanned by u and v, such that |u| = |v| = 1 and Dot(u, v) = 0.
		 (That is, u and v form an orthonormal basis).
		Define A as
		A = [ u v ]
		=>
		[ u_1 v_1 ]
		[ u_2 v_2 ]
		[ u_3 v_3 ]
		Then P = A * Transpose(A) = 
		[ u_1^2+v_1^2  u_12+v_12  u_13+v_13 ]
		[ u_12+v_12  u_2^2+v_2^2  u_23+v_23 ]
		[ u_13+v_13  u_23+v_23  u_3^2+v_3^2 ]
		= Outer(u, u) + Outer(v, v)

		However, in the name of brevity at the expense of completeness, we don't have a mat2x3 or mat3x2 struct to work with.
		(GLM, on the other hand, does, so if you're using GLM you needn't concern yourself with the following.)
		So how can we make this calculation with nothing but square matrices?
		The answer is we pad our matrix with zeroes.
		So for example, in R3, we would instead have A = [ u v 0 ].
		*/

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

		std::cout << "P_nonunit*(0, 0, 1) = " << P_nonunit * Vector3D(0, 0, 1) << "\n\n";

		// A stretching by a factor of 4!
		// Clearly this is undesirable.

		// This can be fixed in one of two ways:
		// 1) instead of P = A*Transpose(A), define P' = A * Inverse(Transpose(A)*A) * Transpose(A)
		//    Clearly this is very computationally intensive and undesirable.
		//    IN ADDITION THIS WILL NOT WORK IF THE ZERO VECTOR IN INCLUDED AS IN THE "TRICK" ABOVE
		// 2) Normalize the vectors *before* constructing the projection matrix.
	}
	std::cout << "Press Enter to continue . . .\n";
	std::cin.get();

	//////////////////////////
	// Rejection revisisted //
	//////////////////////////
	{
		/*
		Similarly to projection as above, recall that we have from vector rejection that
		reject_b(a) = a - proj_b(a).

		We can pull a similar trick by inserting the identity matrix before a (multiplication by the identity doesn't change the result).
		Then
		reject_b(a) = I*a - proj_b(a)
		
		Factoring out the a yields
		reject_b(a) = (I - proj_b) * a

		So just as with projection, we now have a way to represent rejection as a matrix!
		*/

		Vector3D b(randFloat(0, 1), randFloat(0, 1), randFloat(0, 1));
		std::cout << "b = " << b << "\n";

		Matrix3D reject_b = MakeRejection(b);
		std::cout << "reject_b =\n" << reject_b << "\n";

		for (int i = 0; i < 5; i++)
		{
			Vector3D a(randFloat(0, 1), randFloat(0, 1), randFloat(0, 1));
			std::cout << "a = " << a << "\nreject_b(a) (matrix) = " << reject_b * a << "\nreject_b(a) (vector) = " << Reject(a, b) << "\n";
			if (MagSquared((reject_b * a) - Reject(a, b)) < FLT_EPSILON)
			{
				std::cout << "The two formulations are equal!\n\n";
			}
		}
	}

	std::cout << "\nPress Enter to exit . . . ";
	std::cin.get();
	return 0;
}
