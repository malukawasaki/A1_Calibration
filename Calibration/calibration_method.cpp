/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "calibration.h"
#include "matrix_algo.h"


using namespace easy3d;



/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 * @return True on success, otherwise false. On success, the camera parameters are returned by
 */


 /*The camera calibration function calibration() takes in two inputs: points_3d and points_2d, which are arrays of
 3D points and corresponding 2D image points, respectively.

 The function outputs the camera intrinsic parameters
 (focal length fx and fy, principal point coordinates cx and cy, and skew factor skew),
  as well as the camera's rotation matrix R and translation vector t.*/
bool Calibration::calibration(
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx,    /// output: focal length (i.e., K[0][0], which is equal to 'alpha' in our slides).
        double& fy,    /// output: focal length (i.e., K[1][1], which is equal to 'beta/sin(theta)' in our slides).
        double& cx,    /// output: x component of the principal point (i.e., K[0][2], which is 'u0' in our slides).
        double& cy,    /// output: y component of the principal point (i.e., K[1][2], which is 'v0' in our slides).
        double& skew,  /// output: skew factor (i.e., K[0][1], which is equal to '-alpha * cot(theta)' in our slides).
        Matrix33& R,   /// output: the 3x3 rotation matrix encoding camera rotation.
        Vector3D& t)   /// output：a 3D vector encoding camera translation.
{
    std::cout << "\nTODO: I am going to implement the calibration() function in the following file:\n"
                 "\t    - calibration_method.cpp\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tCamera calibration requires computing the SVD and inverse of matrices.\n"
                 "\tIn this assignment, I provide you with a 'Matrix' and a 'Vector' data structures for storing and\n"
                 "\tmanipulating matrices and vectors of arbitrary sizes. I also wrote some code to show you how to:\n"
                 "\t    - compute the SVD of a matrix;\n"
                 "\t    - compute the inverse of a matrix;\n"
                 "\t    - compute the transpose of a matrix.\n\n"
                 "\tFeel free to use any of the provided data structures and functions. The commonly used linear algebra\n"
                 "\tfunctions are provided in the following files:\n"
                 "\t    - Calibration/matrix.h  Matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Calibration/vector.h  Vectors of arbitrary dimensions and related functions.\n"
                 "\t    - Calibration/matrix_algo.h  Determinant, inverse, SVD, linear least-squares...\n"
                 "\tPlease refer to the above files for a complete list of useful functions and their usage.\n\n"
                 "\tIn your final submission, please\n"
                 "\t    - delete ALL unrelated test or debug code and avoid unnecessary output.\n"
                 "\t    - include all the source code (and please do NOT modify the structure of the directories).\n"
                 "\t    - do NOT include the 'build' directory (which contains the intermediate files in a build step).\n"
                 "\t    - make sure your code compiles and can reproduce your results without ANY modification.\n\n" << std::flush;



    std::cout << "\n[Liangliang]:\n"
                 "\tThe input parameters of this function are:\n"
                 "\t\t- points_3d: An array of 3D points (input to this function)\n"
                 "\t\t- points_2d: An array of 2D image points (input to this function)\n"
                 "\tThis function must return either 'true' on success or 'false' otherwise. On success, the camera\n"
                 "\tparameters are returned by the following variables:\n"
                 "\t\t- fx and fy: the focal lengths (in our slides, we use 'alpha' and 'beta')\n"
                 "\t\t- cx and cy: the principal point (in our slides, we use 'u0' and 'v0')\n"
                 "\t\t- skew:      the skew factor ('-alpha * cot_theta')\n"
                 "\t\t- R:         the 3x3 rotation matrix encoding camera orientation\n"
                 "\t\t- t:         a 3D vector encoding camera location.\n"
                 "\tIMPORTANT: don't forget to write your recovered parameters to the above variables." << std::endl;

    // TODO: check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)

    // TODO: construct the P matrix (so P * m = 0).

    // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.

    //Creating an example Matrix called A to work on.
    //TODO: SIMAY: change this matrix with the real one later.
    //Matrix A is going to be size 2Nx12, where N is the number of points that you have

     const int m = 6, n = 12; /*must change "m" with the number of points we have*/
     Matrix A(m, n, array.data());    // 'array.data()' returns a pointer to the array.
     // Matrix CHECK
     //  std::cout << "M: \n" << A << std::endl;

     // matrix-vector product
     Vector3D v = A * Vector4D(1, 2, 3, 4); // A is m to n matrix, v is n-dimensional vector, so the result is m-dimensional vector.
     // "m" is the number of rows of A, "n" is the number of columns of A.
     // U and V are orthogonal matrices, S is a diagonal matrix.
     Matrix U(m, m, 0.0);   // initialized with 0s
     Matrix S(m, n, 0.0);   // initialized with 0s
     Matrix V(n, n, 0.0);   // initialized with 0s

     // Compute the SVD decomposition of A
     // U , S , V^T =SVD(A)
     svd_decompose(A, U, S, V);

    //Now let's CHECK if the SVD result is correct.
         // Check 1: U is orthogonal, so U * U^T must be identity
            std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;

         // Check 2: V is orthogonal, so V * V^T must be identity
            std::cout << "V*V^T: \n" << V * transpose(V) << std::endl;

         // Check 3: S must be a diagonal matrix
            std::cout << "S: \n" << S << std::endl;

         // Check 4: according to the definition, A = U * S * V^T
            std::cout << "M - U * S * V^T: \n" << A - U * S * transpose(V) << std::endl;

     // Compute the INVERSE of Matrix A
     Matrix invA;
     inverse(A, invA); // Inverse of Matrix A is stored in invA.

         // Create an IDENTİTY MATRIX with the same dimensions as A: B
         Matrix B(m, n);
         for (int i = 0; i < m; i++) {
             for (int j = 0; j < n; j++) {
                 B(i, j) = (i == j) ? 1 : 0;
             }
         }

     // Let's CHECK if the inverse is correct
     std::cout << "Check if the inverse is correct: A * invA: \n" << B * invA << std::endl; //

     // Camera Matrix(M) is the last column of V, and you need to reshape to 3x4
        Matrix M(3, 4);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 4; j++) {
                M(i, j) = V(i, 11);
            }
        }

        //Let's CHECK if the camera matrix is correct
        std::cout << "Camera Matrix(M): \n" << M << std::endl;



    // TODO: extract intrinsic parameters from M.

    // TODO: extract extrinsic parameters from M.

    std::cout << "\n\tTODO: After you implement this function, please return 'true' - this will trigger the viewer to\n"
                 "\t\tupdate the rendering using your recovered camera parameters. This can help you to visually check\n"
                 "\t\tif your calibration is successful or not.\n\n" << std::flush;
    return false;
}
