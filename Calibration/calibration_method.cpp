/**
 * GEO10016 Assignment 1: Camera Calibration
 * Group 06: Maria Luisa Tarozzo Kawasaki (5620341), Simay Batum (5715598), Rianne Aalders (4593987)
 *
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * For GNU General Public License, see <http://www.gnu.org/licenses/>.
 */

#include "calibration.h"
#include "matrix_algo.h"

using namespace easy3d;

/**
 * @return True on success, otherwise false. On success, the camera parameters are returned by
 */
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

    /// Check if input is valid (e.g. sizes of 2D/3D points must match and be greater or equal to 6)
    if (points_3d.size() != points_2d.size() || points_2d.size() < 6) {
        std::cout << "Sizes of 2D/3D points do not match or are smaller than 6. This operation is not possible" << std::endl;
        return false;
    }

    else {
        std::cout << "Sizes of 2D/3D points match and are greater or equal to 6. This operation is possible"
                  << std::endl;

        /// Construct the P matrix (so P * m = 0).
        Matrix P (2*points_3d.size(), 12, 0.0);

        for (int i = 0; i < points_2d.size(); i++) {
           P.set_row(2*i,{points_3d[i].x(),points_3d[i].y(),points_3d[i].z(),1.,0.,0.,0.,0.,-(points_2d[i].x()*points_3d[i].x()),-(points_2d[i].x()*points_3d[i].y()),-(points_2d[i].x()*points_3d[i].z()),- points_2d[i].x()});
           P.set_row(2*i+1,{0.,0.,0.,0.,points_3d[i].x(),points_3d[i].y(),points_3d[i].z(),1.,-(points_2d[i].y()*points_3d[i].x()),-(points_2d[i].y()*points_3d[i].y()),-(points_2d[i].y()*points_3d[i].z()),- points_2d[i].y()});
       }
        //Check if the P matrix is correct:
        // std::cout << "P Matrix: \n" << P << std::endl;

       //P * M = 0
       // U and V are orthogonal matrices, S is a diagonal matrix.

        int n_cols = P.cols();
        int n_rows = P.rows();

        /// Solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
        Matrix U(n_rows, n_rows);   // initialized with 0s
        Matrix S(n_rows, n_cols);   // initialized with 0s
        Matrix V(n_cols, n_cols);   // initialized with 0s

        Matrix M (3., 4.);

        svd_decompose(P, U, S, V);

        // Populating M with the last column of V.
        int col_index = 11;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 4; j++) {
                M(i,j) = V(i*4 + j, col_index);
            }
        }

        //Check if the camera matrix is correct:
        // std::cout << "Camera Matrix(M): \n" << M << std::endl;

        /// Extract intrinsic parameters from M.

       // pi = K [R t] Pi = MPi

       // Mu = [A b] where A and b values corresponds to the ones in M.
       // Values in A correspond to the three first columns in M.
        Vector3D a1(M[0][0],M[0][1],M[0][2]);
        Vector3D a2(M[1][0],M[1][1],M[1][2]);
        Vector3D a3(M[2][0],M[2][1],M[2][2]);

        //Calculate intrinsic parameters
        double Ro_positive = 1 / a3.length();
        std::cout << Ro_positive << std::endl;
        double Ro_negative = - 1 / a3.length();
        cx = pow(Ro_positive,2)*(dot(a1, a3));
        cy = pow(Ro_positive, 2)*(dot(a2, a3));
        double cos_theta = dot((cross(a1, a3)),(cross(a2,a3))) / (cross(a1,a3).length() * cross(a2,a3).length());
        double sin_theta = sqrt(1 - pow(cos_theta,2));

        fx = pow(Ro_positive,2) * cross(a1,a3).length() * sin_theta;
        double beta = pow(Ro_positive,2) * cross(a2,a3).length() * sin_theta;

        skew = - fx * (cos_theta/sin_theta);
        fy = beta / sin_theta;

        Matrix K(3, 3);
        K[0][0] = fx;
        K[0][1] = skew;
        K[0][2] = cx;
        K[1][1] = fy;
        K[1][2] = cy;
        K[2][2] = 1.0;

        //Check if the K matrix is correct
        // std::cout << "K matrix is:" << K << std::endl;

        /// Extract extrinsic parameters from M.

        // Values in b correspond to the last column in M.
        Vector3D b = {M[0][3], M[1][3], M[2][3]};

        //Dealing with Ro positive or negative?
        double Ro;
        Matrix invK;
        inverse(K, invK);
        Vector3D t_Ropos = Ro_positive * invK * b;
        Vector3D t_Roneg = Ro_negative * invK * b;
        Vector3D z = {0,0,1};

        if (t_Ropos[2] < 0){
            t = t_Roneg;
            Ro = Ro_negative;
        }
        else {
            t = t_Ropos;
            Ro = Ro_positive;
        }

        Vector3D r1 = cross(a2,a3)/cross(a2,a3).length();
        Vector3D r3 = Ro * a3;
        Vector3D r2 = cross(r3,r1);

        // Populating matrix R
        R.set_row(0,r1);
        R.set_row(1,r2);
        R.set_row(2,r3);

        return true;
        }
}


















