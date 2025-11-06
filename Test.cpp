#include <iostream>
#include <Eigen/Dense>   // Include Eigen library

using namespace std;
using namespace Eigen;

int main() {
    // Create 2D vectors
    Vector2d a(2.0, 3.0);
    Vector2d b(4.0, 1.0);

    // Vector operations
    cout << "a + b = " << a + b << endl;
    cout << "a • b = " << a.dot(b) << endl;
    cout << "‖a‖ = " << a.norm() << endl;

    // Create a 2x2 matrix
    Matrix2d M;
    M << 1, 2,
         3, 4;

    // Matrix-vector multiplication
    Vector2d result = M * a;

    cout << "M * a = " << result.transpose() << endl;
    cout << "Determinant of M = " << M.determinant() << endl;

    return 0;
}
