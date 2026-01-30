#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

/** * Solution of the symmetric band equation system by LDU method. 
 * The matrix is in column format.
 * * Direct port of the provided LUband.java code.
 */
class LUband {
public:
    /** * UtDU decomposition for symmetric band matrix.
     * @param a matrix by columns of constant height
     * @param n number of equations
     * @param h bandwidth (half-bandwidth including diagonal)
     */
    void decompose(std::vector<double>& a, int n, int h) {
        int i, j, k;
        int h1 = h - 1;
        double w;

        for (j = 1; j < n; j++) {
            for (i = std::max(j - h1, 0); i < j; i++) {
                for (k = std::max(j - h1, 0); k < i; k++) {
                    a[i + h1 * (j + 1)] -= a[k + h1 * (i + 1)] * a[k + h1 * (j + 1)];
                }
            }
            for (i = std::max(j - h1, 0); i < j; i++) {
                w = a[i + h1 * (j + 1)];
                a[i + h1 * (j + 1)] /= a[i + h1 * (i + 1)];
                a[j + h1 * (j + 1)] -= a[i + h1 * (j + 1)] * w;
            }
        }
    }

    /** * Reduction and backsubstitution for RHS.
     * @param a decomposed matrix by columns of constant height
     * @param n number of equations
     * @param h bandwidth
     * @param b right-hand side(in)/solution vector(out)
     */
    void solve(const std::vector<double>& a, int n, int h, std::vector<double>& b) {
        int i, j;
        int h1 = h - 1;

        // Forward elimination
        for (j = 1; j < n; j++) {
            for (i = std::max(j - h1, 0); i < j; i++) {
                b[j] -= a[i + h1 * (j + 1)] * b[i];
            }
        }

        // Diagonal scaling
        for (j = 0; j < n; j++) {
            b[j] /= a[j + h1 * (j + 1)];
        }

        // Back substitution
        for (j = n - 1; j >= 0; j--) {
            for (i = std::max(j - h1, 0); i < j; i++) {
                b[i] -= a[i + h1 * (j + 1)] * b[j];
            }
        }
    }
};