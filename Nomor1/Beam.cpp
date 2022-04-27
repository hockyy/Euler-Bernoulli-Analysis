#include <bits/stdc++.h>
using namespace std;
#define fi first
#define se second
#define pb push_back
#define trav(a, x) for(auto& a : x)
const double EPS = 1e-10;
inline bool eq(double x, double y) { return fabs(x-y) < EPS; }
inline bool le(double x, double y) { return x < y + EPS; }
inline bool lt(double x, double y) { return x + EPS < y; }
const int MAGIC = 7;

inline void shiftRight(vector <double> &a, int amount){
    reverse(a.begin(), a.end());
    a.resize(MAGIC+amount);
    reverse(a.begin(), a.end());
    a.resize(MAGIC);
}

inline void shiftLeft(vector <double> &a, int amount){
    reverse(a.begin(), a.end());
    a.resize(MAGIC-amount);
    reverse(a.begin(), a.end());
    a.resize(MAGIC);
}

struct Beam{
    vector<vector<double>> A;
    vector<double> b;
    int n;
    Beam(int _n): n(_n){
        A.pb({0, 0, 16, -9, 8.0/3, -1.0/4, 0});
        A.pb({0, -4, 6, -4, 1, 0, 0});
        for(int i = 2;i < n - 2;i++) A.pb({1, -4, 6, -4, 1, 0, 0});
        A.pb({16.0/17, -60.0/17, 72.0/17, -28.0/17, 0, 0, 0});
        b = getB();
        // Baris terakhir tentunya lebih satu kolom sebelah kiri, mesti dinormalisasi dulu dengan baris ke-(n-1)
        // lakukan A(n) = A(n) - A(n-1) * -0.75;
        // Nanti didapatkan baris terakhir ialah [0, 51/17, -102/17, 51/17].
        // B(n) = B(n) - B(n-1) * -0.75 juga.
        // Zero based indexing
        A.pb({3, -6, 3, 0, 0, 0, 0});
        b[n-1] -= b[n-2] * -0.75;
    }

    vector <double> getB(){
        double L = 2;
        double h = L/n;
        double w = 0.3;
        double d = 0.03;
        double E = 1.3e10;
        double I = w * pow(d, 3) / 12;
        double g = 9.81;
        double f = (-480 * w * d * g);
        f = pow(h, 4)/(E * I) * f;
        return vector<double>(n, f);
    }

    void printA(){
        cout << "____________________________________" << endl;
        trav(cur, A){
            trav(tmp, cur) cout << tmp << " ";
            cout << endl;
        }
    }

    void printB(){
        trav(cur, b){
            cout << cur << endl;
        }
    }

    vector <double> gaussian(){
        // Lakukan gaussian elimination pada compact storage n × 7 ini
        // setiap baris indexing relatif kolomnya ialah -2. Perhatikan juga 0-based
        for(int j = 0;j < n;j++){
            // nolkan semua yang berada tepat dibawahnya, lebar band hanya 2
            // Cari posisi maksimum terlebih dahulu.
            // Simpan pair antara baris dan kolomnya yang relatif
            pair<int, int> maxRow = {j, 2};
            if(j + 1 < n && fabs(A[maxRow.fi][maxRow.se]) < fabs(A[j+1][1])) maxRow = {j+1, 1};
            if(j + 2 < n && fabs(A[maxRow.fi][maxRow.se]) < fabs(A[j+2][0])) maxRow = {j+2, 0};

            // Sesuaikan shift compact storage, yaitu sebesar 2-maxRow.se
            // Swap rownya
            if(maxRow.se != 2){
                shiftLeft(A[j], 2-maxRow.se);
                shiftRight(A[maxRow.fi], 2-maxRow.se);
                swap(A[j], A[maxRow.fi]);
                swap(b[j], b[maxRow.fi]);
            }

            // Eliminasi 2 baris di bawah
            for(int i = 1;i < 3 && j + i < n;i++){
                double m = A[j + i][2-i]/A[j][2];
                for(int k = 0;k+i < MAGIC;k++){
                    A[j+i][k] -= m * A[j][k+i];
                }

                b[j+i] -= m * b[j];
            }
        }

        // Do backward substitution
        for(int i = n - 1;i >= 0;i--){
            for(int j = 3;j < MAGIC && i + j - 2 < n;j++){
                // if(eq(A[i][j], 0)) break; // Pruning jika sudah ketemu 0
                b[i] -= b[i + j - 2] * A[i][j];
                A[i][j] = 0;
            }
            b[i] /= A[i][2];
            A[i][2] = 1;
        }
        return b;
    }
};

int main(){
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    cout.tie(0);
    // cout << scientific << setprecision(10);
    cout << fixed << setprecision(10);
    int n; cin >> n;
    double t0 = clock();
    Beam solve(n);
    // solve.printA();
    // solve.printB();
    auto ret = solve.gaussian();
    // solve.printB();
    cout << fixed << setprecision(10) << (double) (clock() - t0) / CLOCKS_PER_SEC << " Second" << endl;
}