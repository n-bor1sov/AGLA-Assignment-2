//Nikita Borisov DSAI-03
#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;


class Matrix {
public:
    pair<int, int> size;
    vector<vector<double>> matrix;
    Matrix(vector<pair<double, double>> & vec, int degree) {
        size = pair<double, double>(vec.size(), degree + 1);
        matrix.resize(vec.size());
        for(int i = 0; i < matrix.size(); i++) {
            for(double j = 0; j <= degree; j++) {
                matrix[i].emplace_back(pow(vec[i].first, j));
            }
        }
    }
    Matrix(int r, int c) {
        size = pair<int, int>(r, c);
        matrix = vector<vector<double>>(r);
    }
    Matrix() = default;
    void InputSize(){
        cin >> size.first >> size.second;
    }
    friend ostream &operator<<(ostream &output, const Matrix &M ) {
        for(int i = 0; i < M.size.first; i++) {
            for(int j = 0; j < M.size.second; j++) {
                if(abs(M.matrix[i][j]) < 0.00000001) {
                    output << 0.00;
                } else {
                    output << M.matrix[i][j];
                }
                if(j != M.size.second - 1) {
                    cout << " ";
                }
            }
            cout << '\n';
        }
        return output;
    }

    friend istream &operator>>(istream  &input, Matrix &M ) {
        for(int i = 0; i < M.size.first; i++) {
            vector<double> a(M.size.second, 0);
            M.matrix.push_back(a);
            for(int j = 0; j < M.size.second; j++) {
                input >> M.matrix[i][j];
            }
        }
        return input;
    }
    void operator=(const Matrix &M) {
        for(int i = 0; i < M.size.first; i++) {
            for (int j = 0; j < M.size.second; j++) {
                matrix[i].push_back(M.matrix[i][j]);
            }
        }
    }
    Matrix operator+(const Matrix &M) {
        Matrix result(size.first, size.second);
        for (int i = 0; i < size.first; i++) {
            for (int j = 0; j < size.second; j++) {
                result.matrix[i].push_back(matrix[i][j] + M.matrix[i][j]);
            }
        }
        return result;
    }
    Matrix operator-(const Matrix &M) {
        Matrix result(size.first, size.second);
        for(int i = 0; i < size.first; i++) {
            for(int j = 0; j < size.second; j++)  {
                result.matrix[i].push_back(matrix[i][j] - M.matrix[i][j]);
            }
        }
        return result;
    }
    Matrix operator*(const Matrix &M) {
        Matrix result(size.first, M.size.second);
        for(int i = 0; i < size.first; i++) {
            for(int j = 0; j < M.size.second; j++)  {
                double value = 0;
                for(int k = 0; k < M.size.first; k++) {
                    value += matrix[i][k] * M.matrix[k][j];
                }
                result.matrix[i].push_back(value);
            }
        }
        return result;
    }
    Matrix transpose() {
        Matrix result(size.second, size.first);
        for(int i = 0; i < size.second; i++) {
            for(int j = 0; j < size.first; j++) {
                result.matrix[i].push_back(matrix[j][i]);
            }
        }
        return result;
    }
    Matrix inverse() {
        Matrix inversionMatrix(size.first, size.first * 2);
        for(int i = 0; i < size.first; i++) {
            for(int j = 0; j < size.second; j++) {
                inversionMatrix.matrix[i].emplace_back(matrix[i][j]);
            }
        }
        for(int i = 0; i < size.first; i++) {
            for(int j = size.second; j < size.second * 2; j++) {
                if(i != j - size.first) {
                    inversionMatrix.matrix[i].emplace_back(0);
                } else {
                    inversionMatrix.matrix[i].emplace_back(1);
                }
            }
        }
        for (int i = 0; i < inversionMatrix.size.first; i++) {
            double maxPivot = matrix[i][i];
            int newPivotRaw = i;
            for (int j = i; j < inversionMatrix.size.first; j++) {
                if (j != i) {
                    if (inversionMatrix.matrix[j][i] != 0) {
                        double coef = inversionMatrix.matrix[j][i] / inversionMatrix.matrix[i][i];
                        for (int k = 0; k < inversionMatrix.size.second; k++) {
                            inversionMatrix.matrix[j][k] = inversionMatrix.matrix[j][k] - inversionMatrix.matrix[newPivotRaw][k] * coef;
                        }
                    }
                }
            }
        }
        for (int i = inversionMatrix.size.first - 1; i >= 0; i--) {
            for (int j = i; j >= 0; j--) {
                if (j != i) {
                    if (inversionMatrix.matrix[j][i] != 0) {
                        double coef = inversionMatrix.matrix[j][i] / inversionMatrix.matrix[i][i];
                        for (int k = 0; k < inversionMatrix.size.second; k++) {
                            inversionMatrix.matrix[j][k] = inversionMatrix.matrix[j][k] - inversionMatrix.matrix[i][k] * coef;
                        }
                    }
                }
            }
        }
        for (int i = 0; i < inversionMatrix.size.first; i++) {
            if (inversionMatrix.matrix[i][i] != 1) {
                double coef = inversionMatrix.matrix[i][i];
                for(int j = 0; j < inversionMatrix.size.second; j++) {
                    inversionMatrix.matrix[i][j] = inversionMatrix.matrix[i][j] / coef;
                }
            }
        }
        Matrix result(size.first, size.second);
        for(int i = 0; i < size.first; i++) {
            for(int j = size.second; j < size.second * 2; j++) {
                result.matrix[i].emplace_back(inversionMatrix.matrix[i][j]);
            }
        }
        return result;
    }
};


#define GNUPLOT_NAME "C:\\gnuplot\\bin\\gnuplot -persist"

int main() {
    std::cout << std::setprecision(4) << std::fixed;
    int amountOfData;
    int degree;
    cin >> amountOfData;
    vector<pair<double, double>> Data(amountOfData);
    for(int i = 0; i < amountOfData; i++) {
        cin >> Data[i].first >> Data[i].second;
    }
    cin >> degree;
    Matrix A(Data, degree);
    cout << "A:" << '\n';
    cout << A;
    cout << "A_T*A:" << '\n';
    Matrix A_T = A.transpose();
    Matrix transposeTimesOriginal = A_T* A;
    cout << transposeTimesOriginal;
    cout << "(A_T*A)^-1:" << '\n';
    Matrix A_TA_1  = transposeTimesOriginal.inverse();
    cout << A_TA_1;
    Matrix B(A.size.first, 1);
    for(int i = 0; i < A.size.first; i++) {
        B.matrix[i].emplace_back(Data[i].second);
    }
    Matrix A_TB = A_T * B;
    cout << "A_T*b:" << '\n';
    cout << A_TB;
    cout << "x~:" << '\n';
    Matrix x = A_TA_1 * A_TB;
    cout << A_TA_1 * A_TB;

    FILE* pipe = _popen(GNUPLOT_NAME, "w");

    if (pipe == NULL) {
        cout << "Pipe failed" << endl;
        return 0;
    }

    // Setup GNUPlot
    fprintf(pipe, "set xtics 1\n");
    fprintf(pipe, "set ytics 1\n");
    fprintf(pipe, "set grid xtics ytics\n");

    // Function print
    string response = "plot ";

    for(int i = x.matrix.size() - 1; i > 0; i--) {
        response += to_string(x.matrix[i][0]) + "*x**" + to_string(i);
        if(x.matrix[i][0] < 0) {
            response += "+";
        }
    }
    response += to_string(x.matrix[x.matrix.size() - 1][0]);

    cout << response << endl;

    response += " with lines, '-' with points pointtype 6 pointsize 1";

    // Print points
    fprintf(pipe, "%s", response.c_str());
    for (int i = 0; i < Data.size(); i++)
        fprintf(pipe, "%f %f\n", Data[i].first, Data[i].second);

    fflush(pipe);
    pclose(pipe);

}
