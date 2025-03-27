#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

struct Test {
    double (*f)(double);
    string fName;
    double* x;
    double A;
    double B;
    int N;
};

Test TestFunc[];

struct Table {
    int Size;
    double A, B; //граничиные условия f`
    double* x;
    double* y;

    Table(string FileName) {
        ifstream in;
        in.open(FileName);
        if (!in.is_open()) {
            cout << " Cant open file " << endl;
            return;
        }
        bool isCorrectOrder;
        in >> Size >> A >> B;
        if (Size < 3) {
            cout << " IER 1 not enough points " << endl;
            return;
        }
        x = new double[Size];
        y = new double[Size];
        for (int i = 0; i < Size; i++) {
            in >> x[i];
        }
        for (int i = 0; i < Size; i++) {
            in >> y[i];
        }
        in.close();
        isCorrectOrder = true;
        for (int i = 1; i < Size; i++) {
            if (x[i - 1] >= x[i]) {
                isCorrectOrder = false;
                break;
            }
        }
        if (!isCorrectOrder) {
            cout << " IER 2 wrong X order " << endl;
            return;
        }
    }

    Table(double (*f)(double), int N, double A, double B, double* X) {
        Size = N;
        this->A = A;
        this->B = B;
        x = new double[N];
        y = new double[N];
        for (int i = 0; i < N; i++) {
            x[i] = X[i];
            y[i] = f(x[i]);
        }
    }

    Table(int N, double A, double B, double* X, double* Y) {
        Size = N;
        this->A = A;
        this->B = B;
        x = new double[N];
        y = new double[N];
        for (int i = 0; i < N; i++) {
            x[i] = X[i];
            y[i] = Y[i];
        }
    }

    void Display() {
        cout << "A = " << A << " B = " << B << endl;
        cout << " x ";
        for (int i = 0; i < Size; i++) {
            cout << x[i] << " ";
        }
        cout << endl << " y ";
        for (int i = 0; i < Size; i++) {
            cout << y[i] << " ";
        }
        cout << endl;
        cout << endl;
    }

   /* ~Table() {
        delete[] x;
        delete[] y;
    }*/
};

double* SweepMethod(double* A, double* B, double* C, double* vec, int N) {//progonka 
    double* res = new double[N];
    double* mu = new double[N];
    double* nu = new double[N];
    mu[0] = vec[0] / A[0];
    nu[0] = -1 * C[0] / A[0];
    for (int i = 1; i < N - 1; i++) {
        mu[i] = (vec[i] - B[i - 1] * mu[i - 1]) / (A[i] + B[i - 1] * nu[i - 1]);
        nu[i] = (-1 * C[i]) / (A[i] + B[i - 1] * nu[i - 1]);
    }
    mu[N - 1] = vec[N - 1] / A[N - 1];
    nu[N - 1] = -1 * (B[N - 2] / A[N - 1]);
    /*for (int i = 0; i < 3; i++) {
        cout << " Deb mu nu " << mu[i] << " " << nu[i] << " " << i << endl;
    }*/
    res[N - 1] = (mu[N - 1] + nu[N - 1] * mu[N - 2]) / (1 - nu[N - 1] * nu[N - 2]);
    for (int i = N - 2; i >= 0; i--) {
        res[i] = mu[i] + nu[i] * res[i + 1];
    }   
    return res;
}

class CubicSplain {
    double a;
    double b;
    double c;
    double d;
    double left;
    double right;

public:
    CubicSplain() {
        a = b = c = d = left = right = 0;
    }

    CubicSplain(double Left, double Right) {
        a = b = c = d = 0;
        left = Left;
        right = Right;
    }

    void SetBorders(double Left, double Right) {
        left = Left;
        right = Right;
    }

    void SetCoef(double A, double B, double C, double D) {
        a = A;
        b = B;
        c = C;
        d = D;
    }

    string BordersToString() {
        stringstream str;
        str << "[" << left << ":" << right << "]";
        return str.str();
    }

    string SplainToString() {
        stringstream str;
        str << d << "*x**3 + " << c << "*x**2 + " << b << "*x + " << a;
        return str.str();
    }

    friend CubicSplain* CalculateCubicSplainInterpolation(Table data);
    friend ostream& operator<<(ostream& out, CubicSplain const& p);
};

ostream& operator<<(ostream& out, CubicSplain const& p){
    out << p.d << "*x**3 + " << p.c << "*x**2 + " << p.b << "*x + " << p.a;
    return out;
}

CubicSplain* CalculateCubicSplainInterpolation(Table data) {
    int numberOfSplains = data.Size - 1;
    CubicSplain* res = new CubicSplain[numberOfSplains];
    for (int i = 0; i < numberOfSplains; i++) {
        res[i].SetBorders(data.x[i], data.x[i + 1]);
    }

    double* h = new double[numberOfSplains];
    for (int i = 0; i < numberOfSplains; i++) {
        h[i] = data.x[i + 1] - data.x[i];
    }

    double* A = new double[numberOfSplains + 1]; //main
    double* B = new double[numberOfSplains]; //left
    double* C = new double[numberOfSplains]; //right
    double* f = new double[numberOfSplains + 1]; 

    A[0] = 2;
    C[0] = 1;
    f[0] = (3 / h[0]) * ((data.y[1] - data.y[0]) / h[0] - data.A);
    for (int i = 1; i < numberOfSplains ; i++) {
        B[i-1] = h[i - 1];
        A[i] = 2 * (h[i - 1] + h[i]);
        C[i] = h[i];
        f[i] = 3 * ((data.y[i + 1] - data.y[i]) / h[i] - (data.y[i] - data.y[i - 1]) / h[i - 1]);
    }
    A[numberOfSplains] = 2;
    B[numberOfSplains - 1] = 1;
    f[numberOfSplains] = (3 / h[numberOfSplains - 1]) * (data.B - (data.y[data.Size - 1] - data.y[data.Size - 2]) / h[numberOfSplains - 1]);

    double* coefsC = SweepMethod(A, B, C, f, numberOfSplains + 1);

   /* for (int i = 0; i < numberOfSplains + 1; i++) {
        cout << "deb C " << coefsC[i] << endl;
    }*/

    double* coefsA = new double[numberOfSplains];
    double* coefsD = new double[numberOfSplains];
    double* coefsB = new double[numberOfSplains];

    for (int i = 0; i < numberOfSplains; i++) {
        coefsA[i] = data.y[i];
        coefsD[i] = (coefsC[i + 1] - coefsC[i]) / (3 * h[i]);
       // coefsB[i] = 2*h[i] * coefsC[i+1] - 3*h[i] * h[i] * coefsD[i] + (data.y[i + 1] - data.y[i]) / h[i];
        coefsB[i] = (data.y[i + 1] - data.y[i]) / h[i] - h[i] * (coefsC[i + 1] + 2 * coefsC[i]) / 3;
    }

   /* for (int i = 0; i < numberOfSplains; i++) {
        cout << " deb " << coefsA[i] << " " << coefsB[i] << " " << coefsC[i+1] << " " << coefsD[i] << " " << endl;
    }*/

    for (int i = 0; i < numberOfSplains; i++) {
        res[i].SetCoef((coefsA[i] - coefsB[i] * res[i].left + coefsC[i] * pow(res[i].left, 2) - coefsD[i] * pow(res[i].left, 3)), (coefsB[i] - 2 * coefsC[i] * res[i].left + 3 * coefsD[i] * pow(res[i].left, 2)), (coefsC[i] - 3 * coefsD[i] * res[i].left), coefsD[i]); //Ddadaadadadada
    }
    return res;
    
}

//obobshit
void BuildSplainGraph(CubicSplain* splains, int numberOfSplains) {
    ofstream out("temp.dat");
    for (int i = 0; i < numberOfSplains; i++) {
        out << "P" << i << "(x)= " << splains[i].SplainToString() << "\n";
    }
    out << "plot 0 ls 0, ";//costil
    for (int i = 0; i < numberOfSplains - 1; i++) {
        out << splains[i].BordersToString() << " P" << i << "(x) ls " << i + 1 << ", ";
    }
    out << splains[numberOfSplains-1].BordersToString() << " P" << numberOfSplains-1 << "(x) ls " << numberOfSplains;
    out << "\n" << "pause -1" << "\n";
    out.close();
    system("temp.dat");
}

void BuildSplainGraph(CubicSplain* splains, int numberOfSplains, string origFunc) { //для тестов
    ofstream out("temp.dat");
    out << "f(x)= " << origFunc << "\n";
    for (int i = 0; i < numberOfSplains; i++) {
        out << "P" << i << "(x)= " << splains[i].SplainToString() << "\n";
    }
    out << "plot f(x) ls 1, ";//costil
    for (int i = 0; i < numberOfSplains - 1; i++) {
        out << splains[i].BordersToString() << " P" << i << "(x) ls " << i + 2 << ", ";
    }
    out << splains[numberOfSplains - 1].BordersToString() << " P" << numberOfSplains - 1 << "(x) ls " << numberOfSplains;
    out << "\n" << "pause -1" << "\n";
    out.close();
    system("temp.dat");
}

double test1[]{ 0,1,2,3,4,5 };
double test2[]{-3.14, -2.3, -1.57, -0.75, 0, 1, 2, 3, 3.14};
double test3[]{0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 3, 4 };
double test4[]{ 0.1, 0.2, 0.3, 0.5, 1, 2, 4, 12 };

double CalculateFirstDerivative(double (*f)(double), double x, double h) {
    return (f(x + h) - f(x - h)) / (2 * h);
}

Test TestFunc[] = {
        { [](double x) { return pow(x,3); },"x**3", test1, 0, 75, 6 },
        { [](double x) { return pow(x,2); },"x**2", test1, 0, 10, 6 },
        { [](double x) { return x; },"x", test1, 1, 1, 6 },
        { [](double x) { return (double)1; },"1", test1, 0, 0, 6 },
        { [](double x) { return 4 * pow(x,4) + 3 * pow(x,3) + 2 * pow(x,2) + x + 1; },"4*x**4 + 3*x**3 + 2*x**2 + x + 1",test1, 1, 2246, 6},
        { [](double x) { return sin(x); },"sin(x)", test2, cos(test2[0]), cos(test2[7]), 8},
        { [](double x) { return log(x); },"log(x)", test3, 1/test3[0], 1/(test3[8]), 9},
        { [](double x) { return 1/log(x+1)+x; },"1/log(x+1)+x", test4, 1 - (1 / ((test4[0] + 1) * pow(log(test4[0] + 1),2))), 1 - (1 / ((test4[7] + 1) * pow(log(test4[7] + 1),2))), 8}
};

int main()
{
    int numberOfTests = sizeof(TestFunc) / sizeof(Test);
    int chose;
    do {
        cout << " Chose Test 1-" << numberOfTests << endl;
        for (int i = 0; i < numberOfTests; i++) {
            cout << " " << i+1 << " " << TestFunc[i].fName << endl;
        }
        cin >> chose;
        if (chose > 0 && chose <= numberOfTests) {
            Table testTable(TestFunc[chose - 1].f, TestFunc[chose - 1].N, TestFunc[chose - 1].A, TestFunc[chose - 1].B, TestFunc[chose - 1].x);
            testTable.Display();
            CubicSplain* test = CalculateCubicSplainInterpolation(testTable);
            for (int i = 0; i < testTable.Size - 1; i++) {
                cout << " " << test[i].BordersToString() << " " << test[i].SplainToString() << endl;
            }
            BuildSplainGraph(test, TestFunc[chose - 1].N - 1, TestFunc[chose -1 ].fName);
            delete[] test;
        }
    } while (chose != 0);

    string File;
    do {
        cout << " Enter Test file name ";
        cin >> File;
        if (File != "0") {
            Table test(File);

            CubicSplain* splains = CalculateCubicSplainInterpolation(test);

            for (int i = 0; i < test.Size - 1; i++) {
                cout << " " << splains[i].BordersToString() << " " << splains[i].SplainToString() << endl;
            }

            BuildSplainGraph(splains, test.Size - 1);
            delete[] splains;
        }
    } while (File != "0");

    return 0;
}

