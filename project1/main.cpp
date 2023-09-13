#include <iostream>

using namespace std;
using namespace std;


void print_menu();

struct Interval {
private:
    double left_bracket;
    double right_bracket;
    bool square_left;
    bool square_right;
    vector<int> numbers_inside;
public:
    Interval(double left, double right, bool leftSquare, bool rightSquare)
    : left_bracket(left), right_bracket(right), square_left(leftSquare), square_right(rightSquare) {}
    
    bool includes(double number) {
        if (number < left_bracket || number > right_bracket) {
            return false;
        }
        
        if (number == left_bracket && !square_left) {
            return false;
        }
        
        if (number == right_bracket && !square_right) {
            return false;
        }
        this->numbers_inside.push_back(number);
        return true;
    }
    unsigned long returnSize(){
        return numbers_inside.size();
    }
    bool returnSquareLeft(){
        return this->square_left;
    }
    bool returnSquareRight(){
        return this->square_right;
    }
    double returnLeftBracket(){
        return this->left_bracket;
    }
    double returnRightBracket(){
        return this->right_bracket;
    }
};
class RandomGeneratorBase {
private:
    array<Interval, 10> intervals{Interval(0, 0.1, true, true), Interval(0.1, 0.2, false, true), Interval(0.2, 0.3, false, true), Interval(0.3, 0.4, false, true), Interval(0.4, 0.5, false, true), Interval(0.5, 0.6, false, true), Interval(0.6, 0.7, false, true), Interval(0.7, 0.8, false, true), Interval(0.8, 0.9, false, true), Interval(0.9, 1, false, true), };
public:
    void print_result(vector<tuple<long long, float, long long>> result, long long m){
        
        for(auto& tuple: result){
            float number = get<1>(tuple);
            for(auto&  Interval : intervals){
                if (Interval.includes(number)){
                    break;
                }
            }
        }
        cout << "\t"<< "Інтервал"<< "\t" << "Частота"<< endl;
        for(auto&  Interval : intervals){
            char leftBracket = Interval.returnSquareLeft() ? '[' : '(';
            char rightBracket = Interval.returnSquareRight() ? ']' : ')';
            cout << "\t"<< leftBracket << Interval.returnLeftBracket() << ", " << Interval.returnRightBracket() << rightBracket  << "\t"<< Interval.returnSize()/static_cast<float>(m) << endl;
            }
    }
    long long LCM(long long m, long long X0, int a, int c){
        return (a*X0 +c) % m;
    }
    long long QCM(long long m, long long X0, long long X1, int d, int a, int c){
        return (d* X0 * X0 +a *X0 + c ) % m;
    }
    long long FNM(long long m, long long X1, long long X2 ){
        return (X1 + X2)%m;
    }
};


class LinearCongruentialMethod: public RandomGeneratorBase {
private:
    long long m;
    int a = 2, c = 3 ;

public:
    LinearCongruentialMethod(long m)
        : m(m){}

    void linearCongruentialMethod() {
        long X0=1, X1 = 1;
        vector<tuple<long long, float, long long>> linearCongruentialMethodVector;
        linearCongruentialMethodVector.push_back(make_tuple(X0,
                                                      static_cast<float>(X0) / static_cast<float>(m-1), static_cast<float>(X0) / static_cast<float>(m-1) * 100
                                                      ));
        for (int i = 1; i < m; ++i) {
            X1 = LCM(m, X0, a, c);
            X0 = X1;
            linearCongruentialMethodVector.push_back(make_tuple(X1,
                                                          static_cast<float>(X1) / static_cast<float>(m-1),
                                                          static_cast<float>(X1) / static_cast<float>(m-1) * 100
                                                          ));
        }
        RandomGeneratorBase::print_result(linearCongruentialMethodVector, m);
    }
    
};


class QuadraticCongruentialMethod: public RandomGeneratorBase {
private:
    long long m;
    int a = 6  , c = 3 , d = 1;
public:
    QuadraticCongruentialMethod(long m)
        : m(m){}
    void quadraticCongruentialMethod(){
        long long X0 = 13, X1= 0;
        
        vector<tuple<long long, float, long long>> quadraticCongruentialMethodVector;
        quadraticCongruentialMethodVector.push_back(make_tuple(X0,
                                                      static_cast<float>(X0) / static_cast<float>(m-1), static_cast<float>(X0) / static_cast<float>(m-1) * 100
                                                      ));
        for (int i = 1; i < m; ++i) {
            X1 = QCM(m, X0, X1, d, a, c);
            X0 = X1;
            quadraticCongruentialMethodVector.push_back(make_tuple(X1,
                                                          static_cast<float>(X1) / static_cast<float>(m-1),
                                                          static_cast<float>(X1) / static_cast<float>(m-1) * 100
                                                          ));
        }
        RandomGeneratorBase::print_result(quadraticCongruentialMethodVector, m);
    }
    
};


class FibonachiNumbersMethod: public RandomGeneratorBase {
private:
    long long m;
    

public:
    FibonachiNumbersMethod(long m)
        : m(m){}
    void fibonachiNumbersMethod(){
        long long X1= 0 % m, X2 = 1 % m, X3 = 0;
        vector<tuple<long long, float, long long>> fibonachiNumbersMethodVector;
        fibonachiNumbersMethodVector.push_back(make_tuple(X1,
                                                      static_cast<float>(X1) / static_cast<float>(m-1), static_cast<float>(X1) / static_cast<float>(m-1) * 100
                                                      ));
        fibonachiNumbersMethodVector.push_back(make_tuple(X2,
                                                      static_cast<float>(X2) / static_cast<float>(m-1)
                                                      , static_cast<float>(X2) / static_cast<float>(m-1) * 100
                                                      ));
        for (int i = 2; i < m; ++i) {
            X3 = FNM(m, X1, X2);
            X1 = X2;
            X2 = X3;
            fibonachiNumbersMethodVector.push_back(make_tuple(X3,
                                                          static_cast<float>(X3) / static_cast<float>(m-1),
                                                          static_cast<float>(X3) / static_cast<float>(m-1) * 100
                                                          ));
        }
        
        RandomGeneratorBase::print_result(fibonachiNumbersMethodVector, m);
    }
    
};


class InverseCongruentialMethod: public RandomGeneratorBase {
private:
    long long m;
    int a = 2, c = 3, p = 997;

public:
    InverseCongruentialMethod(long m)
        : m(m){}
    void inverseCongruentialMethod(){
        long long X0= 1, X1= 0;
        vector<tuple<long long, float, long long>> inverseCongruentialMethodVector;
        inverseCongruentialMethodVector.push_back(make_tuple(X0,
                                                      static_cast<float>(X0) / static_cast<float>(m-1), static_cast<float>(X0) / static_cast<float>(m-1) * 100
                                                      ));
        for (int i = 1; i < m; ++i) {

            X1 = ICG(p, a, c, X0) ;
            X0 = X1;
                inverseCongruentialMethodVector.push_back(make_tuple(X1,
                                                                     static_cast<float>(X1) / static_cast<float>(m-1),
                                                                     static_cast<float>(X1) / static_cast<float>(m-1) * 100
                                                                     ));

        }
                RandomGeneratorBase::print_result(inverseCongruentialMethodVector, m);
            }
            long long mod_inv(long long a, long long n)
            {
                long long b0 = n, t, q;
                long long x0 = 0, x1 = 1;
                if (n == 1) return 1;
                while (a > 1)
                    {
                    q = a / n;
                    t = n;
                    n = a % n;
                    a = t;
                    t = x0;
                    x0 = x1 - q * x0;
                    x1 = t;
                }
                if (x1 < 0) x1 += b0;
                return x1;
            }

            long long ICG(long long n, int a, int c, long long X0)
            {
                
                return (a * mod_inv(X0, n) + c) % n;
            }

        };

        class UnionMethod: public FibonachiNumbersMethod{
        private:
            long m;
        public:
            UnionMethod(long m)
                    : FibonachiNumbersMethod(static_cast<long long>(m)), m(m) {

                }

            void unionMethod() {
                long long X0 = 6, Y0 = 1, X1, Y1, Z ;
                vector<tuple<long long, float, long long>> unionMethodVector;
                
                for (int i = 0; i < m; ++i) {
                    X1 = LCM(m, X0, 4, 6);
                    X0 = X1;
                    Y1 = LCM(m/1000, Y0, 6, 7);
                    Y0 = Y1;
                    Z = abs(X1 - Y1)%m;
                    unionMethodVector.push_back(make_tuple(Z,
                                                                  static_cast<float>(Z) / static_cast<float>(m-1),
                                                                  static_cast<float>(Z) / static_cast<float>(m-1) * 100
                                                                  ));
                }
                RandomGeneratorBase::print_result(unionMethodVector, m);
            }
        };

        class SigmaMethod{};
        class PolarMethod{};
        class RelationMethod{};
        class LogarithmMethod{};
        class ArensMethod{};
        int main() {
            print_menu();

            int type;
            long long m = 100000;
            cin >> type;

            switch (type) {
                case 1:{
                    LinearCongruentialMethod generator(m);
                    generator.linearCongruentialMethod();
                    break;
                }
                case 2:{
                    QuadraticCongruentialMethod generator(m);
                    generator.quadraticCongruentialMethod();
                    break;
                }
                case 3:
                {
                    FibonachiNumbersMethod generator(m);
                    generator.fibonachiNumbersMethod();
                    break;
                }
                case 4:
                {
                    InverseCongruentialMethod generator(m);
                    generator.inverseCongruentialMethod();
                    break;
                }
                case 5:
                {
                    UnionMethod generator(m);
                    generator.unionMethod();
                    break;
                }
                case 6:
                {
                    //code
                    break;
                }
                case 7:
                {
                    //code
                    break;
                }
                case 8:{
                    //code
                    break;
                    
                }
                case 9:{
                    //code
                    break;
                    
                }
                case 10:{
                    //code
                    break;
                    
                }
                default:{
                    cout << "Wrong input!!!\nOnly 1-10 requires";
                }
            }

            return 0;
        }

        void print_menu() {
            cout << "Enter generator type: " << endl;
            cout << "\t1 - Linear Congruential Generator" << endl;
            cout << "\t2 - Quadratic Congruential Generator" << endl;
            cout << "\t3 - Fibonacci Numbers Generator" << endl;
            cout << "\t4 - Inverse Congruent SequenceGenerator" << endl;
            cout << "\t5 - Union Method Generator" << endl;
            cout << "\t6 - Sigma Method Generator" << endl;
            cout << "\t7 - Polar Coordinate Method Generator" << endl;
            cout << "\t8 - Relation Method Generator" << endl;
            cout << "\t9 - Logarithm Method Generator" << endl;
            cout << "\t10 - Arens Method Generator\n" << endl;
        }
