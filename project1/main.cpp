#include <iostream>

using namespace std;


void print_menu();


struct IntervalDivider {
private:
    double left_edge_, right_edge_;
    int n_;
    
public:
    IntervalDivider(double left_edge, double right_edge, int n): left_edge_(left_edge), right_edge_(right_edge), n_(n){
        
    }
    vector<double> divideIntervals(){
        vector<double> extremes;
        double step = (right_edge_ - left_edge_) / (double)(n_);
        for(int i = 0; i <= n_; i++){
            double value = left_edge_ + i * step;
            extremes.push_back(value);
        }
        return extremes;
    };
    double returnLeftEdge(){
        return this->left_edge_;
    }
    double returnRightEdge(){
        return this->right_edge_;
    }
    int returnNumIntervals(){
        return this->n_;
    }
};


struct Interval {
private:
    double left_bracket_;
    double right_bracket_;
    bool square_left_;
    bool square_right_;
    vector<int> numbers_inside_;
    
public:
    Interval() : left_bracket_(0), right_bracket_(0), square_left_(false),  square_right_(false) {};
    Interval(double left, double right, bool leftSquare, bool rightSquare)
    : left_bracket_(left), right_bracket_(right), square_left_(leftSquare), square_right_(rightSquare){
        
    }
    bool includes(double number) {
        if (number < left_bracket_ || number > right_bracket_) {
            return false;
        }
        
        if (number == left_bracket_ && !square_left_) {
            return false;
        }
        
        if (number == right_bracket_ && !square_right_) {
            return false;
        }
        this->numbers_inside_.push_back(number);
        return true;
    }
    unsigned long returnSize(){
        return this->numbers_inside_.size();
    }
    bool returnSquareLeft(){
        return this->square_left_;
    }
    bool returnSquareRight(){
        return this->square_right_;
    }
    double returnLeftBracket(){
        return this->left_bracket_;
    }
    double returnRightBracket(){
        return this->right_bracket_;
    }
};


class GeneratorBase {
    
public:
    long long LCM(long long m, long long X0, int a, int c){
        return (a*X0 +c) % m;
    }
    long long QCM(long long m, long long X0, long long X1, int d, int a, int c){
        return (d* X0 * X0 +a *X0 + c ) % m;
    }
    long long FNM(long long m, long long X1, long long X2 ){
        return (X1 + X2)%m;
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


class ResultingBase: public GeneratorBase {
private:
    vector<Interval> makeIntervals(const vector<double> extremes){
        vector<Interval> intervals;
        for(int i = 0; i<10; i++){
            intervals.push_back(i == 0 ? Interval(extremes[i], extremes[i+1], true, true): Interval(extremes[i], extremes[i+1], false, true)) ;
        }
        return intervals;
    }
    
    
    
    vector<double> tenths = IntervalDivider(0, 1, 10).divideIntervals();
    vector<double> decades = IntervalDivider(0, 100, 10).divideIntervals();
    vector<Interval> intervalsTenths = makeIntervals(tenths);
    vector<Interval> intervalsDecades = makeIntervals(decades);
    
    
public:
    void print_result(vector<tuple<long long, float, long long>> result, long long m){
        
        for(auto& tuple: result){
            double number = get<2>(tuple);
            for(auto&  Interval : intervalsDecades){
                if (Interval.includes(number)){
                    break;
                }
            }
        }
        cout << "\t"<< "Інтервал"<< "\t" << "Частота"<< endl;
        for(auto&  Interval : intervalsDecades){
            char leftBracket = Interval.returnSquareLeft() ? '[' : '(';
            char rightBracket = Interval.returnSquareRight() ? ']' : ')';
            cout << "\t"<< leftBracket << Interval.returnLeftBracket() << ", " << Interval.returnRightBracket() << rightBracket  << "\t"<< Interval.returnSize()/static_cast<float>(m) << endl;
        }
    }
};


class LinearCongruentialMethod: public ResultingBase {
private:
    long long m_;
    int a_ = 2, c_ = 3 ;
    
public:
    LinearCongruentialMethod(long m)
    : m_(m){}
    
    void linearCongruentialMethod() {
        long X0=1, X1 = 1;
        vector<tuple<long long, float, long long>> linearCongruentialMethodVector;
        linearCongruentialMethodVector.push_back(make_tuple(X0,
                                                            static_cast<float>(X0) / static_cast<float>(m_-1), static_cast<float>(X0) / static_cast<float>(m_-1) * 100
                                                            ));
        for (int i = 1; i < m_; ++i) {
            X1 = LCM(m_, X0, a_, c_);
            X0 = X1;
            linearCongruentialMethodVector.push_back(make_tuple(X1,
                                                                static_cast<float>(X1) / static_cast<float>(m_-1),
                                                                static_cast<float>(X1) / static_cast<float>(m_-1) * 100
                                                                ));
        }
        ResultingBase::print_result(linearCongruentialMethodVector, m_);
    }
};


class QuadraticCongruentialMethod: public ResultingBase {
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
        ResultingBase::print_result(quadraticCongruentialMethodVector, m);
    }
};


class FibonachiNumbersMethod: public ResultingBase {
private:
    long long m_;
    
    
public:
    FibonachiNumbersMethod(long m)
    : m_(m){}
    void fibonachiNumbersMethod(){
        long long X1= 0 % m_, X2 = 1 % m_, X3 = 0;
        vector<tuple<long long, float, long long>> fibonachiNumbersMethodVector;
        fibonachiNumbersMethodVector.push_back(make_tuple(X1,
                                                          static_cast<float>(X1) / static_cast<float>(m_-1), static_cast<float>(X1) / static_cast<float>(m_-1) * 100
                                                          ));
        fibonachiNumbersMethodVector.push_back(make_tuple(X2,
                                                          static_cast<float>(X2) / static_cast<float>(m_-1)
                                                          , static_cast<float>(X2) / static_cast<float>(m_-1) * 100
                                                          ));
        for (int i = 2; i < m_; ++i) {
            X3 = FNM(m_, X1, X2);
            X1 = X2;
            X2 = X3;
            fibonachiNumbersMethodVector.push_back(make_tuple(X3,
                                                              static_cast<float>(X3) / static_cast<float>(m_-1),
                                                              static_cast<float>(X3) / static_cast<float>(m_-1) * 100
                                                              ));
        }
        ResultingBase::print_result(fibonachiNumbersMethodVector, m_);
    }
};


class InverseCongruentialMethod: public ResultingBase {
private:
    long long m_;
    int a_ = 2, c_ = 3, p_ = 997;
    
public:
    InverseCongruentialMethod(long m)
    : m_(m){}
    void inverseCongruentialMethod(){
        long long X0= 1, X1= 0;
        vector<tuple<long long, float, long long>> inverseCongruentialMethodVector;
        inverseCongruentialMethodVector.push_back(make_tuple(X0,
                                                             static_cast<float>(X0) / static_cast<float>(m_-1), static_cast<float>(X0) / static_cast<float>(m_-1) * 100
                                                             ));
        for (int i = 1; i < m_; ++i) {
            X1 = ICG(p_, a_, c_, X0) ;
            X0 = X1;
            inverseCongruentialMethodVector.push_back(make_tuple(X1,
                                                                 static_cast<float>(X1) / static_cast<float>(m_-1),
                                                                 static_cast<float>(X1) / static_cast<float>(m_-1) * 100
                                                                 ));
        }
        ResultingBase::print_result(inverseCongruentialMethodVector, m_);
        
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
        ResultingBase::print_result(unionMethodVector, m);
    }
};


class SigmaMethod{};
class PolarMethod{};
class RelationMethod{};
class LogarithmMethod{};
class ArensMethod{};
int main() {
    print_menu();
    IntervalDivider divider(10.0, 20.0, 10);
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

