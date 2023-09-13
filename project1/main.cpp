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
    Interval(double left_bracket, double right_bracket, bool square_left, bool square_right)
    : left_bracket_(left_bracket), right_bracket_(right_bracket), square_left_(square_left), square_right_(square_right){
        
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


struct IntervalMaker{
private:
    vector<double> extremes_;
public:
    IntervalMaker(vector<double >extremes): extremes_(extremes){
        
    }
    vector<Interval> makeIntervals(){
        vector<Interval> intervals;
        for(int i = 0; i<10; i++){
            intervals.push_back(i == 0 ? Interval(extremes_[i], extremes_[i+1], true, true): Interval(extremes_[i], extremes_[i+1], false, true)) ;
        }
        return intervals;
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
    long long mod_inv(long long param, long long n)
    {
        long long b0 = n, t, q;
        long long x0 = 0, x1 = 1;
        if(param==0) return INT_MAX;
        else if(param == INT_MAX) return 0 ;
        if (n == 1) return 1;
        while (param > 1)
        {
            q = param / n;
            t = n;
            n = param % n;
            param = t;
            t = x0;
            x0 = x1 - q * x0;
            x1 = t;
        }
        if (x1 < 0) x1 += b0;
        
        return x1;
    }
    
    long long ICG(long long n, int a, int c, long long param)
    {
        if (param == 0) return c;
        return (a * mod_inv(param, n) + c) % n;
    }
};


class ResultingEvenlyBase: public GeneratorBase {
private:
    
    vector<double> tenths = IntervalDivider(0, 1, 10).divideIntervals();
    vector<double> decades = IntervalDivider(0, 100, 10).divideIntervals();
    vector<Interval> intervals_tenths = IntervalMaker(tenths).makeIntervals();
    vector<Interval> intervals_decades = IntervalMaker(decades).makeIntervals();
    
    
public:
    void print_result(vector<tuple<long long, float, long long>> result, long long m){
        
        for(auto& tuple: result){
            double number = get<2>(tuple);
            for(auto&  Interval : intervals_decades){
                if (Interval.includes(number)){
                    break;
                }
            }
        }
        cout << "\t"<< "Інтервал"<< "\t" << "Частота"<< endl;
        for(auto&  Interval : intervals_decades){
            char leftBracket = Interval.returnSquareLeft() ? '[' : '(';
            char rightBracket = Interval.returnSquareRight() ? ']' : ')';
            cout << "\t"<< leftBracket << Interval.returnLeftBracket() << ", " << Interval.returnRightBracket() << rightBracket  << "\t"<< Interval.returnSize()/static_cast<float>(m) << endl;
        }
    }
};

class NormalBase: public GeneratorBase{
private:
    long long m_;
public:
    NormalBase(long long m): m_(m){
        cout << "hihi" << endl;
    }
};

class LinearCongruentialMethod: public ResultingEvenlyBase {
private:
    long long m_;
    int a_ = 2, c_ = 3 ;
    
public:
    LinearCongruentialMethod(long m)
    : m_(m){}
    
    void linearCongruentialMethod() {
        long X0=1, X1;
        vector<tuple<long long, float, long long>> linearCongruentialMethod_vector;
        linearCongruentialMethod_vector.push_back(make_tuple(X0,
                                                            static_cast<float>(X0) / static_cast<float>(m_-1), static_cast<float>(X0) / static_cast<float>(m_-1) * 100
                                                            ));
        for (int i = 1; i < m_; ++i) {
            X1 = LCM(m_, X0, a_, c_);
            X0 = X1;
            linearCongruentialMethod_vector.push_back(make_tuple(X1,
                                                                static_cast<float>(X1) / static_cast<float>(m_-1),
                                                                static_cast<float>(X1) / static_cast<float>(m_-1) * 100
                                                                ));
        }
        ResultingEvenlyBase::print_result(linearCongruentialMethod_vector, m_);
    }
};


class QuadraticCongruentialMethod: public ResultingEvenlyBase {
private:
    long long m;
    int a = 6  , c = 3 , d = 1;
    
public:
    QuadraticCongruentialMethod(long m)
    : m(m){}
    void quadraticCongruentialMethod(){
        long long X0 = 13, X1= 0;
        
        vector<tuple<long long, float, long long>> quadraticCongruentialMethod_vector;
        quadraticCongruentialMethod_vector.push_back(make_tuple(X0,
                                                               static_cast<float>(X0) / static_cast<float>(m-1), static_cast<float>(X0) / static_cast<float>(m-1) * 100
                                                               ));
        for (int i = 1; i < m; ++i) {
            X1 = QCM(m, X0, X1, d, a, c);
            X0 = X1;
            quadraticCongruentialMethod_vector.push_back(make_tuple(X1,
                                                                   static_cast<float>(X1) / static_cast<float>(m-1),
                                                                   static_cast<float>(X1) / static_cast<float>(m-1) * 100
                                                                   ));
        }
        ResultingEvenlyBase::print_result(quadraticCongruentialMethod_vector, m);
    }
};


class FibonachiNumbersMethod: public ResultingEvenlyBase {
private:
    long long m_;
    
    
public:
    FibonachiNumbersMethod(long m)
    : m_(m){}
    void fibonachiNumbersMethod(){
        long long X1= 0 % m_, X2 = 1 % m_, X3 = 0;
        vector<tuple<long long, float, long long>> fibonachiNumbersMethod_vector;
        fibonachiNumbersMethod_vector.push_back(make_tuple(X1,
                                                          static_cast<float>(X1) / static_cast<float>(m_-1), static_cast<float>(X1) / static_cast<float>(m_-1) * 100
                                                          ));
        fibonachiNumbersMethod_vector.push_back(make_tuple(X2,
                                                          static_cast<float>(X2) / static_cast<float>(m_-1)
                                                          , static_cast<float>(X2) / static_cast<float>(m_-1) * 100
                                                          ));
        for (int i = 2; i < m_; ++i) {
            X3 = FNM(m_, X1, X2);
            X1 = X2;
            X2 = X3;
            fibonachiNumbersMethod_vector.push_back(make_tuple(X3,
                                                              static_cast<float>(X3) / static_cast<float>(m_-1),
                                                              static_cast<float>(X3) / static_cast<float>(m_-1) * 100
                                                              ));
        }
        ResultingEvenlyBase::print_result(fibonachiNumbersMethod_vector, m_);
    }
};


class InverseCongruentialMethod: public ResultingEvenlyBase {
private:
    long long m_;
    int a_ = 2, c_ = 3, p_ = 997;
    
public:
    InverseCongruentialMethod(long m)
    : m_(m){}
    void inverseCongruentialMethod(){
        long long X0= 1, X1= 0;
        vector<tuple<long long, float, long long>> inverseCongruentialMethod_vector;
        inverseCongruentialMethod_vector.push_back(make_tuple(X0,
                                                             static_cast<float>(X0) / static_cast<float>(m_-1), static_cast<float>(X0) / static_cast<float>(m_-1) * 100
                                                             ));
        for (int i = 1; i < m_; ++i) {
            X1 = ICG(p_, a_, c_, X0) ;
            X0 = X1;
            inverseCongruentialMethod_vector.push_back(make_tuple(X1,
                                                                 static_cast<float>(X1) / static_cast<float>(m_-1),
                                                                 static_cast<float>(X1) / static_cast<float>(m_-1) * 100
                                                                 ));
        }
        ResultingEvenlyBase::print_result(inverseCongruentialMethod_vector, m_);
        
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
        vector<tuple<long long, float, long long>> unionMethod_vector;
        
        for (int i = 0; i < m; ++i) {
            X1 = LCM(m, X0, 4, 6);
            X0 = X1;
            Y1 = LCM(m/1000, Y0, 6, 7);
            Y0 = Y1;
            Z = (X1 - Y1);
            if(Z<0){
                Z+=m;
            }
            Z%=m;
            unionMethod_vector.push_back(make_tuple(Z,
                                                   static_cast<float>(Z) / static_cast<float>(m-1),
                                                   static_cast<float>(Z) / static_cast<float>(m-1) * 100
                                                   ));
        }
        ResultingEvenlyBase::print_result(unionMethod_vector, m);
    }
};

class ResultingNormalBase: public GeneratorBase {
private:
    vector<Interval> makeIntervals(const vector<double> extremes){
        vector<Interval> intervals;
        for(int i = 0; i<10; i++){
            intervals.push_back(i == 0 ? Interval(extremes[i], extremes[i+1], true, true): Interval(extremes[i], extremes[i+1], false, true)) ;
        }
        return intervals;
    }
    
    
    
    vector<double> threes = IntervalDivider(-3, 3, 10).divideIntervals();
    
    vector<Interval> intervals_normal = makeIntervals(threes);
    
    
public:
    void print_result(vector<tuple<long long, float, long long>> result, long long m){
        
        for(auto& tuple: result){
            double number = get<2>(tuple);
            for(auto&  Interval : intervals_normal){
                if (Interval.includes(number)){
                    break;
                }
            }
        }
        cout << "\t"<< "Інтервал"<< "\t" << "Частота"<< endl;
        for(auto&  Interval : intervals_normal){
            char leftBracket = Interval.returnSquareLeft() ? '[' : '(';
            char rightBracket = Interval.returnSquareRight() ? ']' : ')';
            cout << "\t"<< leftBracket << Interval.returnLeftBracket() << ", " << Interval.returnRightBracket() << rightBracket  << "\t"<< Interval.returnSize()/static_cast<float>(m) << endl;
        }
    }
};
class SigmaMethod: public ResultingNormalBase{
private:
    long long m_;
    int a_ = 2, c_ = 3;
public:
    SigmaMethod(long m)
    : m_(m){}
    
    void sigmaMethod() {
        long X0=1, X1;
        vector<tuple<long long, float, long long>> linearMethod_vector, sigmaMethod_vector;
//        array<tuple<long long, float, long long>, 12> selected_numbers;
        linearMethod_vector.push_back(make_tuple(X0,
                                                            static_cast<float>(X0) / static_cast<float>(m_-1), static_cast<float>(X0) / static_cast<float>(m_-1) * 100
                                                            ));
        for (int i = 1; i < m_; ++i) {
            
            X1 = LCM(m_, X0, a_, c_);
            X0 = X1;
            linearMethod_vector.push_back(make_tuple(X1,
                                                                static_cast<float>(X1) / static_cast<float>(m_-1),
                                                                static_cast<float>(X1) / static_cast<float>(m_-1) * 100
                                                                ));
        }
       
        ResultingNormalBase::print_result(sigmaMethod_vector, m_);
    }
};
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
            InverseCongruentialMethod generator(1000);
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
            SigmaMethod generator(m);
            generator.sigmaMethod();
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

