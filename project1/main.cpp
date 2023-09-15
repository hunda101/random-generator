#include <iostream>
#include <math.h>

using namespace std;


void print_menu();


struct IntervalEdges {
private:
    double left_edge_, right_edge_;
    int n_;
    
public:
    IntervalEdges(double left_edge, double right_edge, int n): left_edge_(left_edge), right_edge_(right_edge), n_(n){
        
    }
    vector<float> divideIntervals(){
        vector<float> extremes;
        double step = (right_edge_ - left_edge_) / (double)(n_);
        for(int i = 0; i <= n_; i++){
            double value = left_edge_ + i * step;
            extremes.push_back(value);
        }
        return extremes;
    };
    double returnLeftEdge(){
        return left_edge_;
    };
    double returnRightEdge(){
        return right_edge_;
    };
    int returnSegmentsNumber(){
        return n_;
    };
    
};



struct Interval {
private:
    double left_edge_;
    double right_edge_;
    bool square_left_;
    bool square_right_;
    vector<double> numbers_inside_;
    float frequency_;
public:
    Interval() : left_edge_(0), right_edge_(0), square_left_(false),  square_right_(false) {};
    Interval(double left_bracket, double right_bracket, bool square_left, bool square_right)
    : left_edge_(left_bracket), right_edge_(right_bracket), square_left_(square_left), square_right_(square_right){
        
    }
    bool includes(double number) {
        if (number < left_edge_ || number > right_edge_) {
            return false;
        }
        
        if (number == left_edge_ && !square_left_) {
            return false;
        }
        
        if (number == right_edge_ && !square_right_) {
            return false;
        }
        
        return true;
    }
    void intervalPushNumber(double number){
        if(this->includes(number)){
            this->numbers_inside_.push_back(number);
        }
    }
    vector<double> returnIntervalNumbers(){
        return this->numbers_inside_;
    }
    unsigned long returnSize(){
        return this->numbers_inside_.size();
    }
    void setFrequency(long long m){
        this->frequency_ = this->numbers_inside_.size()/static_cast<float>(m);
    }
    
    double returnFrequency(){
        return this->frequency_;
    }
    
    bool returnSquareLeft(){
        return this->square_left_;
    }
    bool returnSquareRight(){
        return this->square_right_;
    }
    double returnLeftEdge(){
        return this->left_edge_;
    }
    double returnRightEdge(){
        return this->right_edge_;
    }
};


struct IntervalMaker{
private:
    vector<float> extremes_;
public:
    IntervalMaker(vector<float>extremes): extremes_(extremes){
        
    }
    vector<Interval> makeIntervals(){
        vector<Interval> intervals;
        for(int i = 0; i<extremes_.size()-1; i++){
            intervals.push_back(i == 0 ? Interval(extremes_[i], extremes_[i+1], true, true): Interval(extremes_[i], extremes_[i+1], false, true)) ;
        }
        return intervals;
    }
};
struct ValueVector{
private:
    vector<double> vals_;
public:
    ValueVector()  {
        
    }
    void pushEvenlyValue(long long value, long long m){
        this->vals_.push_back(static_cast<double>(value)/static_cast<double>(m));
    }
    void pushValue(double value){
        this->vals_.push_back(value);
    }
    unsigned long returnVectorSize(){
        return vals_.size();
    }
    std::vector<double>::iterator begin() {
            return vals_.begin();
        }

    std::vector<double>::iterator end() {
            return vals_.end();
        }
    void insert(std::size_t position, double value) {
            if (position <= vals_.size()) {
                vals_.insert(vals_.begin() + position, value);
            }
        }
    void insert(std::size_t position, const std::vector<double>& values) {
            if (position <= vals_.size()) {
                vals_.insert(vals_.begin() + position, values.begin(), values.end());
            }
        }
    
};
class GeneratorBase {
    
public:
    long long LCM(long long m, long long X0, int a, int c){
        return (a*X0 +c) % m;
    }
    long long QCM(long long m, long long X0,  int d, int a, int c){
        return (d* X0 * X0 +a *X0 + c ) % m;
    }
    long long FNM(long long m, long long X1, long long X2 ){
        return (X1 + X2)%m;
    }
    long long mod_inv(long long param, long long m)
    {
        long long b0 = m, t, q;
        long long x0 = 0, x1 = 1;
        if(param==0) return INT_MAX;
        else if(param == INT_MAX) return 0 ;
        if (m == 1) return 1;
        while (param > 1)
        {
            q = param / m;
            t = m;
            m = param % m;
            param = t;
            t = x0;
            x0 = x1 - q * x0;
            x1 = t;
        }
        if (x1 < 0) x1 += b0;
        
        return x1;
    }
    
    long long ICG(long long m, int a, int c, long long param)
    {
        if (param == 0) return c;
        return (a * mod_inv(param, m) + c) % m;
    }
    void isIncluded(ValueVector result, long long m, IntervalEdges interval){
        vector<float> decades = interval.divideIntervals();
        vector<Interval> intervals_decades = IntervalMaker(decades).makeIntervals();
        for(auto& num: result){
            for(auto&  Interval : intervals_decades){
                Interval.intervalPushNumber(num);
            }
            for(auto& Interval : intervals_decades){
                Interval.setFrequency(m);
            }
        }
        
        printResult(intervals_decades);
    }
    void printResult(vector<Interval> intervals){
        cout << "\t"<< "Інтервал"<< "\t" << "Частота"<< endl;
        for(auto&  Interval : intervals){
            char leftBracket = Interval.returnSquareLeft() ? '[' : '(';
            char rightBracket = Interval.returnSquareRight() ? ']' : ')';
            cout << "\t"<< leftBracket << Interval.returnLeftEdge() << ", " << Interval.returnRightEdge() << rightBracket  << "\t"<< Interval.returnFrequency() << endl;
        }
    }
};




class LinearCongruentialMethod: public GeneratorBase {
private:
    const long long m_;
    
    
public:
    LinearCongruentialMethod(long m)
    : m_(m){}
    
    void linearCongruentialMethod() {
        long X0=1, X1;
        ValueVector linearCongruentialMethod_vector;
        linearCongruentialMethod_vector.pushEvenlyValue(X0, m_);
        for (int i = 1; i < m_; ++i) {
            X1 = LCM(m_, X0, 2, 3);
            X0 = X1;
            linearCongruentialMethod_vector.pushEvenlyValue(X1, m_);
        }
        GeneratorBase::isIncluded(linearCongruentialMethod_vector, m_, IntervalEdges(0, 1, 10));
    }
};


class QuadraticCongruentialMethod: public GeneratorBase {
private:
    const long long m_;
    
public:
    QuadraticCongruentialMethod(long m)
    : m_(m){}
    void quadraticCongruentialMethod(){
        long long X0 = 13, X1;
        ValueVector quadraticCongruentialMethod_vector;
        quadraticCongruentialMethod_vector.pushEvenlyValue(X0, m_);
        for (int i = 1; i < m_; ++i) {
            X1 = QCM(m_, X0, 1, 6, 3);
            X0 = X1;
            quadraticCongruentialMethod_vector.pushEvenlyValue(X1, m_);
        }
        GeneratorBase::isIncluded(quadraticCongruentialMethod_vector, m_, IntervalEdges(0, 1, 10));
    }
};


class FibonachiNumbersMethod: public GeneratorBase {
private:
    const long long m_;
    
    
public:
    FibonachiNumbersMethod(long m)
    : m_(m){}
    void fibonachiNumbersMethod(){
        long long X0= 0 % m_, X1 = 1 % m_, X2 = 0;
        ValueVector fibonachiNumbersMethod_vector;
        fibonachiNumbersMethod_vector.pushEvenlyValue(X0, m_);
        fibonachiNumbersMethod_vector.pushEvenlyValue(X1, m_);
        for (int i = 2; i < m_; ++i) {
            X2 = FNM(m_, X0, X1);
            X0 = X1;
            X1 = X2;
            fibonachiNumbersMethod_vector.pushEvenlyValue(X2, m_);
        }
        GeneratorBase::isIncluded(fibonachiNumbersMethod_vector, m_, IntervalEdges(0, 1, 10));
    }
};


class InverseCongruentialMethod: public GeneratorBase {
private:
    const long long m_;
    
    
public:
    InverseCongruentialMethod(long m)
    : m_(m){}
    void inverseCongruentialMethod(){
        long long X0= 1, X1= 0;
        ValueVector inverseCongruentialMethod_vector;
        inverseCongruentialMethod_vector.pushEvenlyValue(X1, m_);
        for (int i = 1; i < m_; ++i) {
            X1 = ICG(100003 , 2, 3, X0) ;
            X0 = X1;
            inverseCongruentialMethod_vector.pushEvenlyValue(X1, m_);

        }
        GeneratorBase::isIncluded(inverseCongruentialMethod_vector, m_, IntervalEdges(0, 1, 10));
        
    }
};


class UnionMethod: public GeneratorBase{
private:
    const long long m_;
public:
    UnionMethod(long m)
    : m_(m) {
        
    }
    
    void unionMethod() {
        long long X0 = 6, Y0 = 1, X1, Y1, Z ;
        ValueVector unionMethod_vector;
        unionMethod_vector.pushEvenlyValue(X0-Y0, m_);
        for (int i = 1; i < m_; ++i) {
            X1 = LCM(m_, X0, 4, 6);
            X0 = X1;
            Y1 = LCM(m_/1000, Y0, 6, 7);
            Y0 = Y1;
            Z = (X1 - Y1);
            if(Z<0)
                Z+=m_;
            
            Z%=m_;
            
            unionMethod_vector.pushEvenlyValue(Z, m_);
        }
        GeneratorBase::isIncluded(unionMethod_vector, m_, IntervalEdges(0, 1, 10));
    }
};




class SigmaMethod: public GeneratorBase{
private:
    long long m_;
    
public:
    SigmaMethod(long m)
    : m_(m){}
    
    void sigmaMethod() {
        float X0;
        long Y0 = 1, Y1;
        ValueVector sigmaMethod_vector;
        array<float, 12> selected_numbers;
        float sum=0;
        selected_numbers[0] = Y0/(m_-1);
        for (int i = 0; i < m_; ++i) {
            for(int i = 0; i < 12; ++i){
                Y1 = LCM(m_, Y0, 6, 7);
                Y0 = Y1;
                selected_numbers[i] = static_cast<float>(Y1)/static_cast<float>(m_-1);
                sum += selected_numbers[i];
            }
            
            X0 = sum - 6;
            sigmaMethod_vector.pushValue(X0);
            sum = 0;

        }
        
        GeneratorBase::isIncluded(sigmaMethod_vector, m_, IntervalEdges(-3, 3, 12));
    }
};
class PolarMethod: public GeneratorBase{
private:
    long long m_;
    
public:
    PolarMethod(long m)
    : m_(m){}
    void polarMehod(){
        long Y1 = 1, Y2 = 2, Y3, Y4, numSkipped= 0 ;
        float S,U1 ,U2, V1, V2, X1, X2;
        ValueVector polarMethod_vector;
        for(int i = 0; i < m_/2 + numSkipped; ++i){
            
            Y3 = LCM(m_, Y1, 6, 7);
            Y4 = LCM(m_, Y2, 2, 3);
            Y1 = Y3;
            Y2 = Y4;
            U1 = static_cast<float>(Y3)/static_cast<float>(m_);
            U2 = static_cast<float>(Y4)/static_cast<float>(m_);
            V1 = 2*U1-1;
            V2 = 2*U2-1;
            S = V1*V1 + V2*V2;
            if(S > 1){
                numSkipped+=1;
                continue;
            };
            float ln = -2*log(S);
            float expression1 = (ln)/S;
            float expression2 = (ln)/S;
            X1 = V1*sqrt(expression1);
            X2 = V2*sqrt(expression2);
            polarMethod_vector.insert(polarMethod_vector.returnVectorSize(), {X1, X2});
            polarMethod_vector.pushValue(X1);
            polarMethod_vector.pushValue(X2);
        }
        GeneratorBase::isIncluded(polarMethod_vector, m_, IntervalEdges(-3, 3, 12));

    }
};
class RelationMethod: public GeneratorBase{
private:
    long long m_;
    
public:
    RelationMethod(long m)
    : m_(m){}
    void relationMethod(){
        long Y1 = 1, Y2 = 2, Y3, Y4, numSkipped= 0;
        float U1, V1, X;
        ValueVector relationMethod_vector;
        for(int i = 0; i < m_/2+numSkipped; ++i){
            
            Y3 = LCM(m_, Y1, 6, 7);
            Y4 = LCM(m_, Y2, 2, 3);
            Y1 = Y3;
            Y2 = Y4;
            U1 = static_cast<float>(Y3)/static_cast<float>(m_);
            V1 = static_cast<float>(Y4)/static_cast<float>(m_);
            if(U1 == 0) {
                numSkipped+=1;
                continue;
                
            }
            V1 = static_cast<float>(Y4)/static_cast<float>(m_);
            X = sqrt(8.0/exp(1))*((V1-0.5)/U1);
            
            if(X*X <= 5-4*exp(0.25)*U1){
                relationMethod_vector.pushValue(X);
                continue;
            }
            else if (X*X >= ((4*exp(-1.35)/U1)+1.4)){
                numSkipped+=1;
                continue;
            }else if(X*X <= -4*log(U1)){
                relationMethod_vector.pushValue(X);
                
            }else{
                numSkipped+=1;
            }
            
        }
        GeneratorBase::isIncluded(relationMethod_vector, m_, IntervalEdges(-3, 3, 12));

    }
};



class LogarithmMethod: public GeneratorBase {
private:
    const long long m_;
    
    
public:
    LogarithmMethod(long m)
    : m_(m){}
    
    void logarithmMethod() {
        long Y0=1, Y1;
        float U0, X;
        ValueVector logarithmMethod_vector;
        
        for (int i = 1; i < m_; ++i) {
            Y1 = LCM(m_, Y0, 2, 3);
            Y0 = Y1;
            U0 = static_cast<float>(Y1)/static_cast<float>(m_);
            X = -14*log(U0);
            logarithmMethod_vector.pushValue(X);
        }
        GeneratorBase::isIncluded(logarithmMethod_vector, m_, IntervalEdges(0, 100, 10));
    }
};
class ArensMethod: public GeneratorBase {
private:
    const long long m_;
    
    
public:
    ArensMethod(long m)
    : m_(m){}
    
    void arensMethod() {
        long Z0=1, Z1, numSkipped= 0, H0=3, H1;
        const int a= 70;
        float U0, X,Y, V0 ;
        ValueVector arensMethod_vector;
        
        for (int i = 1; i < m_+numSkipped; ++i) {
            Z1 = LCM(m_, Z0, 2, 3);
            Z0 = Z1;
            U0 = static_cast<float>(Z1)/static_cast<float>(m_);
            Y = tan(M_PI*U0);
            X=sqrt(2*a-1)*Y+a-1;
            if(X<=0){
                numSkipped+=1;
                continue;
            }
            H1 = LCM(m_, H0, 2, 3);
            H0 = H1;
            V0 = static_cast<float>(H1)/static_cast<float>(m_);
            if (V0 > ((1+Y*Y)*exp(((a-1)*log(X/(a-1)-sqrt(2*a-1)*Y))))){
                numSkipped+=1;
                continue;
            }
            arensMethod_vector.pushValue(X);
        }
        GeneratorBase::isIncluded(arensMethod_vector, m_, IntervalEdges(0, 100, 10));
    }
};
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
            SigmaMethod generator(m);
            generator.sigmaMethod();
            break;
        }
        case 7:
        {
            PolarMethod generator(m);
            generator.polarMehod();
            break;
        }
        case 8:{
            RelationMethod generator(m);
            generator.relationMethod();
            break;
            
        }
        case 9:{
            LogarithmMethod generator(m);
            generator.logarithmMethod();
            break;
            
        }
        case 10:{
            ArensMethod generator(m);
            generator.arensMethod();
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

