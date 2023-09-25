

#ifndef RandomGenerators_h
#define RandomGenerators_h

#pragma once 

#include <iostream>
#include <math.h>
#include <algorithm>
#include <limits>
#include <vector>
#include <climits>
#include <array>

using namespace std;
namespace randgenerator {
struct Interval {
private:
    double left_edge_;
    double right_edge_;
    bool square_left_;
    bool square_right_;
    vector<double> numbers_inside_;
    double frequency_;
public:
    Interval() : left_edge_(0), right_edge_(0), square_left_(false),  square_right_(false) {};
    Interval(double left_bracket, double right_bracket, bool square_left, bool square_right)
    : left_edge_(left_bracket), right_edge_(right_bracket), square_left_(square_left), square_right_(square_right){
        
    }
    bool includes(double number) {
        if (number < left_edge_ || number > right_edge_) {
            return false;
        }
        
        else if (number == left_edge_ && !square_left_) {
            return false;
        }
        
        else if (number == right_edge_ && !square_right_) {
            return false;
        }
        
        return true;
    }
    void intervalPushNumber(double number){
        if(includes(number)){
            numbers_inside_.push_back(number);
        }
    }
    vector<double> returnIntervalNumbers(){
        return numbers_inside_;
    }
    unsigned long long returnSize(){
        return numbers_inside_.size();
    }
    void setFrequency(long long m){
        frequency_ = returnSize()/static_cast<double>(m);
        
    }
    
    double returnFrequency(){
        return frequency_;
    }
    
    bool returnSquareLeft(){
        return square_left_;
    }
    bool returnSquareRight(){
        return square_right_;
    }
    double returnLeftEdge(){
        return left_edge_;
    }
    double returnRightEdge(){
        return right_edge_;
    }
};
struct IntervalEdges : public Interval {
private:
    int n_;
public:
    IntervalEdges(double left_edge, double right_edge, int n, bool square_left, bool square_right)
    : Interval(left_edge, right_edge, square_left, square_right), n_(n) {}
    
    vector<double> divideIntervals() {
        vector<double> extremes;
        double step = (returnRightEdge() - returnLeftEdge()) / static_cast<double>(n_);
        for (int i = 0; i <= n_; i++) {
            double value = returnLeftEdge() + i * step;
            extremes.push_back(value);
        }
        return extremes;
    }
    
    int returnSegmentsNumber() {
        return n_;
    }
    vector<Interval> makeIntervals(vector<double>extremes){
        vector<Interval> intervals;
        for(size_t i = 0; i<extremes.size()-1; i++){
            intervals.push_back(Interval(extremes[i], extremes[i+1], i == 0 ? returnSquareLeft(): false, i == extremes.size()-2 ? returnSquareRight(): true));
        }
        return intervals;
    }
    
};


struct NumberVector{
private:
    vector<double> vals_;
public:
    NumberVector()  {}
    void pushEvenlyValue(long long value, long long m){
        vals_.push_back(static_cast<double>(value)/static_cast<double>(m));
    }
    void pushHundredEvenlyValue(long long value, long long m){
        vals_.push_back(static_cast<double>(value)/static_cast<double>(m)*100);
    }
    void pushValue(double value){
        vals_.push_back(value);
    }
    unsigned long long returnVectorSize(){
        return vals_.size();
    }
    vector<double>::iterator begin() {
        return vals_.begin();
    }
    
    vector<double>::iterator end() {
        return vals_.end();
    }
    void insert(size_t position, vector<double> values) {
            if (position <=  this->vals_.size()) {
                this->vals_.insert(this->vals_.begin() + position, values.begin(), values.end());
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
    vector<Interval> makeSubInterval(vector<Interval> interval, double left_edge, double right_edge){
        vector<Interval> sub_interval;
        for(size_t i = 0; i < interval.size(); i++){
            if (interval[i].includes(left_edge+numeric_limits<double>::min())){
                while(!interval[i-1].includes(right_edge-numeric_limits<double>::min())){
                    sub_interval.push_back(interval[i]);
                    i+=1;
                    
                };
                break;
            }
        }
        
        return sub_interval;
    }
    vector<Interval> calcFrequency(NumberVector result, long long m, IntervalEdges edges){
        vector<double> decades = edges.divideIntervals();
        vector<Interval> intervals_decades = edges.makeIntervals(decades);
        for(auto& num: result){
            for(auto&  Interval : intervals_decades){
                Interval.intervalPushNumber(num);
            }
        }
        for(auto& Interval : intervals_decades){
            Interval.setFrequency(m);
        }
        return intervals_decades;
    }
    
    void printResult(vector<vector<Interval>> intervalsVector){
        
        for(auto& intervals: intervalsVector){
            cout << "\t"<< "Інтервал"<< "\t" << "Частота"<< endl;
            for(auto& Interval : intervals){
                char leftBracket = Interval.returnSquareLeft() ? '[' : '(';
                char rightBracket = Interval.returnSquareRight() ? ']' : ')';
                cout << "\t"<< leftBracket << Interval.returnLeftEdge() << ", " << Interval.returnRightEdge() << rightBracket  << "\t"<< Interval.returnFrequency() << endl;
                
            }
            cout << "\n" << endl;
            
        }
    }
    
};




class LinearCongruentialMethod: public GeneratorBase {
private:
    const long long m_;
    
    
public:
    LinearCongruentialMethod(long long m)
    : m_(m){}
    
    void linearCongruentialMethod() {
        long long X0=1, X1;
        NumberVector linearCongruentialMethod_vector;
        linearCongruentialMethod_vector.pushHundredEvenlyValue(X0, m_);
        for (size_t i = 1; i < m_; ++i) {
            X1 = LCM(m_, X0, 2, 3);
            X0 = X1;
            linearCongruentialMethod_vector.pushHundredEvenlyValue(X1, m_);
        }
        vector<Interval> vals = calcFrequency(linearCongruentialMethod_vector, m_, IntervalEdges(0, 100, 10, true, true));
        printResult(vector<vector<Interval>>{vals});
    }
};


class QuadraticCongruentialMethod: public GeneratorBase {
private:
    const long long m_;
    
public:
    QuadraticCongruentialMethod(long long m)
    : m_(m){}
    void quadraticCongruentialMethod(){
        long long X0 = 13, X1;
        NumberVector quadraticCongruentialMethod_vector;
        quadraticCongruentialMethod_vector.pushHundredEvenlyValue(X0, m_);
        for (size_t i = 1; i < m_; ++i) {
            X1 = QCM(m_, X0, 1, 6, 3);
            X0 = X1;
            quadraticCongruentialMethod_vector.pushHundredEvenlyValue(X1, m_);
        }
        vector<Interval> vals = calcFrequency(quadraticCongruentialMethod_vector, m_, IntervalEdges(0, 100, 10, true, true));
        printResult(vector<vector<Interval>>{vals});
    }
};


class FibonachiNumbersMethod: public GeneratorBase {
private:
    const long long m_;
    
    
public:
    FibonachiNumbersMethod(long long m)
    : m_(m){}
    void fibonachiNumbersMethod(){
        long long X0= 0 % m_, X1 = 1 % m_, X2 = 0;
        NumberVector fibonachiNumbersMethod_vector;
        fibonachiNumbersMethod_vector.pushHundredEvenlyValue(X0, m_);
        fibonachiNumbersMethod_vector.pushHundredEvenlyValue(X1, m_);
        for (size_t i = 2; i < m_; ++i) {
            X2 = FNM(m_, X0, X1);
            X0 = X1;
            X1 = X2;
            fibonachiNumbersMethod_vector.pushHundredEvenlyValue(X2, m_);
        }
        vector<Interval> vals = calcFrequency(fibonachiNumbersMethod_vector, m_, IntervalEdges(0, 100, 10, true, true));
        printResult(vector<vector<Interval>>{vals});
    }
};


class InverseCongruentialMethod: public GeneratorBase {
private:
    const long long m_;
    
    
public:
    InverseCongruentialMethod(long long m)
    : m_(m){}
    void inverseCongruentialMethod(){
        long long X0= 1, X1= 0;
        NumberVector inverseCongruentialMethod_vector;
        inverseCongruentialMethod_vector.pushHundredEvenlyValue(X1, m_);
        for (size_t i = 1; i < m_; ++i) {
            X1 = ICG(1000003 , 2, 3, X0) ;
            X0 = X1;
            inverseCongruentialMethod_vector.pushHundredEvenlyValue(X1, m_);
            
        }
        vector<Interval> vals = calcFrequency(inverseCongruentialMethod_vector, m_, IntervalEdges(0, 100, 10, true, true));
        printResult(vector<vector<Interval>>{vals});
    }
};


class UnionMethod: public GeneratorBase{
private:
    const long long m_;
public:
    UnionMethod(long long m)
    : m_(m) {
        
    }
    
    void unionMethod() {
        long long X0 = 6, Y0 = 1, X1, Y1, Z ;
        NumberVector unionMethod_vector;
        unionMethod_vector.pushHundredEvenlyValue(X0-Y0, m_);
        for (size_t i = 1; i < m_; ++i) {
            X1 = LCM(m_, X0, 4, 6);
            X0 = X1;
            Y1 = LCM(m_/1000, Y0, 6, 7);
            Y0 = Y1;
            Z = (X1 - Y1);
            if(Z<0)
                Z+=m_;
            
            Z%=m_;
            
            unionMethod_vector.pushHundredEvenlyValue(Z, m_);
        }
        vector<Interval> vals = calcFrequency(unionMethod_vector, m_, IntervalEdges(0, 100, 10, true, true));
        printResult(vector<vector<Interval>>{vals});
    }
};




class SigmaMethod: public GeneratorBase{
private:
    long long m_;
    
public:
    SigmaMethod(long long m)
    : m_(m){}
    
    void sigmaMethod() {
        double X0;
        long long Y0 = 1, Y1;
        NumberVector sigmaMethod_vector;
        array<double, 12> selected_numbers;
        double sum=0;
        for (size_t i = 0; i < m_; ++i) {
            for(size_t i = 0; i < 12; ++i){
                Y1 = LCM(m_, Y0, 6, 7);
                Y0 = Y1;
                selected_numbers[i] = static_cast<double>(Y1)/static_cast<double>(m_-1);
                sum += selected_numbers[i];
            }
            
            X0 = sum - 6;
            sigmaMethod_vector.pushValue(X0);
            sum = 0;
            
        }
        
        vector<Interval> vals = calcFrequency(sigmaMethod_vector, m_, IntervalEdges(-3, 3, 12, true, true));
        vector<Interval> sub_interval = makeSubInterval(vals, 0, 1);
        printResult(vector<vector<Interval>>{vals, sub_interval});
    }
};
class PolarMethod: public GeneratorBase{
private:
    long long m_;
    
public:
    PolarMethod(long long m)
    : m_(m){}
    void polarMehod(){
        long long Y0 = 1, Y1 = 2, Y2, Y3, numSkipped= 0 ;
        double S,U1 ,U2, V1, V2, X1, X2;
        NumberVector polarMethod_vector;
        for(size_t i = 0; i < m_/2 + numSkipped; ++i){
            
            Y2 = LCM(m_, Y0, 6, 7);
            Y3 = LCM(m_, Y1, 2, 3);
            Y0 = Y2;
            Y1 = Y3;
            U1 = static_cast<double>(Y2)/static_cast<double>(m_);
            U2 = static_cast<double>(Y3)/static_cast<double>(m_);
            V1 = 2*U1-1;
            V2 = 2*U2-1;
            S = V1*V1 + V2*V2;
            if(S > 1){
                numSkipped+=1;
                continue;
            };
            double ln = -2*log(S);
            double expression = (ln)/S;
            X1 = V1*sqrt(expression);
            X2 = V2*sqrt(expression);
            polarMethod_vector.insert(polarMethod_vector.returnVectorSize(), {X1, X2});
            
        }
        vector<Interval> vals = calcFrequency(polarMethod_vector, m_, IntervalEdges(-3, 3, 12, true, true));
        vector <Interval> sub_interval = makeSubInterval(vals, 0, 1);
        printResult(vector<vector<Interval>>{vals, sub_interval});
        
    }
};
class RelationMethod: public GeneratorBase{
private:
    long long m_;
    
public:
    RelationMethod(long long m)
    : m_(m){}
    void relationMethod(){
        long long Y0 = 1, Y1 = 2, Y2, Y3, numSkipped= 0;
        double U1, V1, X;
        NumberVector relationMethod_vector;
        for(size_t i = 0; i < m_+numSkipped; ++i){
            
            Y2 = LCM(m_, Y0, 6, 7);
            Y3 = LCM(m_, Y1, 2, 3);
            Y0 = Y2;
            Y1 = Y3;
            U1 = static_cast<double>(Y2)/static_cast<double>(m_);
            V1 = static_cast<double>(Y3)/static_cast<double>(m_);
            if(U1 == 0) {
                numSkipped+=1;
                continue;
                
            }
            V1 = static_cast<double>(Y3)/static_cast<double>(m_);
            X = sqrt(8.0/exp(1))*((V1-0.5)/U1);
            
            if(X*X <= (5-4*exp(0.25)*U1)){
                relationMethod_vector.pushValue(X);
                continue;
            }
            else if (X*X >= ((4*exp(-1.35)/U1)+1.4)){
                numSkipped+=1;
                continue;
            }else if(X*X <= (-4*log(U1))){
                relationMethod_vector.pushValue(X);
                
            }else{
                numSkipped+=1;
            }
            
        }
        vector<Interval> vals = calcFrequency(relationMethod_vector, m_, IntervalEdges(-3, 3, 12, true, true));
        vector <Interval> sub_interval = makeSubInterval(vals, 0, 1);
        printResult(vector<vector<Interval>>{vals, sub_interval});
    }
};



class LogarithmMethod: public GeneratorBase {
private:
    const long long m_;
    
    
public:
    LogarithmMethod(long long m)
    : m_(m){}
    
    void logarithmMethod() {
        long long Y0=1, Y1;
        double U0, X;
        NumberVector logarithmMethod_vector;
        
        for (size_t i = 1; i < m_; ++i) {
            Y1 = LCM(m_, Y0, 2, 3);
            Y0 = Y1;
            U0 = static_cast<double>(Y1)/static_cast<double>(m_);
            X = -14*log(U0);
            logarithmMethod_vector.pushValue(X);
        }
        vector<Interval> vals = calcFrequency(logarithmMethod_vector, m_, IntervalEdges(0, 100, 10, true, true));
        printResult(vector<vector<Interval>>{vals});
    }
};
class ArensMethod: public GeneratorBase {
private:
    const long long m_;
    
    
public:
    ArensMethod(long long m)
    : m_(m){}
    
    void arensMethod() {
        long long Z0=1, Z1, numSkipped= 0, H0=3, H1;
        const int a= 40;
        double U0, X,Y, V0 ;
        NumberVector arensMethod_vector;
        
        for (size_t i = 1; i < m_+numSkipped; ++i) {
            Z1 = LCM(m_, Z0, 2, 3);
            Z0 = Z1;
            U0 = static_cast<double>(Z1)/static_cast<double>(m_);
            Y = tan(M_PI*U0);
            double extractedExpr = sqrt(2*a-1)*Y;
            X=extractedExpr+a-1;
            if(X<=0){
                numSkipped+=1;
                continue;
            }
            H1 = LCM(m_, H0, 6, 7);
            H0 = H1;
            V0 = static_cast<double>(H1)/static_cast<double>(m_);
            if (V0 > ((1+Y*Y)*exp((a-1)*log(X/(a-1))-extractedExpr*Y))){
                numSkipped+=1;
                continue;
            }
            arensMethod_vector.pushValue(X);
        }
        vector<Interval> vals = calcFrequency(arensMethod_vector, m_, IntervalEdges(0, 100, 10, true, true));
        printResult(vector<vector<Interval>>{vals});
    }
};
}



#endif
