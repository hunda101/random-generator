//
//  RandomGenerators.h
//  project1
//
//  Created by Mark Khomenko on 19.09.2023.
//

#ifndef RandomGenerators_h
#define RandomGenerators_h

#pragma once // Защита от множественного включения

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
struct IntervalEdges : public Interval {
private:
    int n_;
public:
    IntervalEdges(double left_edge, double right_edge, int n, bool square_left, bool square_right)
    : Interval(left_edge, right_edge, square_left, square_right), n_(n) {}
    
    vector<float> divideIntervals() {
        vector<float> extremes;
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
    vector<Interval> makeIntervals(vector<float>extremes){
        vector<Interval> intervals;
        for(int i = 0; i<extremes.size()-1; i++){
            intervals.push_back(Interval(extremes[i], extremes[i+1], i == 0 ? returnSquareLeft(): false, i == extremes.size()-2 ? returnSquareRight(): true));
        }
        return intervals;
    }
    
};


struct NumberVector{
private:
    vector<double> vals_;
public:
    NumberVector()  {
        
    }
    void pushEvenlyValue(long long value, long long m){
        this->vals_.push_back(static_cast<double>(value)/static_cast<double>(m));
    }
    void pushValue(double value){
        this->vals_.push_back(value);
    }
    unsigned long returnVectorSize(){
        return this->vals_.size();
    }
    vector<double>::iterator begin() {
        return this->vals_.begin();
    }
    
    vector<double>::iterator end() {
        return this->vals_.end();
    }
    void insert(size_t position, double value) {
        if (position <= this->vals_.size()) {
            this->vals_.insert(this->vals_.begin() + position, value);
        }
    }
    void insert(size_t position, const vector<double>& values) {
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
    vector<Interval> populateAndCalculateFrequencies(NumberVector result, long long m, IntervalEdges edges){
        vector<float> decades = edges.divideIntervals();
        vector<Interval> intervals_decades = edges.makeIntervals(decades);
        for(auto& num: result){
            for(auto&  Interval : intervals_decades){
                Interval.intervalPushNumber(num);
            }
            for(auto& Interval : intervals_decades){
                Interval.setFrequency(m);
            }
        }
        
        return intervals_decades;
    }
    vector<Interval> makeSubInterval(vector<Interval> parent_interval, double left_edge, double right_edge){
        vector<Interval> sub_interval;
        for(int i = 0; i < parent_interval.size(); i++){
            if (parent_interval[i].includes(left_edge+numeric_limits<float>::min())){
                while(!parent_interval[i-1].includes(right_edge-numeric_limits<float>::min())){
                    sub_interval.push_back(parent_interval[i]);
                    i+=1;
                    
                };
                break;
            }
        }
        
        return sub_interval;
    }
    void printResult(vector<vector<Interval>> intervalsVector){
        
        for(auto& intervals: intervalsVector){
            cout << "\t"<< "Інтервал"<< "\t" << "Частота"<< endl;
            for(size_t i = 0; i < intervals.size(); ++i){
                char leftBracket = intervals[i].returnSquareLeft() ? '[' : '(';
                char rightBracket = intervals[i].returnSquareRight() ? ']' : ')';
                cout << "\t"<< leftBracket << intervals[i].returnLeftEdge() << ", " << intervals[i].returnRightEdge() << rightBracket  << "\t"<< intervals[i].returnFrequency() << endl;
                
            }
            cout << "\n" << endl;
            
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
        NumberVector linearCongruentialMethod_vector;
        linearCongruentialMethod_vector.pushEvenlyValue(X0, m_);
        for (int i = 1; i < m_; ++i) {
            X1 = LCM(m_, X0, 2, 3);
            X0 = X1;
            linearCongruentialMethod_vector.pushEvenlyValue(X1, m_);
        }
        vector<Interval> vals = populateAndCalculateFrequencies(linearCongruentialMethod_vector, m_, IntervalEdges(0, 1, 10, true, true));
        printResult(vector<vector<Interval>>{vals});
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
        NumberVector quadraticCongruentialMethod_vector;
        quadraticCongruentialMethod_vector.pushEvenlyValue(X0, m_);
        for (int i = 1; i < m_; ++i) {
            X1 = QCM(m_, X0, 1, 6, 3);
            X0 = X1;
            quadraticCongruentialMethod_vector.pushEvenlyValue(X1, m_);
        }
        vector<Interval> vals = populateAndCalculateFrequencies(quadraticCongruentialMethod_vector, m_, IntervalEdges(0, 1, 10, true, true));
        printResult(vector<vector<Interval>>{vals});
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
        NumberVector fibonachiNumbersMethod_vector;
        fibonachiNumbersMethod_vector.pushEvenlyValue(X0, m_);
        fibonachiNumbersMethod_vector.pushEvenlyValue(X1, m_);
        for (int i = 2; i < m_; ++i) {
            X2 = FNM(m_, X0, X1);
            X0 = X1;
            X1 = X2;
            fibonachiNumbersMethod_vector.pushEvenlyValue(X2, m_);
        }
        vector<Interval> vals = populateAndCalculateFrequencies(fibonachiNumbersMethod_vector, m_, IntervalEdges(0, 1, 10, true, true));
        printResult(vector<vector<Interval>>{vals});
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
        NumberVector inverseCongruentialMethod_vector;
        inverseCongruentialMethod_vector.pushEvenlyValue(X1, m_);
        for (int i = 1; i < m_; ++i) {
            X1 = ICG(100003 , 2, 3, X0) ;
            X0 = X1;
            inverseCongruentialMethod_vector.pushEvenlyValue(X1, m_);
            
        }
        vector<Interval> vals = populateAndCalculateFrequencies(inverseCongruentialMethod_vector, m_, IntervalEdges(0, 1, 10, true, true));
        printResult(vector<vector<Interval>>{vals});
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
        NumberVector unionMethod_vector;
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
        vector<Interval> vals = populateAndCalculateFrequencies(unionMethod_vector, m_, IntervalEdges(0, 1, 10, true, true));
        printResult(vector<vector<Interval>>{vals});
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
        NumberVector sigmaMethod_vector;
        array<float, 12> selected_numbers;
        float sum=0;
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
        
        vector<Interval> vals = populateAndCalculateFrequencies(sigmaMethod_vector, m_, IntervalEdges(-3, 3, 12, true, false));
        vector<Interval> sub_interval = makeSubInterval(vals, 0, 1);
        printResult(vector<vector<Interval>>{vals, sub_interval});
    }
};
class PolarMethod: public GeneratorBase{
private:
    long long m_;
    
public:
    PolarMethod(long m)
    : m_(m){}
    void polarMehod(){
        long Y0 = 1, Y1 = 2, Y2, Y3, numSkipped= 0 ;
        float S,U1 ,U2, V1, V2, X1, X2;
        NumberVector polarMethod_vector;
        for(int i = 0; i < m_/2 + numSkipped; ++i){
            
            Y2 = LCM(m_, Y0, 6, 7);
            Y3 = LCM(m_, Y1, 2, 3);
            Y0 = Y2;
            Y1 = Y3;
            U1 = static_cast<float>(Y2)/static_cast<float>(m_);
            U2 = static_cast<float>(Y3)/static_cast<float>(m_);
            V1 = 2*U1-1;
            V2 = 2*U2-1;
            S = V1*V1 + V2*V2;
            if(S > 1){
                numSkipped+=1;
                continue;
            };
            float ln = -2*log(S);
            float expression = (ln)/S;
            X1 = V1*sqrt(expression);
            X2 = V2*sqrt(expression);
            polarMethod_vector.insert(polarMethod_vector.returnVectorSize(), {X1, X2});
            
        }
        vector<Interval> vals = populateAndCalculateFrequencies(polarMethod_vector, m_, IntervalEdges(-3, 3, 12, true, true));
        vector<Interval> vals1 = populateAndCalculateFrequencies(polarMethod_vector, m_, IntervalEdges(0, 1, 10, true, true));
        vector <Interval> sub_interval = makeSubInterval(vals, 0, 1);
        printResult(vector<vector<Interval>>{vals,vals1, sub_interval});
        
    }
};
class RelationMethod: public GeneratorBase{
private:
    long long m_;
    
public:
    RelationMethod(long m)
    : m_(m){}
    void relationMethod(){
        long Y0 = 1, Y1 = 2, Y2, Y3, numSkipped= 0;
        float U1, V1, X;
        NumberVector relationMethod_vector;
        for(int i = 0; i < m_+numSkipped; ++i){
            
            Y2 = LCM(m_, Y0, 6, 7);
            Y3 = LCM(m_, Y1, 2, 3);
            Y0 = Y2;
            Y1 = Y3;
            U1 = static_cast<float>(Y2)/static_cast<float>(m_);
            V1 = static_cast<float>(Y3)/static_cast<float>(m_);
            if(U1 == 0) {
                numSkipped+=1;
                continue;
                
            }
            V1 = static_cast<float>(Y3)/static_cast<float>(m_);
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
        vector<Interval> vals = populateAndCalculateFrequencies(relationMethod_vector, m_, IntervalEdges(-3, 3, 12, true, true));
        vector <Interval> sub_interval = makeSubInterval(vals, 0, 1);
        printResult(vector<vector<Interval>>{vals, sub_interval});
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
        NumberVector logarithmMethod_vector;
        
        for (int i = 1; i < m_; ++i) {
            Y1 = LCM(m_, Y0, 2, 3);
            Y0 = Y1;
            U0 = static_cast<float>(Y1)/static_cast<float>(m_);
            X = -14*log(U0);
            logarithmMethod_vector.pushValue(X);
        }
        vector<Interval> vals = populateAndCalculateFrequencies(logarithmMethod_vector, m_, IntervalEdges(0, 100, 10, true, true));
        printResult(vector<vector<Interval>>{vals});
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
        const int a= 50;
        float U0, X,Y, V0 ;
        NumberVector arensMethod_vector;
        
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
        vector<Interval> vals = populateAndCalculateFrequencies(arensMethod_vector, m_, IntervalEdges(0, 100, 10, true, true));
        printResult(vector<vector<Interval>>{vals});
    }
};
}



#endif /* RandomGenerators_h */
