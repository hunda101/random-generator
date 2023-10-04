#include <iostream>
#include <math.h>
#include <algorithm>
#include <limits>
#include <vector>
#include <climits>
#include <array>
#include "RandomGenerators.h"

using namespace std;
using namespace randgenerator;

void print_menu();



int main() {
    print_menu();
    int type;
    long long m = 1000000;
    cin >> type;
    
    switch (type) {
        case 1:{
            LinearCongruentialMethod generator;
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
            UnionMethod generator;
            generator.unionMethod();
            break;
        }
        case 6:
        {
            SigmaMethod generator;
            generator.sigmaMethod();
            break;
        }
        case 7:
        {
            PolarMethod generator;
            generator.polarMehod();
            break;
        }
        case 8:{
            RelationMethod generator;
            generator.relationMethod();
            break;
            
        }
        case 9:{
            LogarithmMethod generator;
            generator.logarithmMethod();
            break;
            
        }
        case 10:{
            ArensMethod generator;
            generator.arensMethod();
            break;
            
        }
        default:{
            cout << "Wrong input!!!\nOnly 1-10 requires";
        }
    }
    
    return 1;
    
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

