// FU11 - Exercise 01 - s1312004

#include<bits/stdc++.h>
#include<iomanip>

using namespace std;

template <typename T>
T getNegativeMachinePrecision(T radix){
    int n = 1;
    // starting with r^-1
    T radix_power = 1.0 / radix;

    while( (1.0 - radix_power ) - 1.0 != 0.0 ){
        radix_power /= radix;
        n++;
    }

    return radix_power;
}

template<typename T>
T getLargestNumber(T radix){
    T epsilon = getNegativeMachinePrecision<double>(radix);
    T f = 1.0 - radix * epsilon;
    T prev_f = f;

    while( !isinf(f) ) {
        prev_f = f;
        f *= radix;
    }

    return prev_f;
}

int main(){
    cout << "Floating-point machine parameters \n" << "--------------------------------------------------" << endl; 

    cout << setprecision(16);

    cout << "Negative machine precision = " << getNegativeMachinePrecision<double>(2) << endl;
    cout << "Largest number = " << getLargestNumber<double>(2) << endl;

    return 0;
}