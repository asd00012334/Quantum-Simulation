#include<bits/stdc++.h>
#include "Quantum.hpp"

using namespace std;


int main(){
    using namespace Quantum;
    srand(time(NULL)); /** Must do this for stochastic algorithm **/

    Qreg hp = Hadamard(Qreg(0,1));
    printf("|+> = ");
    cout << hp <<".\n\n";

    int DJ = DeutschJozsa([](int i){return i&1;}, 16);
    if(DJ == 1) puts("The function given is balanced.\n");
    else puts("The function given is constant.\n");

    vector<int> N = {8,15,23,33,35,39,44};
    vector<int> a = {3,2,7,5,4,11,13};
    /** a, N must be coprime **/
    for(int i=0;i<a.size();++i){
        int t0 = clock();
        int r = ShorPeriod(a[i],N[i]);
        printf("The order of %d mod %d is %d with high probability.\n",a[i],N[i],r);
        printf("Time taken: %d ms.\n\n", (clock()-t0)/(CLOCKS_PER_SEC/1000) );
    }

    vector<int> N2 = {6,12,15,21,28,33,35,39,44};
    /** N2 should not be prime power **/
    for(int i=0;i<N2.size();++i){
        int t0 = clock();
        int divisor = ShorFactor(N2[i]);
        printf("%d divides %d\n",divisor, N2[i]);
        printf("Time taken: %d ms.\n\n", (clock()-t0)/(CLOCKS_PER_SEC/1000) );

    }

}
