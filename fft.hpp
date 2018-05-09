#ifndef __FFT_HPP__
#define __FFT_HPP__

#include <bits/stdc++.h>

using namespace std;

#define BIT_CNT 21

int directReverse(int,int);
int bitReverse(int,int);

int directReverse(int index, int n){
    static int sub[(1<<BIT_CNT)+1]={0,1,1,3};
    if(sub[index] || index<4) return sub[index];
    return sub[index]=bitReverse(index>>1,n-1)|(index&1)<<(n-1);
}

int bitReverse(int index, int n){
    if(!index) return 0;
    for(int m=n;;--m)
        if(index&(1<<(m-1)))
            return directReverse(index,m)<<(n-m);
}

template<typename iter>
void bitShuffle(iter begin, int size){
    /// Assume end-begin is 2 power
    int N = __lg(size);
    for(int i=0;i<size;++i){
        int tmp = bitReverse(i,N);
        if(i<tmp) swap(begin[i],begin[tmp]);
    }
}

template<typename iter>
void fft(iter begin, iter end){
    /// Assume end-begin is 2 power
    typedef typename iterator_traits<iter>::value_type type;
    type sqt2 = sqrt(2);
    int size = end-begin;
    bitShuffle(begin,size);
    for(int hsize=1;hsize<size;hsize<<=1){
        for(iter base=begin;base<end;base+=hsize<<1){
            iter mid = base+hsize;
            static const double pi = acos(-1);
            type w = 1, dw = exp(type(0,-pi/hsize));
            for(int k=0;k<hsize;k++, w*=dw){
                type even = base[k], odd = w*mid[k];
                base[k] = (even+odd)/sqt2;
                mid[k] = (even-odd)/sqt2;
            }
        }
    }
}

template<typename iter>
void ifft(iter begin, iter end){
    /// Assume end-begin is 2 power
    typedef typename iterator_traits<iter>::value_type type;
    type sqt2 = sqrt(2);
    int size = end-begin;
    bitShuffle(begin,size);
    for(int hsize=1;hsize<size;hsize<<=1){
        for(iter base=begin;base<end;base+=hsize<<1){
            iter mid = base+hsize;
            static const double pi = acos(-1);
            type w = 1, dw = exp(type(0,pi/hsize));
            for(int k=0;k<hsize;k++, w*=dw){
                type even = base[k], odd = w*mid[k];
                base[k] = (even+odd)/sqt2;
                mid[k] = (even-odd)/sqt2;
            }
        }
    }
}

#endif // __FFT_HPP__
