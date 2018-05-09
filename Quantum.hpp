#ifndef __QUANTUM_HPP__
#define __QUANTUM_HPP__

#include<bits/stdc++.h>
#include"tensor.hpp"
#include"fft.hpp"

namespace Quantum{

using namespace std;

class Gate;
class Qreg;

class Qreg{
    int const qubit;
    tensor state;
public:

    inline double getL2()const{
        double l2 = 0;
        for(int i=0;i<state.cols();++i)
            l2+=abs(state(0,i))*abs(state(0,i));
        return l2;
    }
    Qreg(int qubit):qubit(qubit), state(1,1<<qubit){assert(qubit>0);}
    Qreg(int init, int qubit):qubit(qubit),state(1,1<<qubit){
        state(0,mask(init,0,qubit)) = 1;
    };
    Qreg(Vec const& state, int qubit):state({state},1,1<<qubit), qubit(qubit){
       checkL2();
    }
    inline cpx& operator[](size_t i){ return state(0,i);}
    inline cpx const& operator[](size_t i)const{ return state(0,i);}
    inline size_t bitNum()const{ return qubit;}
    inline size_t stateNum()const{ return 1<<qubit;}
    static inline int mask(int x,int s, int t){
        return (x&~(-1<<t))>>s;
    }

    vector<pair<Qreg,Qreg> > split(int p)const{
        /// split into |[0,p)> x |[p,end)>
        auto tmp = state.factor({1,1<<(qubit-p)});
        vector<pair<Qreg,Qreg> > out;
        for(int i=0;i<tmp.size();++i)
            out.push_back( {mkQreg(tmp[i].first), mkQreg(tmp[i].second)} );
        return out;
    }
    vector<pair<Qreg,Qreg> > slice(int s, int t)const{
        /// slice out |[s,t)>
        Qreg perm(qubit);
        int totalStat = 1<<qubit;
        for(int i=0;i<totalStat;++i){
            int a = mask(i,t,qubit);
            int b = mask(i,s,t);
            int c = mask(i,0,s);
            int sel = b<<(qubit-t+s) | a<<s | c;
            perm[sel] = state(0,i);
        }
        return perm.split(qubit-t+s);
    }

    friend Qreg operator*(Qreg const& l, Qreg const& r){
        Qreg out(l.bitNum()+r.bitNum());
        out.state = tprod(l.state,r.state);
        return out;
    }

    friend Qreg operator+(Qreg const& l, Qreg const& r){
        assert(l.bitNum()==r.bitNum());
        Qreg out(l.bitNum());
        for(int i=0;i<l.stateNum();++i)
            out[i] = l[i] + r[i];
        return out;
    }

    friend ostream& operator<<(ostream& sout, Qreg const& reg){
        int term = 0;
        double l2 = reg.getL2();
        for(int i=0;i<(1<<reg.qubit);++i){
            if(abs(reg[i])<eps) continue;
            if(term) sout<<" + "; else term = 1;
            sout<<reg[i]<<" ";
            sout<<"|"<<i<<">";
        }
        return sout;
    }

    static Qreg mkQreg(tensor const& t){
        Qreg out(__lg(t.cols()));
        assert(1<<out.bitNum()==t.cols() && t.rows()==1);
        out.state = t;
        return out;
    }

    int measure(int s, int t){
        double p = (double)rand()/RAND_MAX;
        int sel;
        for(sel=0;sel<stateNum()-1;++sel){
            p-=abs(state(0,sel))*abs(state(0,sel));
            if(p<0) break;
        }
        sel = mask(sel,s,t);
        for(int i=0;i<stateNum();++i)
            if(mask(i,s,t)^sel) state(0,i) = 0;
        normalize();
        return sel;
    }

    void checkL2()const{
        double l2 = getL2();
        assert(abs(l2-1)<eps);
    }

    void normalize(){
        double sqtl2 = sqrt(getL2());
        for(int i=0;i<state.cols();++i)
            state(0,i)/=sqtl2;
    }

    Qreg& operator=(Qreg const& r){
        assert(qubit == r.qubit);
        state = r.state;
        return *this;
    }

    friend Qreg applySOP(vector<pair<Qreg,Qreg> > const& sop, function<Qreg(Qreg)>const& l, function<Qreg(Qreg)>const& r){
        assert(sop.size()>0);
        Qreg out(sop[0].first.bitNum()+sop[0].second.bitNum());
        for(int i=0;i<sop.size();++i)
            out = out + l(sop[i].first)*r(sop[i].second);
        return out;
    }

    Qreg applySlice(int s, int t, function<Qreg(Qreg)>const& f)const{
        return applySOP(slice(s,t),f,[](Qreg const& r){return r;});
    }

    friend Gate;
};

class Gate{
    int qubit;
    tensor mat;
public:
    inline size_t bitNum()const{ return qubit;}
    inline size_t stateNum()const{ return 1<<qubit;}
    Gate(int qubit):qubit(qubit), mat(1<<qubit,1<<qubit){}
    Gate(Mat mat, int qubit=0): mat(mat,qubit?1<<qubit:0,qubit?1<<qubit:0){
        assert(Gate::mat.cols()>1 && Gate::mat.rows()>1);
        assert(((Gate::mat.rows()-1)&Gate::mat.rows()) == 0);
        assert(((Gate::mat.cols()-1)&Gate::mat.cols()) == 0);
        assert((Gate::mat.rows() == Gate::mat.cols()));
        Gate::qubit = __lg(Gate::mat.rows());
    }

    friend Gate operator*(Gate const& l, Gate const& r){
        Gate out(l.bitNum()+r.bitNum());
        out.mat = tprod(l.mat,r.mat);
        return out;
    }

    Gate operator()(Gate const& r)const{
        assert(bitNum()==r.bitNum());
        Gate out(bitNum());
        out.mat = tmul(mat,r.mat);
        return out;
    }

    Qreg operator()(Qreg const& r)const{
        auto mT = mat.transpose();
        Qreg out = Qreg::mkQreg(tmul(r.state,mT));
        return out;
    }

    friend ostream& operator<<(ostream& sout, Gate const& r){
        sout<<r.mat;
    }

};

static const Gate
    H({
        {1/sqrt(2),1/sqrt(2)},
        {1/sqrt(2),-1/sqrt(2)}
    }), CNOT({
        {1,0,0,0},
        {0,1,0,0},
        {0,0,0,1},
        {0,0,1,0}
    });

Qreg Hadamard(Qreg const& r){
    static vector<Gate> Hn;
    static int init=1;
    if(init){
        init = 0;
        Hn.push_back(H);
        for(int i=1;i<=10;++i)
            Hn.push_back(Hn.back()*H);
    }
    if(r.bitNum()<=Hn.size())
        return Hn[r.bitNum()-1](r);
    assert(r.bitNum()<=20);
    return applySOP(r.slice(Hn.size(),r.bitNum()),Hn[r.bitNum()-Hn.size()-1],Hn.back());
}

Qreg QFT(Qreg const& r){
    Qreg out(r);
    assert(r.bitNum()<=BIT_CNT);
    fft(&out[0], &out[1<<out.bitNum()]); /// !! notice: access violation?
    return out;
}

Qreg oracle(Qreg const& addrReg, Qreg const& valReg, function<int(int)>const& f){
    Qreg in = addrReg*valReg;
    Qreg out(in.bitNum());
    for(int i=0;i<out.stateNum();++i){
        int addr = i>>valReg.bitNum();
        int val = i & ((1<<valReg.bitNum())-1);
        out[addr<<valReg.bitNum() | f(addr)^val] += in[i];
    }
    out.checkL2();
    out.normalize();
    return out;
}

bool DeutschJozsa(function<int(int)>const& f, int const bitNum){
    /// Given f(x) being either constant of balanced
    /// determine which
    /// return true if balanced
    /// return false if constant
    Qreg init = Hadamard(Qreg(0,bitNum));
    Qreg neg = H(Qreg(1,1));
    Qreg result = Hadamard(oracle(init,neg,f));
    return result.measure(1,bitNum+1)!=0;
}

inline vector<int> ContFrac(int entr, int dntr){
    int g = __gcd(entr,dntr);
    entr /= g, dntr /= g;
    vector<int> out;
    while(dntr){
        out.push_back(entr/dntr);
        entr%=dntr;
        swap(entr,dntr);
    }
    return out;
}

inline int approxPeriod(int a, int N, int b, int Q){
    auto aVec = ContFrac(b,Q);
    int hm1 = 1, hm2 = 0;
    int km1 = 0, km2 = 1;
    for(auto &e: aVec){
        hm2 = e*hm1 + hm2;
        km2 = e*km1 + km2;
        swap(hm2,hm1);
        swap(km2,km1);
        int g = __gcd(km1,hm1);
        km1 /= g, hm1 /=g;
        // cout<<hm1<<"/"<<km1<<", ";
        if(km1>=N) return km2;
    }
    return km1;
}

inline int modPow(int b, int e, int M){
    int val = 1;
    for(;e;b=b*b%M,e>>=1)
        if(e&1) val = val*b%M;
    return val;
}

int ShorPeriod(int a, int N){
    while(1){
        int q = __lg(N*N);
        if(1<<q < N*N) ++q;
        int Q = 1<<q, valQ = __lg(N);
        if(1<<valQ < N) ++valQ;
        Qreg init = Hadamard(Qreg(0,q));
        function<int(int)> f = [=](int x){
            return modPow(a,x,N);
        };
        auto mid = oracle(init, Qreg(0,valQ),f);
        Qreg out = mid.applySlice(valQ,q+valQ,QFT);

        int b = out.measure(valQ,q+valQ);

        int r = approxPeriod(a,N,b,Q);
        if(r>0 && f(r)==1) return r;
    }
}

int ShorFactor(int N){
    /// Assume N is not prime
    while(1){
        int a = rand()%(N-1)+1;
        int g = __gcd(a,N);
        if(g!=1) return g;
        int r = ShorPeriod(a,N);
        if(r&1) continue;
        if(modPow(a,r>>1,N)==N-1) continue;
        return __gcd(modPow(a,r>>1,N)+1,N);
    }
}

}


#endif // __QUANTUM_HPP__
