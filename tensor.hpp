#ifndef __TENSOR_HPP__
#define __TENSOR_HPP__

#include<bits/stdc++.h>
#include<Eigen/SVD>

using namespace std;

#define cpx complex<double>
#define Vec vector<cpx>
#define Mat vector<Vec>
#define size2D pair<int, int>
#define eps (1e-8)

class tensor{
    int m, n;
    Mat mat;
public:
    tensor(int m, int n): m(m), n(n), mat(m,Vec(n,0)){}
    tensor(Mat mat, int m=0, int n=0): mat(mat){
        if(!m && !n){
            tensor::m = mat.size();
            tensor::n = mat[0].size();
            return;
        }
        assert(m&&n);
        tensor::m = m;
        tensor::n = n;
        tensor::mat.resize(m);
        for(int i=0;i<m;++i)
            tensor::mat[i].resize(n,0);
    }
    inline size_t rows()const{ return m;}
    inline size_t cols()const{ return n;}
    inline size2D size()const{ return {m,n};}
    cpx& operator()(int i, int j){ return mat[i][j];}
    cpx const& operator()(int i, int j)const{ return mat[i][j];}

    friend tensor tprod(tensor const& l, tensor const& r){
        tensor out(l.rows()*r.rows(),l.cols()*r.cols());
        for(int i=0;i<out.rows();++i)
        for(int j=0;j<out.cols();++j)
            out(i,j) = l(i/r.rows(), j/r.cols()) * r(i%r.rows(), j%r.cols());
        return out;
    }

    friend tensor tmul(tensor const& l, tensor const& r){
        assert(l.cols()==r.rows());
        int knum = l.cols();
        tensor out(l.rows(), r.cols());
        for(int i=0;i<out.rows();++i)
        for(int j=0;j<out.cols();++j)
            for(int k=0;k<knum;++k)
                out(i,j) += l(i,k)*r(k,j);
        return out;
    }

    friend tensor tadd(tensor const& l, tensor const& r){
        assert(l.size()==r.size());
        tensor out(l);
        for(int i=0;i<out.rows();++i)
        for(int j=0;j<out.cols();++j)
            out(i,j) += r(i,j);
        return out;
    }

    vector<pair<tensor,tensor> > factor(size2D s)const{
        using namespace Eigen;
        assert(rows()%s.first==0 && cols()%s.second==0);
        int r = s.first*s.second, c = rows()*cols()/r;
        typedef Matrix<cpx,Dynamic,Dynamic> MatCpx;
        MatCpx A(r,c);
        for(int i=0, sz = rows()*cols();i<sz;++i)
            A(i/c,i%c) = mat[i/cols()][i%cols()];
        JacobiSVD<MatCpx> svd(A,ComputeThinU | ComputeThinV);
        auto sval = svd.singularValues();
        auto U = svd.matrixU();
        auto Va = svd.matrixV().conjugate();
        vector<pair<tensor, tensor> > out;
        for(int i=0;i<sval.size();++i){
            if(sval(i)<eps) break;
            tensor left(s.first,s.second);
            tensor right(rows()/s.first,cols()/s.second);
            double sqtSval = sqrt(sval(i));
            for(int j=0;j<U.rows();++j)
                left(j/left.cols(),j%left.cols()) = U(j,i)*sqtSval;
            for(int j=0;j<Va.rows();++j)
                right(j/right.cols(),j%right.cols()) = Va(j,i)*sqtSval;
            out.push_back({left,right});
        }
        return out;
    }

    friend ostream& operator<<(ostream& sout, tensor const& m){
        for(int i=0;i<m.rows();++i){
            for(int j=0;j<m.cols();++j)
                sout<<m(i,j)<<" ";
            sout<<"\n";
        }
        return sout;
    }

    tensor transpose()const{
        tensor out(cols(),rows());
        for(int i=0;i<rows();++i)
        for(int j=0;j<cols();++j)
            out(j,i) = mat[i][j];
        return out;
    }
};

#endif // __TENSOR_HPP__
