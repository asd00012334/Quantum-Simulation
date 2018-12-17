#include<bits/stdc++.h>

#define __APPROX__ (1e-3)

#include <Eigen/QR>
#include <Eigen/Eigenvalues>
#include <complex>
#include <cmath>
#include "Quantum.hpp"

using namespace std;
using namespace Quantum;
using namespace Eigen;

Eigen::MatrixXf thinQ;

Qreg UnitaryOP(Qreg const& r, Qreg const& c) {
	int n = c.bitNum(), m = r.bitNum();
	int const unitary_size = 1 << m;
	
	// Convert type from MatrixXf to float array.
	Mat u = Mat(unitary_size, Vec(unitary_size));
	for (int i = 0; i < unitary_size; i++) {
		for (int j = 0; j < unitary_size; j++) {
			u[i][j] = thinQ(i, j);
		}
	}

	// Build unitary quantum gate.
	Qreg out(n + m);
	Gate U(u);
	for(int i=0;i<out.stateNum();++i) assert(abs(out[i])<1e-9);
	auto x = r;

	for(int i = 0; i < c.stateNum(); ++i){
		for(int j = 0; j < x.stateNum(); ++j){
			out[(j<<n)|i] += c[i] * x[j];
		}
		if (i < c.stateNum()-1) x = U(x);
	}

	return out;
}

vector<float> PhaseEstimation(int const& bitNum) {
	vector<float> retval;
	int n = 10, m = bitNum;

	Qreg init = Hadamard(Qreg(0, n));
	Qreg func = UnitaryOP(Qreg(0, m), init);
	Qreg result = func.applySlice(0, n, QFT);
	for(auto e: result.sample(m,n+m,4*m))
		retval.push_back((float)e/(1<<n));
	//retval = ((double) result.sample(m, n+m,4*m)) / (1<<n);

	return retval;
}

void CheckCorrectness(vector<float>const& theta, int unitary_size) {
	complex<float> complex_i;
	complex<float> lambda;
	float pi;

	pi = 2 * asin(1);
	complex_i = complex<float>(0,1);
	for(int i=0;i<theta.size();++i){
		lambda = exp(2*pi*complex_i*theta[i]);
		printf("phase: %lf\nestimate eigenvalue: ", theta[i]);
		cout << lambda << "\n\n";
	}
	

	EigenSolver<MatrixXf> es(thinQ);
	cout << "Correct eigenvalue: \n" << es.eigenvalues() << '\n';
	
	auto x = es.eigenvalues();
	cout << "\nCorrect phase: \n";
	for(int i = 0; i < unitary_size; ++i) {
		float tmp_phase = (log(x[i])/complex<float>(0,2*pi)).real();
		if (tmp_phase < 0) tmp_phase += 1;
		cout << tmp_phase << '\n';
	}
}

int main(){
    srand(time(NULL)); /** Must do this for stochastic algorithm **/

	int t0, unitary_bits = 3;
	int unitary_size = (1<<unitary_bits);

	// Create unitary matrix.
	MatrixXf A(MatrixXf::Random(unitary_size, unitary_size));
	HouseholderQR<MatrixXf> qr(A);
	thinQ = qr.householderQ() * MatrixXf::Identity(unitary_size, unitary_size);
	//cout << thinQ.adjoint()*thinQ << endl<<endl;
	//cout << thinQ*thinQ.adjoint() << endl;
				
	printf("***** Phase estimation *****\nstatus: Start.\n");
	t0 = clock();
		
	auto&& theta = PhaseEstimation(unitary_bits);
	
	printf("status: Complete.\n\n");
	printf("Time taken: %ld ms.\n\n", (clock()-t0)/(CLOCKS_PER_SEC/1000) );	
	CheckCorrectness(theta, unitary_size);
}

