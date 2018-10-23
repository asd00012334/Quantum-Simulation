#include<bits/stdc++.h>

#include "Quantum.hpp"

using namespace std;
using namespace Quantum;

vector<uint32_t> SimonAlgo(function<int(int)>const& f, int const bitNum) {

	Qreg init = Hadamard(Qreg(0, bitNum));
	Qreg func = oracle(init, Qreg(0, bitNum), f);
	func.measure(0, bitNum);
	Qreg result = func.applySlice(bitNum, bitNum<<1, Hadamard);

    vector<uint32_t>retval;
    for(auto e:result.sample(bitNum,bitNum<<1,bitNum<<2) )
        retval.push_back(e);
	//for (int i = 0; i < bitNum<<2; i++) {
	//	auto tmp = result;
	//	retval.push_back(tmp.measure(bitNum, bitNum<<1));
	//}

	return retval;
}

uint32_t GetElement(vector<uint32_t>& eq, int const row, int const col) {
	return eq[row] & (1<<col);
}

vector<uint32_t> GaussianEllimination(vector<uint32_t>& eq, int const bitNum) {
	int rows = eq.size();
	int cols = bitNum;

	assert(rows >= cols);
	for (int i = 0; i < cols; i++) {
		if (GetElement(eq, i, bitNum-i-1) == 0) {
			for (int j = i + 1; j < rows; j++) {
				if (GetElement(eq, j, bitNum-i-1) != 0) {
					swap(eq[i], eq[j]);
					break;
				}
			}
		}

		if (GetElement(eq, i, bitNum-i-1) == 0) continue;

		for (int j = 0; j < rows; j++) {
			if (GetElement(eq, j, bitNum-i-1) && (i != j))
				eq[j] ^= eq[i];
		}
	}

	return eq;
}

int Count1(uint32_t const x, int const bitNum) {
	int retval = 0;

	for (int i = 0; i < bitNum; i++)
		if (x & (1<<i)) retval++;

	return retval;
}

int SolLinearEQ(vector<uint32_t>& y, int const bitNum) {
	uint32_t retval = (1<<bitNum) - 1;

	GaussianEllimination(y, bitNum);

	for (vector<uint32_t>::iterator it = y.begin(); it != y.end(); it++) {
		printf("%d\n", *it);
		if (Count1(*it, bitNum) & 1) {
			for (int i = 0; i < bitNum; i++) {
				retval &= ~ ((*it) & (1<<i));
			}
		}
	}

	return retval;
}

int s = 0;

int OracleFunc(int i) {
	// return i & 0xFFFFFFFE;
	return (i^s)>i? i : (i^s);
}

int main(){
    srand(time(NULL)); /** Must do this for stochastic algorithm **/

	int result_s;

	while (true) {
    	int t0;
		int qubits_num, samples;

		printf("Enter qubit's numbers: ");
		scanf("%d", &qubits_num);
		assert(qubits_num <= 32);
		samples = qubits_num<<2;
		printf("Enter s: ");
		scanf("%d", &s);

		printf("***** Simon's  algorithm *****\nstatus: Start.\n");
		t0 = clock();

		vector<uint32_t> y = SimonAlgo(OracleFunc, qubits_num);
		result_s = SolLinearEQ(y, qubits_num);

		printf("status: Complete.\n\ns = %d\n\n", result_s);
		printf("Time taken: %ld ms.\n\n", (clock()-t0)/(CLOCKS_PER_SEC/1000) );
	}
}

