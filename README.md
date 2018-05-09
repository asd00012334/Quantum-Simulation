# Quantum-Simulation

## Description
Simulating Quantum Algorithms including Shor's Algorithm and Deutsch-Jozsa's Algorithms.

## Dependency
Use C++ [Eigen Library](http://eigen.tuxfamily.org/index.php?title=Main_Page)

## Guide
- First, set up the dependency (C++ Eigen).
- Put all files in the same directory, then include `Quantum.hpp`.
- You might refer to the `test.cpp` for sample code.

## Basic Usage
### Namespace
Our interface is within namespace `Quantum`. You might want to use
```
using namespace Quantum;
```

### Qreg
It represents a quantum register.
- To instantiate a Qreg object, use
```c++
Qreg obj(initialValue, bitNumber);
```
- To measure a Qreg from s'th bit to t'th bit, use
```c++
obj.measure(s,t);
```
- To slice a Qreg, use
```c++
auto stateList = obj.slice(s,t);
```
It represents a sum of product form,
\[(u1,v1),(u2,v2),...(un,vn)\]
where ui, vi are stored as pair within a list.
- Use `a*b` to perform Kronecker product.
### Gate
It represent unitary matrix acting on small number of bits.
- Use `gate(qubit)` to take act on a quantum register and return its output.
- Use `gate1(gate2)` to perform composite of gates.
- Use `gate1*gate2` to perform the Kronecker product of gates.

### Deutsch-Jozsa Algorithm
- Use `DeutschJozsa(f, n)` to perform Deutsch-Jozsa algorithm,
where `f` is of a function handler either constant or balanced
and `n` is the total bit number

### Shor's Algorithm
- Use `ShorFacto(N)` to return a trivial divisor.
- Use `ShorPeriod(a,N)` to return the order of $a \mod N$ while a, N are co-prime.
