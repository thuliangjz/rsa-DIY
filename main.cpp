#include <iostream>
#include "RSA.h"
#include "BigInteger.h"

#include <random>

extern void testDiv();
extern void testMul();
extern void testEuclid();
extern void testRandom();
extern void testExponent();
extern void testPrime();
extern void testShr();
extern void testReverse();
extern void testRightShift();
extern void testExponentNewton();
int main() {
    testPrime();
    return 0;
}

