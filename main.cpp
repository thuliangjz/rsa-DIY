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
int main() {
    testPrime();
    return 0;
}

