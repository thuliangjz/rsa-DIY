#include <iostream>
#include "RSA.h"
#include "BigInteger.h"

#include <random>

extern void testDiv();
extern void testMul();
extern void testEuclid();
extern void testRandom();
extern void testExponent();
int main() {
    testExponent();
    return 0;
}

