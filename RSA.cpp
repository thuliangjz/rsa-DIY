#include <utility>

//
// Created by liangjz on 19-11-2.
//
#include <iostream>
#include "RSA.h"
using std::cout;
using std::endl;
void RSA::keyGen(int length) {
    if(length < MIN_KEY_LENGTH)
        throw RSAException(ERR_KEY_LENGTH_TOO_SHORT);
    p = BigInteger::getPrimeWithin(BigInteger::nBitMin(length / 2), BigInteger::nBitMax(length / 2));
    BigInteger n_max = BigInteger::nBitMax(length), n_min = BigInteger::nBitMin(length);
    BigInteger q_min = BigInteger::div(n_min, p).first, q_max = BigInteger::div(n_max, p).first;
    q = BigInteger::getPrimeWithin(q_min, q_max);
    n = BigInteger::mul(p, q);
    BigInteger q_1 = BigInteger::sub(q, BigInteger::one),
            p_1 = BigInteger::sub(p, BigInteger::one),
            lambda = BigInteger::mul(q_1, p_1);
    std::random_device generator;
    BigInteger gcd_e_n;
    do{
        e = BigInteger::randomWithin(BigInteger::nBitMin(32), BigInteger::nBitMax(32), generator);
        e = BigInteger::add(e, BigInteger::one);
        BigIntegerPair tmp = BigInteger::partialExtendedEuclidean(e, lambda);
        gcd_e_n = tmp.first; d = tmp.second;    //e比较小将使得d的计算也相对快速一些
    }while(BigInteger::compare(gcd_e_n, BigInteger::one) != 0);
    m_is_key_generated = true;
}

void RSA::printKeyInfo() {
    cout << "key length: " << n.getBitCnt() <<endl;
    cout << "n(modulo): "; n.printHex(BigInteger::PRINT_MODE_COMPACT);
    cout << "p: "; p.printHex(BigInteger::PRINT_MODE_COMPACT);
    cout << "q: "; q.printHex(BigInteger::PRINT_MODE_COMPACT);
    cout << "e: "; e.printHex(BigInteger::PRINT_MODE_COMPACT);
    cout << "d: "; d.printHex(BigInteger::PRINT_MODE_COMPACT);
}
