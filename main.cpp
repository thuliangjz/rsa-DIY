#include <iostream>
#include "RSA.h"
#include "BigInteger.h"
#include "BigIntegerTest.h"


int main() {
    RSA rsa;
    auto start = std::chrono::system_clock::now();
    rsa.keyGen(768);
    auto end = std::chrono::system_clock::now();
    rsa.printKeyInfo();
    std::cout << "time used:" << std::dec << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;
    return 0;
}

