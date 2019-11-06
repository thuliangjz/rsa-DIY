#include <utility>

//
// Created by liangjz on 19-11-2.
//
#include <iostream>
#include "RSA.h"


void RSA::keyGen(int length) {
    if(length < MIN_KEY_LENGTH)
        throw RSAException(ERR_KEY_LENGTH_TOO_SHORT);

}
