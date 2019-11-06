//
// Created by liangjz on 19-11-2.
//

#ifndef RSA_RSA_H
#define RSA_RSA_H

#include<vector>
using std::vector;
class RSA {
public:
    class RSAException: public std::exception{
    public:
        explicit RSAException(int code):m_code(code){}
        int getCode(){return m_code;}
    private:
        int m_code;
    };
    static const int ERR_KEY_LENGTH_TOO_SHORT = 0;
    static const int MIN_KEY_LENGTH = 256;
public:
    void keyGen(int length);
};


#endif //RSA_RSA_H
