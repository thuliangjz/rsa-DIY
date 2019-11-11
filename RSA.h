//
// Created by liangjz on 19-11-2.
//

#ifndef RSA_RSA_H
#define RSA_RSA_H

#include<vector>
#include "BigInteger.h"
using std::vector;
class RSA {
public:
    class RSAException: public std::exception{
    public:
        explicit RSAException(int code):m_code(code){}
        int getCode(){return m_code;}
        std::string getMessage(){
            switch (m_code){
                case ERR_KEY_LENGTH_TOO_SHORT:
                    return "key length too short";
                case ERR_KEY_UNKNOWN:
                    return "key not known for encrypt/decrypt";
                case ERR_BAD_PADDING:
                    return "bad padding for decrypted message";
                case ERR_INPUT_FILE_NOT_EXIST:
                    return "input file not exist";
                case ERR_FILE_PARSE_FAILED:
                    return "failed to parse key file";
                default:
                    return "";
            }
        }
    private:
        int m_code;
    };
    static const int ERR_KEY_LENGTH_TOO_SHORT = 0;
    static const int ERR_KEY_UNKNOWN = 1;
    static const int ERR_BAD_PADDING = 2;
    static const int ERR_INPUT_FILE_NOT_EXIST = 3;
    static const int ERR_FILE_PARSE_FAILED = 4;
    static const int MIN_KEY_LENGTH = 256;
public:
    void keyGen(int length);
    vector<char> encrypt(const vector<char>& message);
    vector<char> decrypt(const vector<char>& message);
    RSA():m_is_key_generated(false), m_prime_worker_cnt(6), m_prime_iter_cnt(32),
        m_pub_key_read(false), m_private_key_read(false){}
    void printKeyInfo();
    void setWorkerCntForMR(int worker_cnt);
    void setIterCntForMR(int iter_cnt);
    void saveKey(const std::string& file_public,const std::string& file_private);
    void readPubKey(const std::string& file);
    void readPrivateKey(const std::string& file);
private:
    static void writeBigIntegerToFile(std::ofstream& file_out, const BigInteger& big_int);
    static BigInteger readBigIntegerFromFile(std::ifstream& file_in);
private:
    bool m_is_key_generated, m_private_key_read, m_pub_key_read;
    BigInteger p, q, n, e, d;
    int m_prime_worker_cnt, m_prime_iter_cnt;
};


#endif //RSA_RSA_H
