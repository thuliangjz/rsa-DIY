#include <utility>

//
// Created by liangjz on 19-11-2.
//
#include <iostream>
#include "RSA.h"
#include <thread>
#include <string.h>
#include <fstream>
#include <iomanip>

using std::cout;
using std::endl;
void RSA::keyGen(int length) {
    if(length < MIN_KEY_LENGTH)
        throw RSAException(ERR_KEY_LENGTH_TOO_SHORT);
    p = BigInteger::getPrimeWithin(BigInteger::nBitMin(length / 2), BigInteger::nBitMax(length / 2), m_prime_worker_cnt, m_prime_iter_cnt);
    BigInteger n_max = BigInteger::nBitMax(length), n_min = BigInteger::nBitMin(length);
    BigInteger q_min = BigInteger::div(n_min, p).first, q_max = BigInteger::div(n_max, p).first;
    q = BigInteger::getPrimeWithin(q_min, q_max, m_prime_worker_cnt, m_prime_iter_cnt);
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
    cout << "n(modulo): 0x"; n.printHex(BigInteger::PRINT_MODE_COMPACT);
    cout << "p: 0x"; p.printHex(BigInteger::PRINT_MODE_COMPACT);
    cout << "q: 0x"; q.printHex(BigInteger::PRINT_MODE_COMPACT);
    cout << "e: 0x"; e.printHex(BigInteger::PRINT_MODE_COMPACT);
    cout << "d: 0x"; d.printHex(BigInteger::PRINT_MODE_COMPACT);
}

void RSA::setWorkerCntForMR(int worker_cnt) {
    int max_thread = std::thread::hardware_concurrency();
    max_thread = max_thread > 0 ? max_thread : 1;
    worker_cnt = worker_cnt <= 0 ? max_thread : worker_cnt;
    worker_cnt = worker_cnt > max_thread ? max_thread : worker_cnt;
    cout << "max thread supported: " << max_thread << ", work count set to: " << worker_cnt << endl;
    m_prime_worker_cnt = worker_cnt;
}

void RSA::setIterCntForMR(int iter_cnt) {
    iter_cnt = iter_cnt < 16 ? 16 : iter_cnt;
    iter_cnt = iter_cnt > 50 ? 50 : iter_cnt;
    cout << "iteration count set to: " << iter_cnt << endl;
    m_prime_iter_cnt = iter_cnt;
}

vector<char> RSA::encrypt(const vector<char> &message) {
    if(!m_pub_key_read)
        throw RSAException(ERR_KEY_UNKNOWN);
    int size_buffer = n.getConstVector().size();
    int size_msg = message.size();
    int key_length = n.getBitCnt(), padding_msg = key_length / 8 - 1;
    int padding_result = (key_length - 1) / 8 + 1;
    u64vec buffer(size_buffer, 0);
    vector<char> result(padding_result * ((size_msg - 1) / padding_msg + 1));
    for(int i = 0, j = 0; i < size_msg; i += padding_msg, ++j){
        int copy_len = size_msg - i >= padding_msg ? padding_msg : size_msg - i;
        memcpy(&buffer[0], &message[0] + i, copy_len);
        BigInteger num(buffer);
        num = BigInteger::add(num, BigInteger::two);
        num = n.fastExponentNewton(num, e);     //使用公钥进行加密
        u64vec buffer_result = num.getBits();
        buffer_result.resize(size_buffer, 0);   //剩下的位必须用0补全
        memcpy(&result[padding_result * j], &buffer_result[0], padding_result);
        std::fill(buffer.begin(), buffer.end(), 0);
    }
    return result;
}

vector<char> RSA::decrypt(const vector<char> &message) {
    if(!m_private_key_read)
        throw RSAException(ERR_KEY_UNKNOWN);
    int size_buffer = n.getConstVector().size();
    int size_encrypted = message.size();
    int key_length = n.getBitCnt();
    int padding_encrypted = (key_length - 1) / 8 + 1;
    int padding_decrypted = key_length / 8 - 1; //相当于encrypt中的padding_msg
    if(size_encrypted % padding_encrypted != 0)
        throw RSAException(ERR_BAD_PADDING);
    u64vec buffer(size_buffer, 0);
    vector<char> decrypted(padding_decrypted * size_encrypted / padding_encrypted);
    for(int i = 0, j = 0; i < size_encrypted; i += padding_encrypted, ++j){
        memcpy(&buffer[0], &message[0] + i, padding_encrypted);
        BigInteger num(buffer);
        num = n.fastExponentNewton(num, d);     //解密用的是私钥
        num = BigInteger::sub(num, BigInteger::two);
        u64vec buffer_decrypted = num.getBits();
        buffer_decrypted.resize(size_buffer, 0);
        memcpy(&decrypted[j * padding_decrypted], &buffer_decrypted[0], padding_decrypted);
        std::fill(buffer.begin(), buffer.end(), 0);
    }
    return decrypted;
}

void RSA::writeBigIntegerToFile(std::ofstream& f,const BigInteger& big_int){
    const u64vec& vec = big_int.getConstVector();
    f << std::hex << vec.size() << " ";
    for(auto &e:vec){
        f << std::setw(16) << std::setfill('0') << std::hex << e << " ";
    }
    f<<endl;
}

BigInteger RSA::readBigIntegerFromFile(std::ifstream &file_in) {
    int size_vec;
    file_in >> std::hex >> size_vec;
    u64vec vec;
    for(int i = 0; i < size_vec; ++i){
        uint64_t tmp;
        file_in >> std::hex >> tmp;
        if(file_in.eof())
            throw RSAException(ERR_FILE_PARSE_FAILED);
        vec.push_back(tmp);
    }
    return BigInteger(vec);
}

void RSA::readPubKey(const std::string &file) {
    std::ifstream file_in(file);
    if(!file_in)
        throw RSAException(ERR_INPUT_FILE_NOT_EXIST);
    n = readBigIntegerFromFile(file_in);
    e = readBigIntegerFromFile(file_in);
    m_pub_key_read = true;
}

void RSA::readPrivateKey(const std::string &file) {
    std::ifstream file_in(file);
    if(!file_in)
        throw RSAException(ERR_INPUT_FILE_NOT_EXIST);
    n = readBigIntegerFromFile(file_in);
    d = readBigIntegerFromFile(file_in);
    p = readBigIntegerFromFile(file_in);
    q = readBigIntegerFromFile(file_in);
    m_private_key_read = true;
}

void RSA::saveKey(const std::string& file_public, const std::string& file_private) {
    if(!m_is_key_generated)
        throw RSAException(ERR_KEY_UNKNOWN);
    std::ofstream f_pub(file_public, std::ios::out), f_prvt(file_private, std::ios::out);
    writeBigIntegerToFile(f_pub, n);
    writeBigIntegerToFile(f_pub, e);

    writeBigIntegerToFile(f_prvt, n);
    writeBigIntegerToFile(f_prvt, d);
    writeBigIntegerToFile(f_prvt, p);
    writeBigIntegerToFile(f_prvt, q);
}