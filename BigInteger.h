//
// Created by liangjz on 19-11-2.
//

#ifndef RSA_BIGINTEGER_H
#define RSA_BIGINTEGER_H

#include <vector>
#include <cstdint>
#include <tuple>
#include <random>
using std::vector;
using std::pair;
using std::tuple;
class BigInteger;

typedef  vector<uint64_t> u64vec;
typedef pair<BigInteger, BigInteger> BigIntegerPair;
typedef tuple<BigInteger, BigInteger, BigInteger> BigIntegerTriplet;
class BigInteger {
    //vec[0]的LSB表示个位，以此类推
    public:
    class BigIntegerException: public std::exception{
    public:
        explicit BigIntegerException(int code):code(code){}
        int getCode()const{return code;}
    private:
        int code;
    };
    public:
        static BigInteger nBitMin(int n);
        static BigInteger nBitMax(int n);
        static BigInteger add(const BigInteger& a, const BigInteger& b);
        static BigInteger sub(const BigInteger& a, const BigInteger& b);    //compute a - b
        static BigInteger mul(const BigInteger& a, const BigInteger& b);
        inline static int shrink(vector<uint64_t>& vec);   //将末尾0移除向量,返回表示大数的位数
        static pair<BigInteger, BigInteger> div(const BigInteger& a, const BigInteger& b); //compute a/b, return(q, r)
        static BigInteger euclidean(const BigInteger& a, const BigInteger& b);
        static BigIntegerPair partialExtendedEuclidean(const BigInteger& a, const BigInteger& n); //return (gcd, v) where av \equiv gcd \pmod n
        static BigInteger getPrimeWithin(const BigInteger& min, const BigInteger& max);
        static BigInteger randomWithin(const BigInteger& min,const BigInteger& max, std::random_device& generator); //thread safe so long as generator are thread-local
        static BigInteger fastExponent(const BigInteger& a, const BigInteger& e, const BigInteger& n); //compute a^e (mode n)
        static int compare(const BigInteger& left,const BigInteger& right); //return -1 if left < right, 1 if left > right, else 0
        static const BigInteger zero;
        static const BigInteger one;
        static const BigInteger two;
        int getBitCnt()const{return m_cnt_bits;}
        vector<uint64_t> getBits()const{return m_vec_bits;}
        explicit BigInteger(const vector<uint64_t>&   vec);
        void printHex(int mode = PRINT_MODE_SPACED);
        inline bool isEven()const;
    public:
        static const int ERR_EMPTY_VEC = 0;
        static const int ERR_DIV_ZERO = 1;      //zero divisor
        static const int ERR_MAX_SMALLER_THAN_RANGE = 2;
        static const int PRINT_MODE_SPACED = 0;
        static const int PRINT_MODE_COMPACT = 1;
    private:
        vector<uint64_t> m_vec_bits;
        int m_cnt_bits;
};


//entry_max = (bit_count - 1) / 64;
//local_max = (bit_count - 1) % 64;
#endif //RSA_BIGINTEGER_H
