#include <utility>

#include <utility>
#include <iostream>
#include <iomanip>
#include <random>
//
// Created by liangjz on 19-11-2.
//

#include "BigInteger.h"

#define ASM_LOAD "lahf; shr $8, %%ax; mov %%al, %%r15b;"
#define ASM_RESTORE "mov %%r15b, %%al; sal $8, %%ax; sahf;"

using std::cout;
using std::endl;

const BigInteger BigInteger::zero(u64vec(1, 0));
const BigInteger BigInteger::one(u64vec(1, 1));
const BigInteger BigInteger::two(u64vec(1, 2));

inline int getBitAt(const vector<uint64_t>& v, int m){
    return (v[m/64] >> static_cast<uint64_t>(m%64)) & UINT64_C(1);
}
//借助汇编实现逻辑右移
inline uint64_t __asm_shr(uint64_t a, int s){
    uint64_t r;
    asm volatile (
            "mov %1, %0;"
            "shr %%cl, %0;"
            :"=r"(r)
            :"r"(a), "c"(s)
            :
            );
    return s == 64 ? 0 : r;
}


int BigInteger::shrink(vector<uint64_t> &vec) {
    int cnt_bits = 64 * vec.size();
    while(cnt_bits > 0) {
        uint64_t entry = vec[(cnt_bits - 1) / 64];
        uint64_t mask = static_cast<uint64_t>(0x1) << static_cast<uint64_t>((cnt_bits - 1) % 64);
        if((entry & mask) != 0)
            break;
        --cnt_bits;
        if(cnt_bits % 64 == 0 && cnt_bits != 0)
            vec.pop_back();
    }
    cnt_bits = cnt_bits == 0 ? 1 : cnt_bits;
    return cnt_bits;
}

BigInteger::BigInteger(const vector<uint64_t>& vec):m_vec_bits(vec) {
    if(m_vec_bits.empty())
        throw BigIntegerException(ERR_EMPTY_VEC);
    m_cnt_bits = shrink(m_vec_bits);
}

inline bool BigInteger::isEven() const {
    return getBitAt(m_vec_bits, 0) == 0;
}

BigInteger BigInteger::nBitMax(int n) {
    vector<uint64_t> bits;
    if(n <= 0)
        return BigInteger(bits);
    int cnt_entry = (n-1) / 64 + 1;
    for(int i = 0; i < cnt_entry; ++i){
        bits.push_back(0);
    }
    for(int i = 0; i < cnt_entry; ++i){
        for(int j = 0; 64 * i + j < n; ++j){
            bits[i] = bits[i] | (UINT64_C(1) << (uint8_t)j);
        }
    }
    return BigInteger(bits);
}

BigInteger BigInteger::nBitMin(int n) {
    vector<uint64_t> bits;
    if(n <= 0)
        return BigInteger(bits);
    int cnt_entry = (n-1) / 64 + 1;
    for(int i = 0; i < cnt_entry; ++i){
        bits.push_back(0);
    }
    bits.back() |= UINT64_C(1) << (uint8_t)((n - 1) % 64);
    return BigInteger(bits);
}

void BigInteger::printHex(int mode) {
    std::string sep = mode == PRINT_MODE_SPACED ? " " : "";
    for(int i = m_vec_bits.size() - 1; i >= 0; i--){
        cout << std::setw(16) << std::setfill('0') << std::hex << m_vec_bits[i] << sep;
    }
    cout<<endl;
}

int BigInteger::compare(const BigInteger& left,const BigInteger& right) {
    if(left.getBitCnt() < right.getBitCnt())
        return -1;
    if(left.getBitCnt() > right.getBitCnt())
        return 1;
    int entry_cnt = (left.getBitCnt() - 1) / 64;
    vector<uint64_t> v1 = left.getBits(), v2 = right.getBits();
    for(int i = entry_cnt; i >= 0; --i){
        if(v1[i] > v2[i])
            return 1;
        if(v1[i] < v2[i])
            return -1;
    }
    return 0;
}

inline void __asm_add(vector<uint64_t>& a, vector<uint64_t>& b, vector<uint64_t>& c){
    //默认a, b长度相等, c.size() = a.size() + 1
    uint64_t *pt_a = &a[0],
        *pt_b = &b[0],
        *pt_c = &c[0];
    uint64_t len = a.size();
        asm volatile (
            "movq %0, %%r8;"
            "movq %1, %%r9;"
            "movq %2, %%r10;"
            "movzx %3, %%rdx;"
            "movq $0, %%rcx;"
            "clc;"
            "lahf;"
            "loop_add:"
                "sahf;" //从寄存器中读取flag
                "movq (%%r8, %%rcx, 8), %%rax;"
                "movq (%%r9, %%rcx, 8), %%rbx;"
                "adc %%rax, %%rbx;"
                "lahf;" //保存flag,后面的inc可能会破坏carry
                "movq %%rbx, (%%r10, %%rcx, 8);"
                "inc %%rcx;"
                "cmp %%rcx, %%rdx;"
                "jne loop_add;"
                "movq $0, (%%r10, %%rcx, 8);"
                "sahf;"
                "adc $0, (%%r10, %%rcx, 8);"
            :
            :"m"(pt_a), "m"(pt_b), "m"(pt_c), "m"(len)
            :"memory", "cc", "r8", "r9", "r10", "rdx", "rcx", "rax", "rbx");
}

inline void __asm_sub_from(vector<uint64_t>& a, vector<uint64_t> &b){
    //默认a和b一样长,将a-b的结果保存在a中
    uint64_t *pt_a = &a[0], *pt_b = &b[0], length = a.size();
    asm volatile (
        "movq %0, %%r8;"    //减数放在r8
        "movq %1, %%r9;"    //被减数r9
        "mov %2, %%rdx;"    //保存最大计数
        "movq $0, %%rcx;"
        "clc;"
        "lahf;"
        "1:"
            "sahf;" //从寄存器中读取flag
            "movq (%%r8, %%rcx, 8), %%rax;"
            "sbb %%rax, (%%r9, %%rcx, 8);"
            "lahf;" //保存flag,后面的inc可能会破坏carry
            "inc %%rcx;"
            "cmp %%rcx, %%rdx;"
            "jne 1b;"
        :
        :"m"(pt_b), "m"(pt_a), "m"(length)
        :"memory", "cc", "r8", "r9", "rdx", "rcx", "rax", "rbx");
}

BigInteger BigInteger::add(const BigInteger &a, const BigInteger& b) {
    vector<uint64_t> a_bits = a.getBits(), b_bits = b.getBits();
    int max_len = std::max(a_bits.size(), b_bits.size());
    a_bits.resize(max_len, 0);
    b_bits.resize(max_len, 0);
    vector<uint64_t> result(max_len + 1, 0);
    __asm_add(a_bits, b_bits, result);
    return BigInteger(result);
}

BigInteger BigInteger::sub(const BigInteger &a, const BigInteger &b) {
 vector<uint64_t> a_bits = a.getBits(), b_bits = b.getBits();
    int max_len = std::max(a_bits.size(), b_bits.size());
    a_bits.resize(max_len, 0);
    b_bits.resize(max_len, 0);
    vector<uint64_t> result(max_len, 0);
    uint64_t *pt_a_bits = &a_bits[0],
        *pt_b_bits = &b_bits[0],
        *pt_result = &result[0];
    asm volatile (
            "movq %0, %%r8;"    //减数放在r8
            "movq %1, %%r9;"    //被减数r9
            "movq %2, %%r10;"
            "movzx %3, %%rdx;"
            "movq $0, %%rcx;"
            "clc;"
            "lahf;"
            "1:"
                "sahf;" //从寄存器中读取flag
                "movq (%%r8, %%rcx, 8), %%rax;"
                "movq (%%r9, %%rcx, 8), %%rbx;"
                "sbb %%rax, %%rbx;"
                "lahf;" //保存flag,后面的inc可能会破坏carry
                "movq %%rbx, (%%r10, %%rcx, 8);"
                "inc %%rcx;"
                "cmp %%rcx, %%rdx;"
                "jne 1b;"
            :
            :"m"(pt_b_bits), "m"(pt_a_bits), "m"(pt_result), "m"(max_len)
            :"memory", "cc", "r8", "r9", "r10", "rdx", "rcx", "rax", "rbx");
    return BigInteger(result);
}

BigInteger BigInteger::mul(const BigInteger &a, const BigInteger &b) {
    vector<uint64_t> a_bits = a.getBits(), b_bits = b.getBits();
    vector<uint64_t> result(a_bits.size() + b_bits.size(), 0);
    vector<uint64_t> tmp(a_bits.size() + b_bits.size(), 0); //tmp仅仅保存每次乘法的中间值
    uint64_t *pt_a = &a_bits[0], *pt_b = &b_bits[0], *pt_result = &result[0], *pt_tmp = &tmp[0];
    uint64_t length_a = a_bits.size(), length_tmp = tmp.size();
    for(int i = 0; i < b_bits.size(); ++i){
        uint64_t* pt_tmp_i = pt_tmp + i;
        asm volatile (
                "movq $0, %%rdx;"
                "movq $0, %%r14;"   //r14保存上一次乘法的溢出结果
                "movq $0, %%rcx;"   //rcx用于计数
                "movq %1, %%rsi;"   //rsi保存大数的首地址
                "movq %3, %%rdi;"   //rdi保存中间结果（就是tmp）的首地址
                "movq %0, %%r8;"    //r8保存小乘数,维持不变(寄存器访问速度更快)
                "clc;"
                ASM_LOAD    //carry保存在r15中
                "1:"
                    "mov %%r8, %%rax;"
                    "mov (%%rsi, %%rcx, 8), %%rbx;"
                    "mul %%rbx;"    //乘法完成,结果保留在rdx:rax中
                    "mov %%rax, %%rbx;"
                    ASM_RESTORE
                    "adc %%r14, %%rbx;" //连带carry外加上一次的乘法溢出一起加到本次结果上
                    ASM_LOAD
                    "mov %%rdx, %%r14;"  //保存本次溢出
                    "mov %%rbx, (%%rdi, %%rcx, 8);" //本次计算结果保存到tmp数组
                    "inc %%rcx;"
                    "cmp %%rcx, %2;"    //是否已经执行了length_a次
                    "jne 1b;"
                    ASM_RESTORE
                    "mov $0, %%rax;"
                    "adc %%r14, %%rax;"
                    "mov %%rax, (%%rdi, %%rcx, 8);"
                :
                :"m"(pt_b[i]), "m"(pt_a), "m"(length_a), "m"(pt_tmp_i)
                :"memory", "cc", "rax", "rbx", "rcx", "rdx", "rsi", "rdi", "r8", "r14", "r15"
                );
        //tmp更后面的结果赋值为0,避免每次对tmp重新赋值
        std::fill(tmp.begin(), tmp.begin() + i, 0);
        //这段汇编完成的是将tmp的值加到result上面,注意两个向量长度相等
        asm volatile (
            "movq %0, %%r8;"
            "movq %1, %%r9;"
            "mov %2, %%rdx;"
            "movq $0, %%rcx;"
            "clc;"
            "lahf;"
            "1:"
                "sahf;"
                "movq (%%r8, %%rcx, 8), %%rax;"
                "adc %%rax, (%%r9, %%rcx, 8);"  //将r8指向的数加到r9上
                "lahf;"
                "inc %%rcx;"
                "cmp %%rcx, %%rdx;"
                "jne 1b;"
            :
            :"m"(pt_tmp), "m"(pt_result), "m"(length_tmp)
            :"memory", "cc", "r8", "r9", "rdx", "rcx", "rax", "rbx");
    }
    return BigInteger(result);
}

pair<BigInteger, BigInteger> BigInteger::div(const BigInteger &a, const BigInteger &b) {
    vector<uint64_t> a_bits = a.getBits(), b_bits = b.getBits();
    if(b.getBitCnt() == 0 && b_bits[0] == UINT64_C(0))
        throw BigIntegerException(ERR_DIV_ZERO);
    int l1 = a.getBitCnt(), l2 = b.getBitCnt();
    vector<uint64_t> q(std::max((l1 - l2) / 64, 0) + 1, 0), r(l2, 0);
    vector<uint64_t> tmp(a_bits.size(), 0); //保存中间移位得到的结果
    while(l1 >= l2){
        int m = l2 - 1;
        int pos_q;  //首先确定本轮除法商的位置
        while(m >= 0 && getBitAt(b_bits, m) == getBitAt(a_bits, l1 - l2 + m))
            --m;
        if(m < 0 || getBitAt(b_bits, m) < getBitAt(a_bits, l1 - l2 + m))
            pos_q = l1 - l2;
        else if(l1 > l2)    //此时已经有m >= 0 && getBitAt(b_bits, m) > getBitAt(a_bits, l1 - l2 + m)
            pos_q = l1 - l2 - 1;
        else
            break;
        int displace_entry = pos_q / 64, displace_bit = pos_q % 64;
        q[displace_entry] |= static_cast<uint64_t >(1) << static_cast<uint64_t >(displace_bit); //设置商
        //计算移位,这里合并了两个循环可以加快减少数据依赖
        tmp[displace_entry] = b_bits[0] << static_cast<uint64_t >(displace_bit);
        for(int i = 1; i < b_bits.size(); ++i){
            tmp[i + displace_entry] = b_bits[i] << static_cast<uint64_t >(displace_bit);
            tmp[i + displace_entry] |= __asm_shr(b_bits[i - 1], 64 - displace_bit);
        }
        tmp[b_bits.size() + displace_entry] |= __asm_shr(b_bits.back(), 64 - displace_bit);
        //将tmp从a_bits中减去
        __asm_sub_from(a_bits, tmp);
        //更新l1的值,缩小a_bits
        l1 = shrink(a_bits);
        tmp.resize(a_bits.size());
        std::fill(tmp.begin(), tmp.end(), 0);
    }
    r = a_bits;
    return std::make_pair(BigInteger(q), BigInteger(r));
}

BigInteger BigInteger::euclidean(const BigInteger &a, const BigInteger &b) {
    BigInteger small(a), big(b);
    if(compare(a, b) == 1){
        small = b;
        big = a;
    }
    while(compare(small, zero) == 1){
        BigIntegerPair qr = div(big, small);
        big = small; small = qr.second;
    }
    return big;
}

BigIntegerPair BigInteger::partialExtendedEuclidean(const BigInteger &a, const BigInteger &n) {
    BigInteger small = div(a, n).second, big = n, n_square = mul(n, n);
    BigInteger v(u64vec(1, 0)), f(u64vec(1, 1));
    while(compare(small, zero) == 1){
        BigIntegerPair qr = div(big, small);
        big = small; small = qr.second;
        BigInteger tmp = f;
        f = div(sub(add(v, n_square), mul(qr.first, f)), n).second; //f = v + n^2 - qf \pmod n
        v = tmp;
    }
    return std::make_pair(big, v);
}
BigInteger BigInteger::randomWithin(const BigInteger& min, const BigInteger &max, std::random_device &generator) {
    int comp_result = compare(min, max);
    if(comp_result == 1)
        throw BigIntegerException(ERR_MAX_SMALLER_THAN_RANGE);
    if(comp_result == 0)
        return BigInteger(min);
    u64vec range_bits = sub(max, min).getBits();
    u64vec result_bits(range_bits.size(), 0);
    uint64_t mask = 0xffff;
    bool is_smaller = false;
    std::uniform_int_distribution<int> distribution_default(0, 0xffff);
    for(int i = range_bits.size() - 1; i >= 0; --i){
        for(int j = 3; j >= 0; --j){
            int tmp_range = (range_bits[i] >> static_cast<uint64_t >(j * 16)) & mask;
            uint64_t tmp_result = 0;
            if(!is_smaller && tmp_range != 0){
                std::uniform_int_distribution<int> distribution_tmp(0, tmp_range);
                tmp_result = distribution_tmp(generator);
                if(tmp_result == tmp_range)
                    is_smaller = true;
            }
            else if(!is_smaller)
                continue;
            else
                tmp_result = distribution_default(generator);
            result_bits[i] |= tmp_result << static_cast<uint64_t >(j * 16);
        }
    }
    return add(BigInteger(result_bits), min);
}

BigInteger BigInteger::fastExponent(const BigInteger &a, const BigInteger &e, const BigInteger &n) {
    u64vec e_bits = e.getBits();
    int e_len = e.getBitCnt();
    BigInteger basis = a, result(u64vec(1, 1));
    for(int i = 0; i < e_len; ++i){
        int bit = getBitAt(e_bits, i);
        if(bit){
            result = div(mul(result, basis), n).second;
        }
        basis = div(mul(basis, basis), n).second;
    }
    return result;
}

BigInteger BigInteger::getPrimeWithin(const BigInteger &min, const BigInteger &max) {
    std::random_device generator;
    BigInteger n(u64vec(1, 0));
    vector<BigInteger> n_sampled;
    int tested = 0;
    while(true){
        while (true) {
            n = randomWithin(min, max, generator);
            if (n.isEven())
                n = add(n, one);
            if (compare(n, max) == 1)
                continue;
            bool sampled = false;
            for(auto &b : n_sampled){
                if(compare(b, n) == 0){
                    sampled = true;
                    break;
                }
            }
            if(!sampled)
                break;
        }
        int s = 0;

//        u64vec tmp;
//        tmp.push_back(0x4b0ceacd9b81bb47);
//        tmp.push_back(0xb);
//        n = BigInteger(tmp);

        BigInteger n_minus_1 = sub(n, one);
        u64vec n_minus_1_bits = n_minus_1.getBits();
        while(getBitAt(n_minus_1_bits, s) == 0)
            ++s;
        u64vec d_bits(n_minus_1_bits.size(), 0);
        for(int i = s; i < n.getBitCnt(); ++i){
            d_bits[(i - s) / 64] |= static_cast<uint64_t >(getBitAt(n_minus_1_bits, i)) << static_cast<uint64_t>((i - s) % 64);
        }
        BigInteger d(d_bits);
        bool is_prime = true;
        vector<BigInteger> checked;
        int upperbound = compare(n, BigInteger(u64vec(1, 50))) == -1 ? n.getBits()[0] - 2 : 50;
        while(checked.size() < upperbound){
            BigInteger a = zero;
            //有一些数被抽取到的概率大一些,避免这些数被重复检测到
            while(true){
                a = randomWithin(two, n_minus_1, generator);
                bool is_in_checked = false;
                for(auto &b:checked){
                    if(compare(a, b) == 0){
                        is_in_checked = true;
                        break;
                    }
                }
                if(!is_in_checked){
                    checked.push_back(a);
                    break;
                }
            }

//            u64vec fake_a;
//            fake_a.push_back(0x468c654b82f33f45);
//            fake_a.push_back(1);
//            a = BigInteger(fake_a);

            BigInteger basis = fastExponent(a, d, n);
            if(compare(basis, one) == 0){
                checked.push_back(a);
                continue;
            }
            else{
                bool neg_one_found = false;
                for(int i = 0; i < s; ++i){
                    if(compare(n_minus_1, basis) == 0){
                        neg_one_found = true;
                        break;
                    }
                    basis = div(mul(basis, basis), n).second;
                }
                if(!neg_one_found){
                    is_prime = false;
                    break;
                }
            }
        }
        if(is_prime){
            cout<< "number tested:" <<tested << endl;
            return n;
        }
        ++tested;
    }
}

void testShr(){
    uint64_t a = 0xffffffffffffffff;
    cout << std::setw(16) << std::setfill('0') << std::hex << __asm_shr(a, 4) << endl;
    cout << std::setw(16) << std::setfill('0') << std::hex << __asm_shr(a, 1) <<endl;
    cout << std::setw(16) << std::setfill('0') << std::hex << __asm_shr(a, 0) <<endl;
    cout << std::setw(16) << std::setfill('0') << std::hex << __asm_shr(a, 64) <<endl;
    cout << std::setw(16) << std::setfill('0') << std::hex << __asm_shr(a, 65) << endl;
    cout << std::setw(16) << std::setfill('0') << std::hex << __asm_shr(a, 63) << endl;
}

void test_asm_sub(){
    vector<uint64_t> a, b;
    a.push_back(0);
    a.push_back(1);
    a.push_back(1);
    b.push_back(1);
    b.push_back(0);
    b.push_back(0);
    __asm_sub_from(a, b);
    BigInteger(a).printHex();
}

void testConstructor(){
    vector<uint64_t> bits1;
    bits1.push_back(2);
    bits1.push_back(4);
    bits1.push_back(0);
    BigInteger a(bits1);
    a.printHex();
}

void testAdd(){
    vector<uint64_t> bits1, bits2;
    bits1.push_back(0xffffffffffffffff);
    //bits1.push_back();
    bits2.push_back(1);
    BigInteger a(bits1), b(bits2);
    BigInteger::add(a, b).printHex();
}

void testSub(){
    vector<uint64_t> bits1, bits2;
    bits1.push_back(0);
    bits1.push_back(0);
    bits1.push_back(1);
    bits2.push_back(2);
    BigInteger::sub(BigInteger(bits1), BigInteger(bits2)).printHex();
}

void testMul(){
//    vector<uint64_t> bits1, bits2;
//    bits1.push_back(0xffffffffffffffff);
//    bits1.push_back(0xffff);
//    bits2.push_back(4);
//    bits2.push_back(2);
//    BigInteger a = BigInteger(bits1), b(bits2);
//    a.printHex();
//    b.printHex();
//    BigInteger::mul(a, b).printHex();


    vector<uint64_t> bits1, bits2;
    bits1.push_back(7);
    bits2.push_back(7);
    BigInteger a = BigInteger(bits1), b(bits2);
    BigInteger::mul(a, b).printHex();
}

void testDiv(){
//    int compact_mode = BigInteger::PRINT_MODE_COMPACT;
//    u64vec a_bits, b_bits;
//    a_bits.push_back(5);
//    a_bits.push_back(0xf);
//    b_bits.push_back(0x1);
//    b_bits.push_back(0x1);
//    BigInteger a(a_bits), b(b_bits);
//    BigIntegerPair p = BigInteger::div(a, b);
//    cout << "divident: ";
//    a.printHex(compact_mode);
//    cout<< "divisor: ";
//    b.printHex(compact_mode);
//
//    cout<< "q: ";
//    p.first.printHex(compact_mode); //print q
//    cout << "r: ";
//    p.second.printHex(compact_mode);    //print r
    u64vec a, b;
    a.push_back(0x9f357a2209ab0161);
    a.push_back(0x595266a26e2d89a5);
    a.push_back(0x2);
    b.push_back(0x4b0ceacd9b81bb47);
    b.push_back(0xb);
    BigIntegerPair qr = BigInteger::div(BigInteger(a), BigInteger(b));
    BigInteger(a).printHex(BigInteger::PRINT_MODE_COMPACT);
    BigInteger(b).printHex(BigInteger::PRINT_MODE_COMPACT);
    qr.first.printHex(BigInteger::PRINT_MODE_COMPACT);
    qr.second.printHex(BigInteger::PRINT_MODE_COMPACT);

}

void testEuclid(){
    u64vec a_bits, b_bits;
    a_bits.push_back(128);
    b_bits.push_back(81);
    BigInteger a(a_bits), b(b_bits);
    BigIntegerPair gcd_v = BigInteger::partialExtendedEuclidean(b, a);
    gcd_v.first.printHex();
    gcd_v.second.printHex();
    BigInteger::euclidean(a, b).printHex();
}

void testRandom(){
    u64vec min_bits, max_bits;
    min_bits.push_back(5);
    max_bits.push_back(10);
    max_bits.push_back(1);
    std::random_device generator;
    for(int i = 0; i < 100; ++i)
    BigInteger::randomWithin(BigInteger(min_bits), BigInteger(max_bits), generator).printHex();
}

void testExponent(){
    u64vec a, e, n;
    a.push_back(7);
    e.push_back(0);
    e.push_back(0);
    e.push_back(1);
    n.push_back(10);
    BigInteger::fastExponent(BigInteger(a), BigInteger(e), BigInteger(n)).printHex();
}

void testPrime(){
    BigInteger::getPrimeWithin(BigInteger::nBitMin(768),BigInteger::nBitMax(768)).printHex(BigInteger::PRINT_MODE_COMPACT);
}