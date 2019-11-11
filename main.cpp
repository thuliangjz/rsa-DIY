#include <iostream>
#include <fstream>
#include "RSA.h"
#include "BigInteger.h"
#include "BigIntegerTest.h"

using std::cout;
using std::endl;
using std::string;

const int CMD_UNKNOWN = -1;
const int CMD_KEYGEN = 0;
const int CMD_ENCRYPT = 1;
const int CMD_DECRYPT = 2;

int strToCmdInt(const std::string& cmd){
    if(cmd == "keygen")
        return CMD_KEYGEN;
    if(cmd == "encrypt")
        return CMD_ENCRYPT;
    if(cmd == "decrypt")
        return CMD_DECRYPT;
    return CMD_UNKNOWN;
}

int safeReadIntFromCmd(const string& s){
    int n;
    try{
        n = std::stoi(s);

    }catch(...){
        std::cerr << "erro parsing int in cmd: "<< s << endl;
        n = -1;
    }
    return n;
}

bool getCharVectorFromFile(string& file, vector<char>& vec){
    std::ifstream file_in(file, std::ios::in|std::ios::binary);
    if(!file_in)
        return false;
    file_in.seekg(0, std::ios::end);
    int length = file_in.tellg();
    file_in.seekg(0, std::ios::beg);
    vec.clear();
    vec.resize(length);
    file_in.read(&vec[0], length);
    return true;
}

bool charVectorToFile(string& file, const vector<char>& vec){
    std::ofstream file_out(file, std::ios::out|std::ios::binary);
    if(!file_out)
        return false;
    file_out.write(&vec[0], vec.size());
    return true;
}

int main(int argc, char** argv) {
//    RSA rsa;
//    rsa.setIterCntForMR(16);
//    rsa.setWorkerCntForMR(8);
//    auto start = std::chrono::system_clock::now();
//    rsa.keyGen(768);
//    rsa.readPrivateKey("id_rsa");
//    rsa.readPubKey("id_rsa.pub");
//    auto end = std::chrono::system_clock::now();
//    rsa.printKeyInfo();
//    std::cout << "time used:" << std::dec << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
//    rsa.saveKey("id_rsa.pub", "id_rsa");
//    std::string msg = "hello from Jianzhe Liang HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHLHHS";
//    vector<char> encrypted = rsa.encrypt(vector<char>(msg.c_str(), msg.c_str() + msg.size()));
//    vector<char> decrypted = rsa.decrypt(encrypted);
//    std::string msg_1(&decrypted[0]);
//    std::cout << msg_1 << std::endl;

    //返回值含义:-1表示不识别的指令(包括空),-2表示参数列表错误, -3表示RSAException, -4为BigIntegerException
    //-5表示要加密的文件打不开, -6表示结果无法写回
    const string help_info = "usage: \n"
                                  "rsa keygen length [out_pub] [out_private] [worker] [iter]\n"
                                  "rsa encrypt pub_key_file msg_file out\n"
                                  "rsa decrypt private_key_file msg_file out\n";
    if(argc < 2){
        cout << help_info;
        return -1;
    }
    string private_key_file, public_key_file, msg_file, out;
    int key_length = 0;
    int worker = 0, iter = 32;
    int cmd = strToCmdInt(argv[1]);
    //先确定参数类型和初始值,再调用RSA
    switch(cmd){
        case CMD_KEYGEN:
            if(argc > 2){
                key_length = safeReadIntFromCmd(argv[2]);
            }
            else{
                cout << "key length need to be specified" << endl;
                return -2;
            }
            public_key_file = argc > 3 ? argv[3] : "id_rsa.pub";
            private_key_file = argc > 4 ? argv[4] : "id_rsa";
            worker = argc > 5 ? safeReadIntFromCmd(argv[5]) : worker;
            iter = argc > 6 ? safeReadIntFromCmd(argv[6]) : iter;
            break;
        case CMD_ENCRYPT:
            if(argc >= 5){
                public_key_file = argv[2];
                msg_file = argv[3];
                out = argv[4];
            }
            else{
                cout << "pub_key_file and msg_file shall be provided"<<endl;
                return -2;
            }
            break;
        case CMD_DECRYPT:
            if(argc >= 5){
                private_key_file = argv[2];
                msg_file = argv[3];
                out = argv[4];
            }
            else{
                cout << "private_key_file and msg_file shall be provided"<<endl;
                return -2;
            }
            break;
        default:
            cout << "unknown command: " << argv[1] <<endl;
            cout << help_info;
            return -1;
    }
    RSA rsa;
    try{
        std::chrono::time_point<std::chrono::system_clock> start, end;
        vector<char> msg;
        switch(cmd){
            case CMD_KEYGEN:
                rsa.setWorkerCntForMR(worker);
                rsa.setIterCntForMR(iter);
                start = std::chrono::system_clock::now();
                rsa.keyGen(key_length);
                end = std::chrono::system_clock::now();
                rsa.printKeyInfo();
                std::cout << "generation time: " << std::dec << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
                rsa.saveKey(public_key_file, private_key_file);
                break;
            case CMD_ENCRYPT:
                rsa.readPubKey(public_key_file);
                if(!getCharVectorFromFile(msg_file, msg)){
                    std::cerr << "error opening msg file" << endl;
                    return -5;
                }
                msg = rsa.encrypt(msg);
                if(!charVectorToFile(out, msg)){
                    std::cerr << "error writing results"<<endl;
                    return -6;
                }
                break;
            case CMD_DECRYPT:
                rsa.readPrivateKey(private_key_file);
                if(!getCharVectorFromFile(msg_file, msg)){
                    std::cerr << "error opening msg file" << endl;
                    return -5;
                }
                msg = rsa.decrypt(msg);
                while(msg.back() == '\0')
                    msg.pop_back();
                if(!charVectorToFile(out, msg)){
                    std::cerr << "error writing results"<<endl;
                    return -6;
                }
                break;
        }
    }catch(RSA::RSAException &e){
        std::cerr << "rsa algorithm exception occurred: " << e.getMessage() << endl;
        return -3;
    }
    catch(BigInteger::BigIntegerException &e){
        std::cerr << "BigInteger runtime exception, err_code: " << e.getCode() << endl;
        return -4;
    }
    return 0;
}

