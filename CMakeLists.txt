cmake_minimum_required(VERSION 3.14)
project(rsa)

set(CMAKE_CXX_STANDARD 11)

add_executable(rsa main.cpp RSA.cpp RSA.h BigInteger.cpp BigInteger.h)