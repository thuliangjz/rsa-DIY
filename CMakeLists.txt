cmake_minimum_required(VERSION 3.14)
project(rsa)

set(CMAKE_CXX_STANDARD 11)
find_package(Threads)

add_executable(rsa main.cpp RSA.cpp RSA.h BigInteger.cpp BigInteger.h BigIntegerTest.h)
target_link_libraries(rsa ${CMAKE_THREAD_LIBS_INIT})