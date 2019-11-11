# RSA-DIY

一个只借助stl实现的简单且相对高性能的c++ rsa加密算法。1024位以内秘钥可以实现1s内生成，2048位5s内生成。

## 程序构建和使用说明

### 硬件要求

由于使用了gcc的内联汇编特性，最终的可执行文件中包含x64汇编指令，请确保使用的cpu能够运行x64指令集。

### 操作系统要求

开发操作系统为ubuntu 18.04，主要依赖gcc，g++和cmake，绝大多数Linux平台都应该可以支持。

### 构建

程序的构建采用标准的`cmake`构建流程：

~~~bash
mkdir build
cd build
cmake ../ -DCMAKE_BUILD_TYPE=Release
make -j
~~~

执行完成之后目录下的`rsa`即为可运行文件，该文件包含一个简单的命令行接口。

### 使用说明

__密钥生成__

~~~bash
./rsa keygen length [out_pub] [out_private] [worker] [iter]
~~~

其中`length`表示生成的密钥长度，有`[]`的表示可选参数，其中，`out_pub`表示生成的公钥保存的文件名称(默认为`id_rsa.pub`)，`out_private`表示生成的私钥保存的文件名称（默认为`id_rsa`），`worker`表示在用[Miller-Rabin](https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test)检测生成素数时使用线程的个数（默认为0表示使用最大并行线程数），`iter`表示在Miller-Rabin检测中生成素数时检验数的个数，默认为32。

下面的指令：

~~~bash
./rsa keygen 768 rsa_768.pub rsa_768 8 16
~~~

即生成一个768位的密钥，将公钥保存在当前目录下名为`rsa_768.pub`的文件中，将私钥保存在`rsa_768`中，检测素数时使用8个线程同时对每个待检验的素数测试16次。

参考的运行结果如下：

~~~
max thread supported: 12, work count set to: 8
iteration count set to: 50
key length: 768
n(modulo): 0xdc...
p: 0x83...
q: 0x1a...
d: 0x7a...
generation time: 342ms

~~~

运行结果中包含了机器支持的最大运行线程数以及实际使用的线程数(第一行)，输出的n, d, e, p, q的含义都和[标准RSA算法](https://en.wikipedia.org/wiki/RSA_(cryptosystem)#Key_generation)一致，输出的值均为16进制。

程序同样还测试了运行的时间(time cost行)，对于768位及以下的密钥，绝大多数情况下生成时间都远远小于1s

__加密__

~~~bash
./rsa encrypt pub_key_file msg_file out
~~~

后面的三个参数分别表示`keygen`指令生成的`out_pub`文件，要加密的文件和输出文件

~~~bash
./rsa keygen 1024
wget www.baidu.com	-O index.html
./rsa encrypt id_rsa.pub index.html encrypted
~~~

上面三条指令，第一条按照默认配置生成了一个长度为1024位的秘钥，第二条获取了`www.baidu.com`的首页并将结果存储在`index.html`中。第三条指令读取生成的1024位公钥，将`index.html`加密之后将结果保存在名为`encrypted`的文件中

__解密__

~~~bash
./rsa  decrypt private_key_file msg_file out
~~~

后面的三个参数分别表示`keygen`指令生成的`out_private`文件，要解密的文件和输出文件

~~~bash
./rsa keygen 2048
wget www.baidu.com	 -O  index.html
./rsa encrypt id_rsa.pub index.html encrypted
./rsa decrypt id_rsa encrypted decrypted
diff decrypted index.html
~~~

前三句和加密的过程相同。第四句表示将加密后的`encrypted`文件用私钥`id_rsa`解密之后输出到decrypted。最后一句用`diff`比较了两个文件的`decrypted`和`index.html`的区别，没有任何输出，说明decrypted和index之间没有区别。

## 加速方案介绍

项目中涉及到的加速方案包括：基于汇编的大数加减乘法，基于定点逆的快速幂，多线程素性测试

### 基于汇编的大数加减乘法

项目中所有的大数均采用二进制进行表示。一方面，这样可以方便地得到大数本身的位数信息；另一方面，借助于底层的汇编指令，我们可以一次性完成多位的加减乘法运算同时获取该次运算的进位/溢出。

__加法：__

实现在`BigInteger.cpp`的`__asm_add`函数中。核心的代码为

~~~asm
sahf;
movq (%%r8, %%rcx, 8), %%rax;
movq (%%r9, %%rcx, 8), %%rbx;
adc %%rax, %%rbx;
lahf;
movq %%rbx, (%%r10, %%rcx, 8);
~~~

在两个大数按位相加的过程中，每一次相加之前先将上一次相加的进位恢复（`sahf`），然后借助`adc`指令将该位上的两个数连同上一次的进位(`CARRY FLAG`)一起加起来，同时保存本次进位，保存结果之后增加索引进入下一轮相加即可

__减法__

对于x64指令集，减法借位的效果和加法进位的效果一样，都是设置`CARRY FLAG`。把上面的`adc`改为`sbb`即可，实现在`BigInteger::sub`

__乘法__

总体上乘法分成两步。

+ 用一个乘数中的某一个整型去乘以另一个大整数得到一个（移位的）结果
+ 将该结果加到最终的结果上

第一步的核心代码如下：

~~~asm
mov %%r8, %%rax;	//小乘数保存在r8中
mov (%%rsi, %%rcx, 8), %%rbx;	//取出大乘数
mul %%rbx;
mov %%rax, %%rbx;
ASM_RESTORE
adc %%r14, %%rbx;	//本轮结果即本次乘法+上次乘法溢出+上次加法溢出
ASM_LOAD
mov %%rdx, %%r14;	//mul指令溢出保存在r14中
mov %%rbx, (%%rdi, %%rcx, 8);
~~~

这一步结果中的每一个整型都是小乘数乘以大乘数中的一个整型加上上一次相乘的溢出（包括乘法和加法一起）得到的。每次`mul`指令的乘法溢出被保留到`r14`寄存器中。`ASM_RESTORE`和`ASM_LOAD`是一对宏，其效果和`sahf`以及`lahf`相同，只不过是将`CARRY FLAG`的值保存到`r15`而非`ah`中，因为`rax`要保存`mul`的乘数。

第二步的基本思路仍然是借助于`adc`指令。但是每次相加的长度不需要等于最终结果的长度，只要等于第一步长度的最大值加1即可。

实现请参看`BigInteger::mul`

__除法__

项目中实现的除法为长除法(`BigInteger::div`)，注意到在确定商的时候由于除数可能存在很多（> 64）位，所以不能使用汇编的整型除法指令。故这里实现的实际上是基于二进制的整形除法，每次确定商的一个比特位。长除法中的乘法全部被移位所替代，因为每次除数需要乘的都是2的幂次，这在一定程度上加快了除法的速度。然而，考虑到在RSA生成密钥以及加密的过程中取模运算是不可避免的，完全借助长除法来实现取模仍然会带来比较大的开销。所以项目中还实现了更高效的快速幂算法。

### 基于定点逆的快速幂

基本的想法是借助[Newton-Raphson](https://en.wikipedia.org/wiki/Division_algorithm#Newton%E2%80%93Raphson_division)算法，用乘法和加减法来代替除法。在计算$x \pmod{n}$时先计算并保留$\frac 1 n$，并考虑将$\lfloor \frac 1n \times x \rfloor$作为x除以n的商。这里将$\frac 1 n$以定点数的形式保存，这样可以调用之前整数的运算来实现小数和小数，小数和整数之间的加法，减法和乘法。

关于定点数的精度，对于计算 $ a^b \pmod{n}$，在使用[快速幂算法]( https://en.wikipedia.org/wiki/Modular_exponentiation#Right-to-left_binary_method )进行迭代计算的过程中，每一步实际上都是对于某个$x \le n$计算$x^2 \pmod n$，所以实际上被除数不大于$n^2$。如果$n$有$p$位，则$n^2$不超过$2p$位，从而对于这样的$x^2$，逆的精度只要精确到小数点后$2p$位即误差不超过$2^{-2p}$即可保证与被除数相乘与商的误差不超过1。

注意到如果将$\lfloor \frac 1n \times x \rfloor $作为商，即使误差本身不超过1，向下取整之后可能会达到1（向下取整在这里实现为移位运算），从而按上述方法计算的结果可能需要减去n才能得到真正比n小的余数。

实现参看`BigInteger::fastExponentNewton`，另外一个`BigInteger::fastExponent`是用长除法完成的快速幂，速度慢不知道到哪去了...

### 多线程素性测试

在使用Miller-Rabin检测生成随机素数时，项目中开启多个线程并行地生成随机数并进行检测，多个线程之间共享一个标志变量，一旦有一个线程检查到标志变量则修改该变量，其他线程检测到修改之后停止查找并退出。考虑到当位数比较多的时候，不同线程之间找到相同数的概率非常小，故没有设置一个互斥的容器来存储所有线程检查过的数，能够有效地增加并行度。

参看`BigInteger.cpp`中的`getPrimeWorker`

## 收获

+ 复习并进一步学习了汇编语言。本科的时候学习的是基于x86的MSAM，这次作业使用的是x64的汇编，且使用AT&T的语法，对gcc的内联汇编特性也有了一定的了解。
+ 课外自主学习了Miller-Rabin检测以及Newton-Raphson除法，自主推导了逆的定点数精度，通过移位和反转等操作借助已有的整数加减乘法操作实现定点小数和整数的加减乘法
+ 在调试一些非常深的bug时增加了经验
  + 长除法一开始会在结尾莫名其妙的出错，用gdb断点一步一步跟进之后最后发现是逻辑右移64位会导致结果不变而不是全零
  + 大量的向量角标索引操作难免会出现越界访问和修改。很多时候结果都正确但是退出的时候报错。有一次是`free`函数出了问题，最后是在gdb中获取了free的`chunk`地址，在`chunk`的头部打了硬件断点之后监视每次的修改最后定位到具体的越界访问位置。