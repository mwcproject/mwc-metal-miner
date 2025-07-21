#include "test_mem_pool.h"
#include "../miner/mem_pool.h"
#include "../miner/metal.h"

void test_mem_pool() {
    MemPool pool;
    pool.reset(4, 1024*4*32);

    std::string d1 = pool.getAvailableBlockDump();
    assert("Range(0,0, 8), Range(1,0, 8), Range(2,0, 8), Range(3,0, 8)" == d1);

    MetalEncoderManager enc_man(nullptr);

    MemRange ranges[10];
    ranges[0] = pool.allocate(10, enc_man);
    ranges[1] = pool.allocate(4096*3, enc_man); // 4
    ranges[2] = pool.allocate(4096*2, enc_man); // 5
    ranges[3] = pool.allocate(4096*7, enc_man); // 12
    ranges[4] = pool.allocate(4096*5, enc_man); // 17
    ranges[5] = pool.allocate(4096*1, enc_man); // 18

    std::string d2 = pool.getAvailableBlockDump();
    assert("Range(0,6, 8), Range(2,5, 8), Range(3,0, 8)" == d2);

    /*pool.release(ranges[0]);
    pool.release(ranges[2]);
    pool.release(ranges[5]);

    ranges[0] = pool.allocate(4096*8); // 26
    ranges[2] = pool.allocate(4096*4); // 30
    ranges[5] = pool.allocate(4096*3); // 33

    std::string d3 = pool.getAvailableBlockDump();
    assert("Range(0,0, 1), Range(1,7, 8)" == d3); */

}
