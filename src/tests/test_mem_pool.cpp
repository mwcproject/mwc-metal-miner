#include "test_mem_pool.h"
#include "../miner/mem_pool.h"

void test_mem_pool() {
    MemPool pool;
    pool.reset(4, 1024*4*32);

    std::string d1 = pool.getAvailableBlockDump();
    assert("Range(0,0, 8), Range(1,0, 8), Range(2,0, 8), Range(3,0, 8)" == d1);

    MemRange ranges[10];
    ranges[0] = pool.allocate(10);
    ranges[1] = pool.allocate(4096*3); // 4
    ranges[2] = pool.allocate(4096*2); // 5
    ranges[3] = pool.allocate(4096*7); // 12
    ranges[4] = pool.allocate(4096*5); // 17
    ranges[5] = pool.allocate(4096*1); // 18

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
