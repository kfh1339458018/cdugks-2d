# include "cdugks_2d.h"
# include "mesh_read.h"
# include "mpi_process.h"
# include "output.h"
# include "config.h"
# include <iostream>
# include <ostream>

//===测试算例
#include "test_cylinder_half.h"
#include "test_cylinder.h"
#include "test_cavity.h"
#include "test_riemann_2d.h"
#include "test_riemann_1d.h"


int main(int argc, char** argv)
{

    
    cylinder_half();  //高超声速半圆柱绕流
    //cylinder();  //高超声速圆柱扰流
    //cavity();  //顶盖驱动方腔流动
    //riemann_2d();  //二维黎曼问题
    //riemann_1d();  //一维黎曼问题

    
    return 0;
}

