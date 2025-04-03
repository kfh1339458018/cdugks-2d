#ifndef MPI_PROCESS_H
#define MPI_PROCESS_H


# pragma once
#include <mpi.h>
#include <iostream>  
#include <string>
#include <vector>
#include "config.h"

class MPIProcess {  
private:  
    int defaultpid=0;
public:  

    int pid;   //每个进程编号
    int numprocs;  //总共的进程数（不包括0）
    int pid_to_whole;  //速度空间偏移量
    MPIProcess(); // 构造函数，用于初始化MPI和world_rank  
    void endMPIProcess(); // 析构函数，用于清理MPI环境  
    void printRank() const; // 成员函数，用于打印world_rank  

    void print(std::string info) const;
    void printnum(double info) const;

    void velAllocate(Block2d& block,std::vector<Velocity2d>& vel,std::vector<Velocity2d>& velq);  //速度空间分配
  
    // 可以添加其他需要world_rank的成员函数  
};  


#endif



  
