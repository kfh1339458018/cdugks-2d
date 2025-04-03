#include "mpi_process.h"  


  
MPIProcess::MPIProcess() 
{  
    MPI_Init(nullptr, nullptr); // 初始化MPI环境  
    MPI_Comm_rank(MPI_COMM_WORLD, &pid); // 获取每个进程的编号  
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);  //获取总的进程数
    // 可以在这里添加其他初始化代码  
}  
  
void MPIProcess::endMPIProcess() 
{  
    MPI_Finalize(); // 清理MPI环境  
}  
  
void MPIProcess::printRank() const 
{  
    std::cout << "Process " << pid << " is running," ;
    std::cout << "Total Process Num is " << numprocs << std::endl; 
    MPI_Barrier(MPI_COMM_WORLD);  
}


void MPIProcess::print(std::string info) const
{
    if (pid == defaultpid) {  
            std::cout  << info << std::endl;  
        }
}

void MPIProcess::printnum(double info) const
{
    if (pid == defaultpid) {  
            std::cout  << info << std::endl;  
        }
}


//速度空间并行，给每个线程分配任务
void MPIProcess::velAllocate(Block2d& block,std::vector<Velocity2d>& vel,std::vector<Velocity2d>& velq)
{
    int qi_start=(pid*block.q)/numprocs;
    int qi_end=((pid+1)*block.q)/numprocs;
    int qi_local=qi_end-qi_start;  //不同进程分配的长度可能不同

    MPI_Barrier(MPI_COMM_WORLD);
    //输出变量值并标记进程数
    printf("Process：%d, element in vel space：%d, Process num：%d, qi_start：%d, qi_end(do not included)：%d, qi_local：%d\n",
           pid, block.q, numprocs, qi_start, qi_end, qi_local);


    //汇总数据
    int EachNum[numprocs];
    MPI_Allgather(&qi_local, 1, MPI_INT, EachNum, 1, MPI_INT, MPI_COMM_WORLD);


    //偏移量
    pid_to_whole=0;
    for(int procs=1;procs<=pid;procs++)
    {
        pid_to_whole+=EachNum[procs-1];
    }
    //std::cout << pid_to_whole << std::endl;
    
    block.q_local=qi_local;
    block.pid_to_whole=pid_to_whole;
    
    //将vel的数据导入到velq
    for(int direc=0;direc<block.q_local;direc++)
    {
        //相对坐标->绝对坐标
        int direc_ref=direc+block.pid_to_whole;
        velq.push_back(vel[direc_ref]);
    }


    print("vel allocation compeleted.");

    
}
