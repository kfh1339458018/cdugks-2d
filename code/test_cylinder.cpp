#include "test_cylinder.h"
#include "cdugks_2d.h"
#include "mesh_read.h"


//初始化参数
void cylinder_param(Block2d& block)
{
    
    block.grid_path=add_mesh_directory_modify_for_linux()+"/../mesh_file/cylinder_kn1_2D_GSI.neu";  //物理文件（.neu）的路径
    block.velspace_type=VELSPACE_TYPE::Unstruct;  
    block.vel_path=add_mesh_directory_modify_for_linux()+"/../mesh_file/unstruct_4521.neu";  //非结构速度文件（.neu）的路径   unstruct_1192

    block.endtime_type=ENDTIME_TYPE::Set;
    block.EndTime=10000000;  //迭代步数
    block.WriteInterval=100;  //输出间隔
    block.timemax=100.0;  

    
    block.tau_type=TAU_TYPE::NS;
    block.l_ref=1.0;  //参考长度
    block.alpha=0.8;  //CFL数
    block.R=0.5;  //气体常数
    block.Pr=2.0/3.0;  //普朗特数
    block.K=0.0;  //分子的振动自由度，单原子为0，双原子为2，多原子为3
    block.gama=(block.K+5.0)/(block.K+3.0);  //比热比
    block.rho_init=1.0;  //参考密度
    block.temp_init=1.0;  //参考温度
    block.Ma=5.0;  //参考马赫数
    block.Kn=1.0;  //初始努森数
    //block.Re;  //初始雷诺数
    block.omega=0.81;  //omega,计算粘温变化时用，硬球分子模型omega=0.5
    block.K_limit=0.0;  //限制器中的可调参数K_limit
    block.error_limit=1e-9;

    block.Ma_inlet=block.Ma;  //inlet界面来流马赫数
    
}

//初场条件
void cylinder_init_param(Block2d& block,std::vector<Fluid2d>& fluids)
{

    for (auto& fluid : fluids) {
    
        fluid.rho = block.rho_init;
        fluid.temp = block.temp_init;
        fluid.pressure = fluid.rho*block.R*block.temp_init; // 计算压强
        
        //计算声速
        double sound_speed=sqrt(block.gama*block.R*block.temp_init);
        if(sqrt(pow(fluid.x,2)+pow(fluid.y,2))>=5.0)
        {
            fluid.vx = block.Ma*sound_speed;
        }
        else
        {
            fluid.vx = block.Ma*sound_speed*(sqrt(pow(fluid.x,2)+pow(fluid.y,2))-1)/(5.0-1.0);
        }
        fluid.vy = 0.0;
    } 
}



void cylinder()
{
    MPIProcess process; // 创建MPIProcess实例
    process.printRank(); //每个进程print


    make_directory_for_result();  //创建结果文件夹
    
    Block2d block;
    
    std::vector<Point2d> points;
    std::vector<Point2d> points_v;
    std::vector<Fluid2d> fluids;
    std::vector<Velocity2d> vel;  //总的速度空间
    std::vector<Velocity2d> velq;  //相对速度空间
    
    
    block.BeginSimulTime = MPI_Wtime();  //记录开始时间
    

    //config
    cylinder_param(block);
    
    
    
    Read_mesh_neu_2d(process,fluids,block,points);
    if(block.velspace_type==VELSPACE_TYPE::Unstruct)
    {
        Unstruct_velocity(process,vel,block,points_v);
    }
    else if(block.velspace_type==VELSPACE_TYPE::D2Q12)
    {
        Struct_velocity_D2Q12(process,vel,block,points_v);
    }
    else if(block.velspace_type==VELSPACE_TYPE::D2Q37)
    {
        Struct_velocity_D2Q37(process,vel,block,points_v);
    }
    	

    process.velAllocate(block,vel,velq);  //速度空间分配
    

    //流场初始量输入
    cylinder_init_param(block,fluids);
    //流场初始化
    init(process,block,fluids,points,vel,velq); 

    //判断程序是否应该退出
    bool stop=false;

    //迭代循环
    for(block.step=1;block.step<=block.EndTime;block.step++)
    {
        //辅助分布函数
        AuxDistribution(process,block,fluids,vel,velq,DistributionType::G);
        AuxDistribution(process,block,fluids,vel,velq,DistributionType::H);
        
        //梯度插值
        InterPolation(process,block,fluids,vel,velq,DistributionType::G);
        InterPolation(process,block,fluids,vel,velq,DistributionType::H);

        //界面宏观量
        BoundaryMacroscopicVariables(process,block,fluids,vel,velq);
        
        //界面平衡态与dt/2时刻界面真实分布
        EquilibriumAndNext(process,block,fluids,vel,velq,DistributionType::G);
        EquilibriumAndNext(process,block,fluids,vel,velq,DistributionType::H);

        //边界条件
        BoundaryCondition(process,block,fluids,vel,velq);
        //微观通量
        CalcFmeso(process,block,fluids,vel,velq);
        //宏观通量
        CalcFmacro(process,block,fluids,vel,velq);
        //格心宏观量更新
        UpdateMacroscopicVariables(process,block,fluids,vel,velq);
        //格心微观分布更新
        UpdataDistribution(process,block,fluids,vel,velq);

        MPI_Barrier(MPI_COMM_WORLD);

        stop=SteadyJudge(process,block,fluids,vel,velq);

        if (stop) { 
            if (block.endtime_type  == ENDTIME_TYPE::Converge) 
            { 
                outPutFile_dat(process, "steady", block, fluids, points); 
            } else if (block.endtime_type  == ENDTIME_TYPE::Set) 
            { 
                outPutFile_dat(process, "final", block, fluids, points); 
            } 
            break; 
        }
        
        if(block.step!=0 && block.step%block.WriteInterval==0)
        {
            //输出图
            outPutFile_dat(process,"result_",block,fluids,points);
        }
    }
    
    /*输出运行的时长*/
    double EndSimulTime = MPI_Wtime();
    if (process.pid == 0) 
    {
        std::cout << "CDUGKS end" << std::endl;
        std::cout << "Running Time: " << EndSimulTime - block.BeginSimulTime << 's'
                  << std::endl;
    }

    process.endMPIProcess();
}

