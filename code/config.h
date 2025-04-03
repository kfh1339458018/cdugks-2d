#ifndef CONFIG_H
#define CONFIG_H

#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <memory>

using namespace std;

enum class DistributionType {G,H}; //计算g还是h
enum class TAU_TYPE {Euler, NS}; //松弛时间类型
enum class ENDTIME_TYPE {Converge, Set}; // Converge达到收敛条件停止，Set达到设定时间停止
enum class VELSPACE_TYPE {Unstruct, D2Q12, D2Q37};  //速度空间类型 非结构, D2Q12, D2Q37


class Block2d
{
public:

    int step;  //迭代次数
    double time_calculation=0.0;  //模拟时间
    double timemax;  //退出时间

    //===计算可选项
    TAU_TYPE tau_type;
    ENDTIME_TYPE endtime_type;
    VELSPACE_TYPE velspace_type;
    
    //物理格子
    int TotalCells;  //总的格子数
    int TotalPoints;  //总的格点数
    int BoundaryType;  //边界类型的个数

    double min_length;  //物理网格中，最小的边长
    double max_vel;  //速度空间中，最大的速度绝对值
    double BeginSimulTime;  //开始计算的时间

    //速度格子
    int q_node;  //速度空间的格点数
    int q;  //速度空间的控制体数
    int q_local;  //当强进程中，速度空间的格点数
    int pid_to_whole;  //从进程中的序号转换到整体
    int EndTime;  //迭代步数
    int WriteInterval;  //输出间隔

    //===初始化信息
    string grid_path;  //物理文件（.neu）的路径
    string vel_path;  //速度文件（.neu）的路径
    double alpha;  //CFL数
    double R;  //气体常数
    double Pr;  //普朗特数
    int K;  //分子的振动自由度，单原子为0，双原子为2，多原子为3
    double gama;  //比热比
    double rho_init;  //初始密度
    double temp_init;  //初始温度
    double speed_init;  //初始速度
    double Ma;  //初始马赫数
    double Kn;  //初始努森数
    double Re;  //初始雷诺数
    double Ma_inlet;  //来流马赫数
    double miu_init;  //初始miu0
    double omega;  //omega,计算粘温变化时用，硬球分子模型omega=0.5
    double l_ref;  //参考长度
    double K_limit;  //Venkatakrishnan限制器
    double error_limit;  //判断收敛条件


};

// point coordinate 3d
class Point2d
{
public:
    double x;
    double y;
};


// to remember the cell avg value and geometric information
class Fluid2d
{
public:

	double x;  //格心坐标
	double y;
	
    int edge_num;  //triangle or quadrilateral
    double area;  //格子面积
    int v[5];  //格点编号
    double xv[5];  //格点坐标
    double yv[5];
    double edgelength[5];  //各边边长


    int adj[5];//相邻控制体的编号，如果为边界控制体，则赋值为0
    int adj_xb[5];  //相邻控制体，对边的编号

    //边界上的坐标点（中点）
    double adj_line_x[5];  //中点x坐标
    double adj_line_y[5];  //中点y坐标

    //边界上的坐标点（垂点）
    double adj_vert_x[5];  
    double adj_vert_y[5];

    //垂直于界面的单位向量
    double nx[5];
    double ny[5];

    //判断特殊位置
    //第一维度代表界面的序号位置，第二维度代表边界类型
	int boundary[5][10];  //该界面是否在边界


    //===定义迭代循环中需要的量
    double rho;
    double vx;
    double vy;
    double temp;
    double pressure;
    double e;
    double ev;
    double miu;
    double dt;
    double ds;  //ds=dt/2
    double tau;
    double _tau;
    double qx_flux;
    double qy_flux;
    double vx_history[1000];
    double vy_history[1000];

    double* g;
    double* h;
    double* gs;
    double* hs;
    double* ga;
    double* ha;
    double* gap;
    double* hap;
    double* g_gradient_x;
    double* g_gradient_y;
    double* h_gradient_x;
    double* h_gradient_y;
    double** gapb;
    double** hapb;
    double** Gapb;
    double** Hapb;
    double** gppb;
    double** hppb;

    double* rhob;
    double* vxb;
    double* vyb;
    double* eb;
    double* evb;
    double* tempb;
    double* qxb_flux;
    double* qyb_flux;
    double* taub;  //界面上的松弛时间

    double* Gmeso;
    double* Hmeso;

    double Fmacro_rho;
    double Fmacro_vx;
    double Fmacro_vy;
    double Fmacro_e;

    
};

//velocity space
class Velocity2d
{
public:
    int v[5];  //速度空间格点的编号（针对非结构速度空间）
    
    double vx;  //格心的x坐标
    double vy;  //格心的y坐标

    //面积(权重)
    double area;  //控制体的面积

};



/*分配一个三维数组*/
template <class T>
T*** Allocate3D(int M, int N, int K) {
    T*** tmp = nullptr;
    T* Data = new T[M * N * K];  //分配连续内存块
    tmp = new T**[M];
    for (int i = 0; i < M; i++) {
        tmp[i] = new T*[N];
        for (int j = 0; j < N; j++)
            tmp[i][j] = Data + i * N * K + j * K;
    }
    return tmp;
}

/*删除一个三维数组*/
template <class T>
void Delte3D(T*** buf, int M, int N) {
    if (buf[0][0]) delete[] buf[0][0];
    for (int i = 0; i < M; i++) {
        if (buf[i]) delete[] buf[i];
    }
    if (buf) delete[] buf;
}

/*分配一个二维数组*/
template <class T>
T** Allocate2D(int M, int N) {
    T** tmp = nullptr;
    tmp = new T*[M];
    T* Data = new T[M * N];
    for (int i = 0; i < M * N; i++)
        Data[i] = 0;
    for (int i = 0; i < M; i++)
        tmp[i] = Data + i * N;

    return tmp;
}

/*删除一个二维数组*/
template <class T>
void Delte2D(T** buf, int M) {
    if (buf[0]) delete[] buf[0];
    if (buf) delete[] buf;
}

/*分配一个一维数组*/
template <class T>
T* Allocate1D(int M) {
    T* tmp = new T[M];
    return tmp;
}

/*删除一个一维数组*/
template <class T>
void Delte1D(T* buf) {
    if (buf) delete[] buf;
}

template<typename T>  
void resize2DVector(std::vector<std::vector<T>>& vec, int rows, int cols) {  
    vec.resize(rows);  
    for (auto& row : vec) {  
        row.resize(cols);  
    }  
} 


#endif
