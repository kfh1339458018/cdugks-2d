# include "cdugks_2d.h"


//平衡态
double Calc_s(Block2d* block, Fluid2d* fluids, Velocity2d* velq,DistributionType type)
{
    switch (type) {  
        case DistributionType::G:  
            return CalcGs(block, fluids, velq);  
        case DistributionType::H:  
            return CalcHs(block, fluids, velq);  
        default:  
            // 处理错误情况，比如未知的DistributionType  
            throw std::invalid_argument("Unknown DistributionType");  
    } 
}


  
double CalcGs(Block2d* block, Fluid2d* fluids, Velocity2d* velq) {  
    double geq = 0.0, gpr = 0.0, gs = 0.0;  
    double cx = 0.0, cy = 0.0;  // 特征速度  
    double _pressure = fluids->rho * block->R * fluids->temp;  
  
    // 计算特征速度  
    cx = velq->vx - fluids->vx;  
    cy = velq->vy - fluids->vy;  
  
    // 计算geq  
    geq = (fluids->rho / (2.0 * M_PI * block->R * fluids->temp)) *  
          exp(-(cx * cx + cy * cy) / (2.0 * block->R * fluids->temp));  
  
    // 计算gpr  
    gpr = (1.0 - block->Pr) * (cx * fluids->qx_flux + cy * fluids->qy_flux) /  
          (5.0 * _pressure * block->R * fluids->temp) *  
          (((cx * cx + cy * cy) / (block->R * fluids->temp)) - 2.0 - 2.0) * geq;  
  
    // 计算gs  
    gs = geq + gpr;  
    gs = gs * velq->area;  
  
    return gs;  
}

double CalcHs(Block2d* block, Fluid2d* fluids, Velocity2d* velq) {  
    double geq = 0.0, heq = 0.0, hpr = 0.0, hs = 0.0;  
    double cx = 0.0, cy = 0.0;  // 特征速度  
    double _pressure = fluids->rho * block->R * fluids->temp;  
  
    // 计算特征速度  
    cx = velq->vx - fluids->vx;  
    cy = velq->vy - fluids->vy;  
  
    // 计算geq  
    geq = (fluids->rho / (2.0 * M_PI * block->R * fluids->temp)) *  
          exp(-(cx * cx + cy * cy) / (2.0 * block->R * fluids->temp));  
  
    // 计算heq  
    heq = (block->K + 3.0 - 2.0) * block->R * fluids->temp * geq;  
  
    // 计算hpr  
    hpr = (1.0 - block->Pr) * (cx * fluids->qx_flux + cy * fluids->qy_flux) /  
          (5.0 * _pressure * block->R * fluids->temp) *  
          (((cx * cx + cy * cy) / (block->R * fluids->temp) - 2.0) * (block->K + 3.0 - 2.0) - 2.0 * block->K) *  
          block->R * fluids->temp * geq;  
  
    // 计算hs  
    hs = heq + hpr;  
    hs = hs * velq->area;  
  
    return hs;  
}


//界面上平衡态
double bCalc_s(Block2d* block, Fluid2d* fluids, Velocity2d* velq,int adj,DistributionType type)
{
    switch (type) {  
        case DistributionType::G:  
            return bCalcGs(block, fluids, velq, adj);  
        case DistributionType::H:  
            return bCalcHs(block, fluids, velq, adj);  
        default:  
            // 处理错误情况，比如未知的DistributionType  
            throw std::invalid_argument("Unknown DistributionType");  
    } 
}

double bCalcGs(Block2d* block, Fluid2d* fluids, Velocity2d* velq, int adj) {  
    double geq = 0.0, gpr = 0.0, gs = 0.0;  
    double cx = 0.0, cy = 0.0;  // 特征速度  
    double _pressure = fluids->rhob[adj] * block->R * fluids->tempb[adj];  
  
    // 计算特征速度  
    cx = velq->vx - fluids->vxb[adj];  
    cy = velq->vy - fluids->vyb[adj];  
  
    geq = (fluids->rhob[adj] / (2.0 * M_PI * block->R * fluids->tempb[adj])) *  
          exp(-(cx * cx + cy * cy) / (2.0 * block->R * fluids->tempb[adj]));  
  
    gpr = (1.0 - block->Pr) * (cx * fluids->qxb_flux[adj] + cy * fluids->qyb_flux[adj]) /  
          (5.0 * _pressure * block->R * fluids->tempb[adj]) *  
          ((cx * cx + cy * cy) / (block->R * fluids->tempb[adj]) - 2.0 - 2.0) * geq;  
  
    gs = geq + gpr;  
    gs = gs * velq->area;  
    return gs;  
}

double bCalcHs(Block2d* block, Fluid2d* fluids, Velocity2d* velq, int adj) {  
    double geq = 0.0, heq = 0.0, hpr = 0.0, hs = 0.0;  
    double cx = 0.0, cy = 0.0;  // 特征速度  
    double _pressure = fluids->rhob[adj] * block->R * fluids->tempb[adj];  
  
    // 计算特征速度  
    cx = velq->vx - fluids->vxb[adj];  
    cy = velq->vy - fluids->vyb[adj];  
  
    // geq  
    geq = (fluids->rhob[adj] / (2.0 * M_PI * block->R * fluids->tempb[adj])) *  
          exp(-(cx * cx + cy * cy) / (2.0 * block->R * fluids->tempb[adj]));  
    
    // heq  
    heq = (block->K + 3.0 - 2.0) * block->R * fluids->tempb[adj] * geq;  
  
    // hpr  
    hpr = (1.0 - block->Pr) * (cx * fluids->qxb_flux[adj] + cy * fluids->qyb_flux[adj]) /  
          (5.0 * _pressure * block->R * fluids->tempb[adj]) *  
          (((cx * cx + cy * cy) / (block->R * fluids->tempb[adj]) - 2.0) * (block->K + 3.0 - 2.0) - 2.0 * block->K) *  
          block->R * fluids->tempb[adj] * geq;  
  
    hs = heq + hpr;  
    hs = hs * velq->area;  
    return hs;  
}

// 判断点是否在三角形内，返回值：0-边界上，1-内部，-1-外部
int pointLocationInTriangle(double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Px, double Py) 
{
    // 向量AB、BC、CA
    double AB_x = Bx - Ax;
    double AB_y = By - Ay;
    double BC_x = Cx - Bx;
    double BC_y = Cy - By;
    double CA_x = Ax - Cx;
    double CA_y = Ay - Cy;

    // 向量AP、BP、CP
    double AP_x = Px - Ax;
    double AP_y = Py - Ay;
    double BP_x = Px - Bx;
    double BP_y = Py - By;
    double CP_x = Px - Cx;
    double CP_y = Py - Cy;

    // 计算点与每条边的外积
    double cross1 = AB_x * AP_y - AB_y * AP_x;
    double cross2 = BC_x * BP_y - BC_y * BP_x;
    double cross3 = CA_x * CP_y - CA_y * CP_x;

    // 误差范围
    double epsilon = 1e-9;

    // 如果点在边上，返回0；如果在内部，返回1；如果在外部，返回-1
    if (fabs(cross1) < epsilon || fabs(cross2) < epsilon || fabs(cross3) < epsilon) {
        return 0;  // 在边界上(cross1 cross2 cross3为同一符号)
    } else if ((cross1 >= -epsilon && cross2 >= -epsilon && cross3 >= -epsilon) ||
               (cross1 <= epsilon && cross2 <= epsilon && cross3 <= epsilon)) {
        return 1;  // 在内部
    } else {
        return -1; // 在外部
    }
}



// 判断点是否在四边形内，返回值：0-边界上，1-内部，-1-外部
int pointLocationInQuadrilateral(double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Dx, double Dy, double Px, double Py)
{

    int result1 = pointLocationInTriangle(Ax, Ay, Bx, By, Cx, Cy, Px, Py);
    int result2 = pointLocationInTriangle(Ax, Ay, Cx, Cy, Dx, Dy, Px, Py);

    // 如果点在边上，返回0；如果在内部，返回1；如果在外部，返回-1
    if ((result1 == 0 && result2!=0) || (result2 == 0 && result1!=0)) {
        return 0;  // 在边界上
    } else if ((result1 == 1 || result2 == 1) || (result1==0 && result2==0)) {
        return 1;  // 在内部
    } else {
        return -1; // 在外部
    }
}


//流场初始化
void init(MPIProcess& process,Block2d& block,std::vector<Fluid2d>& fluids,std::vector<Point2d>& points,std::vector<Velocity2d>& vel,std::vector<Velocity2d>& velq)
{
    block.step=1;
    if(block.tau_type==TAU_TYPE::NS)
    {
        //雷诺数
        block.Re=sqrt(block.gama/M_PI)*(sqrt(2.0)*(5.0-2.0*block.omega)*(7.0-2.0*block.omega)/15.0)*block.Ma/block.Kn;
        //process.printnum(block.Re);  //输出雷诺数
        block.speed_init=block.Ma*sqrt(block.gama*block.R*block.temp_init);        
        
        //粘度
        block.miu_init=block.rho_init*block.speed_init*block.l_ref/block.Re;
        for(auto&fluid : fluids)
        {
            fluid.miu=block.miu_init*pow(fluid.temp/block.temp_init,block.omega);
            
        }
        //printf("%lf",fluids[1].miu);

        //最大温度
        double temp_max = fluids[0].temp;
        // 遍历数组，找到最大值
        for (int tt = 1; tt <block.TotalCells; tt++) {
            if (fluids[tt].temp > temp_max) {
                temp_max = fluids[tt].temp;
            }
        }
        //计算dt,ds
        for(auto&fluid : fluids)
        {
            fluid.dt=block.alpha*block.min_length/((block.max_vel+block.Ma*sqrt(block.gama*block.R*temp_max)));
            fluid.ds=fluid.dt/2.0;
        }
        //tau
        for(auto&fluid : fluids)
        {
            fluid.tau=fluid.miu/fluid.pressure;
        }
    }

    if(block.tau_type==TAU_TYPE::Euler)
    {       

        //最大温度
        double temp_max = fluids[0].temp;
        // 遍历数组，找到最大值
        for (int tt = 1; tt <block.TotalCells; tt++) {
            if (fluids[tt].temp > temp_max) {
                temp_max = fluids[tt].temp;
            }
        }

        double C1=0.05;
        double C2=1.0;

        for(auto&fluid : fluids)
        {
            fluid.dt=block.alpha*block.min_length/((block.max_vel+block.Ma*sqrt(block.gama*block.R*temp_max)));
            fluid.ds=fluid.dt/2.0;
            fluid.miu=C1*fluid.dt;
        }
        //printf("%lf",fluids[1].miu);

       
        //tau
        for(auto&fluid : fluids)
        {
            fluid.tau=fluid.miu/fluid.pressure;
        }
    }

    //内存分配
    for(auto&fluid : fluids)
    {
        fluid.g=Allocate1D<double>(block.q_local);
        fluid.h=Allocate1D<double>(block.q_local);
        fluid.gs=Allocate1D<double>(block.q_local);
        fluid.hs=Allocate1D<double>(block.q_local);
        fluid.ga=Allocate1D<double>(block.q_local);
        fluid.ha=Allocate1D<double>(block.q_local);
        fluid.gap=Allocate1D<double>(block.q_local);
        fluid.hap=Allocate1D<double>(block.q_local);
        fluid.g_gradient_x=Allocate1D<double>(block.q_local);
        fluid.g_gradient_y=Allocate1D<double>(block.q_local);
        fluid.h_gradient_x=Allocate1D<double>(block.q_local);
        fluid.h_gradient_y=Allocate1D<double>(block.q_local);
        fluid.Gmeso=Allocate1D<double>(block.q_local);
        fluid.Hmeso=Allocate1D<double>(block.q_local);

        
        int row=5;
        int cols=block.q_local;
        fluid.gapb=Allocate2D<double>(row,cols);
        fluid.hapb=Allocate2D<double>(row,cols);
        fluid.Gapb=Allocate2D<double>(row,cols);
        fluid.Hapb=Allocate2D<double>(row,cols);
        fluid.gppb=Allocate2D<double>(row,cols);
        fluid.hppb=Allocate2D<double>(row,cols);

        fluid.rhob=Allocate1D<double>(row);
        fluid.vxb=Allocate1D<double>(row);
        fluid.vyb=Allocate1D<double>(row);
        fluid.eb=Allocate1D<double>(row);
        fluid.evb=Allocate1D<double>(row);
        fluid.tempb=Allocate1D<double>(row);
        fluid.qxb_flux=Allocate1D<double>(row);
        fluid.qyb_flux=Allocate1D<double>(row);
        fluid.taub=Allocate1D<double>(row);

        
    } 

    //初始分布函数
    for(auto&fluid : fluids)
    {
        fluid.qx_flux=0.0;
        fluid.qy_flux=0.0;
        double add=0.0;

        for(int direc=0;direc<block.q_local;direc++)
        {
            fluid.g[direc]=CalcGs(&block,&fluid,&velq[direc]);
            fluid.h[direc]=CalcHs(&block,&fluid,&velq[direc]);   
            //add+=fluid.g[direc];  
        }
        
             
    }

    //初始能量
    for(auto&fluid : fluids)
    {
        fluid.e=0.0;
        fluid.ev=0.0;
        for(int direc=0;direc<block.q_local;direc++)
        {
            fluid.e+=0.5*((velq[direc].vx*velq[direc].vx+velq[direc].vy*velq[direc].vy)*fluid.g[direc]+fluid.h[direc]);
        }
        MPI_Allreduce(MPI_IN_PLACE, &fluid.e, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    for(auto&fluid : fluids)
    {
        fluid.e/=fluid.rho;
        fluid.ev=fluid.e-0.5*(fluid.vx*fluid.vx+fluid.vy*fluid.vy);
    }

    //output init
    outPutFile_dat(process,"init",block,fluids,points);
    
}





//辅助分布函数
void AuxDistribution(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq, DistributionType type) 
{  
    for (auto& fluid : fluids) {  
        for (int direc = 0; direc < block.q_local; direc++) {  

            
            // 注意：这里的 ap 计算可能需要根据实际情况调整，上面的公式可能是一个简化的例子  
            if (type == DistributionType::G) {  
                fluid.gs[direc] = CalcGs(&block,&fluid,&velq[direc]); 
                fluid.ga[direc]=(2.0 * fluid.tau + fluid.dt) / (2.0 * fluid.tau) * fluid.g[direc] - fluid.dt / (2.0 * fluid.tau) * fluid.gs[direc];  
                fluid.gap[direc]=(2.0 * fluid.tau - fluid.ds) / (2.0 * fluid.tau + fluid.dt) * fluid.ga[direc] + (3.0 * fluid.ds) / (2.0 * fluid.tau + fluid.dt) * fluid.gs[direc];
            } else if (type == DistributionType::H) {  
                fluid.hs[direc] = CalcHs(&block,&fluid,&velq[direc]); 
                fluid.ha[direc]=(2.0 * fluid.tau + fluid.dt) / (2.0 * fluid.tau) * fluid.h[direc] - fluid.dt / (2.0 * fluid.tau) * fluid.hs[direc];  
                fluid.hap[direc]=(2.0 * fluid.tau - fluid.ds) / (2.0 * fluid.tau + fluid.dt) * fluid.ha[direc] + (3.0 * fluid.ds) / (2.0 * fluid.tau + fluid.dt) * fluid.hs[direc];
            }
            
        }  
    }  
}



//Venkatakrishnan 限制器
double Venkatakrishnan(double max,double min,double delta,double h,double K_limit)
{
    double error_=1e-30;
    double omega_=pow(K_limit*h,3);

    double a=0.0,b=0.0;

    double fai=0.0;

    //delta>0
    if(delta>error_)
    {
        a=max;
        b=delta;

        fai=(a*a+2.0*a*b+omega_)/(a*a+2.0*b*b+a*b+omega_);

        if(fai>1.0)
        {
            fai=1.0;
        }
    }

    //delta<0
    if(delta< -1.0*error_)
    {
        a=min;
        b=delta;

        fai=(a*a+2.0*a*b+omega_)/(a*a+2.0*b*b+a*b+omega_);

        if(fai>1.0)
        {
            fai=1.0;
        }
    }

    return fai;

}


//界面插值
void InterPolation(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq, DistributionType type) 
{
    //WM
    double** max=Allocate2D<double>(block.TotalCells,block.q_local);
    //Wm
    double** min=Allocate2D<double>(block.TotalCells,block.q_local);

    //最小二乘法求梯度
    for(int i=0;i<block.TotalCells;i++)
    {
        //相邻的格心坐标
        double centerx_adj[5],centery_adj[5];
        for(int adj=1;adj<=fluids[i].edge_num;adj++)
        {
            //如果相邻的控制体为0,跳过
            if(fluids[i].adj[adj]==0)continue;
            centerx_adj[adj]=fluids[fluids[i].adj[adj]].x;
            centery_adj[adj]=fluids[fluids[i].adj[adj]].y;
        }
        
        for (int direc = 0; direc <block.q_local; direc++) 
        {
            
            // 初始化变量
            double gradientx = 0.0; // 初始化x方向梯度
            double gradienty = 0.0; // 初始化y方向梯度

            if(type==DistributionType::G)
            {
                max[i][direc]=fluids[i].gap[direc];
                min[i][direc]=fluids[i].gap[direc];
            }
            else if(type==DistributionType::H)
            {
                max[i][direc]=fluids[i].hap[direc];
                min[i][direc]=fluids[i].hap[direc];
            }
            
            
            double x11 = 0.0, x12 = 0.0, x21 = 0.0, x22 = 0.0, y1 = 0.0, y2 = 0.0;
            
            for (int adj = 1; adj <= fluids[i].edge_num; adj++) 
            {
                x11 += (centerx_adj[adj] - fluids[i].x) * (centerx_adj[adj] - fluids[i].x);  
                x12 += (centerx_adj[adj] - fluids[i].x) * (centery_adj[adj] - fluids[i].y);  
                x21 += (centerx_adj[adj] - fluids[i].x) * (centery_adj[adj] - fluids[i].y);  
                x22 += (centery_adj[adj] - fluids[i].y) * (centery_adj[adj] - fluids[i].y);  
        
                double diff=0.0;  
                if (type == DistributionType::G) {  
                    diff = fluids[fluids[i].adj[adj]].gap[direc] - fluids[i].gap[direc];  
                } else if (type == DistributionType::H) {  
                    diff = fluids[fluids[i].adj[adj]].hap[direc] - fluids[i].hap[direc];  
                }  
        
                y1 += (centerx_adj[adj] - fluids[i].x) * diff;  
                y2 += (centery_adj[adj] - fluids[i].y) * diff;  
        
                if (diff > 0) 
                {   
                    auto value = (type == DistributionType::G) ? fluids[i].gap[direc] : fluids[i].hap[direc];  
                    max[i][direc] = std::max(max[i][direc], value);  
                }
                if (diff < 0) 
                {  
                    auto value = (type == DistributionType::G) ? fluids[i].gap[direc] : fluids[i].hap[direc];   
                    min[i][direc] = std::min(min[i][direc], value);  
                }
                
            }
            
            double determinant = x11 * x22 - x12 * x21;
            
            if (determinant == 0) {
                process.print("error! gradient can not be solve!\n"); //分母为零，梯度方程无解
            }
                        
            //解方程求梯度
            //gradientx = (y1 * x22 - y2 * x12) / determinant;
            //gradienty = (x11 * y2 - x21 * y1) / determinant;
            

            if (type == DistributionType::G) {  
                fluids[i].g_gradient_x[direc] = gradientx;  
                fluids[i].g_gradient_y[direc] = gradienty;  
            } else if (type == DistributionType::H) {  
                fluids[i].h_gradient_x[direc] = gradientx;  
                fluids[i].h_gradient_y[direc] = gradienty;  
            } 
            
        }
        
    }

    //delta_ij
    double*** delta=Allocate3D<double>(block.TotalCells,5,block.q_local);

    
    //限制器求值
    double* fai=Allocate1D<double>(block.TotalCells);
    for(int i=0;i<block.TotalCells;i++)
    {
        //赋初始值
        fai[i]=1.0;

        //对不同边界处理
        for(int adj=1;adj<=fluids[i].edge_num;adj++)
        {
            for(int direc=0;direc<block.q_local;direc++)
            {
                double x_correct=0.0,y_correct=0.0;//修正后的值
                x_correct=(fluids[i].adj_line_x[adj]-fluids[i].x)-fluids[i].ds*velq[direc].vx;
                y_correct=(fluids[i].adj_line_y[adj]-fluids[i].y)-fluids[i].ds*velq[direc].vy;


                // 待检查的点在内部/边界上,按内部点处理
                if ((fluids[i].nx[adj]*velq[direc].vx+fluids[i].ny[adj]*velq[direc].vy)>=0.0 || fluids[i].adj[adj]==0)
                {
                    if (type == DistributionType::G)
                    {
                        fluids[i].Gapb[adj][direc]= fluids[i].gap[direc] +  (x_correct * fluids[i].g_gradient_x[direc] + y_correct * fluids[i].g_gradient_y[direc]);

                        //Venkatakrishnan limiter
                        delta[i][adj][direc]=fluids[i].Gapb[adj][direc]-fluids[i].gap[direc];

                        double fai_=Venkatakrishnan(max[i][direc]-fluids[i].gap[direc],min[i][direc]-fluids[i].gap[direc],delta[i][adj][direc],sqrt(fluids[i].area),block.K_limit);
                        
                        fai[i]=std::min(fai_,fai[i]);
                        
                    }

                    if (type == DistributionType::H)
                    {
                        
                        fluids[i].Hapb[adj][direc]= fluids[i].hap[direc] +  (x_correct * fluids[i].h_gradient_x[direc] + y_correct * fluids[i].h_gradient_y[direc]);

                        //Venkatakrishnan limiter
                        delta[i][adj][direc]=fluids[i].Hapb[adj][direc]-fluids[i].hap[direc];

                        double fai_=Venkatakrishnan(max[i][direc]-fluids[i].hap[direc],min[i][direc]-fluids[i].hap[direc],delta[i][adj][direc],sqrt(fluids[i].area),block.K_limit);
                        
                        fai[i]=std::min(fai_,fai[i]);
                    }
                }
            }
        }
        
    }

    MPI_Allreduce(MPI_IN_PLACE, &fai[0],block.TotalCells, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    //利用梯度插值
    for(int i=0;i<block.TotalCells;i++)
    {
        //对不同边界处理
        for(int adj=1;adj<=fluids[i].edge_num;adj++)
        {
            for(int direc=0;direc<block.q_local;direc++)
            {
                double x_correct=0.0,y_correct=0.0;//修正后的值
                x_correct=(fluids[i].adj_line_x[adj]-fluids[i].x)-fluids[i].ds*velq[direc].vx;
                y_correct=(fluids[i].adj_line_y[adj]-fluids[i].y)-fluids[i].ds*velq[direc].vy;

                //printf("%lf %lf\n",x_correct,y_correct);

                
                // 待检查的点在内部/边界上,按内部点处理
                if ((fluids[i].nx[adj]*velq[direc].vx+fluids[i].ny[adj]*velq[direc].vy)>=0.0 || fluids[i].adj[adj]==0)
                {
                    if (type == DistributionType::G)
                    {
                        fluids[i].Gapb[adj][direc]= fluids[i].gap[direc] +  fai[i]*(x_correct * fluids[i].g_gradient_x[direc] + y_correct * fluids[i].g_gradient_y[direc]);
                        //fluids[i].Gapb[adj][direc]= fluids[i].gap[direc] +  (x_correct * fluids[i].g_gradient_x[direc] + y_correct * fluids[i].g_gradient_y[direc]);
                    }
                    if (type == DistributionType::H)
                    {
                        fluids[i].Hapb[adj][direc]= fluids[i].hap[direc] +  fai[i]*(x_correct * fluids[i].h_gradient_x[direc] + y_correct * fluids[i].h_gradient_y[direc]);
                        //fluids[i].Hapb[adj][direc]= fluids[i].hap[direc] +  (x_correct * fluids[i].h_gradient_x[direc] + y_correct * fluids[i].h_gradient_y[direc]);
                    }
                } 
            }
        }
        
    }
    
    Delte3D(delta,block.TotalCells,5);
    Delte2D(max,block.TotalCells);
    Delte2D(min,block.TotalCells);
    Delte1D(fai);

    
    //相邻控制体的辅助分布函数合并
    for(int i=0;i<block.TotalCells;i++)
    {
        //对不同的边界进行处理
        for(int adj=1;adj<=fluids[i].edge_num;adj++)
        {

            for(int direc=0;direc<block.q_local;direc++)
            {
                if (type == DistributionType::G)
                {
                    fluids[i].gapb[adj][direc]=0.0;
                    
                    //如果自己的界面处fb已求，则使用本体积内计算的值
                    if((fluids[i].nx[adj]*velq[direc].vx+fluids[i].ny[adj]*velq[direc].vy)>=0.0 || fluids[i].adj[adj]==0)
                    {
                        fluids[i].gapb[adj][direc]=fluids[i].Gapb[adj][direc];

                    }
                    //如果自己的界面处fb未求，则使用相邻界面同方向的fb
                    else
                    {
                        fluids[i].gapb[adj][direc]=fluids[fluids[i].adj[adj]].Gapb[fluids[i].adj_xb[adj]][direc];
                    }
                }
                if (type == DistributionType::H)
                {
                    fluids[i].hapb[adj][direc]=0.0;
                    //如果自己的界面处fb已求，则使用本体积内计算的值
                    if((fluids[i].nx[adj]*velq[direc].vx+fluids[i].ny[adj]*velq[direc].vy)>=0.0 || fluids[i].adj[adj]==0)
                    {
                        fluids[i].hapb[adj][direc]=fluids[i].Hapb[adj][direc];

                    }
                    //如果自己的界面处fb未求，则使用相邻界面同方向的fb
                    else
                    {
                        fluids[i].hapb[adj][direc]=fluids[fluids[i].adj[adj]].Hapb[fluids[i].adj_xb[adj]][direc];
                    }
                }

            }
        }
    }

}



//计算界面宏观量
void BoundaryMacroscopicVariables(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq)
{
    //界面上的宏观量
    for(auto& fluid : fluids)
    {
        //密度
        for (int adj = 1; adj <= fluid.edge_num; adj++)
        {
            //置零
            fluid.rhob[adj] = 0.0;
            fluid.vxb[adj] = 0.0;
            fluid.vyb[adj] = 0.0;
            fluid.eb[adj]=0.0;
            fluid.evb[adj]=0.0;
            fluid.tempb[adj]=0.0;
            fluid.qxb_flux[adj]=0.0;
            fluid.qyb_flux[adj]=0.0;
            fluid.taub[adj]=0.0;


            //计算密度
            for (int direc = 0; direc < block.q_local; direc++) {
                fluid.rhob[adj] += fluid.gapb[adj][direc];
            }
            //printf("%lf--%d\n",rhob[i][adj],pid);

            MPI_Allreduce(MPI_IN_PLACE, &fluid.rhob[adj], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            //process.printnum(fluid.rhob[adj]);
        }


    }
    
    for(auto& fluid : fluids)
    {
        //速度
        for (int adj = 1; adj <= fluid.edge_num; adj++) {
            for (int direc=0; direc < block.q_local; direc++) {
                
                fluid.vxb[adj] += velq[direc].vx * fluid.gapb[adj][direc];
                fluid.vyb[adj] += velq[direc].vy * fluid.gapb[adj][direc];
            }

            fluid.vxb[adj] /= fluid.rhob[adj];
            fluid.vyb[adj] /= fluid.rhob[adj];

            MPI_Allreduce(MPI_IN_PLACE, &fluid.vxb[adj], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &fluid.vyb[adj], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            //process.printnum(fluid.vyb[adj]);
        }   
    }
    
    for(auto& fluid : fluids)
    {
        for (int adj = 1; adj <= fluid.edge_num; adj++) 
        {
            //能量
            for (int direc = 0; direc < block.q_local; direc++)
            {
                fluid.eb[adj]+=0.5*((pow(velq[direc].vx,2)+pow(velq[direc].vy,2))*fluid.gapb[adj][direc]+fluid.hapb[adj][direc]);
            }
            MPI_Allreduce(MPI_IN_PLACE, &fluid.eb[adj], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            fluid.eb[adj] /= fluid.rhob[adj];
        }
    }
    
    
    for(auto& fluid : fluids)
    {
        for (int adj = 1; adj <= fluid.edge_num; adj++) 
        {    

            fluid.evb[adj]=fluid.eb[adj]-0.5*(pow(fluid.vxb[adj],2)+pow(fluid.vyb[adj],2));
            //温度=内能÷定压比热
            fluid.tempb[adj]=fluid.evb[adj]/(block.R/(block.gama-1.0));
            //process.printnum(fluid.tempb[adj]);

            fluid.qxb_flux[adj]=0.0;
            fluid.qyb_flux[adj]=0.0;
            //热流
            for(int direc=0;direc<block.q_local;direc++)
            {
                double cx=0.0,cy=0.0;
                cx=velq[direc].vx-fluid.vxb[adj];
                cy=velq[direc].vy-fluid.vyb[adj];

                //辅助分布的热流
                fluid.qxb_flux[adj]+=0.5*cx*((cx*cx+cy*cy)*fluid.gapb[adj][direc]+fluid.hapb[adj][direc]);
                fluid.qyb_flux[adj]+=0.5*cy*((cx*cx+cy*cy)*fluid.gapb[adj][direc]+fluid.hapb[adj][direc]);
            }
            //辅助分布的热流->实际的界面热流
            fluid.qxb_flux[adj]=2.0*fluid.tau/(2.0*fluid.tau+fluid.ds*block.Pr)*fluid.qxb_flux[adj];
            fluid.qyb_flux[adj]=2.0*fluid.tau/(2.0*fluid.tau+fluid.ds*block.Pr)*fluid.qyb_flux[adj];

            MPI_Allreduce(MPI_IN_PLACE, &fluid.qxb_flux[adj], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &fluid.qyb_flux[adj], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            //process.printnum(fluid.qyb_flux[adj]);
        }
    }

    //计算界面的taub
    if(block.tau_type==TAU_TYPE::NS)
    {
    for(auto& fluid : fluids)
    {
        for (int adj = 1; adj <= fluid.edge_num; adj++) 
        {
            fluid.taub[adj]=(block.miu_init*pow(fluid.tempb[adj]/block.temp_init,block.omega))/(fluid.rhob[adj]*block.R*fluid.tempb[adj]);
        }
    }
    }
    
}


//界面平衡态与dt/2时刻界面真实分布
void EquilibriumAndNext(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq, DistributionType type) 
{
    
    for(auto& fluid : fluids)
    {
        for (int adj = 1; adj <= fluid.edge_num; adj++)
        {
            for(int direc=0;direc<block.q_local;direc++)
            {
                if (type == DistributionType::G) 
                {  
                    double bs=0.0;
                    //界面平衡态
                    bs= bCalcGs(&block,&fluid,&velq[direc],adj);
                    fluid.gppb[adj][direc]=2.0*fluid.taub[adj]/(2.0*fluid.taub[adj]+fluid.ds)*fluid.gapb[adj][direc]+fluid.ds/(2.0*fluid.taub[adj]+fluid.ds)*bs;
                } 
                else if (type == DistributionType::H) 
                {  
                    double bs=0.0;
                    //界面平衡态
                    bs= bCalcHs(&block,&fluid,&velq[direc],adj);
                    fluid.hppb[adj][direc]=2.0*fluid.taub[adj]/(2.0*fluid.taub[adj]+fluid.ds)*fluid.hapb[adj][direc]+fluid.ds/(2.0*fluid.taub[adj]+fluid.ds)*bs;
                } 
                
            }
        }
    }
}




//边界条件
void BoundaryCondition(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq) 
{
    double error_bound=1e-30;
    
    for(auto& fluid : fluids)
    {   
        for (int adj = 1; adj <= fluid.edge_num; adj++)
        {
            //Wall 固定边界   type1(fluids[i].boundary[adj][1]==1)
            if(fluid.boundary[adj][1]==1)
            {
                fluid.rhob[adj]=1.0;
                fluid.vxb[adj]=0.0;
                fluid.vyb[adj]=0.0;
                //fluid.tempb[adj]=block.temp_init;
                fluid.tempb[adj]=block.temp_init;
                fluid.evb[adj]=0.0;
                fluid.eb[adj]=0.0;
                fluid.qxb_flux[adj]=0.0;
                fluid.qyb_flux[adj]=0.0;
                
                //process.printnum(fluid.rhob[adj]);
                
                //计算边界上的密度
                double result_out=0.0,result_in=0.0;
                for(int direc=0;direc<block.q_local;direc++)
                {
                    double calc_condition=0.0;
                    calc_condition=(velq[direc].vx*fluid.nx[adj]+velq[direc].vy*fluid.ny[adj]);

                    //朝外
                    if(calc_condition>error_bound)
                    {
                        result_out+=fluid.gppb[adj][direc]*calc_condition;
                    }
                    //朝内
                    if(calc_condition<-error_bound)
                    {
                        result_in+= bCalcGs(&block,&fluid,&velq[direc],adj)*calc_condition;
                        
                    }
                }
                
                MPI_Allreduce(MPI_IN_PLACE, &result_out, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, &result_in, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
                
                //密度
                fluid.rhob[adj]=-result_out/result_in;

                
                
                //process.printnum(result_out);
                //process.printnum(result_in);
                //process.print("#######");
                
                fluid.evb[adj]=block.R*fluid.tempb[adj]/(block.gama-1.0);
                fluid.eb[adj]=0.5*(pow(fluid.vxb[adj],2)+pow(fluid.vyb[adj],2))+fluid.evb[adj];

                //根据得到的密度，补齐朝内的分布
                for(int direc=0;direc<block.q_local;direc++)
                {
                    double calc_condition=0.0;
                    //解析该方向上的单位向量(流出为正)
                    calc_condition=(velq[direc].vx*fluid.nx[adj]+velq[direc].vy*fluid.ny[adj]);

                    //朝内
                    if(calc_condition<-error_bound)
                    {
                        fluid.gppb[adj][direc]= bCalcGs(&block,&fluid,&velq[direc],adj);
                        fluid.hppb[adj][direc]= bCalcHs(&block,&fluid,&velq[direc],adj);
                    }

                }
        

            }

            //MoveWall  移动边界   type2(fluids[i].boundary[adj][2]==1)
            if(fluid.boundary[adj][2]==1)
            {
                fluid.rhob[adj]=block.rho_init;
                fluid.vxb[adj]=block.Ma*sqrt(block.gama*block.R*block.temp_init);
                //fluid.vxb[adj]=0.0;
                fluid.vyb[adj]=0.0;
                fluid.tempb[adj]=block.temp_init;
                fluid.evb[adj]=0.0;
                fluid.eb[adj]=0.0;
                fluid.qxb_flux[adj]=0.0;
                fluid.qyb_flux[adj]=0.0;

                //计算边界上的密度
                double result_out=0.0,result_in=0.0;
                for(int direc=0;direc<block.q_local;direc++)
                {
                    double calc_condition=0.0;
                    calc_condition=(velq[direc].vx*fluid.nx[adj]+velq[direc].vy*fluid.ny[adj]);

                    //朝外
                    if(calc_condition>error_bound)
                    {
                        result_out+=fluid.gppb[adj][direc]*calc_condition;
                    }
                    //朝内
                    if(calc_condition<-error_bound)
                    {
                        result_in+= bCalcGs(&block,&fluid,&velq[direc],adj)*calc_condition;
                        
                        //process.printnum(bCalcGs(&block,&fluid,&velq[direc],adj));
                        /*
                        if(process.pid==0)
                        {
                            std::cout<<velq[direc].vx<<"   "<<velq[direc].vy<<"   "<<adj<<"   "<<direc<<std::endl;
                        }
                        */


                    }
                }
                
                MPI_Allreduce(MPI_IN_PLACE, &result_out, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Allreduce(MPI_IN_PLACE, &result_in, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
                
                //密度
                fluid.rhob[adj]=-result_out/result_in;
                //process.printnum(fluid.rhob[adj]);

                fluid.evb[adj]=block.R*fluid.tempb[adj]/(block.gama-1.0);
                fluid.eb[adj]=0.5*(pow(fluid.vxb[adj],2)+pow(fluid.vyb[adj],2))+fluid.evb[adj];

                //根据得到的密度，补齐朝内的分布
                for(int direc=0;direc<block.q_local;direc++)
                {
                    double calc_condition=0.0;
                    //解析该方向上的单位向量(流出为正)
                    calc_condition=(velq[direc].vx*fluid.nx[adj]+velq[direc].vy*fluid.ny[adj]);

                    //朝内
                    if(calc_condition<=0.0)
                    {
                        fluid.gppb[adj][direc]= bCalcGs(&block,&fluid,&velq[direc],adj);
                        fluid.hppb[adj][direc]= bCalcHs(&block,&fluid,&velq[direc],adj);
                    }

                }

            }

            //inLet  入口   type3(fluids[i].boundary[adj][3]==1)
            if(fluid.boundary[adj][3]==1)
            {
                fluid.rhob[adj]=block.rho_init;
                fluid.vxb[adj]=block.Ma_inlet*sqrt(block.gama*block.R*block.temp_init);
                fluid.vyb[adj]=0.0;
                fluid.tempb[adj]=block.temp_init;
                fluid.evb[adj]=block.R*block.temp_init/(block.alpha-1.0);
                fluid.eb[adj]=block.R*block.temp_init/(block.alpha-1.0)+0.5*(pow(fluid.vxb[adj],2)+pow(fluid.vyb[adj],2));
                fluid.qxb_flux[adj]=0.0;
                fluid.qyb_flux[adj]=0.0;


                
                //补齐朝内的分布
                for(int direc=0;direc<block.q_local;direc++)
                {
                    double calc_condition=0.0;
                    //解析该方向上的单位向量(流出为正)
                    calc_condition=(velq[direc].vx*fluid.nx[adj]+velq[direc].vy*fluid.ny[adj]);

                    //朝内
                    if(calc_condition<=0.0)
                    {
                        fluid.gppb[adj][direc]= bCalcGs(&block,&fluid,&velq[direc],adj);
                        fluid.hppb[adj][direc]= bCalcHs(&block,&fluid,&velq[direc],adj);
                    }

                }

            }
            //OutLet  出口   type4(fluids[i].boundary[adj][4]==1)
            if(fluid.boundary[adj][4]==1)
            {
                //不处理

            }

            //Sym 对称 type5(fluids[i].boundary[adj][5])
            if(fluid.boundary[adj][5]==1)
            {
                /*
                for(int direc=0;direc<block.q_local;direc++)
                {
                    fluid.gppb[(adj+2)%4][direc]= bCalcGs(&block,&fluid,&velq[direc],adj);
                    fluid.hppb[(adj+2)%4][direc]= bCalcHs(&block,&fluid,&velq[direc],adj);
                }
                */
                
            }
        }
    }

    
}


//微观通量
void CalcFmeso(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq)
{
    //控制体的微观通量
    for(auto& fluid : fluids)
    {        
        //控制体的微观通量
        for(int direc=0;direc<block.q_local;direc++)
        {
            
            //置0
            fluid.Gmeso[direc]=0.0;
            fluid.Hmeso[direc]=0.0;

            //为防止出现±90°的情况，采用向量相乘的方法
            //对于不同的边界
            for(int adj=1;adj<= fluid.edge_num;adj++)
            {
                fluid.Gmeso[direc]+=(fluid.nx[adj]*velq[direc].vx+fluid.ny[adj]*velq[direc].vy)*fluid.gppb[adj][direc]*fluid.edgelength[adj];
                fluid.Hmeso[direc]+=(fluid.nx[adj]*velq[direc].vx+fluid.ny[adj]*velq[direc].vy)*fluid.hppb[adj][direc]*fluid.edgelength[adj];
            }
        }
    }
}


//宏观通量
void CalcFmacro(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq)
{
    for(int i=0;i<block.TotalCells;i++)
    {   
        //置0
        fluids[i].Fmacro_rho= 0.0;
        fluids[i].Fmacro_vx = 0.0;
        fluids[i].Fmacro_vy = 0.0;
        fluids[i].Fmacro_e=0.0;
        //rho
        //vx,vy,energy
        for(int adj=1;adj <= fluids[i].edge_num;adj++)
        {
            for(int direc=0;direc<block.q_local;direc++)
            {                
                fluids[i].Fmacro_rho+=(fluids[i].nx[adj]*velq[direc].vx+fluids[i].ny[adj]*velq[direc].vy)*fluids[i].gppb[adj][direc]*fluids[i].edgelength[adj];
                fluids[i].Fmacro_vx+=velq[direc].vx*(fluids[i].nx[adj]*velq[direc].vx+fluids[i].ny[adj]*velq[direc].vy)*fluids[i].gppb[adj][direc]*fluids[i].edgelength[adj];
                fluids[i].Fmacro_vy+=velq[direc].vy*(fluids[i].nx[adj]*velq[direc].vx+fluids[i].ny[adj]*velq[direc].vy)*fluids[i].gppb[adj][direc]*fluids[i].edgelength[adj];
                fluids[i].Fmacro_e+=0.5*(fluids[i].nx[adj]*velq[direc].vx+fluids[i].ny[adj]*velq[direc].vy)*((pow(velq[direc].vx,2)+pow(velq[direc].vy,2))*fluids[i].gppb[adj][direc]+fluids[i].hppb[adj][direc])*fluids[i].edgelength[adj];
     
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, &fluids[i].Fmacro_rho, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &fluids[i].Fmacro_vx, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &fluids[i].Fmacro_vy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &fluids[i].Fmacro_e, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
    }
}


//格心宏观量更新
void UpdateMacroscopicVariables(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq)
{
    //计算格心的宏观量（n+1）(宏观量更新)
    double* rho_next=Allocate1D<double>(block.TotalCells);
    
    for(int i=0;i<block.TotalCells;i++)
    {

        //密度
        rho_next[i]=fluids[i].rho-fluids[i].dt*fluids[i].Fmacro_rho/fluids[i].area;

        fluids[i].vx_history[block.step%1000]=fluids[i].vx;
        fluids[i].vy_history[block.step%1000]=fluids[i].vy;
        
        //速度
        fluids[i].vx = (fluids[i].vx * fluids[i].rho - fluids[i].dt * fluids[i].Fmacro_vx/ fluids[i].area) / (rho_next[i]);
        fluids[i].vy = (fluids[i].vy * fluids[i].rho - fluids[i].dt * fluids[i].Fmacro_vy/ fluids[i].area) / (rho_next[i]);
        //printf("e=%lf\n",e[i]);
        //单位质量的总能量
        fluids[i].e= (fluids[i].e * fluids[i].rho - fluids[i].dt * fluids[i].Fmacro_e / fluids[i].area) / (rho_next[i]);

        fluids[i].rho = rho_next[i];

        fluids[i].ev=fluids[i].e-0.5*(pow(fluids[i].vx,2)+pow(fluids[i].vy,2));
        //温度
        fluids[i].temp=fluids[i].ev/(block.R/(block.gama-1.0));

    }
    
    
    //找到temp_max;
    double temp_max = fluids[0].temp;
    // 遍历数组，找到最大值
    for (int tt = 1; tt <block.TotalCells; tt++) {
        if (fluids[tt].temp > temp_max) {
            temp_max = fluids[tt].temp;
        }
    }

    //更新tau,pressure,miu,dt,ds
    for(auto& fluid : fluids)
    {
        fluid.pressure=fluid.rho*block.R*fluid.temp;
        //粘度与温度，omega有关
        fluid.miu=block.miu_init*pow(fluid.temp/block.temp_init,block.omega);
        fluid._tau=fluid.tau;  //上一歩时间的tau
        if(block.tau_type==TAU_TYPE::NS)
        {
            fluid.tau=fluid.miu/fluid.pressure;
        }
        fluid.dt=block.alpha*block.min_length/(block.max_vel+block.Ma*sqrt(block.gama*block.R*temp_max));
        fluid.ds=fluid.dt/2.0;
    }
    
    Delte1D(rho_next);
}


//格心微观分布更新
void UpdataDistribution(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq)
{
    for(auto& fluid : fluids)
    {
        for(int direc=0;direc<block.q_local;direc++)
        {

            fluid.g[direc]=(1.0/(1.0+fluid.dt/(2.0*fluid.tau)))*(fluid.g[direc]-fluid.dt/fluid.area*fluid.Gmeso[direc]+fluid.dt/2.0*(CalcGs(&block,&fluid,&velq[direc])/fluid.tau+(fluid.gs[direc]-fluid.g[direc])/(fluid._tau)));

            fluid.h[direc]=(1.0/(1.0+fluid.dt/(2.0*fluid.tau)))*(fluid.h[direc]-fluid.dt/fluid.area*fluid.Hmeso[direc]+fluid.dt/2.0*(CalcHs(&block,&fluid,&velq[direc])/fluid.tau+(fluid.hs[direc]-fluid.h[direc])/(fluid._tau)));

        }
    }

    //更新热流（n+1时刻）
    for(auto& fluid : fluids)
    {
        //热流
        fluid.qx_flux=0.0;
        fluid.qy_flux=0.0;

        for(int direc=0;direc<block.q_local;direc++)
        {
            double cx=0.0,cy=0.0;
            cx=velq[direc].vx-fluid.vx;
            cy=velq[direc].vy-fluid.vy;
            fluid.qx_flux+=0.5*cx*((cx*cx+cy*cy)*fluid.g[direc]+fluid.h[direc]);
            fluid.qy_flux+=0.5*cy*((cx*cx+cy*cy)*fluid.g[direc]+fluid.h[direc]);
        }

        MPI_Allreduce(MPI_IN_PLACE, &fluid.qx_flux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &fluid.qy_flux, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
}

//判断收敛
bool SteadyJudge(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq)
{
    bool stop_=false;

    block.time_calculation+=fluids[0].dt;

    double L2_error=0.0,s1=0.0,s2=0.0;
    for(int i=0;i<block.TotalCells;i++)
    {
        //blow up
        if (isnan(fluids[i].vx) || isnan(fluids[i].vy) || isnan(fluids[i].rho) || isnan(fluids[i].pressure) || isnan(fluids[i].temp)) 
        {
            if(process.pid==0)
            {
                printf("blow up point: (%f, %f)\n", fluids[i].x, fluids[i].y);
            }
            
            stop_=true;
        }
        
        if(block.step<=1000)
        {
            s1+=(sqrt(pow((fluids[i].vx-0.0),2)+pow((fluids[i].vy-0.0),2)));
            s2+=(sqrt(pow(fluids[i].vx,2)+pow(fluids[i].vy,2)));   
        }
        else
        {
            s1+=(sqrt(pow((fluids[i].vx-fluids[i].vx_history[block.step%1000]),2)+pow((fluids[i].vy-fluids[i].vy_history[block.step%1000]),2)));
            s2+=(sqrt(pow(fluids[i].vx,2)+pow(fluids[i].vy,2)));
        }
    }

    L2_error=s1/s2;


    //达到收敛条件或达到设定时长
    if(L2_error<=block.error_limit || ((block.time_calculation>block.timemax) && block.endtime_type==ENDTIME_TYPE::Set))
    {
        stop_=true;
    }
    
    
    double SimulTime = MPI_Wtime()-block.BeginSimulTime;

    if (process.pid == 0) std::cout << "Step=" << block.step << "  L2_error:"<<L2_error<<"  Time_simulated:"<<block.time_calculation<<"  Time_consumption:"<<SimulTime<<std::endl;

    return stop_;
}
