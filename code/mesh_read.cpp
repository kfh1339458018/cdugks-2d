#include "mesh_read.h"
#include "cdugks_2d.h"




//兼容linux
string add_mesh_directory_modify_for_linux()
{
#if defined(_WIN32)
	return "";
#else
	return "../";
#endif
}



void Read_mesh_neu_2d(MPIProcess& process,std::vector<Fluid2d>& fluids,Block2d& block,std::vector<Point2d>& points)
{



    //===read grid file   
    ifstream in;
    std::istringstream line_data;
    std::string str;
    std::vector<std::string> TempVec;

    in.open(block.grid_path, ifstream::binary);


    process.print("grid mesh path:  "+block.grid_path);
    
    //如果找不到文件
    if(!in.is_open())
    {   
        process.print("cannot find input geometry grid file(.neu)");
    }
    process.print("reading geometry mesh......");
    

    std::string line;
    while (std::getline(in, line)) { // 逐行读取文件


        //find number of cells,nodes
        if (line.find("NUMNP") != std::string::npos)
        {
            getline(in, line);
            TempVec.resize(0);
            line_data = std::istringstream(line);  //read from line

            while (line_data >> str)
            {
                TempVec.push_back(str);
            }
            block.TotalPoints = stoi(TempVec[0]);  //stoi()字符串转整数
            block.TotalCells = stoi(TempVec[1]);
            block.BoundaryType = stoi(TempVec[3]);
        }

        //extract node information
        if(line.find("NODAL COORDINATES")!=std::string::npos)
        {
            Point2d temppoint;
            //point
            for(int count=0;count<block.TotalPoints;count++)
            {
                getline(in,line);
                TempVec.resize(0);
                line_data=std::istringstream(line);
                while(line_data>>str)
                {
                    TempVec.push_back(str);
                }
                temppoint.x=stod(TempVec[1]);
                temppoint.y=stod(TempVec[2]);
                points.push_back(temppoint);
                //process.printnum(points[count].y);
                //process.printnum(points[count].x);
            }
        }

        //extract elements/cells information
        if(line.find("ELEMENTS/CELLS")!=std::string::npos)
        {
            Fluid2d tempfluid;
            //elements/cells
            for(int count=0;count<block.TotalCells;count++)
            {
                getline(in,line);
                TempVec.resize(0);
                line_data=std::istringstream(line);
                while(line_data>>str)
                {
                    TempVec.push_back(str);
                }
                if(stoi(TempVec[1])==2) //quadrilateral
                {
                    tempfluid.edge_num=4;
                    tempfluid.v[1]=(stoi(TempVec[3])-1);//neu文件从1开始计数，转成从0开始计数
                    tempfluid.v[2]=(stoi(TempVec[4])-1);
                    tempfluid.v[3]=(stoi(TempVec[5])-1);
                    tempfluid.v[4]=(stoi(TempVec[6])-1);
                }
                if(stoi(TempVec[1])==3) //triangle
                {
                    tempfluid.edge_num=4;
                    tempfluid.v[1]=(stoi(TempVec[3])-1);
                    tempfluid.v[2]=(stoi(TempVec[4])-1);
                    tempfluid.v[3]=(stoi(TempVec[5])-1);
                    tempfluid.v[4]=(stoi(TempVec[5])-1); //v3=v4
                }
                fluids.push_back(tempfluid);
            }  
        }

        //extract boundary conditon
        if(line.find("BOUNDARY CONDITIONS")!=std::string::npos)
        {
            getline(in,line);  //再读一行
            TempVec.resize(0);
            line_data=std::istringstream(line);  //read from "line"
            while(line_data>>str)
            {
                TempVec.push_back(str);
            }
            //process.print(line);
            //Wall 固定边界   type1(fluids[i].boundary[adj][1])
            if(TempVec[0]=="Wall")
            {
                for(int count=0;;count++)
                {
                    getline(in,line);
                    if(line.find("ENDOFSECTION")!=std::string::npos)
                    {
                        break;  //如果读到该边界条件的结尾，退出循环
                    }
                    TempVec.resize(0);
                    line_data=std::istringstream(line);
                    while(line_data>>str)
                    {
                        TempVec.push_back(str);
                    }
                    //录入边界条件
                    fluids[stoi(TempVec[0])-1].boundary[stoi(TempVec[2])][1]=1;
                }
            }
            //MoveWall  移动边界   type2(fluids[i].boundary[adj][2])
            if(TempVec[0]=="MoveWall")
            {
                for(int count=0;;count++)
                {
                    getline(in,line);
                    if(line.find("ENDOFSECTION")!=std::string::npos)
                    {
                        break;  //如果读到该边界条件的结尾，退出循环
                    }
                    TempVec.resize(0);
                    line_data=std::istringstream(line);
                    while(line_data>>str)
                    {
                        TempVec.push_back(str);
                    }
                    //录入边界条件
                    fluids[stoi(TempVec[0])-1].boundary[stoi(TempVec[2])][2]=1;
                } 
            }
            //InLet  入口   type3(fluids[i].boundary[adj][3])
            if(TempVec[0]=="InLet")
            {
                for(int count=0;;count++)
                {
                    getline(in,line);
                    if(line.find("ENDOFSECTION")!=std::string::npos)
                    {
                        break;  //如果读到该边界条件的结尾，退出循环
                    }
                    TempVec.resize(0);
                    line_data=std::istringstream(line);
                    while(line_data>>str)
                    {
                        TempVec.push_back(str);
                    }
                    //录入边界条件
                    fluids[stoi(TempVec[0])-1].boundary[stoi(TempVec[2])][3]=1;
                } 
            }
            //OutLet  出口   type4(fluids[i].boundary[adj][4])
            if(TempVec[0]=="OutLet")
            {
                for(int count=0;;count++)
                {
                    getline(in,line);
                    if(line.find("ENDOFSECTION")!=std::string::npos)
                    {
                        break;  //如果读到该边界条件的结尾，退出循环
                    }
                    TempVec.resize(0);
                    line_data=std::istringstream(line);
                    while(line_data>>str)
                    {
                        TempVec.push_back(str);
                    }
                    //录入边界条件
                    fluids[stoi(TempVec[0])-1].boundary[stoi(TempVec[2])][4]=1;
                } 
            }

            //Sym 对称 type5(fluids[i].boundary[adj][5])
            if(TempVec[0]=="Sym")
            {
                for(int count=0;;count++)
                {
                    getline(in,line);
                    if(line.find("ENDOFSECTION")!=std::string::npos)
                    {
                        break;  //如果读到该边界条件的结尾，退出循环
                    }
                    TempVec.resize(0);
                    line_data=std::istringstream(line);
                    while(line_data>>str)
                    {
                        TempVec.push_back(str);
                    }
                    //录入边界条件
                    fluids[stoi(TempVec[0])-1].boundary[stoi(TempVec[2])][5]=1;
                } 
            }
        }      
    }
    process.print("mesh reading completed.");
    in.close(); // 关闭文件 

    //===analyze geometric relationship
    //每个控制体的边长
    for(int i=0;i<block.TotalCells;i++)
    {
        //录入每个点的坐标
        for(int r=1;r<=fluids[i].edge_num;r++)
        {
            fluids[i].xv[r]=points[fluids[i].v[r]].x;
            fluids[i].yv[r]=points[fluids[i].v[r]].y;
        }

        //计算边长
        for(int r=1;r<=fluids[i].edge_num;r++)
        {
            int r_next=(r%fluids[i].edge_num+1);
            fluids[i].edgelength[r]=sqrt(pow(fluids[i].xv[r]-fluids[i].xv[r_next],2)+pow(fluids[i].yv[r]-fluids[i].yv[r_next],2));
        }
        
    }

    //找到最短的边长
    block.min_length=fluids[0].edgelength[1];  //假设第一个控制体的第一条边最短
    for(int i=0;i<block.TotalCells;i++)
    {
        for(int r=1;r<=fluids[i].edge_num;r++)
        {
            if(fluids[i].edgelength[r]<block.min_length)
            {
                block.min_length=fluids[i].edgelength[r];
            }
        }
    }
    process.printnum(block.min_length);
    
    //寻找相邻控制体编号
    for(int i=0;i<block.TotalCells;i++)
    {
        for(int j=0;j<block.TotalCells;j++)
        {
            if(j==i)continue;  //跳过当前的控制体
            
            //检查控制体j的顶点是否与控制体i的顶点重合
            for(int k=1;k<=fluids[j].edge_num;k++)
            {
                //v1 v2 轮换格点
                int v1=fluids[j].v[k];
                int v2=fluids[j].v[k%fluids[j].edge_num+1];

                for(int adj=1;adj<=fluids[i].edge_num;adj++)
                {
                    //如果两个边的两个顶点重合
                    if((fluids[i].v[adj]==v1||fluids[i].v[adj]==v2) && (fluids[i].v[adj%fluids[i].edge_num+1]==v1||fluids[i].v[adj%fluids[i].edge_num+1]==v2))
                    {
                        fluids[i].adj[adj]=j;  //相邻的控制体j
                        fluids[i].adj_xb[adj]=k;  //相邻控制体，对应接触边的编号
                    }
                }
            }
            
        }
    }

    //格心坐标，控制体面积，界面中点坐标，界面上垂点的坐标，单位向量
    for(int i=0;i<block.TotalCells;i++)
    {
        //triangle
        if(fluids[i].edge_num==3)
        {
            //格心坐标
            fluids[i].x=(fluids[i].xv[1]+fluids[i].xv[2]+fluids[i].xv[3])/3.0;
            fluids[i].y=(fluids[i].yv[1]+fluids[i].yv[2]+fluids[i].yv[3])/3.0;

            double a = sqrt((fluids[i].xv[2] - fluids[i].xv[3]) * (fluids[i].xv[2] - fluids[i].xv[3]) + (fluids[i].yv[2] - fluids[i].yv[3]) * (fluids[i].yv[2] - fluids[i].yv[3]));
            double b = sqrt((fluids[i].xv[1] - fluids[i].xv[3]) * (fluids[i].xv[1] - fluids[i].xv[3]) + (fluids[i].yv[1] - fluids[i].yv[3]) * (fluids[i].yv[1] - fluids[i].yv[3]));
            double c = sqrt((fluids[i].xv[1] - fluids[i].xv[2]) * (fluids[i].xv[1] - fluids[i].xv[2]) + (fluids[i].yv[1] - fluids[i].yv[2]) * (fluids[i].yv[1] - fluids[i].yv[2]));

            // 计算三角形的半周长
            double s = (a + b + c) / 2.0;
            // 计算三角形的面积
            fluids[i].area = sqrt(s * (s - a) * (s - b) * (s - c));

        }

        //quadrilateral
        if(fluids[i].edge_num==4)
        {
            //格心坐标
            fluids[i].x=(fluids[i].xv[1]+fluids[i].xv[2]+fluids[i].xv[3]+fluids[i].xv[4])/4.0;
            fluids[i].y=(fluids[i].yv[1]+fluids[i].yv[2]+fluids[i].yv[3]+fluids[i].yv[4])/4.0;
            
            double cross1=(fluids[i].xv[1]*fluids[i].yv[2]-fluids[i].xv[2]*fluids[i].yv[1]);
            double cross2=(fluids[i].xv[2]*fluids[i].yv[3]-fluids[i].xv[3]*fluids[i].yv[2]);
            double cross3=(fluids[i].xv[3]*fluids[i].yv[4]-fluids[i].xv[4]*fluids[i].yv[3]);
            double cross4=(fluids[i].xv[4]*fluids[i].yv[1]-fluids[i].xv[1]*fluids[i].yv[4]);

            fluids[i].area=0.5*fabs(cross1 + cross2 + cross3 + cross4);
        }

        //界面中点坐标
        for(int adj=1;adj<=fluids[i].edge_num;adj++)
        {
            fluids[i].adj_line_x[adj]=(fluids[i].xv[adj]+fluids[i].xv[adj%fluids[i].edge_num+1])/2.0;
            fluids[i].adj_line_y[adj]=(fluids[i].yv[adj]+fluids[i].yv[adj%fluids[i].edge_num+1])/2.0;
            //process.printnum(fluids[i].adj_line_x[adj]);
        }

        //界面上垂点的坐标
        for(int adj=1;adj<=fluids[i].edge_num;adj++)
        {
            double right=1.0,left=0.0;
            double dx,dy,f_right,f_left;
            double tarx,tary,mid,f_mid;


            dx=fluids[i].xv[adj%fluids[i].edge_num+1]-fluids[i].xv[adj];
            dy=fluids[i].yv[adj%fluids[i].edge_num+1]-fluids[i].yv[adj];

            f_left = (fluids[i].xv[adj]- fluids[i].x-10.0*dx) * dx + (fluids[i].yv[adj]- fluids[i].y-10.0*dy) * dy;
            f_right = (fluids[i].xv[adj%fluids[i].edge_num+1]- fluids[i].x+10.0*dx) * dx + (fluids[i].yv[adj%fluids[i].edge_num+1] - fluids[i].y+10.0*dy) * dy;
            
            double epsilon=1e-9;
            if (f_left * f_right > 0) {
                // 如果左右边界都在同一侧（上穿或下穿），则无法找到根
                process.print("Both left and right boundaries are on the same side, cannot find a root.\n");
            }

            mid = (right + left) / 2.0;
            while (fabs(right - left) > epsilon) {
                // 计算tarx, tary
                tarx = fluids[i].xv[adj] + mid * dx;
                tary = fluids[i].yv[adj] + mid * dy;
                f_mid = (tarx - fluids[i].x) * dx + (tary - fluids[i].y) * dy;

                if (f_mid == 0) {
                    break;
                } else if (f_mid * f_left < 0) {
                    right = mid;
                } else {
                    left = mid;
                }
                mid = (right + left) / 2.0;
            }

            // 将边上切点的坐标写入结构体中
            fluids[i].adj_vert_x[adj] = fluids[i].xv[adj] + mid * dx;
            fluids[i].adj_vert_y[adj] = fluids[i].yv[adj] + mid * dy;
            //process.printnum(mid);
        }

        //垂直于界面的单位向量(朝外为正)
        for(int adj=1;adj<=fluids[i].edge_num;adj++)
        {
            fluids[i].nx[adj]=(fluids[i].adj_vert_x[adj]-fluids[i].x)/sqrt(pow(fluids[i].adj_vert_x[adj]-fluids[i].x,2)+pow(fluids[i].adj_vert_y[adj]-fluids[i].y,2));
            fluids[i].ny[adj]=(fluids[i].adj_vert_y[adj]-fluids[i].y)/sqrt(pow(fluids[i].adj_vert_x[adj]-fluids[i].x,2)+pow(fluids[i].adj_vert_y[adj]-fluids[i].y,2));

            //process.printnum(fluids[i].nx[adj]);
            //process.printnum(fluids[i].ny[adj]);
        }


    }

    process.print("mesh analysis completed.");
    
}


void Unstruct_velocity(MPIProcess& process,std::vector<Velocity2d>& vel,Block2d& block,std::vector<Point2d>& points)
{
    //===read velocity file   
    ifstream in;
    std::istringstream line_data;
    std::string str;
    std::vector<std::string> TempVec;

    in.open(block.vel_path, ifstream::binary);

    process.print("velocity mesh path:  "+block.vel_path);
    //如果找不到文件
    if(!in.is_open())
    {   
        process.print("cannot find input velocity grid file(.neu)");
    }
    process.print("reading velocity mesh......");
    
    std::string line;

    while (std::getline(in, line)) // 逐行读取文件
    { 
        
        //find number of cells,nodes
        if (line.find("NUMNP") != std::string::npos)
        {
            getline(in, line);
            TempVec.resize(0);
            line_data = std::istringstream(line);  //read from line

            while (line_data >> str)
            {
                TempVec.push_back(str);
            }
            block.q_node = stoi(TempVec[0]);  
            block.q = stoi(TempVec[1]);
        }

        //extract node information
        if(line.find("NODAL COORDINATES")!=std::string::npos)
        {
            Point2d temppoint;
            //point
            for(int count=0;count<block.q_node;count++)
            {
                getline(in,line);
                TempVec.resize(0);
                line_data=std::istringstream(line);
                while(line_data>>str)
                {
                    TempVec.push_back(str);
                }
                temppoint.x=stod(TempVec[1]);
                temppoint.y=stod(TempVec[2]);
                points.push_back(temppoint);
                //process.printnum(points[count].y);
                //process.printnum(points[count].x);
            }
        }

        //extract area velocity information
        if(line.find("ELEMENTS/CELLS")!=std::string::npos)
        {
            Velocity2d tempvel;
            //elements/cells
            for(int count=0;count<block.q;count++)
            {
                getline(in,line);
                TempVec.resize(0);
                line_data=std::istringstream(line);
                while(line_data>>str)
                {
                    TempVec.push_back(str);
                }
                if(stoi(TempVec[1])==2) //quadrilateral
                {

                    double x[5];
                    double y[5];
                    for(int adj=1;adj<=4;adj++)
                    {
                        x[adj]=points[stoi(TempVec[adj+2])-1].x;
                        y[adj]=points[stoi(TempVec[adj+2])-1].y;
                    }
                    //neu文件从1开始计数，转成从0开始计数
                    tempvel.vx=(x[1]+x[2]+x[3]+x[4])/4.0;
                    tempvel.vy=(y[1]+y[2]+y[3]+y[4])/4.0;

                    double cross1=(x[1]*y[2]-x[2]*y[1]);
                    double cross2=(x[2]*y[3]-x[3]*y[2]);
                    double cross3=(x[3]*y[4]-x[4]*y[3]);
                    double cross4=(x[4]*y[1]-x[1]*y[4]);

                    tempvel.area=0.5*fabs(cross1 + cross2 + cross3 + cross4);
                    
                }
                if(stoi(TempVec[1])==3) //triangle
                {
                    double x[4];
                    double y[4];
                    for(int adj=1;adj<=3;adj++)
                    {
                        x[adj]=points[stoi(TempVec[adj+2])-1].x;
                        y[adj]=points[stoi(TempVec[adj+2])-1].y;
                    }

                    tempvel.vx=(x[1]+x[2]+x[3])/3.0;
                    tempvel.vy=(y[1]+y[2]+y[3])/3.0;

                    double a = sqrt(pow((x[2]-x[3]),2)+pow((y[2]-y[3]),2));
                    double b = sqrt(pow((x[1]-x[3]),2)+pow((y[1]-y[3]),2));
                    double c = sqrt(pow((x[1]-x[2]),2)+pow((y[1]-y[2]),2));

                    // 计算三角形的半周长
                    double s = (a + b + c) / 2.0;
                    // 计算三角形的面积
                    tempvel.area = sqrt(s * (s - a) * (s - b) * (s - c));

                }
                vel.push_back(tempvel);
            }  
        }
    }

    //找到速度的最大绝对值（用于确定dt）
    //假定第一个是最大的速度
    block.max_vel=sqrt(pow(vel[0].vx,2)+pow(vel[0].vy,2));
    for(auto& v : vel)
    {
        if(sqrt(pow(v.vx,2)+pow(v.vy,2))>block.max_vel)
        {
            block.max_vel=sqrt(pow(v.vx,2)+pow(v.vy,2));
        }
    }


    process.print("unstructed velocity analysis completed.");
}


void Struct_velocity_D2Q12(MPIProcess& process,std::vector<Velocity2d>& vel,Block2d& block,std::vector<Point2d>& points)
{
    block.q=12;
    block.q_node=12;
    double r=sqrt(6.0*block.R*block.temp_init);
    double s=sqrt((9.0-3.0*sqrt(5.0))*block.R*block.temp_init/4.0);
    double t=sqrt((9.0+3.0*sqrt(5.0))*block.R*block.temp_init/4.0);
    double aa=1.0/36.0;
    double bb=((5.0+2.0*sqrt(5.0))/45.0);
    double cc=((5.0-2.0*sqrt(5.0))/45.0);


    double cx[12]={r,0,-r,0,s,-s,-s,s,t,-t,-t,t};  
    double cy[12]={0,r,0,-r,s,s,-s,-s,t,t,-t,-t};
    double weight[12]={aa,aa,aa,aa,bb,bb,bb,bb,cc,cc,cc,cc};

    for(int i=0;i<block.q;i++)
    {
        Velocity2d tempvel;
        tempvel.vx=cx[i];
        tempvel.vy=cy[i];


        tempvel.area=weight[i]*2.0*M_PI*block.R*block.temp_init*exp((cx[i]*cx[i]+cy[i]*cy[i])/(2.0*block.R*block.temp_init));

        vel.push_back(tempvel);
    }       
       process.print("structed velocity D2Q12 analysis completed.");

}


void Struct_velocity_D2Q37(MPIProcess& process,std::vector<Velocity2d>& vel,Block2d& block,std::vector<Point2d>& points)
{
    block.q=37;
    block.q_node=37;

    double r=1.19697977039307435897239*sqrt(block.R*block.temp_init);

    double aa=0.23315066913235250228650;
    double bb=0.10730609154221900241246;
    double cc=0.05766785988879488203006;
    double dd=0.01420821615845075026469;
    double ee=0.00535304900051377523273;
    double ff=0.00101193759267357547541;
    double gg=0.00024530102775771734547;
    double hh=0.00028341425299419821740;

    double cx[37]={
    0.0,
    1.0,0.0,-1.0,0.0,
    1.0,-1.0,-1.0,1.0,
    2.0,0.0,-2.0,0.0,
    2.0,1.0,-1.0,-2.0,-2.0,-1.0,1.0,2.0,
    2.0,-2.0,-2.0,2.0,
    3.0,0.0,-3.0,0.0,
    3.0,1.0,-1.0,-3.0,-3.0,-1.0,1.0,3.0};

    double cy[37]={
    0.0,
    0.0,1.0,0.0,-1.0,
    1.0,1.0,-1.0,-1.0,
    0.0,2.0,0.0,-2.0,
    1.0,2.0,2.0,1.0,-1.0,-2.0,-2.0,-1.0,
    2.0,2.0,-2.0,-2.0,
    0.0,3.0,0.0,-3.0,
    1.0,3.0,3.0,1.0,-1.0,-3.0,-3.0,-1.0};


    double weight[37]={
    aa,
    bb,bb,bb,bb,
    cc,cc,cc,cc,
    dd,dd,dd,dd,
    ee,ee,ee,ee,ee,ee,ee,ee,
    ff,ff,ff,ff,
    gg,gg,gg,gg,
    hh,hh,hh,hh,hh,hh,hh,hh
    };

    for(int i=0;i<block.q;i++)
    {
        cx[i]*=r;
        cy[i]*=r;

        Velocity2d tempvel;
        tempvel.vx=cx[i];
        tempvel.vy=cy[i];

        tempvel.area=weight[i]*2.0*M_PI*block.R*block.temp_init*exp((cx[i]*cx[i]+cy[i]*cy[i])/(2.0*block.R*block.temp_init));
        vel.push_back(tempvel);
    }
    process.print("structed velocity D2Q37 analysis completed.");



}