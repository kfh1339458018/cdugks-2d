#include "output.h"



//创建结果文件夹
void make_directory_for_result()
{
	string newfolder = add_mesh_directory_modify_for_linux()+"/../result";

#if defined(_WIN32)
	//一个简单的宏定义判断，如果是win系统 则运行以下代码
	if (_mkdir(newfolder.c_str()) != 0)
	{
		_mkdir(newfolder.c_str());
		cout << "creating result folder" << endl;
	}
	else
	{
		cout << "fail to create new folder for result file" << endl;
	}
	//end
#else 
	//如果是linux系统 则运行以下代码
	mkdir(newfolder.c_str(), 0777);
	if (mkdir(newfolder.c_str(), 0777) != 0)
	{
		mkdir(newfolder.c_str(), 0777);
	}
	else
	{
		cout << "fail to create new folder for result file" << endl;
	}
	//end
#endif
}



//输出.dat文件
void outPutFile_dat(MPIProcess& process,const std::string& outputname,Block2d block,std::vector<Fluid2d>& fluids,std::vector<Point2d>& points)
{

	MPI_Barrier(MPI_COMM_WORLD);
	if(process.pid==0)
	{
		std::string filepath = "../../result/" + outputname + std::to_string(block.step)+ ".dat";  

		//打开文件
		std::ofstream outputFile(filepath);
		std::string datapacking = "block";
		std::string varlocation[2] = {"[1-2]=NODAL", "[3-14]=CELLCENTERED"};
		std::string zonetype = "fequadrilateral";

		if (outputFile.is_open()) 
		{
			outputFile << "TITLE = "<< outputname << std::endl;
			outputFile << R"(VARIABLES = "X", "Y", "Rho", "temp", "pressure","qx_flux","qy_flux","vx","vy","Kn-loc","Kn-ref","Kn-gll","vx-ref","vy-ref" )" << std::endl;

			// 添加TECPLOT Grid文件头部信息
			outputFile << "zone N=" << block.TotalPoints << ", E=" << block.TotalCells << ", datapacking=" << datapacking << "\n";
			outputFile << "varlocation=(" << varlocation[0] << "," << varlocation[1] << ")\n";
			outputFile << "zonetype=" << zonetype << "\n";
		} 
		else 
		{
			std::cerr << "Error opening the file!" << std::endl;
		}
		
		//写入数据
		for(auto&point : points)
		{
			outputFile << point.x << "\n";
		}
		for(auto&point : points)
		{
			outputFile << point.y << "\n";
		}
		for(auto&fluid : fluids)
		{
			outputFile << fluid.rho << "\n";
		}
		for(auto&fluid : fluids)
		{
			outputFile << fluid.temp << "\n";
		}
		for(auto&fluid : fluids)
		{
			outputFile << fluid.pressure << "\n";
		}
		for(auto&fluid : fluids)
		{
			outputFile << fluid.qx_flux << "\n";
		}
		for(auto&fluid : fluids)
		{
			outputFile << fluid.qy_flux << "\n";
		}
		for(auto&fluid : fluids)
		{
			outputFile << fluid.vx << "\n";
		}
		for(auto&fluid : fluids)
		{
			outputFile << fluid.vy << "\n";
		}
		for(auto&fluid : fluids)
		{
			outputFile << (2.0*fluid.miu*(5.0-2.0*block.omega)*(7.0-2.0*block.omega)/(15.0*fluid.rho*block.l_ref*sqrt(2.0*M_PI*block.R*fluid.temp)))<< "\n";
		}
		for(auto&fluid : fluids)
		{
			outputFile << (2.0*fluid.miu*(5.0-2.0*block.omega)*(7.0-2.0*block.omega)/(15.0*fluid.rho*block.l_ref*sqrt(2.0*M_PI*block.R*fluid.temp)))/block.Kn << "\n";
		}
		//最小二乘，计算Kn-gll
        for(int i=0;i<block.TotalCells;i++)
        {

            //相邻格心的坐标
            double centerx_adj[5], centery_adj[5];

            for (int adj = 1; adj <= fluids[i].edge_num; adj++) {
                //如果相邻的控制体为0，跳过
                if (fluids[i].adj[adj] == 0) continue;
                centerx_adj[adj] = fluids[fluids[i].adj[adj]].x;
                centery_adj[adj] = fluids[fluids[i].adj[adj]].y;
            }

            //格心的坐标
            double xc, yc;
            xc = fluids[i].x;
            yc = fluids[i].y;
            double x11 = 0.0, x12 = 0.0, x21 = 0.0, x22 = 0.0, y1 = 0.0, y2 = 0.0;

            for (int adj=1; adj <= fluids[i].edge_num; adj++) {
                //如果相邻的控制体为0，跳过
                if (fluids[i].adj[adj] == 0) continue;
                x11 += (centerx_adj[adj] - xc) * (centerx_adj[adj] - xc);
                x12 += (centerx_adj[adj] - xc) * (centery_adj[adj] - yc);
                x21 += (centerx_adj[adj] - xc) * (centery_adj[adj] - yc);
                x22 += (centery_adj[adj] - yc) * (centery_adj[adj] - yc);
                y1 += (centerx_adj[adj] - xc) * (fluids[fluids[i].adj[adj]].rho - fluids[i].rho);
                y2 += (centery_adj[adj] - yc) * (fluids[fluids[i].adj[adj]].rho - fluids[i].rho);

            }

            double determinant = x11 * x22 - x12 * x21;
            if (determinant == 0) {
                printf("error! gradient can not be solve!\n"); //分母为零，梯度方程无解
            }
            //解方程求梯度
            double gradientx = (y1 * x22 - y2 * x12) / determinant;
            double gradienty = (x11 * y2 - x21 * y1) / determinant;


            //自由程
            double lambda_=(2.0*fluids[i].miu*(5.0-2.0*block.omega)*(7.0-2.0*block.omega)/(15.0*fluids[i].rho*sqrt(2.0*M_PI*block.R*fluids[i].temp)));

            double Kn_gll=lambda_/fluids[i].rho*(abs(gradientx)+abs(gradienty));

            outputFile << Kn_gll <<"\n";

            

        }



		for(auto&fluid : fluids)
		{
			outputFile << fluid.vx/block.Ma << "\n";
		}
		for(auto&fluid : fluids)
		{
			outputFile << fluid.vy/block.Ma << "\n";
		}

		
		for(auto&fluid : fluids)
		{
			//格点存储从0开始，需要+1以满足文件要求
			outputFile << " " << fluid.v[1]+1 << " " << fluid.v[2]+1 << " " << fluid.v[3]+1 << " " << fluid.v[4]+1 << std::endl;
		}

		outputFile.close();
	}
    
}

