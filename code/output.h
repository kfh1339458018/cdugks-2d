#ifndef OUTPUT_H
#define OUTPUT_H

#include<string>
#include<cmath>
#include<iostream>
#include<fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include "config.h"
#include "mpi_process.h"
#include "mesh_read.h"

using namespace std;

void make_directory_for_result();

void outPutFile_dat(MPIProcess& process,const std::string& outputname,Block2d block,std::vector<Fluid2d>& fluids,std::vector<Point2d>& points);




#endif
