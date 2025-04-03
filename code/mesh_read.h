#ifndef MESH_RDAD_H
#define MESH_RDAD_h


#include "cdugks_2d.h"
#include "mpi_process.h"
#include "config.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream> 
#include <vector>



using namespace std;

string add_mesh_directory_modify_for_linux();

void Read_mesh_neu_2d(MPIProcess& process,vector<Fluid2d>& fluids,Block2d& block,std::vector<Point2d>& points);
void Unstruct_velocity(MPIProcess& process,std::vector<Velocity2d>& vel,Block2d& block,std::vector<Point2d>& points);
void Struct_velocity_D2Q12(MPIProcess& process,std::vector<Velocity2d>& vel,Block2d& block,std::vector<Point2d>& points);
void Struct_velocity_D2Q37(MPIProcess& process,std::vector<Velocity2d>& vel,Block2d& block,std::vector<Point2d>& points);


#endif
