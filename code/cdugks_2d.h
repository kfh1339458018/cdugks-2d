#ifndef CDUGKS_2D_H
#define CDUGKS_2D_H

# pragma once
# include <cmath>
# include <string>
# include "mpi_process.h"
#include "config.h"
# include "output.h"


inline double CalcGs(Block2d* block, Fluid2d* fluids, Velocity2d* velq);
inline double CalcHs(Block2d* block, Fluid2d* fluids, Velocity2d* velq);
inline double Calc_s(Block2d* block, Fluid2d* fluids, Velocity2d* velq,DistributionType type);
inline double bCalc_s(Block2d* block, Fluid2d* fluids, Velocity2d* velq,int adj,DistributionType type);
inline double bCalcGs(Block2d* block, Fluid2d* fluids, Velocity2d* velq, int adj);
inline double bCalcHs(Block2d* block, Fluid2d* fluids, Velocity2d* velq, int adj);
inline int pointLocationInTriangle(double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Px, double Py);
inline int pointLocationInQuadrilateral(double Ax, double Ay, double Bx, double By, double Cx, double Cy, double Dx, double Dy, double Px, double Py);
void init(MPIProcess& process,Block2d& block,std::vector<Fluid2d>& fluids,std::vector<Point2d>& points,std::vector<Velocity2d>& vel,std::vector<Velocity2d>& velq);
void AuxDistribution(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq, DistributionType type);
void InterPolation(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq, DistributionType type);
void BoundaryMacroscopicVariables(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq);
void EquilibriumAndNext(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq, DistributionType type);
void BoundaryCondition(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq);
void CalcFmeso(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq);
void CalcFmacro(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq);
void UpdateMacroscopicVariables(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq);
void UpdataDistribution(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq);
bool SteadyJudge(MPIProcess& process, Block2d& block, std::vector<Fluid2d>& fluids, std::vector<Velocity2d>& vel, std::vector<Velocity2d>& velq);


#endif