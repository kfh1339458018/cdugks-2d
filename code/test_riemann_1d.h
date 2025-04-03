#ifndef TEST_RIEMANN_1D_H
#define TEST_RIEMANN_1D_H

#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include "config.h"

using namespace std;

void riemann_1d();

void riemann_1d_param(Block2d& block);

void riemann_1d_init_param(Block2d& block,Fluid2d& initfluid);



#endif
