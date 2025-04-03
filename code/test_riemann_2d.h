#ifndef TEST_RIEMANN_2D_H
#define TEST_RIEMANN_2D_H

#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include "config.h"

using namespace std;

void riemann_2d();

void riemann_2d_param(Block2d& block);

void riemann_2d_init_param(Block2d& block,Fluid2d& initfluid);



#endif
