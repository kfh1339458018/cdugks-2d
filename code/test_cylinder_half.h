#ifndef TEST_CYLINDER_HALF_H
#define TEST_CYLINDER_HALF_H

#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include "config.h"

using namespace std;

void cylinder_half();

void cylinder_half_param(Block2d& block);

void cylinder_half_init_param(Block2d& block,Fluid2d& initfluid);



#endif
