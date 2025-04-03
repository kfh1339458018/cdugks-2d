#ifndef TEST_CYLINDER_H
#define TEST_CYLINDER_H

#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include "config.h"

using namespace std;

void cylinder();

void cylinder_param(Block2d& block);

void cylinder_init_param(Block2d& block,Fluid2d& initfluid);



#endif
