#ifndef TEST_CAVITY_H
#define TEST_CAVITY_H

#pragma once
#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include "config.h"

using namespace std;

void cavity();

void cavity_param(Block2d& block);

void cavity_init_param(Block2d& block,Fluid2d& initfluid);



#endif
