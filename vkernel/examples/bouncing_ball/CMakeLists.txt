# SPDX-FileCopyrightText: 2026 INRIA
#
# SPDX-License-Identifier: Apache-2.0
#
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.

cmake_minimum_required(VERSION 3.20)
project(vkernel_bouncing_ball CXX)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(Siconos REQUIRED)
add_executable(bouncing_ball main.cpp)
target_link_libraries(bouncing_ball PRIVATE Siconos::vkernel)
