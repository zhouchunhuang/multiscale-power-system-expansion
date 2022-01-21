#pragma once
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <set>
#include <ctime>
#include "mpi.h"
#include <chrono>

#include "gurobi_c++.h"

//#define LPMODEL
//#define FILEOUTON
#define SolverStreamOff
#define StreamAlgInfo
//#define RestrictSameSiteExpansion
//#define PhaseLinearMP
//#define MULTI_COLS
#define PhaseII
//#define TEST_MODE_2
//#define TEST_MODE
#ifdef TEST_MODE_2
#ifndef LPMODEL
#define LPMODEL
#endif
#endif
//#define GetModelSize

typedef std::vector< std::vector<double> > vector2double;
typedef std::vector<vector2double> vector3double;
typedef std::vector<vector3double> vector4double;
typedef std::vector< std::vector<int> > vector2int;

constexpr auto HOURSPERWEEK = 168;
constexpr auto HOURSPERDAY = 24;
constexpr auto DAYSPERSEASON = 45;
constexpr auto WKDYSPERSEASON = 65;
constexpr auto WKNDSPERSEASON = 26;
constexpr auto WEEKSPERSEASON = 13;
constexpr auto BIGM = 1e20;
constexpr auto inf_rate = 0.1;
constexpr auto stddev = 0.1;			// the standard dev for generating load data and solar energy supply equals stddev multiplied by mean
constexpr auto POWER_WIND_MIN = 3.0;
constexpr auto POWER_WIND_MAX = 12.0;
constexpr auto WkndPeakLoadRatio = 0.8;
constexpr auto NUM_COLUMNS = 3;
constexpr auto MAX_COLUMNS = 50000;
constexpr auto YEARS_PER_STAGE = 5;
constexpr auto depression_rate = 0.3;

//#define WEEK_AHEAD_UC
#ifndef WEEK_AHEAD_UC
#define DAY_AHEAD_UC
#endif

#define PRINT_SECTION(log) {std::cout << "=============" << log << "==============" << std::endl;}
#define PRINT_SUBSECTION(log) {std::cout << "--" << log << "--" << std::endl;}

enum Algorithm
{
   DIRECT = 1,
   EXPECT = 2,
   COLGEN = 3,
   BENDERS = 4,
   NESTED = 5,
   TEST = 6
};

enum Result {
   OPTIMAL = GRB_OPTIMAL,
   INFEASIBLE = GRB_INFEASIBLE,
   UNBOUNDED = GRB_UNBOUNDED
};

enum SpStatus {
   Optimal = 1,
   LRoptimal = 2,
   LRinfeasible = 3,
   MIPinfeasible = 4
};

enum Season {
   SpringWkdy = 0,
   SpringWknd = 1,
   SummerWkdy = 2,
   SummerWknd = 3,
   FallWkdy = 4,
   FallWknd = 5,
   WinterWkdy = 6,
   WinterWknd = 7,
   numSeasons = 8
};

enum FacilityType {
   GENERATOR = 0,
   LINE = 1,
   STORAGE = 2
};
