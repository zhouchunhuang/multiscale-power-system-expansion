// Created by Zhouchun Huang on 2021/9/11.
#pragma once

#include "util.h"

enum RenewableType {
   NONE = 0,
   WIND = 1,
   SOLAR = 2
};

struct Generator
{
   int id;
   std::string name;
   int loc;
   int type;
   double startUpCost;
   double shutDownCost;
   int minDownTime;
   int minUpTime;
   double powerMax;
   double powerMin;
   double rampUp;
   double rampDown;
   double maxSpinRsv;
   double capacity;
   double currentCapacity;
   double varCost;
   double fuelCostA;
   double fuelCostB;
   double fuelCostC;

   int* cost;					// upper level uncertain costs: cost[n] for each scenario-tree node n
   double* fuelPriceRate;

   bool is_renewable;
   vector2int supply_mean;				// mean for uncertain supply for renewable energy
   // int** supply_std;				// std for uncertain supply for renewable energy
   vector2double distn_params;	// statistical distribution parameters. for wind it is weibull distribution shape and scale parameters; for solar it is normal distribution mean and std
   std::vector<int> supply_scenarios;

   Generator(int _numNodes)
   {
      cost = new int[_numNodes];
      fuelPriceRate = new double[_numNodes];
      is_renewable = false;
   }

   bool exists() const { return (currentCapacity > 0); }

   RenewableType get_rn_type() const {
      if (name.find("WIND") != std::string::npos) {
         return RenewableType::WIND;
      }
      if (name.find("SOLAR") != std::string::npos) {
         return RenewableType::SOLAR;
      }
      return RenewableType::NONE;
   }

   static bool isSameGen(const Generator &gen1, const Generator &gen2){
      if((!gen1.is_renewable) || (!gen2.is_renewable)) return false;
      if(gen1.loc != gen2.loc) return false;
      if(gen1.name != gen2.name) return false;
      return true;
   }

   //~Generator()
   //{
   //	delete[] cost;
   //}
};
