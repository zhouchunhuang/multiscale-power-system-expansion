#pragma once

enum Direction{
   Forward = 0,
   Backward = 1,
   numDirections = 2
};

struct Line
{
   int id;
   int from;
   int to;
   double imped;
   double volt;
   double capacity;
   double currentCapacity;
   double varCost;         // capex cost per MW
   double length;          // length of the line
   double lossRatio;       // loss per MW

   int* cost;

   Line(int numScenarios, int numPeriods): id(-1), from(-1), to(-1), imped(0), volt(0), capacity(0),
   currentCapacity(0), varCost(0), length(0), lossRatio(0)
   {
      cost = new int[numScenarios * numPeriods];
   }

   bool exists() const { return (currentCapacity > 0); }

   /*~Line()
   {
      delete[] cost;
   }*/
};
