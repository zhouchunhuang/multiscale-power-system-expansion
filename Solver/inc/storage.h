#pragma once

struct Storage
{
   int id;
   int loc;
   double capacity;
   double currentCapacity;
   double varCost;

   int* cost;

   Storage(int numScenarios, int numPeriods)
   {
      cost = new int [numScenarios * numPeriods];
   }

   bool exists() const { return (currentCapacity > 0); }

   static bool isSameStorage(const Storage &stg1, const Storage &stg2){
      return (stg1.loc == stg2.loc);
   }

   /*~Storage()
   {
      delete[]cost;
   }*/
};
