#pragma once
#include <vector>

class ParamRegistry{
private:
   static ParamRegistry* paramInstance;
   ParamRegistry();

public:

   static ParamRegistry* instance(){
      if(!paramInstance){
         paramInstance = new ParamRegistry();
      }
      return paramInstance;
   }

   int numScenarios;
   int numPeriods;
   int numRealizations;
   double prob_scenario;

   std::vector<double> fossilCapex;
   std::vector<double> gasCapex;
   std::vector<double> renewableCapex;
   std::vector<double> storageCapex;
   std::vector<double> lineCapex;
   std::vector<double> demandChange;
   std::vector<double> fuelPriceChange;

   int numDmdScenarios;
   int numSPsPerSeason;			// how many subproblems per season week
   double prob_dmd_scenario;

   double time_lmt;
   double smp_timlmt;
   double sp_tmlmt;
   double remove_cuts_threshold;
   double sp_gap;
   double mp_gap;
   double smp_gap;
   double max_smp_gap;

   double lossPower;				// 1% power will loss from power flow
   double stgCeff;					//certain percentage of power can be used from that withdrawn from storage facility
   double unmetDmdPnlty;

   double bigM;
   bool keep_benders_cuts;
   double reserveMarginRate;
};

