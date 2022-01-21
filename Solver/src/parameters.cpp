#include "../inc/parameter.h"

ParamRegistry* ParamRegistry::paramInstance = nullptr;


ParamRegistry::ParamRegistry()
{
   numPeriods = 3;
   numRealizations = 2;
   numScenarios = 4;

   fossilCapex = {0.95, 1.0};
   renewableCapex = {0.9, 1.0};
   storageCapex = {0.9, 1.0};
   lineCapex = {1.0, 1.0};
   demandChange = {1.03, 1.10};
   fuelPriceChange = {0.95, 1.03};

   numDmdScenarios = 20;
   numSPsPerSeason = 20;
   reserveMarginRate = 1.135;

   time_lmt = 3600 * 48;
   smp_timlmt = 600;
   sp_tmlmt = 300;
   remove_cuts_threshold = 600;
   mp_gap = 0.01;
   smp_gap = 0.005;
   sp_gap = 0.005;
   max_smp_gap = 0.2;

   lossPower = 0.0001;    // 1% per 100 miles = 0.01% per mile
   stgCeff = 0.95;
   unmetDmdPnlty = 1e10;

   bigM = 1e20;
   keep_benders_cuts = true;
}

