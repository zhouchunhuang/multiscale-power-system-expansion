#pragma once

class Bus
{
private:
   double peak_load;			// peak load in first stage. next stages the peak load will be uncertain

public:
   friend bool operator==(const Bus &lhs, const Bus &rhs) { return lhs.id == rhs.id; }
   friend bool operator<(const Bus &lhs, const Bus &rhs) { return lhs.id < rhs.id; }
   friend bool operator>(const Bus &lhs, const Bus &rhs) { return lhs.id > rhs.id; }
   friend bool operator!=(const Bus &lhs, const Bus &rhs) { return lhs.id != rhs.id; }

   int id;

   std::vector<int> peak_loads_in_scenario_nodes;																// peak load[node]
   std::vector<std::vector<std::vector<int> > > season_hourly_loads;											// expected hourly load[node][season][hour]
   std::vector<int> season_hourly_load_scenarios;					// load scenarios[get2ndStageIndex(n, season, d, h)]

   explicit Bus(int _numNodes): id(-1), peak_load(0)
   {
      peak_loads_in_scenario_nodes.clear();
      season_hourly_loads.clear();
      season_hourly_load_scenarios.clear();
      peak_loads_in_scenario_nodes = std::vector<int>(_numNodes);
      // season_hourly_loads = std::vector<std::vector<std::vector<int> > >(numPeriods, std::vector<std::vector<int> >(Season::numSeasons, std::vector<int>(ucPeriods)));
   }

   void set_peak_load(double _pk) { peak_load = _pk; }
   double get_peak_load() const { return peak_load; }
};
