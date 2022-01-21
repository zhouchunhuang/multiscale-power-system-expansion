#pragma once

class Model;

class ReadWrite
{
private:
   std::string configFile, busFile, genFile, lineFile, stgFile;
   std::string weeklyLoadFile, dailyLoadFile, hourlyLoadFile;
   std::string boundsFile;
   std::string inputDirectory, outputDirectory;
   std::string paramFile;

public:
   ReadWrite() {
      configFile = "/Configuration.dat";
      paramFile = "/Parameters.dat";
      genFile = "/Generators.csv";
      lineFile = "/Lines.csv";
      stgFile = "/Storages.csv";
      busFile = "/Bus.csv";

      weeklyLoadFile = "/weekly_load.txt";
      dailyLoadFile = "/daily_load.txt";
      hourlyLoadFile = "/hourly_load.txt";

      boundsFile = "/bounds.txt";
   };

   void readParameters();
   void readSystemData(Model &_model);
   void readBusData(Model &_model);
   void readLoadData(Model &_model);
   void readGenerators(Model &_model);
   void readRenewable(Model &_model);		// renewable energy data
   void readLines(Model &_model);
   void readStorages(Model &_model);

   void writeGeneratedData(Model &_model);
   void printSolution(Model &_model);
   static void dataGeneration(Model &_model);		// generate upper level cost scenario and lower level demand scenario data
   void printBounds(int iteration, double ub, double* reducedCost, int numNodes, double lb, double gap);

   std::pair<double, double> split_range(const std::string & s);
   template <typename T> T NumGen(T minNum, T maxNum);
   std::vector<std::string> split(const std::string &s) const;

   void set_inputDirectory(const std::string &_inputDir) { inputDirectory = _inputDir; }
   void set_outputDirectory(const std::string &_outputDir) { outputDirectory = _outputDir; }
   std::string get_outputDirectory() const { return outputDirectory; }
};
