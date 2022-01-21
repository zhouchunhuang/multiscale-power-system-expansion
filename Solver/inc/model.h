#pragma once

#include "util.h"
#include "parameter.h"
#include "generator.h"
#include "transmission.h"
#include "storage.h"
#include "bus.h"
#include "rw.h"

class Node
{
public:
	int id;
	int start;		// start scenario
	int end;		// end scenario
	int parent;		// parent node
	int period;		// time period of the node
	std::vector<int> nodesInPath;	// the nodes from root node the the current node in the scenario path
	double probability;

	Node(int _id, int _t, int _start, int _end): id(_id), period(_t), start(_start), end(_end) {}
};

class Model
{
private:
	Algorithm alg;
	ParamRegistry* paramReg;

	// ======================= data transferred among computer nodes ================ //
	double** ColExpPlanGen;										// coefficients of lambda in the constraints, obtained from secondary MP
	double** ColExpPlanLine;
	double** ColExpPlanStg;
	double** ColExpPlanGenTemp;									// temporarily store expansion column
	double** ColExpPlanLineTemp;
	double** ColExpPlanStgTemp;
	double** DualGen;											// dual solution of the primary master problem
	double** DualLine;
	double** DualStg;
	double* DualConvex;
	double* ReducedCost;
	double* LambdaCoef;
	double* SmpOFV;
	double* SmpGap;
	double* PspTime;        // time spent in solving psp
	double* PspPhase2Time;  // time spent in solving psp phase 2, only applicable for NCD
	double* PspSmpTime;     // time spent in solving smp, only applicable for NCD

	double rankSolutionTime;
	double totalSolutionTime;     // used to compute the parallel computing efficiency

	// multiple columns
    double*** col_gen_expansion;
    double*** col_line_expansion;
    double*** col_stg_expansion;
    double** col_obj;

	std::vector<Node> scenarioNodes;							// upper level scenario nodes
	int numNodes;

	std::string output_directory;

	MPI_Comm sub_comm;
	int numDmdScenariosPerSubproblem;
	int numSPsPerNode, numSPsPerCore, numCoresPerNode;
	std::chrono::steady_clock::time_point smp_start;

protected:
	GRBEnv* env = new GRBEnv();
	// ======================= Master problem ====================== //
	GRBModel MPmodel = GRBModel(*env);
	GRBLinExpr objMP;

	//Decision variables
	GRBVar** var_xg;								//variables indicate the capacity made for each facility at each node
	GRBVar** var_xl;
	GRBVar** var_xs;
	std::vector< std::vector<GRBVar> > columns;	// generated columns
	GRBVar** var_y;									// lower bound for operational cost
	// GRBVar** var_lambda;

	//Constraints
	GRBConstr** SplitGen;					// multi-dimensional constraints for splitting variables [n][tau][g]
	GRBConstr** SplitLine;
	GRBConstr** SplitStg;
	GRBConstr* Convex;											// convexity constraints in primary master problem
	double MpUB, MpLB, MpNewLB, MpNewUB;
	double MpLpObj;

	// =================== Secondary Master Problem ================//
	GRBModel SMPmodel = GRBModel(*env);
	GRBLinExpr objSMP;
	double smp_rdCost, smp_opnCost;
	//Decision Variables
	GRBVar* var_zg;													// var_zg[tau][g]; whether a facility has been built until the current scenario node
	GRBVar* var_zl;													// var_zl[tau][l]
	GRBVar* var_zs;													// var_zs[tau][k]
	GRBVar* var_theta;
	std::vector<GRBConstr> BendersInfCuts;
	std::vector<GRBConstr> BendersOptCuts;
	std::vector<GRBConstr> LshapedInfCuts;
	std::vector<GRBConstr> LshapedOptCuts;
	double* sol_theta;
	double** SecDualGen;												//store dual of lower level subproblem
	double** SecDualLine;												//store dual of lower level subproblem
	double** SecDualStg;												//store dual of lower level subproblem
	double** SecDualDmd;
	double** SecDualRampUp;
	double** SecDualRampDown;
	double** SecDualMinDown;
	double* CapGen;														// cumulative capacity of each generator, obtained from solving SMP
	double** GenStatus;
	double* CapLine;
	double* CapStg;
	double* recFuncVal;													// used in combinatorial cut
	double* recFuncLB;
	int* SPstatus;														// optimization status of solving subproblem
	double* SPofv;															// objective of subproblem
	double SmpUB, SmpLB, SmpNewUB;
	int sub_message;
	int sub_iteration;
	int zEqOne;					// used in combinatorial cut
	GRBLinExpr combnXpr;		// used in combinatorial cut

	// =========================== Subproblem ======================//
	GRBModel SPmodel = GRBModel(*env);
	GRBLinExpr objSP;
	//Decision Variables
	GRBVar** varGenStatus;														//binary variable for on/off status of plants alpha[season][d][k][g]
	GRBVar** varStartUp;														//binary variable for start-up/shut-down of plants alpha[season][d][k][g]
	GRBVar** varPowerGen;														//amount of power generated from generator p[season][d][k][g]
	GRBVar*** varPowerFlow;														//power flow between bus i and bus j f[season][d][k][l]
	GRBVar** varPowerWithdrawal;												//total power withdrawn from storage facilities of bus u[season][d][k][sg]
	GRBVar** varPowerInject;													//total power injected into storage facilities of bus v[season][d][k][sg]
	GRBVar** varPowerRemaining;													//total remaining power in storage facilities of bus at the beginning of k v[season][d][k][sg]
	GRBVar** varUnmetDemand;
	GRBVar** varLinearize;

	GRBConstr* GenLimitConstrs;													// generation limit constraints in SP model	
	GRBConstr* GenLimit2Constrs;
	GRBConstr* LineLimitConstrs;
	GRBConstr* StgLimitConstrs;
	GRBConstr* DemandConstrs;
	GRBConstr* RampUpConstrs;
	GRBConstr* RampDownConstrs;
	std::vector<GRBConstr> MinDownConstrs;

	/* Store solution information in Single-level Benders Decomposition */
	double** CapGenInNodes;														// cumulative capacity of each generator, obtained from solving SMP
	double** GenStatusInNodes;
	double** CapLineInNodes;
	double** CapStgInNodes;
	double** recFuncValInNodes;													// used in combinatorial cut
	double** recFuncLBInNodes;
	int** SPstatusInNodes;														// optimization status of solving subproblem
	double** SPofvInNodes;
	double*** SecDualGenInNodes;												//store dual of lower level subproblem[n][sp][constraint] - benders decomposition												
	double*** SecDualLineInNodes;												//store dual of lower level subproblem
	double*** SecDualStgInNodes;												//store dual of lower level subproblem
	double*** SecDualDmdInNodes;
	double*** SecDualRampUpInNodes;
	double*** SecDualRampDownInNodes;
	double*** SecDualMinDownInNodes;

	// Solution
	double** sol_xg;
	double** sol_xl;
	double** sol_xs;
	double cmp_time;
	double investmentCost;
	double* operationalCost;
	
	//Parallel computing 
	int iteration;	
	double gap;
	double ofv;
	int message;
	int rank, master, num_ranks;
	int color, sub_rank, sub_master;
	
	// System data
	int numGen, numLine, numStg, numBus;
	int numRnGen;		// number of renewable energy generators
   //	int numScenarios, numPeriods;
   //	int numDmdScenarios;
   //	int numRealizations;
	int ucPeriods;		// number of unit commitment periods
	std::vector<Generator>	generators;
	std::vector<Generator*> rn_generators;		// renewable energy generators
	std::vector<Line>	lines;
	std::vector<Storage>	storages;
	std::vector<Bus> buses;
	
public:
	Model() :ucPeriods(HOURSPERWEEK){
#ifdef WEEK_AHEAD_UC
		ucPeriods = HOURSPERWEEK;
#else
		ucPeriods = HOURSPERDAY;
	   paramReg = ParamRegistry::instance();

	   rankSolutionTime = 0;
	   totalSolutionTime = 0;
#endif
	}

	// solution algorithms
	int optimize();
	
	Result columnGeneration();
	Result direct();
	Result expectedValForm();
	Result bendersDecomposition();
	Result nestedDecomposition();

	void computeStochasticValue();
	void getExpansionSolution();
	
	// column generation
	void initialize(bool nested = false);
	void createMasterProblem();										//create model of (primary) master problem
	void createUpperSubproblem();										//create model of subproblem
	void getMpDuals();
	Result solveUpperSubproblem();
	void get_initial_columns();
	void get_second_columns();
	void addColumns(bool nested = false);
	void getMoreColumns(int phase = 1);
	void addMoreColumns();				// add multiple columns per iteration
	void getMPbounds();
	Result solveMasterProblem(bool nested = false);										// solve the master problem after the last CG iteration

	// nested benders decomposition for solving smp
	void createSecondaryMP();
	void createSubproblem(bool nested=false);
	Result solvePrimarySubproblem(int phase = 1);
	void updateSecondaryMP(int phase = 1);
	Result solveSecondaryMP();
	void getSecondaryMpBounds();
	void addNestBendersCuts(int phase = 1);
	Result getRecFuncLB();			// get lower bound of recourse function by solving every subproblem
	Result solveSubproblem(int sp/*no.SubproblemInCore*/, int phase = 1);		// solve the lower-level subproblem at node n in _season, _s'th subproblem

	// single-level benders decomposition - subproblems will be dealt with in parallel
	void bendersLoop(int phase=1);
	void createBendersMasterProblem();
	void initializeBenders();
	void dataTransfer();
	Result solveBendersMP();
	Result solveBendersSubproblem(int sp/*no.SubproblemInCore*/, int phase = 1);
	void getRecFuncLBInNodes();
	void getBendersMpBounds();
	void addBendersCuts(int phase=1);

	void test();
	void testSolveSPsPhaseI();
	void testSolveSPsPhaseII();
	
	void set_rank(int _r) { rank = _r; }
	int get_rank() const { return rank; }
	void set_num_ranks(int _c) { num_ranks = _c; }
	int get_num_ranks() const { return num_ranks; }
	void set_master(int _m) { master = _m; }
	int get_master() const { return master; }
	bool is_master_proc() const { return rank == master; }
	void set_sub_rank(int _r) { sub_rank = _r; }
	int get_sub_rank() const { return sub_rank; }
	void set_sub_master(int _m) { sub_master = _m; }
	int get_sub_master() const { return sub_master; }
	bool is_sub_master_proc() const { return sub_rank == sub_master; }
	void set_algorithm(Algorithm _alg) { alg = _alg; }
	int get_algorithm() const { return alg; }
	void set_output_directory(std::string _s) { output_directory = _s; }
	std::string get_output_directory() const { return output_directory; }

	/* ********lower level two stages index************ */
	int get1stStageIndex(int n, int season, int h) const {
		return n * static_cast<int>(Season::numSeasons) * ucPeriods + season * ucPeriods + h;
	}
	int get2ndStageIndex(int n, int season, int d, int h) const { 
	   return n * static_cast<int>(Season::numSeasons) * paramReg->numDmdScenarios * ucPeriods + season * paramReg->numDmdScenarios * ucPeriods + d * ucPeriods + h;
	}
	int getLowerIndex(int n, int season, int d, int h) const {
	   return n * static_cast<int>(Season::numSeasons) * paramReg->numDmdScenarios * ucPeriods + season * paramReg->numDmdScenarios * ucPeriods + d * ucPeriods + h;
	}
	int getLowerIndexInSP(int d, int h) const {
		return d * ucPeriods + h;
	}
	int getLowerIndexExpectedForm(int n, int season, int h) const {
		return n * static_cast<int>(Season::numSeasons) * ucPeriods + season * ucPeriods + h;
	}
	int get1stStageIndexInNode(int season, int h) const {			// get index in lower level subproblem in single-layer column generation
		return season * ucPeriods + h;
	}
	int get2ndStageIndexInNode(int season, int d, int h) const {			// get index in lower level subproblem in single-layer column generation
	   return season * paramReg->numDmdScenarios * ucPeriods + d * ucPeriods + h;
	}	
	/* ******************************************************** */

	std::pair<int, int> getScenarioStage(int index) const { return std::make_pair(int(index / paramReg->numPeriods), index % paramReg->numPeriods); }
	std::vector<int> getPeriodDmdscnHour(int index) const {
	   return { int(index / (paramReg->numDmdScenarios * ucPeriods)), int(int(index % (paramReg->numDmdScenarios * ucPeriods)) / ucPeriods), int(index % (paramReg->numDmdScenarios * ucPeriods)) % ucPeriods };
	}
	Node* getNode(int s, int t)											// get the scenario node index by scenario and period
	{
		auto it = std::find_if(scenarioNodes.begin(), scenarioNodes.end(), [s, t](auto & node) {
			return (node.period == t && node.start <= s && s < node.end);
		});
		return &(*it);
	}
	
	int getMpConstrIndex(int tau, int facility_no, FacilityType _type = FacilityType::GENERATOR) const {
		int numFacility = 0;
		switch (_type)
		{
		case FacilityType::GENERATOR:
			numFacility = numGen;
			break;
		case FacilityType::LINE:
			numFacility = numLine;
			break;
		case FacilityType::STORAGE:
			numFacility = numStg;
			break;
		default:
			break;
		}
		return tau * numFacility + facility_no;
	}
	int getNumMpConstrs(const Node &node, FacilityType _type = FacilityType::GENERATOR) const { 
		int numFacility = 0;
		switch (_type)
		{
		case FacilityType::GENERATOR:
			numFacility = numGen;
			break;
		case FacilityType::LINE:
			numFacility = numLine;
			break;
		case FacilityType::STORAGE:
			numFacility = numStg;
			break;
		default:
			break;
		}
		return (node.period + 1) * numFacility;
	}

	int get_subproblem_no_in_node(int sp/*no.subproblem per core*/) const { return sub_rank * numSPsPerCore + sp; }

	int get_dmd_scenario(int d, int noSpInNode/*no.subproblem per node*/) const { return (noSpInNode % paramReg->numSPsPerSeason) * numDmdScenariosPerSubproblem + d; }		// given the scenario number "d" in subproblem i, get the scenario number in a season

	int get_season_of_subproblem(int noSpInNode) const { return int(noSpInNode / paramReg->numSPsPerSeason); }		// in which season the subproblem solved in the subproblem i belongs to

	double get_opn_cost_factor(int _season) const {
		double cost_factor = double(WEEKSPERSEASON);
#ifdef DAY_AHEAD_UC
		if ((_season == SpringWkdy) || (_season == SummerWkdy) || (_season == FallWkdy) || (_season == WinterWkdy)) {
			cost_factor = double(WKDYSPERSEASON);
		}
		else {
			cost_factor = double(WKNDSPERSEASON);
		}
#endif
		return cost_factor * YEARS_PER_STAGE;
	}

	void constructScnStructure();

	friend class ReadWrite;
};




