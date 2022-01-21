#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include "model.h"
#include <mpi.h>

//#define DIRECT_SOLVE_MULTIPLE_COLUMNS

using namespace std;

Result Model::nestedDecomposition()
{
	if (is_master_proc())
	{
		PRINT_SECTION("Solving the Problem with Nested Column Generation");
		PRINT_SUBSECTION("Initializing the algorithm");
		std::cout.precision(4);		
	}

	double mpRunTime = 0;
	ReadWrite rw;
	rw.set_outputDirectory(output_directory);

	initialize(true);

	createMasterProblem();

	createSecondaryMP();

	createSubproblem(true);

	get_initial_columns();

	//get_second_columns();
	
	if (is_master_proc())
	{																	// create master problem model
		rw.printBounds(0, MpUB, nullptr, numNodes, MpLB, MpUB);
#ifdef StreamAlgInfo
		std::cout << "UB\t\t";
		for (int n = 0; n < numNodes; n++)  std::cout << n << "\t\t";
		std::cout << "LB" << std::endl;
#endif
	}
														
	/* *********************************************************** Nested Decomposition Algorithm********************************************************** */
	auto _start = std::chrono::steady_clock::now();
	iteration = 0;
	rankSolutionTime = 0;
	while (true) {
		iteration++;
		//................Master processor..................//
		if (is_master_proc())
		{
		   auto mp_start_time = std::chrono::steady_clock::now();

			MPmodel.optimize();																		//solve master problem
			if (MPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
#ifdef LPMODEL
				char buf[100];
				std::sprintf(buf, "%s/mp_infeasible.lp", output_directory.c_str());
				MPmodel.write(buf);
#endif // LPMODEL
				std::cerr << "Primary Master Problem is infeasible!" << std::endl;
				exit(EXIT_FAILURE);
			}
			else if (MPmodel.get(GRB_IntAttr_Status) == GRB_NUMERIC) {
				MPmodel.set(GRB_IntParam_Crossover, 1);
				MPmodel.optimize();
				MPmodel.set(GRB_IntParam_Crossover, 0);
			}
			mpRunTime = MPmodel.get(GRB_DoubleAttr_Runtime);
			MpUB = MPmodel.get(GRB_DoubleAttr_ObjVal);
			getMpDuals();

			auto mp_end_time = std::chrono::steady_clock::now();
			rankSolutionTime += double(std::chrono::duration_cast<std::chrono::seconds>(mp_end_time - mp_start_time).count());;
		}

		for (int n = 0; n < numNodes; n++) {
			MPI_Bcast(DualGen[n], numGen, MPI_DOUBLE, master, MPI_COMM_WORLD);
			MPI_Bcast(DualLine[n], numLine, MPI_DOUBLE, master, MPI_COMM_WORLD);
			MPI_Bcast(DualStg[n], numStg, MPI_DOUBLE, master, MPI_COMM_WORLD);
			MPI_Bcast(DualConvex, numNodes, MPI_DOUBLE, master, MPI_COMM_WORLD);
		}

		auto psp_start_time = std::chrono::steady_clock::now();
		PspSmpTime[color] = 0;

		solvePrimarySubproblem(0);

		solvePrimarySubproblem(1);																			/* ******** Solve subproblem and Get a feasible expansion plan as new columns ******** */

#ifdef PhaseII
      auto psp_phase2_start_time = std::chrono::steady_clock::now();
		solvePrimarySubproblem(2);
		auto psp_phase2_end_time = std::chrono::steady_clock::now();
		PspPhase2Time[color] = double(std::chrono::duration_cast<std::chrono::seconds>(psp_phase2_end_time - psp_phase2_start_time).count());
#endif

		auto psp_end_time = std::chrono::steady_clock::now();
		PspTime[color] = double(std::chrono::duration_cast<std::chrono::seconds>(psp_end_time - psp_start_time).count());

		addColumns(true);																					/* ******** Add columns to master problem ******** */

		getMPbounds();

		if(is_sub_master_proc()){
		   if(PspSmpTime[color] > paramReg->remove_cuts_threshold){
		      for (auto &it : BendersInfCuts) {
		         SMPmodel.remove(it);
		      }
		      BendersInfCuts.clear();
		      for (auto &it : BendersOptCuts) {
		         SMPmodel.remove(it);
		      }
		      BendersOptCuts.clear();
		   }
		}

#ifdef MULTI_COLS
        getMoreColumns(1);
        addMoreColumns();
#endif
		if (is_master_proc()) {
		   if (gap <= paramReg->mp_gap) message = 1;															// generate columns until the gap is relatively small
			auto _end = std::chrono::steady_clock::now();
			cmp_time = double(std::chrono::duration_cast<std::chrono::seconds>(_end - _start).count());
			if (cmp_time >= paramReg->time_lmt)
				message = 1;

#ifdef StreamAlgInfo
			std::cout << "-" << iteration << "-(" << cmp_time << ")\t";
			for (int n = 0; n < numNodes; n++) {
			   std::cout << SmpGap[n] << "(";
			   std::cout <<  PspSmpTime[n] << "+" << PspTime[n] - PspPhase2Time[n] << "+" << PspPhase2Time[n];
			   std::cout << ")\t\t";
			}
			std::cout << std::endl;
			std::cout << MpUB << "(" << mpRunTime << ")\t\t";
			for (int n = 0; n < numNodes; n++)  std::cout << ReducedCost[n] << "\t\t";
			std::cout << MpLB << std::endl;
#endif

			rw.printBounds(iteration, MpUB, ReducedCost, numNodes, MpLB, gap);
		}
		MPI_Bcast(&gap, 1, MPI_DOUBLE, master, MPI_COMM_WORLD);
		MPI_Bcast(&message, 1, MPI_INT, master, MPI_COMM_WORLD);
		if (message) break;
	}

	// solve the primary master problem last time
	solveMasterProblem(true);

	MPI_Reduce(&rankSolutionTime, &totalSolutionTime, 1, MPI_DOUBLE, MPI_SUM, master,
              MPI_COMM_WORLD);

	return Result::OPTIMAL;
}

/* ************************************ Benders Decomposition in Lower Level: Solve the MIP at each node n within N ************************************ */
Result Model::solvePrimarySubproblem(int phase)
{
#ifdef TEST_MODE
	if (is_master_proc()) PRINT_SUBSECTION("Solve the CG subproblem with Benders Decomposition");
#endif

	updateSecondaryMP(phase);					// update secondary MP objective; constraints will not be changed; benders cuts are valid in future upper-level iterations. 

	smp_start = std::chrono::steady_clock::now();
	sub_iteration = 0;
	while (true) {
		sub_iteration++;

		if(is_sub_master_proc()){
		   if(phase == 0){
		      for (int g = 0; g < numGen; g++)
		         var_zg[g].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
		      for (int l = 0; l < numLine; l++)
		         var_zl[l] .set(GRB_CharAttr_VType, GRB_CONTINUOUS);
		      for (int k = 0; k < numStg; k++)
		         var_zs[k].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
		   }
		   else{
		      for (int g = 0; g < numGen; g++)
		         var_zg[g].set(GRB_CharAttr_VType, GRB_BINARY);
		      for (int l = 0; l < numLine; l++)
		         var_zl[l] .set(GRB_CharAttr_VType, GRB_BINARY);
		      for (int k = 0; k < numStg; k++)
		         var_zs[k].set(GRB_CharAttr_VType, GRB_BINARY);
		   }
		}

		auto psp_smp_start = std::chrono::steady_clock::now();
		solveSecondaryMP();
		auto psp_smp_end = std::chrono::steady_clock::now();
		PspSmpTime[color] += double(std::chrono::duration_cast<std::chrono::seconds>(psp_smp_end - psp_smp_start).count());


		MPI_Bcast(CapGen, numGen, MPI_DOUBLE, sub_master, sub_comm);
		MPI_Bcast(CapLine, numLine, MPI_DOUBLE, sub_master, sub_comm);
		MPI_Bcast(CapStg, numStg, MPI_DOUBLE, sub_master, sub_comm);

		// all cores are called in parallel solving all the subproblems
		auto sp_start = std::chrono::steady_clock::now();
		for (int sp = 0; sp < numSPsPerCore; sp++) {	// every CPU core solves a set of subproblems in series;
			solveSubproblem(sp, phase);
		}
		auto sp_end = std::chrono::steady_clock::now();
		rankSolutionTime += double(std::chrono::duration_cast<std::chrono::seconds>(sp_end - sp_start).count());

		if (!is_sub_master_proc()) {
			for (int sp = 0; sp < numSPsPerCore; sp++) {
				const int noSpInNode = get_subproblem_no_in_node(sp);
            MPI_Send(SecDualGen[noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 0, sub_comm);
            MPI_Send(SecDualLine[noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numLine, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 1, sub_comm);
            MPI_Send(SecDualStg[noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numStg, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 2, sub_comm);
            MPI_Send(SecDualDmd[noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numBus, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 3, sub_comm);
            MPI_Send(SecDualRampUp[noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 4, sub_comm);
            MPI_Send(SecDualRampDown[noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 5, sub_comm);
            MPI_Send(SecDualMinDown[noSpInNode], static_cast<int>(MinDownConstrs.size()), MPI_DOUBLE, sub_master, 10000 * noSpInNode + 6, sub_comm);
				MPI_Send(&SPstatus[noSpInNode], 1, MPI_INT, sub_master, 10000 * noSpInNode + 7, sub_comm);
				if ((SPstatus[noSpInNode] == SpStatus::Optimal) || (SPstatus[noSpInNode] == SpStatus::LRoptimal)) {
					MPI_Send(&recFuncVal[noSpInNode], 1, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 8, sub_comm);
					MPI_Send(&SPofv[noSpInNode], 1, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 9, sub_comm);
				}
			}
		}
		else {
			for (int i = 0; i < numSPsPerNode; i++) {
				const int sender = int(i / numSPsPerCore);
				if (sender == sub_master) {
					continue;
				}
            MPI_Recv(SecDualGen[i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, sender, 10000 * i + 0, sub_comm, MPI_STATUS_IGNORE);
            MPI_Recv(SecDualLine[i], numDmdScenariosPerSubproblem * ucPeriods * numLine, MPI_DOUBLE, sender, 10000 * i + 1, sub_comm, MPI_STATUS_IGNORE);
            MPI_Recv(SecDualStg[i], numDmdScenariosPerSubproblem * ucPeriods * numStg, MPI_DOUBLE, sender, 10000 * i + 2, sub_comm, MPI_STATUS_IGNORE);
            MPI_Recv(SecDualDmd[i], numDmdScenariosPerSubproblem * ucPeriods * numBus, MPI_DOUBLE, sender, 10000 * i + 3, sub_comm, MPI_STATUS_IGNORE);
            MPI_Recv(SecDualRampUp[i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, sender, 10000 * i + 4, sub_comm, MPI_STATUS_IGNORE);
            MPI_Recv(SecDualRampDown[i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, sender, 10000 * i + 5, sub_comm, MPI_STATUS_IGNORE);
            MPI_Recv(SecDualMinDown[i], static_cast<int>(MinDownConstrs.size()), MPI_DOUBLE, sender, 10000 * i + 6, sub_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&SPstatus[i], 1, MPI_INT, sender, 10000 * i + 7, sub_comm, MPI_STATUS_IGNORE);
				if ((SPstatus[i] == SpStatus::Optimal) || (SPstatus[i] == SpStatus::LRoptimal)) {
					MPI_Recv(&recFuncVal[i], 1, MPI_DOUBLE, sender, 10000 * i + 8, sub_comm, MPI_STATUS_IGNORE);
					MPI_Recv(&SPofv[i], 1, MPI_DOUBLE, sender, 10000 * i + 9, sub_comm, MPI_STATUS_IGNORE);
				}
			}

			auto smp_temp_start = std::chrono::steady_clock::now();

			getSecondaryMpBounds();

			addNestBendersCuts(phase);

			auto _end = std::chrono::steady_clock::now();
			const double smp_gap = fabs((SmpUB - SmpLB) / SmpUB);
			const double test_time = double(std::chrono::duration_cast<std::chrono::seconds>(_end - smp_start).count());
			rankSolutionTime += double(std::chrono::duration_cast<std::chrono::seconds>(_end - smp_temp_start).count());

#ifdef TEST_MODE
			if (color == 0) {
				std::cout << "Node " << color << " gap = (" << SmpUB << "-" << SmpLB << ")/" << SmpUB << "=" << smp_gap << std::endl;
				char buf[100];
				std::sprintf(buf, "%s/smp-%d-%d.lp", output_directory.c_str(), color, sub_iteration);
				SMPmodel.write(buf);
			}
#endif
			bool break_loop = false;
			if ((smp_gap <= paramReg->max_smp_gap) && ((smp_gap <= paramReg->smp_gap) || (test_time >= paramReg->smp_timlmt))) break_loop = true;
			if ((smp_gap <= paramReg->max_smp_gap) && (phase > 1) && (smp_gap < gap / 5.0)) break_loop = true;

			if (break_loop)
			{
				sub_message = 1;
				SmpGap[color] = smp_gap;
            // PspTime[color] = test_time;
			}
		}
/* ******************************* */
		MPI_Bcast(&sub_message, 1, MPI_INT, sub_master, sub_comm);
		if (sub_message == 1) break;
	}

	return Result::OPTIMAL;
}

void Model::createSecondaryMP()
{
	if (!is_sub_master_proc()) return;
	if(is_master_proc()) PRINT_SUBSECTION("create Secondary Master Problem");

	char buf[100];
	const auto &node = scenarioNodes[color];
	SMPmodel.set(GRB_DoubleParam_MIPGap, paramReg->smp_gap);	// the MIP gap tolerance will be reduced once SP's are feasible
	SMPmodel.set(GRB_DoubleParam_TimeLimit, paramReg->smp_timlmt);
#ifdef SolverStreamOff
	SMPmodel.set(GRB_IntParam_OutputFlag, 0);
#else
	if (!is_master_proc()) {
		SMPmodel.set(GRB_IntParam_OutputFlag, 0);
	}
#endif

	/* ********************** Decision variables ************************** */
	// capacity expansion
    var_zg = SMPmodel.addVars(numGen, GRB_BINARY);
    var_zl = SMPmodel.addVars(numLine, GRB_BINARY);
    var_zs = SMPmodel.addVars(numStg, GRB_BINARY);
#ifdef LPMODEL
    for (int g = 0; g < numGen; g++) {
        std::sprintf(buf, "zg(%d,%d)", color, g);
        var_zg[g].set(GRB_StringAttr_VarName, buf);
    }
    for (int l = 0; l < numLine; l++) {
        std::sprintf(buf, "zl(%d,%d)", color, l);
        var_zl[l].set(GRB_StringAttr_VarName, buf);
    }
    for (int k = 0; k < numStg; k++) {
        std::sprintf(buf, "zs(%d,%d)", color, k);
        var_zs[k].set(GRB_StringAttr_VarName, buf);
    }
#endif // LPMODEL
	// recourse value used in benders decomposition
	var_theta = SMPmodel.addVars(numSPsPerNode);
	sol_theta = new double[numSPsPerNode];
#ifdef LPMODEL	
	for (int i = 0; i < numSPsPerNode; i++) {
		std::sprintf(buf, "theta(%d)", i);
		var_theta[i].set(GRB_StringAttr_VarName, buf);
	}
#endif // LPMODEL
	
	/* ************************* Objective **************************** */
	GRBLinExpr xpr = 0;
    for (int g = 0; g < numGen; g++)
        xpr += var_zg[g] * DualGen[color][g];
    for (int l = 0; l < numLine; l++)
        xpr += var_zl[l] * DualLine[color][l];
    for (int k = 0; k < numStg; k++)
        xpr += var_zs[k] * DualStg[color][k];

	for (int s = 0; s < numSPsPerNode; s++) {
		xpr += var_theta[s];
	}
	SMPmodel.setObjective(xpr, GRB_MINIMIZE);

	/* ************************* Constraints ************************* */
	for (int g = 0; g < numGen; g++) {
		SMPmodel.addConstr(var_zg[g], GRB_LESS_EQUAL, 1 - int(generators[g].exists()));				// at most one expansion along the scenario path
	}
	for (int l = 0; l < numLine; l++) {
		SMPmodel.addConstr(var_zl[l], GRB_LESS_EQUAL, 1 - int(lines[l].exists()));
	}
	for (int k = 0; k < numStg; k++) {
		SMPmodel.addConstr(var_zs[k], GRB_LESS_EQUAL, 1 - int(storages[k].exists()));

		// build energy storage only if there is a renewable energy plant in the same site
		GRBLinExpr capGenInSite = 0;
		for (int g = 0; g < numGen; g++) {
			if ((generators[g].loc != storages[k].loc) || (!generators[g].is_renewable)) continue;
			capGenInSite += int(generators[g].exists()) + var_zg[g];
		}
		SMPmodel.addConstr(var_zs[k], GRB_LESS_EQUAL, capGenInSite);
	}

	// total capacity after expansion must exceed total demand
	GRBLinExpr totalCapacity = 0;
	for (int g = 0; g < numGen; g++) {
		const auto &gen = generators[g];
		totalCapacity += gen.capacity * (int(gen.exists()) + var_zg[g]);
	}
	for (int k = 0; k < numStg; k++) {
		const auto &stg = storages[k];
        totalCapacity += stg.capacity * (int(stg.exists()) + var_zs[k]);
	}

	double totalDemand = 0;
	for (int _s = 0; _s != Season::numSeasons; _s++) {
	   for (int d = 0; d < paramReg->numDmdScenarios; d++) {
			for (int h = 0; h < ucPeriods; h++) {			
				double temp = 0;
				for (int b = 0; b < numBus; b++) {
					temp += buses[b].season_hourly_load_scenarios[get2ndStageIndex(color, _s, d, h)];
				}
				totalDemand = (temp > totalDemand) ? temp : totalDemand;
			}
		}
	}

	SMPmodel.addConstr(totalCapacity, GRB_GREATER_EQUAL, totalDemand);

	// for each location, the total capacity should exceed the demand
	for(int b = 0; b < numBus; b++){
	    const auto &bus = buses[b];
	    double bus_max_hrly_dmd = 0;
	    GRBLinExpr bus_capacity = 0;
        for (int _s = 0; _s != Season::numSeasons; _s++) {
           for (int d = 0; d < paramReg->numDmdScenarios; d++) {
                for (int h = 0; h < ucPeriods; h++) {
                    double temp = bus.season_hourly_load_scenarios[get2ndStageIndex(color, _s, d, h)];
                    if(temp > bus_max_hrly_dmd){
                        bus_max_hrly_dmd = temp;
                    }
                }
            }
        }
        for(int g = 0; g < numGen; g++){
            const auto &gen = generators[g];
            if(gen.loc != b) continue;
            bus_capacity += gen.capacity * (int(gen.exists()) + var_zg[g]);
        }
        for(int l = 0; l < numLine; l++){
            const auto &line = lines[l];
            if((line.from != b) && (line.to != b)) continue;
            bus_capacity += line.capacity * (int(line.exists()) + var_zl[l]);
        }
        for(int k = 0; k < numStg; k++){
            const auto& stg = storages[k];
            if(stg.loc == b){
                bus_capacity += stg.capacity * (int(stg.exists()) + var_zs[k]);
            }
        }
        SMPmodel.addConstr(bus_capacity, GRB_GREATER_EQUAL, bus_max_hrly_dmd);
	}

	// same location only one potential generator to be built
#ifdef RestrictSameSiteExpansion
	for (int b = 0; b < numBus; b++) {
		GRBLinExpr numGensInSite = 0;
		bool hasCandidateGens = false;
		for (int g = 0; g < numGen; g++) {
			const auto &gen = generators[g];
			if (gen.exists() || (gen.loc != b)) continue;
			for (int tau = 0; tau <= node.period; tau++) {
				numGensInSite += var_zg[tau][g];
			}
			hasCandidateGens = true;
		}
		if (hasCandidateGens) {
			SMPmodel.addConstr(numGensInSite, GRB_LESS_EQUAL, 1);
		}
	}
#endif

	// reserve margin requirement
	GRBLinExpr totalGenCap = 0;
	for (int g = 0; g < numGen; g++) {
	   const auto &gen = generators[g];
	   totalGenCap += gen.capacity * (int(gen.exists()) + var_zg[g]);
	}
	double rmReq = 0;
	for (int b = 0; b < numBus; b++) {
	   rmReq += buses[b].peak_loads_in_scenario_nodes[color];
	}
	SMPmodel.addConstr(totalGenCap, GRB_GREATER_EQUAL, rmReq * paramReg->reserveMarginRate);

	// same site renewable and storage facilities
   for(int i = 0; i < numGen - 1; i++){
      const auto &gen1 = generators[i];
      const auto &gen2 = generators[i + 1];
      if((!gen1.is_renewable) || (!gen2.is_renewable)) continue;
      if(Generator::isSameGen(gen1, gen2) && (gen1.capacity == gen2.capacity)){
         SMPmodel.addConstr(var_zg[i + 1], GRB_LESS_EQUAL, var_zg[i]);
      }
   }
   for(int i = 0; i < numStg - 1; i++){
      const auto &stg1 = storages[i];
      const auto &stg2 = storages[i + 1];
      if(Storage::isSameStorage(stg1, stg2)){
         SMPmodel.addConstr(var_zs[i + 1], GRB_LESS_EQUAL, var_zs[i]);
      }
   }
}

void Model::createSubproblem(bool nested)
{
	if (is_master_proc()) PRINT_SUBSECTION("create Subproblem");

	char buf[100];
	const auto &node = scenarioNodes[color];
	SPmodel.set(GRB_IntParam_InfUnbdInfo, 1);		// Additional info for infeasible/unbounded models
	//SPmodel.set(GRB_IntParam_Method, GRB_METHOD_BARRIER);
	//SPmodel.set(GRB_IntParam_Presolve, 0);		// disable presolving in subproblem 	
	//SPmodel.set(GRB_IntParam_Method, 1);		// use dual simplex 
	SPmodel.set(GRB_DoubleParam_MIPGap, paramReg->sp_gap);
	//SPmodel.set(GRB_DoubleParam_TimeLimit, paramReg->sp_tmlmt);

#ifdef SolverStreamOff
	SPmodel.set(GRB_IntParam_OutputFlag, 0);
#else
	if (!is_master_proc()) {
		SPmodel.set(GRB_IntParam_OutputFlag, 0);
	}
#endif

	int index = 0;

	/* ********************** Decision variables ************************** */
	const int numLowerLevelIndices = numDmdScenariosPerSubproblem * ucPeriods;
	varGenStatus = new GRBVar *[numLowerLevelIndices];
	varStartUp = new GRBVar *[numLowerLevelIndices];
	varPowerGen = new GRBVar *[numLowerLevelIndices];
	varPowerFlow = new GRBVar **[numLowerLevelIndices];
	varPowerWithdrawal = new GRBVar *[numLowerLevelIndices];
	varPowerInject = new GRBVar *[numLowerLevelIndices];
	varPowerRemaining = new GRBVar *[numLowerLevelIndices];
	varUnmetDemand = new GRBVar *[numLowerLevelIndices];

	index = 0;
	for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
		for (int h = 0; h < ucPeriods; h++, index++) {
			varGenStatus[index] = SPmodel.addVars(numGen, GRB_BINARY);
			varStartUp[index] = SPmodel.addVars(numGen);
			varPowerGen[index] = SPmodel.addVars(numGen);
			varPowerFlow[index] = new GRBVar *[numLine];
			varPowerWithdrawal[index] = SPmodel.addVars(numStg);
			varPowerInject[index] = SPmodel.addVars(numStg);
			varPowerRemaining[index] = SPmodel.addVars(numStg);
			varUnmetDemand[index] = SPmodel.addVars(numBus);

#ifdef LPMODEL
			for (int g = 0; g < numGen; g++) {
				std::sprintf(buf, "GenState(%d,%d,%d)", d, h, g);
				varGenStatus[index][g].set(GRB_StringAttr_VarName, buf);
				std::sprintf(buf, "startUp(%d,%d,%d)", d, h, g);
				varStartUp[index][g].set(GRB_StringAttr_VarName, buf);
				std::sprintf(buf, "powerGen(%d,%d,%d)", d, h, g);
				varPowerGen[index][g].set(GRB_StringAttr_VarName, buf);
			}
			for (int k = 0; k < numStg; k++) {
				std::sprintf(buf, "powerWithdrawl(%d,%d,%d)", d, h, k);
				varPowerWithdrawal[index][k].set(GRB_StringAttr_VarName, buf);
				std::sprintf(buf, "powerInject(%d,%d,%d)", d, h, k);
				varPowerInject[index][k].set(GRB_StringAttr_VarName, buf);
				std::sprintf(buf, "powerRemain(%d,%d,%d)", d, h, k);
				varPowerRemaining[index][k].set(GRB_StringAttr_VarName, buf);
			}
#endif // LPMODEL

			for (int l = 0; l < numLine; l++) {
				varPowerFlow[index][l] = SPmodel.addVars(Direction::numDirections);
#ifdef LPMODEL	
				std::sprintf(buf, "powerFlow(%d,%d,%d,+)", d, h, l);
				varPowerFlow[index][l][0].set(GRB_StringAttr_VarName, buf);
				std::sprintf(buf, "powerFlow(%d,%d,%d,-)", d, h, l);
				varPowerFlow[index][l][1].set(GRB_StringAttr_VarName, buf);
#endif // LPMODEL				
			}
		}
	}

	/* ************************* Objective **************************** */
	GRBLinExpr spobj_exp = 0;
	index = 0;
	for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
		for (int h = 0; h < ucPeriods; h++, index++) {
			for (int g = 0; g < numGen; g++) {
				const auto &gen = generators[g];
				spobj_exp += node.probability * paramReg->prob_dmd_scenario * (gen.startUpCost * varStartUp[index][g] + gen.fuelCostA * varGenStatus[index][g]) +
				      node.probability * paramReg->prob_dmd_scenario * gen.fuelCostB * gen.fuelPriceRate[color] * varPowerGen[index][g];
			}
			for (int b = 0; b < numBus; b++) {
			   spobj_exp += node.probability * paramReg->prob_dmd_scenario * paramReg->unmetDmdPnlty * varUnmetDemand[index][b];
			}
		}
	}
	SPmodel.setObjective(get_opn_cost_factor(0) * spobj_exp, GRB_MINIMIZE);     // will update in each iteration

	/* ************************* Constraints **************************** */
	GenLimitConstrs = new GRBConstr[numDmdScenariosPerSubproblem * ucPeriods * numGen];
	GenLimit2Constrs = new GRBConstr[numDmdScenariosPerSubproblem * ucPeriods * numGen];
	RampUpConstrs = new GRBConstr[numDmdScenariosPerSubproblem * ucPeriods * numGen];
	RampDownConstrs = new GRBConstr[numDmdScenariosPerSubproblem * ucPeriods * numGen];
	LineLimitConstrs = new GRBConstr[numDmdScenariosPerSubproblem * ucPeriods * numLine];
	StgLimitConstrs = new GRBConstr[numDmdScenariosPerSubproblem * ucPeriods * numStg];
	DemandConstrs = new GRBConstr[numDmdScenariosPerSubproblem * ucPeriods * numBus];
	MinDownConstrs.clear();

	index = 0;
	for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
		for (int h = 0; h < ucPeriods; h++, index++) {
			for (int g = 0; g < numGen; g++) {
				const auto &gen = generators[g];
				GenLimitConstrs[d * ucPeriods * numGen + h * numGen + g] =
					SPmodel.addConstr(gen.capacity * varGenStatus[index][g], GRB_LESS_EQUAL, 0);	//RHS will be updated in benders algorithm	

				GenLimit2Constrs[d * ucPeriods * numGen + h * numGen + g] =
					SPmodel.addConstr(varPowerGen[index][g] - gen.capacity * varGenStatus[index][g], GRB_LESS_EQUAL, 0);	// update varGenStatus coefficient while solving subproblem
				// ramping constraints
            if (h == 0) {
               RampUpConstrs[d * ucPeriods * numGen + h * numGen + g] =
                     SPmodel.addConstr(varPowerGen[index][g], GRB_LESS_EQUAL, generators[g].capacity);
               RampDownConstrs[d * ucPeriods * numGen + h * numGen + g] =
                     SPmodel.addConstr(varPowerGen[index][g], GRB_LESS_EQUAL, generators[g].capacity);
            }
            else {
               RampUpConstrs[d * ucPeriods * numGen + h * numGen + g] =
                     SPmodel.addConstr(varPowerGen[index][g] - varPowerGen[index - 1][g], GRB_LESS_EQUAL, generators[g].rampUp);
               RampDownConstrs[d * ucPeriods * numGen + h * numGen + g] =
                     SPmodel.addConstr(varPowerGen[index - 1][g] - varPowerGen[index][g], GRB_LESS_EQUAL, generators[g].rampDown);
            }

				// action and min up/down constraints
            if (h == 0) {
               SPmodel.addConstr(varStartUp[index][g] - varGenStatus[index][g], GRB_GREATER_EQUAL, 0);
            }
            else
            {
               SPmodel.addConstr(varStartUp[index][g] - varGenStatus[index][g] + varGenStatus[getLowerIndexInSP(d, h - 1)][g], GRB_GREATER_EQUAL, 0);
               for (int tau = h; tau < ucPeriods && tau < h + gen.minUpTime; tau++) {
                  SPmodel.addConstr(varGenStatus[index][g] - varGenStatus[getLowerIndexInSP(d, h - 1)][g], GRB_LESS_EQUAL, varGenStatus[getLowerIndexInSP(d, tau)][g]);
               }
               for (int tau = h; tau < ucPeriods && tau < h + gen.minDownTime; tau++) {
                  MinDownConstrs.push_back(SPmodel.addConstr(varGenStatus[getLowerIndexInSP(d, h - 1)][g] - varGenStatus[index][g] + varGenStatus[getLowerIndexInSP(d, tau)][g],
                                                             GRB_LESS_EQUAL, 1));
               }
            }
			}
			for (int l = 0; l < numLine; l++) {
				const auto &line = lines[l];
				LineLimitConstrs[d * ucPeriods * numLine + h * numLine + l] =
					SPmodel.addConstr(varPowerFlow[index][l][Direction::Forward] + varPowerFlow[index][l][Direction::Backward], GRB_LESS_EQUAL, line.currentCapacity);						//RHS will be updated in benders algorithm											
			}
			for (int k = 0; k < numStg; k++) {
				const auto &storage = storages[k];
				StgLimitConstrs[d * ucPeriods * numStg + h * numStg + k] =
					SPmodel.addConstr(varPowerRemaining[index][k], GRB_LESS_EQUAL, storage.currentCapacity);		//RHS will be updated in benders algorithm

				SPmodel.addConstr(varPowerWithdrawal[index][k] - varPowerRemaining[index][k], GRB_LESS_EQUAL, 0);
				if (h == 0) {
					SPmodel.addConstr(varPowerRemaining[index][k], GRB_EQUAL, 0);
				}
				else {
					SPmodel.addConstr(varPowerRemaining[index][k] - varPowerRemaining[index - 1][k]
						- varPowerInject[index - 1][k] + varPowerWithdrawal[index - 1][k], GRB_LESS_EQUAL, 0);
				}
			}
			for (int b = 0; b < numBus; b++) {
				GRBLinExpr tempExpr = 0;
				for (int l = 0; l < numLine; l++) {
				   if (lines[l].from == b)	tempExpr -= varPowerFlow[index][l][Direction::Forward] - varPowerFlow[index][l][Direction::Backward] * (1 - lines[l].lossRatio);
				   if (lines[l].to == b) tempExpr += varPowerFlow[index][l][Direction::Forward] * (1 - lines[l].lossRatio) - varPowerFlow[index][l][Direction::Backward];
				}
				for (int g = 0; g < numGen; g++) {
					if (generators[g].loc == b) tempExpr += varPowerGen[index][g];
				}
				for (int k = 0; k < numStg; k++) {
					if (storages[k].loc == b) tempExpr += varPowerWithdrawal[index][k] * paramReg->stgCeff - varPowerInject[index][k];
				}
				tempExpr += varUnmetDemand[index][b];
				DemandConstrs[d * ucPeriods * numBus + h * numBus + b] = SPmodel.addConstr(tempExpr, GRB_EQUAL, 0);
			}
		}
	}

	// SolutionInformation
	if (nested) {
		SecDualGen = new double *[numSPsPerNode];
		SecDualRampUp = new double *[numSPsPerNode];
		SecDualRampDown = new double *[numSPsPerNode];
		SecDualLine = new double *[numSPsPerNode];
		SecDualStg = new double *[numSPsPerNode];
		SecDualDmd = new double *[numSPsPerNode];
		SecDualMinDown = new double *[numSPsPerNode];

		for (int i = 0; i < numSPsPerNode; i++) {
			SecDualGen[i] = new double[numDmdScenariosPerSubproblem * ucPeriods * numGen];
			SecDualRampUp[i] = new double[numDmdScenariosPerSubproblem * ucPeriods * numGen];
			SecDualRampDown[i] = new double[numDmdScenariosPerSubproblem * ucPeriods * numGen];
			SecDualLine[i] = new double[numDmdScenariosPerSubproblem * ucPeriods * numLine];
			SecDualStg[i] = new double[numDmdScenariosPerSubproblem * ucPeriods * numStg];
			SecDualDmd[i] = new double[numDmdScenariosPerSubproblem * ucPeriods * numBus];
			SecDualMinDown[i] = new double[static_cast<int>(MinDownConstrs.size())];
		}

		CapGen = new double[numGen];
		CapLine = new double[numLine];
		CapStg = new double[numStg];
		recFuncVal = new double[numSPsPerNode];
		recFuncLB = new double[numSPsPerNode];
		SPstatus = new int[numSPsPerNode];
		SPofv = new double[numSPsPerNode];
	}
	else {
		SecDualGenInNodes = new double **[numNodes];
		SecDualRampUpInNodes = new double **[numNodes];
		SecDualRampDownInNodes = new double **[numNodes];
		SecDualLineInNodes = new double **[numNodes];
		SecDualStgInNodes = new double **[numNodes];
		SecDualDmdInNodes = new double **[numNodes];
		SecDualMinDownInNodes = new double **[numNodes];

		for (int n = 0; n < numNodes; n++) {
			SecDualGenInNodes[n] = new double *[numSPsPerNode];
			SecDualRampUpInNodes[n] = new double *[numSPsPerNode];
			SecDualRampDownInNodes[n] = new double *[numSPsPerNode];
			SecDualLineInNodes[n] = new double *[numSPsPerNode];
			SecDualStgInNodes[n] = new double *[numSPsPerNode];
			SecDualDmdInNodes[n] = new double *[numSPsPerNode];
			SecDualMinDownInNodes[n] = new double *[numSPsPerNode];

			for (int i = 0; i < numSPsPerNode; i++) {
				SecDualGenInNodes[n][i] = new double[numDmdScenariosPerSubproblem * ucPeriods * numGen];
				SecDualRampUpInNodes[n][i] = new double[numDmdScenariosPerSubproblem * ucPeriods * numGen];
				SecDualRampDownInNodes[n][i] = new double[numDmdScenariosPerSubproblem * ucPeriods * numGen];
				SecDualLineInNodes[n][i] = new double[numDmdScenariosPerSubproblem * ucPeriods * numLine];
				SecDualStgInNodes[n][i] = new double[numDmdScenariosPerSubproblem * ucPeriods * numStg];
				SecDualDmdInNodes[n][i] = new double[numDmdScenariosPerSubproblem * ucPeriods * numBus];
				SecDualMinDownInNodes[n][i] = new double[static_cast<int>(MinDownConstrs.size())];
			}
		}

		CapGenInNodes = new double *[numNodes];
		GenStatusInNodes = new double*[numNodes * Season::numSeasons * ucPeriods];
		CapLineInNodes = new double *[numNodes];
		CapStgInNodes = new double *[numNodes];
		recFuncValInNodes = new double *[numNodes];
		recFuncLBInNodes = new double *[numNodes];
		SPstatusInNodes = new int *[numNodes];
		SPofvInNodes = new double *[numNodes];
		
		for (int n = 0; n < numNodes; n++) {
			CapGenInNodes[n] = new double[numGen];
			for (int s = 0; s < Season::numSeasons; s++) {
				for (int h = 0; h < ucPeriods; h++) {
					GenStatusInNodes[get1stStageIndex(n, s, h)] = new double[numGen];
				}
			}
			CapLineInNodes[n] = new double[numLine];
			CapStgInNodes[n] = new double[numStg];
			recFuncValInNodes[n] = new double[numSPsPerNode];
			recFuncLBInNodes[n] = new double[numSPsPerNode];
			SPstatusInNodes[n] = new int[numSPsPerNode];
			SPofvInNodes[n] = new double[numSPsPerNode];
		}
	}
}

Result Model::solveSecondaryMP()
{
	if (!is_sub_master_proc())
		return Result::OPTIMAL;
	auto smp_start_time = std::chrono::steady_clock::now();
	try {
		const auto &node = scenarioNodes[color];
		SMPmodel.optimize();

		if (SMPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
			std::cerr << "The Secondary Master Problem " << color << " is infeasible!" << std::endl;
			exit(EXIT_FAILURE);
		}

		SmpLB = SMPmodel.get(GRB_DoubleAttr_ObjBound) > SmpLB ? SMPmodel.get(GRB_DoubleAttr_ObjBound) : SmpLB;
		SmpNewUB = SMPmodel.get(GRB_DoubleAttr_ObjVal);
		rankSolutionTime += SMPmodel.get(GRB_DoubleAttr_Runtime);

		// get solution from Secondary master problem
		for (int i = 0; i < numSPsPerNode; i++) {
			sol_theta[i] = var_theta[i].get(GRB_DoubleAttr_X);
			SmpNewUB -= sol_theta[i];
		}
        // get columns from subproblems
        for (int g = 0; g < numGen; g++)
            ColExpPlanGenTemp[color][g] = abs(var_zg[g].get(GRB_DoubleAttr_X));
        for (int l = 0; l < numLine; l++)
            ColExpPlanLineTemp[color][l] = abs(var_zl[l].get(GRB_DoubleAttr_X));
        for (int k = 0; k < numStg; k++)
            ColExpPlanStgTemp[color][k] = abs(var_zs[k].get(GRB_DoubleAttr_X));

		for (int g = 0; g < numGen; g++) {
			CapGen[g] = generators[g].currentCapacity + abs(var_zg[g].get(GRB_DoubleAttr_X) * generators[g].capacity);
		}
		for (int l = 0; l < numLine; l++) {
			CapLine[l] = lines[l].currentCapacity + abs(var_zl[l].get(GRB_DoubleAttr_X) * lines[l].capacity);
		}
		for (int k = 0; k < numStg; k++) {
			CapStg[k] = storages[k].currentCapacity + abs(var_zs[k].get(GRB_DoubleAttr_X) * storages[k].capacity);
		}
	}
	catch (GRBException &e) {
		std::cerr << "Error while solving SMP " << color << std::endl;
		exit(EXIT_FAILURE);
	}

	auto smp_end_time = std::chrono::steady_clock::now();
	rankSolutionTime += double(std::chrono::duration_cast<std::chrono::seconds>(smp_end_time - smp_start_time).count());;
	
	return Result::OPTIMAL;
}

void Model::updateSecondaryMP(int phase)
{
	// initialize the SMP with the dual obtained from MP
	if (is_sub_master_proc()) {
        for (int g = 0; g < numGen; g++)
            var_zg[g].set(GRB_DoubleAttr_Obj, DualGen[color][g]);
        for (int l = 0; l < numLine; l++)
            var_zl[l].set(GRB_DoubleAttr_Obj, DualLine[color][l]);
        for (int k = 0; k < numStg; k++)
            var_zs[k].set(GRB_DoubleAttr_Obj, DualStg[color][k]);
	
		if (!paramReg->keep_benders_cuts) {
			for (auto &it : BendersInfCuts) {
				SMPmodel.remove(it);
			}
			BendersInfCuts.clear();
			for (auto &it : BendersOptCuts) {
			   SMPmodel.remove(it);
			}
			BendersOptCuts.clear();
			for (auto &it : LshapedInfCuts) {
				SMPmodel.remove(it);
			}
			LshapedInfCuts.clear();
		}

		for (auto &it : LshapedOptCuts) {
			SMPmodel.remove(it);
		}
		LshapedOptCuts.clear();
	}

	SmpUB = +paramReg->bigM; SmpLB = -paramReg->bigM; SmpNewUB = paramReg->bigM;
	sub_message = 0;
	sub_iteration = 0;
}

Result Model::solveSubproblem(int sp, int phase)
{
	const int noSpInNode = get_subproblem_no_in_node(sp);		// no. subproblem in a node
	const int season = get_season_of_subproblem(noSpInNode);
	const auto& node = scenarioNodes[color];
    /* ************************* Update Subproblem **************************** */
	for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
		for (int h = 0; h < ucPeriods; h++) {
			for (int g = 0; g < numGen; g++) {
				const auto &gen = generators[g];
				GenLimitConstrs[d * ucPeriods * numGen + h * numGen + g].set(GRB_DoubleAttr_RHS, CapGen[g]);
				SPmodel.chgCoeff(GenLimit2Constrs[d * ucPeriods * numGen + h * numGen + g], varGenStatus[getLowerIndexInSP(d, h)][g],
					gen.is_renewable ? (-gen.supply_scenarios[getLowerIndex(color, season, get_dmd_scenario(d, noSpInNode), h)]) : (-gen.capacity));
			}
			for (int l = 0; l < numLine; l++) {
				LineLimitConstrs[d * ucPeriods * numLine + h * numLine + l].set(GRB_DoubleAttr_RHS, CapLine[l]);
			}
			for (int k = 0; k < numStg; k++) {
				StgLimitConstrs[d * ucPeriods * numStg + h * numStg + k].set(GRB_DoubleAttr_RHS, CapStg[k]);
			}
			for (int b = 0; b < numBus; b++) {
				DemandConstrs[d * ucPeriods * numBus + h * numBus + b].set(GRB_DoubleAttr_RHS, 
					buses[b].season_hourly_load_scenarios[get2ndStageIndex(color, season, get_dmd_scenario(d, noSpInNode), h)]);
			}
		}
	}
    // objective
    GRBLinExpr spobj_exp = 0;
    int index = 0;
    for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
        for (int h = 0; h < ucPeriods; h++, index++) {
            for (int g = 0; g < numGen; g++) {
                const auto &gen = generators[g];
                spobj_exp += node.probability * paramReg->prob_dmd_scenario * (gen.startUpCost * varStartUp[index][g] + gen.fuelCostA * varGenStatus[index][g]) +
                      node.probability * paramReg->prob_dmd_scenario * gen.fuelCostB * gen.fuelPriceRate[color] * varPowerGen[index][g];
            }
            for (int b = 0; b < numBus; b++) {
               spobj_exp += node.probability * paramReg->prob_dmd_scenario * paramReg->unmetDmdPnlty * varUnmetDemand[index][b];
            }
        }
    }
    SPmodel.setObjective(get_opn_cost_factor(season) * spobj_exp, GRB_MINIMIZE);

	try {
      /* ************************** Solve SP linear relaxation ************************************** */
      for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
         for (int h = 0; h < ucPeriods; h++) {
            for (int g = 0; g < numGen; g++) {
               varGenStatus[getLowerIndexInSP(d, h)][g].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
               varGenStatus[getLowerIndexInSP(d, h)][g].set(GRB_DoubleAttr_UB, GRB_INFINITY);
            }
         }
      }

      SPmodel.optimize();
      if (SPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
         exit(EXIT_FAILURE);
      }
      for (int i = 0; i < numDmdScenariosPerSubproblem * ucPeriods * numGen; i++) {
         SecDualGen[noSpInNode][i] = GenLimitConstrs[i].get(GRB_DoubleAttr_Pi);
         SecDualRampUp[noSpInNode][i] = RampUpConstrs[i].get(GRB_DoubleAttr_Pi);
         SecDualRampDown[noSpInNode][i] = RampDownConstrs[i].get(GRB_DoubleAttr_Pi);
      }
      for (int i = 0; i < numDmdScenariosPerSubproblem * ucPeriods * numLine; i++) {
         SecDualLine[noSpInNode][i] = LineLimitConstrs[i].get(GRB_DoubleAttr_Pi);
      }
      for (int i = 0; i < numDmdScenariosPerSubproblem * ucPeriods * numStg; i++) {
         SecDualStg[noSpInNode][i] = StgLimitConstrs[i].get(GRB_DoubleAttr_Pi);
      }
      for (int i = 0; i < numDmdScenariosPerSubproblem * ucPeriods * numBus; i++) {
         SecDualDmd[noSpInNode][i] = DemandConstrs[i].get(GRB_DoubleAttr_Pi);
      }
      int constr_index = 0;
      for (const auto &it : MinDownConstrs) {
         SecDualMinDown[noSpInNode][constr_index++] = it.get(GRB_DoubleAttr_Pi);
      }

      SPstatus[noSpInNode] = SpStatus::LRoptimal;
      SPofv[noSpInNode] = SPmodel.get(GRB_DoubleAttr_ObjVal);
      recFuncVal[noSpInNode] = SPmodel.get(GRB_DoubleAttr_ObjVal);

		if(phase > 1) {
			/* ************************** Solve SP MIP model ************************************** */
			for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
				for (int h = 0; h < ucPeriods; h++) {
					for (int g = 0; g < numGen; g++) {
					   if(generators[g].is_renewable) continue;
						varGenStatus[getLowerIndexInSP(d, h)][g].set(GRB_CharAttr_VType, GRB_INTEGER);
					}
				}
			}
			SPmodel.optimize();
			if (SPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
			   std::cerr << "Infeasible subproblem " << color << "-" << sub_rank << std::endl;
			   exit(EXIT_FAILURE);
			}

			SPstatus[noSpInNode] = SpStatus::Optimal;
			SPofv[noSpInNode] = SPmodel.get(GRB_DoubleAttr_ObjVal);
			recFuncVal[noSpInNode] = SPmodel.get(GRB_DoubleAttr_ObjVal);
			if (iteration == 0) {
				recFuncLB[noSpInNode] = SPmodel.get(GRB_DoubleAttr_ObjVal);
			}
		}
	}
	catch (GRBException &e) {
		std::cerr << "Error while solving subproblem at node " << color << " season " << season << " #" << noSpInNode << ": " << e.getMessage() << std::endl;
		std::cerr << "The obj of the subproblem = " << SPofv[noSpInNode] << std::endl;
		std::cerr << "The status of the subproblem = " << SPmodel.get(GRB_IntAttr_Status) << std::endl;
		char buf[100];
		std::sprintf(buf, "%s/sp-%d-%d.lp", output_directory.c_str(), season, noSpInNode);
		SPmodel.write(buf);
		exit(EXIT_FAILURE);
	}

	return Result::OPTIMAL;
}

Result Model::getRecFuncLB()
{
	if (is_master_proc()) PRINT_SUBSECTION("Calcuate the lower bound of recourse functions");
	// use maximum possible capacity
	for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
		for (int h = 0; h < ucPeriods; h++) {
			for (int g = 0; g < numGen; g++) {
				GenLimitConstrs[d * ucPeriods * numGen + h * numGen + g].set(GRB_DoubleAttr_RHS, generators[g].capacity);
			}
			for (int l = 0; l < numLine; l++) {
				LineLimitConstrs[d * ucPeriods * numLine + h * numLine + l].set(GRB_DoubleAttr_RHS, lines[l].capacity);
			}
			for (int k = 0; k < numStg; k++) {
				StgLimitConstrs[d * ucPeriods * numStg + h * numStg + k].set(GRB_DoubleAttr_RHS, storages[k].capacity);
			}
		}
	}

	SPmodel.optimize();
	if (SPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
		std::cerr << "The subproblem " << sub_rank << " in scenario node " << color << "is infeasible!" << std::endl;
		return Result::INFEASIBLE;
	}

	recFuncLB[sub_rank] = SPmodel.get(GRB_DoubleAttr_ObjBound);

	MPI_Bcast(&sub_message, 1, MPI_INT, sub_master, sub_comm);
	if (!is_sub_master_proc()) {
		MPI_Send(&recFuncLB[sub_rank], 1, MPI_DOUBLE, sub_master, 100, sub_comm);
	}
	else {
		for (int i = 0; i < numSPsPerNode; i++) {
			if (i == sub_master) continue;
			MPI_Recv(&recFuncLB[i], 1, MPI_DOUBLE, i, 100, sub_comm, MPI_STATUS_IGNORE);
		}
	}

	return Result::OPTIMAL;
}

void Model::getSecondaryMpBounds()
{
	for (int i = 0; i < numSPsPerNode; i++) {
		SmpNewUB += SPofv[i];
	}

	ReducedCost[color] = SmpLB - DualConvex[color];											// reduced cost
	//update the expansion plan column if a better UB is found
	if (SmpNewUB < SmpUB) {
		SmpUB = SmpNewUB;
		SmpOFV[color] = SmpLB;		
		LambdaCoef[color] = SmpUB;			// get columns from subproblems
        for (int g = 0; g < numGen; g++)
        {
            ColExpPlanGen[color][g] = ColExpPlanGenTemp[color][g];
            LambdaCoef[color] -= ColExpPlanGen[color][g] * DualGen[color][g];
        }
        for (int l = 0; l < numLine; l++)
        {
            ColExpPlanLine[color][l] = ColExpPlanLineTemp[color][l];
            LambdaCoef[color] -= ColExpPlanLine[color][l] * DualLine[color][l];
        }
        for (int k = 0; k < numStg; k++)
        {
            ColExpPlanStg[color][k] = ColExpPlanStgTemp[color][k];
            LambdaCoef[color] -= ColExpPlanStg[color][k] * DualStg[color][k];
        }
	}
}

void Model::addNestBendersCuts(int phase)
{
    zEqOne = 0;
    combnXpr = 0;
    for (int g = 0; g < numGen; g++) {
        if (var_zg[g].get(GRB_DoubleAttr_X) > 0.5) {
            zEqOne++;
            combnXpr += var_zg[g];
        }
        else {
            combnXpr -= var_zg[g];
        }
    }
    for (int l = 0; l < numLine; l++) {
        if (var_zl[l].get(GRB_DoubleAttr_X) > 0.5) {
            zEqOne++;
            combnXpr += var_zl[l];
        }
        else {
            combnXpr -= var_zl[l];
        }
    }
    for (int k = 0; k < numStg; k++) {
        if (var_zs[k].get(GRB_DoubleAttr_X) > 0.5) {
            zEqOne++;
            combnXpr += var_zs[k];
        }
        else {
            combnXpr -= var_zs[k];
        }
    }

	//add benders cuts to SMP
   for (int i = 0; i < numSPsPerNode; i++) {
      const int season = get_season_of_subproblem(i);
      GRBLinExpr cutRhs = 0;

      try {
         for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
            for (int h = 0; h < ucPeriods; h++) {
               for (int g = 0; g < numGen; g++) {
                  const auto &gen = generators[g];
                  cutRhs += SecDualGen[i][d * ucPeriods * numGen + h * numGen + g] * (gen.currentCapacity + var_zg[g] * gen.capacity);

                  cutRhs += SecDualRampUp[i][d * ucPeriods * numGen + h * numGen + g] * ((h == 0) ? gen.capacity : gen.rampUp);
                  cutRhs += SecDualRampDown[i][d * ucPeriods * numGen + h * numGen + g] * ((h == 0) ? gen.capacity : gen.rampDown);
               }
               for (int l = 0; l < numLine; l++) {
                  cutRhs += SecDualLine[i][d * ucPeriods * numLine + h * numLine + l] * (lines[l].currentCapacity + var_zl[l] * lines[l].capacity);
               }
               for (int k = 0; k < numStg; k++) {
                  cutRhs += SecDualStg[i][d * ucPeriods * numStg + h * numStg + k] * (storages[k].currentCapacity + var_zs[k] * storages[k].capacity);
               }
               for (int b = 0; b < numBus; b++) {
                  cutRhs += SecDualDmd[i][d * ucPeriods * numBus + h * numBus + b] *
                     buses[b].season_hourly_load_scenarios[get2ndStageIndex(color, season, get_dmd_scenario(d, i), h)];
               }
            }
         }

         int constr_index = 0;
         for (const auto &it : MinDownConstrs) {
            cutRhs += SecDualMinDown[i][constr_index++];	// RHS of min down constraints are always 1
         }

         BendersOptCuts.push_back(SMPmodel.addConstr(var_theta[i], GRB_GREATER_EQUAL, cutRhs, "BendersOptCut"));
         if(phase > 1){
            LshapedOptCuts.push_back(SMPmodel.addConstr(var_theta[i], GRB_GREATER_EQUAL,
                                                        (recFuncVal[i] - recFuncLB[i]) * (combnXpr - zEqOne + 1) +
                                                        recFuncLB[i], "IntegerLshapedOptCut"));
         }
      }
      catch (GRBException &e) {
         std::cerr << e.getMessage() << "---" << e.getErrorCode() << std::endl;
      }
      catch (const std::out_of_range& oor) {
         std::cerr << "Out of Range error: " << oor.what() << '\n';
      }
   }
}

void Model::get_initial_columns()
{
	if (is_master_proc())
	{
		PRINT_SUBSECTION("get initial capacity expansion columns");
	}

	for (int n = 0; n < numNodes; n++) {
		for (int g = 0; g < numGen; g++) DualGen[n][g] = 0;
		for (int l = 0; l < numLine; l++) DualLine[n][l] = 0;
		for (int k = 0; k < numStg; k++) DualStg[n][k] = 0;
		DualConvex[n] = 0;
	}

	iteration = 0;

	updateSecondaryMP();

	/* Fix expansion variables by adding temporary constraints */
	std::vector<GRBConstr> tempConstrs;

	if (is_sub_master_proc()) {
		const auto &node = scenarioNodes[color];

#ifdef RestrictSameSiteExpansion
		// same location only one potential generator to be built
		for (int b = 0; b < numBus; b++) {
			GRBLinExpr numGensInSite = 0;
			int max_cap_gen = -1;
			double max_cap = 0;
			for (int g = 0; g < numGen; g++) {
				if (generators[g].exists() || (generators[g].loc != b)) continue;
				if (generators[g].capacity > max_cap) {
					max_cap = generators[g].capacity;
					max_cap_gen = g;
				}
			}

			for (int g = 0; g < numGen; g++) {
				if (generators[g].loc != b) continue;
				if (!generators[g].exists()) {
					for (int tau = 0; tau <= node.period; tau++) {
						tempConstrs.push_back(SMPmodel.addConstr(var_zg[tau][g], GRB_EQUAL, int((tau == 0) && (g == max_cap_gen))));
					}
				}
			}
		}
#else
        for (int g = 0; g < numGen; g++) {
            tempConstrs.push_back(SMPmodel.addConstr(var_zg[g], GRB_EQUAL, 1 - int(generators[g].exists())));
        }
#endif

        for (int l = 0; l < numLine; l++) {
            tempConstrs.push_back(SMPmodel.addConstr(var_zl[l], GRB_EQUAL, 1 - int(lines[l].exists())));
        }
        for (int k = 0; k < numStg; k++) {
            tempConstrs.push_back(SMPmodel.addConstr(var_zs[k], GRB_EQUAL, 1 - int(storages[k].exists())));
        }
	}

	solveSecondaryMP();

	MPI_Bcast(CapGen, numGen, MPI_DOUBLE, sub_master, sub_comm);
	MPI_Bcast(CapLine, numLine, MPI_DOUBLE, sub_master, sub_comm);
	MPI_Bcast(CapStg, numStg, MPI_DOUBLE, sub_master, sub_comm);

	if (is_master_proc()) PRINT_SUBSECTION("solve the problem with fixed capacity expansion");
	/* ******** Solve subproblem and Get a feasible expansion plan as new columns ******** */
	auto sp_start = std::chrono::steady_clock::now();
    for (int sp = 0; sp < numSPsPerCore; sp++) {	// every CPU core solves a set of subproblems in series;
		if ((solveSubproblem(sp) == Result::INFEASIBLE))
		{
		   std::string modelFile = output_directory + "/SP-" + std::to_string(color) + "-" + std::to_string(sp) + ".lp";
		   SPmodel.write(modelFile);

			std::cerr << "Subproblem in node " << color << " No." << get_subproblem_no_in_node(sp) << " is infeasible with max capacity expansion in phase 1" << std::endl;
			exit(EXIT_FAILURE);
		}

#ifdef PhaseII
        if (solveSubproblem(sp, 2) == Result::INFEASIBLE) {
		   std::string modelFile = output_directory + "/SP-" + std::to_string(color) + "-" + std::to_string(sp) + ".lp";
		   SPmodel.write(modelFile);

			std::cerr << "Subproblem in node " << color << " No." << get_subproblem_no_in_node(sp) << " is infeasible with max capacity expansion in phase 2" << std::endl;
			exit(EXIT_FAILURE);
		}
#endif		
    }
    auto sp_end = std::chrono::steady_clock::now();
    rankSolutionTime += double(std::chrono::duration_cast<std::chrono::seconds>(sp_end - sp_start).count());

    if(!is_sub_master_proc()) {
        for (int sp = 0; sp < numSPsPerCore; sp++) {
            const int noSpInNode = get_subproblem_no_in_node(sp);
            MPI_Send(&SPstatus[noSpInNode], 1, MPI_INT, sub_master, 10000 * noSpInNode + 7, sub_comm);
            if ((SPstatus[noSpInNode] == SpStatus::Optimal) || (SPstatus[noSpInNode] == SpStatus::LRoptimal)) {
                MPI_Send(&recFuncVal[noSpInNode], 1, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 8, sub_comm);
                MPI_Send(&SPofv[noSpInNode], 1, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 9, sub_comm);
				MPI_Send(&recFuncLB[noSpInNode], 1, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 6, sub_comm);
            }
        }
    }
    else {
        for (int i = 0; i < numSPsPerNode; i++) {
            const int sender = int(i / numSPsPerCore);
            if (sender == sub_master) {
                continue;
            }
            MPI_Recv(&SPstatus[i], 1, MPI_INT, sender, 10000 * i + 7, sub_comm, MPI_STATUS_IGNORE);
            if ((SPstatus[i] == SpStatus::Optimal) || (SPstatus[i] == SpStatus::LRoptimal)) {
                MPI_Recv(&recFuncVal[i], 1, MPI_DOUBLE, sender, 10000 * i + 8, sub_comm, MPI_STATUS_IGNORE);
                MPI_Recv(&SPofv[i], 1, MPI_DOUBLE, sender, 10000 * i + 9, sub_comm, MPI_STATUS_IGNORE);
				MPI_Recv(&recFuncLB[i], 1, MPI_DOUBLE, sender, 10000 * i + 6, sub_comm, MPI_STATUS_IGNORE);
            }
        }

        getSecondaryMpBounds();
    }

	addColumns(true);

	/* Remove temp constraints */
	if (is_sub_master_proc()){
		for (auto it : tempConstrs) {
			SMPmodel.remove(it);
		}
	}
	tempConstrs.clear();
}

void Model::get_second_columns()
{
   if (is_master_proc())
   {
      PRINT_SUBSECTION("get second capacity expansion columns");
   }
   iteration = 1;
   /* ********************** set z coefficient as expansion cost ************************** */
   for (int n = 0; n < numNodes; n++) {
      for (int g = 0; g < numGen; g++) DualGen[n][g] = generators[g].cost[n];
      for (int l = 0; l < numLine; l++) DualLine[n][l] = lines[l].cost[n];
      for (int k = 0; k < numStg; k++) DualStg[n][k] = storages[k].cost[n];
      DualConvex[n] = 0;
   }
   /* ********************** solve pricing sub problems with benders decomposition ************************** */
   solvePrimarySubproblem(1);

#ifdef PhaseII
   solvePrimarySubproblem(2);
#endif

   addColumns(true);
}

void Model::getMoreColumns(int phase)
{
    int column = 0;
    const auto &node = scenarioNodes[color];
    std::vector<GRBConstr> tempConstrs;
    // initialize and reset column objectives
    if(is_sub_master_proc()) {
        for (int n = 0; n < numNodes; n++) {
            for (int i = 0; i < NUM_COLUMNS; i++) {
                col_obj[n][i] = -1; // obj is negative value if fewer than 10 solutions are found
            }
        }
    }
    try {
#ifndef DIRECT_SOLVE_MULTIPLE_COLUMNS
        sub_message = 0;
        while (true) {
            // solve smp to get a column
            if (is_sub_master_proc()) {
                SMPmodel.optimize();

                if (SMPmodel.get(GRB_IntAttr_Status) != GRB_OPTIMAL) {
                    sub_message = 1;
                }
                else {
                    for (int g = 0; g < numGen; g++) {
                        CapGen[g] = generators[g].currentCapacity + abs(var_zg[g].get(GRB_DoubleAttr_X) * generators[g].capacity);
                    }
                    for (int l = 0; l < numLine; l++) {
                        CapLine[l] = lines[l].currentCapacity + abs(var_zl[l].get(GRB_DoubleAttr_X) * lines[l].capacity);
                    }
                    for (int k = 0; k < numStg; k++) {
                        CapStg[k] = storages[k].currentCapacity + abs(var_zs[k].get(GRB_DoubleAttr_X) * storages[k].capacity);
                    }

                    // get columns
                    col_obj[color][column] = SMPmodel.get(GRB_DoubleAttr_ObjVal);
                    // get solution from Secondary master problem
                    for (int i = 0; i < numSPsPerNode; i++) {
                        col_obj[color][column] -= var_theta[i].get(GRB_DoubleAttr_X);
                    }
                    for (int g = 0; g < numGen; g++) {
                        col_gen_expansion[color][column][g] = double(static_cast<int>(var_zg[g].get(GRB_DoubleAttr_X) > 0.5));
                    }
                    for (int l = 0; l < numLine; l++) {
                        col_line_expansion[color][column][l] = double(static_cast<int>(var_zl[l].get(GRB_DoubleAttr_X) > 0.5));
                    }
                    for (int k = 0; k < numStg; k++) {
                        col_stg_expansion[color][column][k] = double(static_cast<int>(var_zs[k].get(GRB_DoubleAttr_X) > 0.5));
                    }

                    // get combinatorial cut info
                    {
                        zEqOne = 0;
                        combnXpr = 0;
                        for (int g = 0; g < numGen; g++) {
                            if (var_zg[g].get(GRB_DoubleAttr_X) > 0.5) {
                                zEqOne++;
                                combnXpr += var_zg[g];
                            } else {
                                combnXpr -= var_zg[g];
                            }
                        }
                        for (int l = 0; l < numLine; l++) {
                            if (var_zl[l].get(GRB_DoubleAttr_X) > 0.5) {
                                zEqOne++;
                                combnXpr += var_zl[l];
                            } else {
                                combnXpr -= var_zl[l];
                            }
                        }
                        for (int k = 0; k < numStg; k++) {
                            if (var_zs[k].get(GRB_DoubleAttr_X) > 0.5) {
                                zEqOne++;
                                combnXpr += var_zs[k];
                            } else {
                                combnXpr -= var_zs[k];
                            }
                        }
                    }
                }
            }
            MPI_Bcast(&sub_message, 1, MPI_INT, sub_master, sub_comm);
            if (sub_message == 1) break;

            MPI_Bcast(CapGen, numGen, MPI_DOUBLE, sub_master, sub_comm);
            MPI_Bcast(CapLine, numLine, MPI_DOUBLE, sub_master, sub_comm);
            MPI_Bcast(CapStg, numStg, MPI_DOUBLE, sub_master, sub_comm);

            // all cores are called in parallel solving each subproblem
            for (int sp = 0; sp < numSPsPerCore; sp++) {    // every CPU core solves a set of subproblems in series;
                solveSubproblem(sp, phase);
            }

            if (!is_sub_master_proc()) {
                for (int sp = 0; sp < numSPsPerCore; sp++) {
                    const int noSpInNode = get_subproblem_no_in_node(sp);
                    MPI_Send(SecDualGen[noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE,
                             sub_master, 10000 * noSpInNode + 0, sub_comm);
                    MPI_Send(SecDualLine[noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numLine, MPI_DOUBLE,
                             sub_master, 10000 * noSpInNode + 1, sub_comm);
                    MPI_Send(SecDualStg[noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numStg, MPI_DOUBLE,
                             sub_master, 10000 * noSpInNode + 2, sub_comm);
                    MPI_Send(SecDualDmd[noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numBus, MPI_DOUBLE,
                             sub_master, 10000 * noSpInNode + 3, sub_comm);
                    MPI_Send(SecDualRampUp[noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE,
                             sub_master, 10000 * noSpInNode + 4, sub_comm);
                    MPI_Send(SecDualRampDown[noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE,
                             sub_master, 10000 * noSpInNode + 5, sub_comm);
                    MPI_Send(SecDualMinDown[noSpInNode], static_cast<int>(MinDownConstrs.size()), MPI_DOUBLE,
                             sub_master, 10000 * noSpInNode + 6, sub_comm);
                    MPI_Send(&SPstatus[noSpInNode], 1, MPI_INT, sub_master, 10000 * noSpInNode + 7, sub_comm);
                    if ((SPstatus[noSpInNode] == SpStatus::Optimal) || (SPstatus[noSpInNode] == SpStatus::LRoptimal)) {
                        MPI_Send(&recFuncVal[noSpInNode], 1, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 8, sub_comm);
                        MPI_Send(&SPofv[noSpInNode], 1, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 9, sub_comm);
                    }
                }
            }
            else {
                for (int i = 0; i < numSPsPerNode; i++) {
                    const int sender = int(i / numSPsPerCore);
                    if (sender == sub_master) {
                        continue;
                    }
                    MPI_Recv(SecDualGen[i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, sender,
                             10000 * i + 0, sub_comm, MPI_STATUS_IGNORE);
                    MPI_Recv(SecDualLine[i], numDmdScenariosPerSubproblem * ucPeriods * numLine, MPI_DOUBLE, sender,
                             10000 * i + 1, sub_comm, MPI_STATUS_IGNORE);
                    MPI_Recv(SecDualStg[i], numDmdScenariosPerSubproblem * ucPeriods * numStg, MPI_DOUBLE, sender,
                             10000 * i + 2, sub_comm, MPI_STATUS_IGNORE);
                    MPI_Recv(SecDualDmd[i], numDmdScenariosPerSubproblem * ucPeriods * numBus, MPI_DOUBLE, sender,
                             10000 * i + 3, sub_comm, MPI_STATUS_IGNORE);
                    MPI_Recv(SecDualRampUp[i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, sender,
                             10000 * i + 4, sub_comm, MPI_STATUS_IGNORE);
                    MPI_Recv(SecDualRampDown[i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, sender,
                             10000 * i + 5, sub_comm, MPI_STATUS_IGNORE);
                    MPI_Recv(SecDualMinDown[i], static_cast<int>(MinDownConstrs.size()), MPI_DOUBLE, sender,
                             10000 * i + 6, sub_comm, MPI_STATUS_IGNORE);
                    MPI_Recv(&SPstatus[i], 1, MPI_INT, sender, 10000 * i + 7, sub_comm, MPI_STATUS_IGNORE);
                    if ((SPstatus[i] == SpStatus::Optimal) || (SPstatus[i] == SpStatus::LRoptimal)) {
                        MPI_Recv(&recFuncVal[i], 1, MPI_DOUBLE, sender, 10000 * i + 8, sub_comm, MPI_STATUS_IGNORE);
                        MPI_Recv(&SPofv[i], 1, MPI_DOUBLE, sender, 10000 * i + 9, sub_comm, MPI_STATUS_IGNORE);
                    }
                }

                for (int i = 0; i < numSPsPerNode; i++) {
                    if ((SPstatus[i] == SpStatus::LRinfeasible) || (SPstatus[i] == SpStatus::MIPinfeasible)) {
                        col_obj[color][column] = -1;
                        break;
                    }
                    col_obj[color][column] += SPofv[i];
                }
                if (col_obj[color][column] != -1) {
                    column++;
                }

                addNestBendersCuts(phase);

                // add combinatorial cut to remove the current solution
                tempConstrs.push_back(SMPmodel.addConstr(combnXpr, GRB_LESS_EQUAL, zEqOne - 1));

                auto _end = std::chrono::steady_clock::now();
                const double test_time = double(
                        std::chrono::duration_cast<std::chrono::seconds>(_end - smp_start).count());

                if ((column >= NUM_COLUMNS) || (test_time >= 2 * paramReg->smp_timlmt)) {
                    sub_message = 1;
                   PspTime[color] = test_time;
                }
            }
            MPI_Bcast(&sub_message, 1, MPI_INT, sub_master, sub_comm);
            if (sub_message == 1) break;
        }
#else
        // change parameters in smp so multiple solutions can be retrieved
        SMPmodel.set(GRB_IntParam_PoolSearchMode, 2);
        SMPmodel.set(GRB_IntParam_PoolSolutions, NUM_COLUMNS);

        // solve the smp and get solution pool
        SMPmodel.optimize();

        // get columns by retrieving solutions to SMP
        const auto &node = scenarioNodes[color];
        for(int i = 0; i < NUM_COLUMNS && i < SMPmodel.get(GRB_IntAttr_SolCount); i++) {
            SMPmodel.set(GRB_IntParam_SolutionNumber, i);
            col_obj[color][i] = SMPmodel.get(GRB_DoubleAttr_PoolObjVal);
            for (int tau = 0; tau <= node.period; tau++) {
                for (int g = 0; g < numGen; g++) {
                    col_gen_expansion[color][i][getMpConstrIndex(tau, g)] =
                            static_cast<int>(var_zs[tau][k].get(var_zg[tau][g].get(GRB_DoubleAttr_Xn) > 0.5);
                }
                for (int l = 0; l < numLine; l++) {
                    col_line_expansion[color][i][getMpConstrIndex(tau, l, FacilityType::LINE)] =
                            static_cast<int>(var_zs[tau][k].get(var_zl[tau][l].get(GRB_DoubleAttr_Xn) > 0.5);
                }
                for (int k = 0; k < numStg; k++) {
                    col_stg_expansion[color][i][getMpConstrIndex(tau, k, FacilityType::STORAGE)] =
                            static_cast<int>(var_zs[tau][k].get(var_zs[tau][k].get(GRB_DoubleAttr_Xn) > 0.5);
                }
            }
        }
        // reset parameter settings
        SMPmodel.set(GRB_IntParam_PoolSearchMode, 0);
#endif
    }catch (GRBException &e) {
        std::cerr << "Error while obtaining multiple columns:" << e.getMessage() << "---" << e.getErrorCode() << std::endl;
        exit(EXIT_FAILURE);
    }catch (const std::out_of_range& oor) {
        std::cerr << "Out of Range error: " << oor.what() << '\n';
    }catch(const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    // send columns to master node and add to RMP
    if ((!is_master_proc()) && (is_sub_master_proc())) {
        MPI_Send(col_obj[color], NUM_COLUMNS, MPI_DOUBLE, master, 0, MPI_COMM_WORLD);
        for(int i = 0; i < NUM_COLUMNS; i++) {
            MPI_Send(col_gen_expansion[color][i], numGen, MPI_DOUBLE, master, 100 + i, MPI_COMM_WORLD);
            MPI_Send(col_line_expansion[color][i], numLine, MPI_DOUBLE, master, 200 + i,
                     MPI_COMM_WORLD);
            MPI_Send(col_stg_expansion[color][i], numStg, MPI_DOUBLE, master, 300 + i,
                     MPI_COMM_WORLD);
        }
    }
    else if(is_master_proc()){
        for(int n = 0; n < numNodes; n++){
            if(n == master) continue;
            const int sender = n * numCoresPerNode;
            MPI_Recv(col_obj[n], NUM_COLUMNS, MPI_DOUBLE, sender, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
            for(int i = 0; i < NUM_COLUMNS; i++) {
                MPI_Recv(col_gen_expansion[n][i], numGen, MPI_DOUBLE, sender, 100 + i, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                MPI_Recv(col_line_expansion[n][i], numLine, MPI_DOUBLE,
                         sender, 200 + i, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
                MPI_Recv(col_stg_expansion[n][i], numStg, MPI_DOUBLE,
                         sender, 300 + i, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
            }
        }
    }

    // remove temp constraints
    if(is_sub_master_proc()) {
        for (auto it : tempConstrs) {
            SMPmodel.remove(it);
        }
        tempConstrs.clear();
    }
}

void Model::addMoreColumns()
{
	if (!is_master_proc()) return;
	// add columns
    for (int n = 0; n < numNodes; n++) {
        for(int i = 0; i < NUM_COLUMNS; i++) {
            if(col_obj[n][i] < 0) continue;
            columns[n].push_back(MPmodel.addVar(0, GRB_INFINITY, col_obj[n][i], GRB_CONTINUOUS));
            auto &newVar = columns[n].back();
#ifdef LPMODEL
            char name_buf[30];
            std::sprintf(name_buf, "column(%d,%d)", n, static_cast<int>(columns[n].size() - 1));
            newVar.set(GRB_StringAttr_VarName, name_buf);
#endif
            MPmodel.chgCoeff(Convex[n], newVar, 1);
            for (int g = 0; g < numGen; g++) {
                MPmodel.chgCoeff(SplitGen[n][g], newVar, -col_gen_expansion[n][i][g]);
            }
            for (int l = 0; l < numLine; l++) {
                MPmodel.chgCoeff(SplitLine[n][l], newVar, -col_line_expansion[n][i][l]);
            }
            for (int k = 0; k < numStg; k++) {
                MPmodel.chgCoeff(SplitStg[n][k], newVar, -col_stg_expansion[n][i][k]);
            }
        }
    }
    MPmodel.update();

#ifdef TEST_MODE_2
    if(is_master_proc()){
	    MPmodel.update();
        try {
            std::cout << "Multiple Columns" << std::endl;
            /*for(int n = 0; n < numNodes; n++){
                for(int i = 0; i < NUM_COLUMNS; i++) {
                    std::cout << " -- Solution #" << i << std::endl;
                    std::cout << " ---- obj: " << col_obj[n][i] << "\n";
                    for(int g = 0; g < numGen; g++){
                        std::cout << " --- G " << g << ":";
                        std::cout << "\t" << col_gen_expansion[n][i][g];
                        std::cout << std::endl;
                    }
                }
            }*/

            for(int n = 0; n < numNodes; n++) {
                std::cout << "\t Node " << n << " has " << columns[n].size() << " columns: ";
                for(auto column : columns[n]){
                    //std::cout << column.get(GRB_StringAttr_VarName) << "-" << column.get(GRB_DoubleAttr_Obj) << "\t";
                    std::cout << column.get(GRB_StringAttr_VarName) << "\t";
                }
                std::cout << std::endl;
            }
        }
        catch (GRBException &e) {
            std::cerr << e.getMessage() << "---" << e.getErrorCode() << std::endl;
            char buf[100];
            std::sprintf(buf, "%s/MP.lp", output_directory.c_str());
            MPmodel.write(buf);
            exit(EXIT_FAILURE);
        }
    }
#endif
}

