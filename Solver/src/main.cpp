// powinv_revised.cpp : Defines the entry point for the console application.

#include "../inc/model.h"

#define PRINT_SECTION(log) {std::cout << "=============" << log << "==============" << std::endl;}
#define PRINT_SUBSECTION(log) {std::cout << "--" << log << "--" << std::endl;}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	ReadWrite rw;
	Model cepModel;

	/* ******************* Parallel Setup ******************** */
	int _rank = 0, count = 1;
	MPI_Comm_rank(MPI_COMM_WORLD, &_rank);							// get rank: processor id
	cepModel.set_rank(_rank);
	cepModel.set_master(0);
	MPI_Comm_size(MPI_COMM_WORLD, &count);
	cepModel.set_num_ranks(count);

	if(cepModel.is_master_proc()) PRINT_SECTION("A Multistage and Multiscale Approach for Capacity Expansion in Power System");

	/* ******************* Process arguments ***************** */
	if (argc >= 3) {
		const std::string directory = argv[1];
		rw.set_inputDirectory(directory + "/input");
		rw.set_outputDirectory(directory + "/output");
		cepModel.set_output_directory(directory + "/output");
		const std::string alg = argv[2];
		if (alg == "-d") {
			cepModel.set_algorithm(Algorithm::DIRECT);
		}
		else if (alg == "-e") {
			cepModel.set_algorithm(Algorithm::EXPECT);
		}
		else if (alg == "-c") {
			cepModel.set_algorithm(Algorithm::COLGEN);
		}
		else if (alg == "-b") {
			cepModel.set_algorithm(Algorithm::BENDERS);
		}
		else if (alg == "-n") {
			cepModel.set_algorithm(Algorithm::NESTED);
		}
		else if (alg == "-t") {
			cepModel.set_algorithm(Algorithm::TEST);
		}
	}

	/* ***************** Read Parameters Data ************ */
	if (cepModel.is_master_proc()) PRINT_SECTION("Reading Parameters");
   rw.readParameters();

   if (cepModel.is_master_proc()) PRINT_SECTION("Constructing Scenario Tree");
	cepModel.constructScnStructure();		// this has to be done before data generation

	/* ******************* Read Data of Bus System ************************ */
	if (cepModel.is_master_proc()) PRINT_SECTION("Reading System Data");
	rw.readSystemData(cepModel);

	/* ******************* Generate Scenario Data and Broadcast ********** */
	rw.dataGeneration(cepModel);

#ifdef FILEOUTON
	rw.writeGeneratedData(cepModel);
#endif

	/* ******************* Solve the model ****************** */
	auto status = static_cast<Result>(cepModel.optimize());

	/* ******************* Write solution ****************** */
	if (status == Result::OPTIMAL) {
		rw.printSolution(cepModel);
	}

	MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return (status);
}
