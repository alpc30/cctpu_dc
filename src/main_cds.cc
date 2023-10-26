/*
 ============================================================================
 Name        : cds.cc
 Author      : Jing Cui
 Version     : 0.0.0.1
 Copyright   : Your copyright notice
 Functions   : 1. Re-implement the conflict directed search algorithm
 ============================================================================
 */

#include "stnu.h"
#include "cds.h"
#include "reader.h"
//#include "bisearch.h"
//#include "mip_baseline.h"
#include "controllability.h"

CDS::Cds cdru;
CDS::S_PARAM input_param; //input parameters
//CDS::MipBaseline solver_mip;

void print_usage()
{
//	printf("Usage: cds_stnu <problem_file>|<problem_set_path>\n");
	printf("\nUsage: check [OPTIONS]... FILE\n");
	printf("A consistency checking implementation for \n"
			"CCTPU (STNU plus controllable options)"
			"\n\n");

	//Define the objective function
	printf("-o Obj:\n	0 - Maximum delay\n");
//	printf("	20 - Maximum Earliness\n");
//	printf("	2 - Relaxing over-constrained problem\n");
//	printf("	12- Relaxing over-constrained problem without rewards\n");
	printf("	5 - DC of CSTNU\n");


	//Define the constraints
	printf("-v Variable:\n	0 - DC \n");
	printf("	1 - SC\n");

	//Define solver
//	printf("-s Solver:\n	0 - Default(CDS) \n");
//	printf("	5 - Baseline\n");
//	printf("	6 - Binary Search + Checking\n");

	//Choose the way to present congtingent links
	printf("-O Output Path\n");
	printf("-T Runtime Limitation\n");
	printf("-e Extra Parameter\n");
	printf("-c Constraint: (0,1,2,3) - DC, SC, WC, consistent\n");

//	printf("-e extra parameter\n");
	printf("-debug Debug:\n");
	printf("	0 - No Debug Output\n");
	printf("	1 - Key Debug Output\n");
	printf("	2 - All Debug Output\n");

	exit(0);
}

void signal_callback_handler(int signum)
{
	printf("Caught signal %d\n", signum);
	if(input_param.id_solver == 0) // cds
	{
	 // cdru.generateResult(-1);
	 // printf("%s\n", cdru.result.c_str());
	}
	else if(input_param.id_solver == S_BASELINE)
	{
//	  solver_mip.generateResult(input_param, 1);
//	  printf("%s\n",solver_mip.result.c_str());
	}
	  

	exit(signum);
}

int main(int argc, char* argv[])
{
	string problem_file; //input file
	string o_file;  //output file
	string output_path = "./";

	signal( SIGTERM, signal_callback_handler);
	
	//signal( 9, signal_callback_handler);

	if (argc >= 2)
	{
		//read
		//printf("%s\n",argv[argc-1]);
		problem_file = (string) argv[argc - 1];

		bool b_disvar = 0;
		for (int i = 1; i < argc - 1; i++)
			if (strcmp(argv[i], "-c") == 0)
				input_param.id_constrs = atoi(argv[++i]);
			else if (strcmp(argv[i], "-o") == 0)
				input_param.id_obj = atoi(argv[++i]);
			else if (strcmp(argv[i], "-e") == 0)
				input_param.extra = atof(argv[++i]);
			else if (strcmp(argv[i], "-s") == 0)
				input_param.id_solver = atoi(argv[++i]);
			else if (strcmp(argv[i], "-O") == 0)
			{
				o_file = (string) argv[++i];
			}
			else if (strcmp(argv[i], "-d") == 0)
			{
				b_disvar = 1;
			}
			else if (strcmp(argv[i], "-debug") == 0)
			{
				input_param.id_debug = atoi(argv[++i]);
			}
			else if (strcmp(argv[i], "-op") == 0)
			{
				output_path = argv[++i];
			}
			else if (strcmp(argv[i], "-v") == 0)
			{
				input_param.id_variable = atoi(argv[++i]);
			}
			else if (strcmp(argv[i], "-T") == 0)
				input_param.f_runtime = atof(argv[++i]);

		if (o_file[0] == 0)
		{
			o_file = (string) "cds_tmp_out.txt";
		}

		//printf("%s\n",problem_file.c_str());
		CDS::stnu current_stnu;		//problem_file, input_param);
		//current_stnu.readSTNU(problem_file);
		CDS::Reader current_reader(&current_stnu);
		current_stnu.s_info = input_param;
		current_reader.read(problem_file, current_stnu);
		current_stnu.initSTNU(input_param);
		if (input_param.id_debug)
			current_stnu.print_dot("input.dot");

		
		string result;

		if (input_param.id_obj == O_DC_CHECK) // DC check
		{
			if (input_param.id_variable == V_SC) // fixed options
			{
				if(input_param.id_solver == S_CDS_E)
					cdru.CdsCCTPU(current_stnu); //envelope
				else
					cdru.CDRU1(current_stnu); //find a conflict each iteration
			}
			else
				cdru.DCCheckCSTNU1(current_stnu); // dynamic options
			result = cdru.result;
		}
		else //optimisation
		{
			/*if (input_param.id_variable == V_SC) // fixed options
			{
				if (input_param.id_solver == S_BASELINE) //constraint
				{
					
					result = solver_mip.solve(current_stnu);
				}
				else if (input_param.id_solver == S_BISEARCH) //binary search + checking
				{
					BiSearch current_search(current_stnu);
					current_search.init();
					current_search.bisearch();
					result = current_search.result;
				}
				else //conflict-directed
				{
					if(input_param.id_solver == S_CDS_E) // envelope
						cdru.CdsCCTPU(current_stnu) ;
					else
						cdru.CDRU1(current_stnu); //one conflict
					result = cdru.result;
				}
			}
			else // dynamic options
			{
				if (input_param.id_solver == S_BASELINE)
				{
       	  current_stnu.print_dot("before.dot");
				  current_stnu.exp4mip();
				  current_stnu.print_dot("test.dot");
				  //CDS::MipBaseline solver_mip;
				  result = solver_mip.solve(current_stnu);
				  //	printf("baseline solver cannot deal with dynamic options!\n");
				}
				else if (input_param.id_solver == S_BISEARCH)
				{
					BiSearch current_search(current_stnu);
					current_search.init();
					current_search.bisearch();
					result = current_search.result;
				}
				else
				{
					printf(
							"conflict directed solve cannot optimise wiht dynamic options.\n");
				}
			}
*/
		}

		printf("%s\n", result.c_str());

		if (input_param.id_debug)
		{
			FILE *fout;
			//if(o_file)
			fout = fopen(o_file.c_str(), "a");
			fprintf(fout, "%s\n", result.c_str());
			fclose(fout);
		}

	}
	else
	{
		print_usage();
	}

	return 0;
}
