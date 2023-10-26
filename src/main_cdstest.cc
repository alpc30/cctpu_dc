/*
============================================================================
 Name        : cds.cc
 Author      : Jing Cui
 Version     : 0.0.0.1
 Copyright   : Your copyright notice
 Functions   : 1. Re-implement the conflict directed search algorithm
============================================================================
*/

#include <stdio.h>
#include <string>
#include <string.h>
#include <stdlib.h>
#include "stnu.h"
#include "cds.h"
#include "reader.h"
#include "controllability.h"
#include "testgenerator.h"
using namespace std;


void print_usage1()
{
//	printf("Usage: cds_stnu <problem_file>|<problem_set_path>\n");
	printf("\nUsage: cds [OPTION]... FILE\n");
	printf("A conflict directed search algorithm for solving optimization\n"
			"problem(s) of Simple Temporal Network with "
			"Uncertainty.\n\n");

	//Define the objective function
	printf("-o Obj:\n	0 - Maximum delay\n");
//	printf("	20 - Maximum Earliness\n");
	printf("	2 - Relaxing over-constrained problem\n");
	printf("	5 - DC of CSTNU\n");

	//Define the constraints
	printf("-v Variable:\n	0 - DC \n");
	printf("	1 - SC\n");
	//printf("	5 - DC of CSTNU\n");

	//Choose the way to present congtingent links

	printf("-op Output Path\n");
	printf("-T Runtime Limitation\n");
	printf("-e Extra Parameter\n");
	//printf("-c Constraint: (0,1,2,3) - DC, SC, WC, consistent\n");

//	printf("-e extra parameter\n");
	printf("-debug Debug\n");

	exit(0);
}

int main(int argc, char* argv[])
{
	string problem_file; //input file
	string o_file;  //output file
	string output_path ="./";
	CDS::S_PARAM input_param; //input parameters

	if( argc >= 2 )
	{
		//read
		//printf("%s\n",argv[argc-1]);
		problem_file = (string)argv[argc-1];

		bool b_disvar = 0;
		for(int i = 1;i < argc-1;i++)
		if(strcmp(argv[i],"-c")==0)
		  input_param.id_constrs = atoi(argv[++i]);
		else if(strcmp(argv[i],"-o")==0)
		  input_param.id_obj = atoi(argv[++i]);
		else if(strcmp(argv[i],"-e")==0)
		  input_param.extra = atof(argv[++i]);
		else if(strcmp(argv[i],"-O")==0)
		{
		  o_file = (string) argv[++i];
		}
		else if(strcmp(argv[i],"-d") == 0)
		{
		  b_disvar = 1;
		}
		else if(strcmp(argv[i],"-debug") == 0)
		{
		  input_param.id_debug = atoi(argv[++i]);
		}
		else if(strcmp(argv[i],"-op")==0)
		{
		  output_path = argv[++i];
		}
		else if(strcmp(argv[i],"-v") == 0)
		{
		  input_param.id_variable = atoi(argv[++i]);
		}
		else if(strcmp(argv[i],"-T") == 0)
			input_param.f_runtime = atof(argv[++i]);

		if(o_file[0] == 0)
		{
		  o_file = (string)"cds_tmp_out.txt";
		}

		if(input_param.extra >= 1.0)
		{
			input_param.f_runtime/=(input_param.extra+1);
		}

 		for(int i_plan = 0; i_plan < input_param.extra+1; i_plan++)
		{
			  CDS::stnu current_stnu(problem_file, input_param);
			  //current_stnu.readSTNU(problem_file);
			  CDS::Reader current_reader(&current_stnu);

			  current_reader.read(problem_file, current_stnu);
			  current_stnu.initSTNU(input_param);
			  if(input_param.id_debug == DEBUG_T)
				  current_stnu.print_dot("input2.dot");
			  if(input_param.extra)
			  {
				  current_stnu.relaxMission(i_plan);

				  CSTNU::TestGenerator tg = CSTNU::TestGenerator();
				  string res = tg.createTest(current_stnu);
			  }

			  char buffer[10];
			  if(input_param.extra)
				  sprintf(buffer, "_%d_%d.cst",i_plan/3,i_plan%3);
			  else
				  sprintf(buffer, ".cst");
			  string output_file = output_path +"/"+ current_reader.getName(problem_file)
					  +(string)buffer;
			  printf("%s\n",output_file.c_str());
			  current_stnu.print_cstnu(output_file);

			  FILE *fout;
			 // if(o_file)
			  fout = fopen(o_file.c_str(), "a");
			  fclose(fout);
		  }

	}
	else
	{
		print_usage1();
	}


	return 0;
}
