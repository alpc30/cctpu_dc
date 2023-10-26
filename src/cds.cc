/*
 * cds.cpp
 *
 *  Created on: 21/09/2015
 *      Author: jing
 */

#include "cds.h"

namespace CDS
{
long long getSystemTime();

struct T
{
	double lb, ub;
	int id;
};

vector<T> intervals;
bool cmpl(T a, T b)
{
	return a.lb < b.lb;
}

bool cmpu(T a, T b)
{
	return a.ub > b.ub;
}
Cds::Cds()
{
	// TODO Auto-generated constructor stub
	this->v_conflicts.clear();
	n_combines = 0;
	n_leaves = 0;
}

void Cds::setWholeSTNU(const stnu & t)
{
	this->v_conflicts.clear();
	whole_stnu = t;
}

int Cds::DCCheckCSTNU1(const stnu & a)
{
	//usage: the top function of DC checking process for CCTPU
	/*
	 * process: 1. generate STNUs by expanding on variables
	 * 			2. run the CDRU to return conflicts from the candidate
	 * 			3. combine the conflicts together, check DC
	 * 			*/
	whole_stnu = a;
	mlayer = 0;
	tlayer = a.n_vars + 1;
	long long s_time = getSystemTime();

	// the first candidate has assignment '0' for each variable.
	std::vector<pair<int, int> > v_value_dis_var;
	for (int i = 0; i < a.n_vars; i++)
		v_value_dis_var.push_back(make_pair(a.v_dis_var[i].id, 0));

	bool fl = 0; // the flag of exhausting combinations
	v_conflicts.clear();

	// store the sub_solutions with assignment of different length
	std::vector<SubSolution> sol;
	for (int i = 0; i < a.n_vars; i++)
		sol.push_back(SubSolution());

	do
	{
		//freopen("output.log","w",stdout);
		run_time = getSystemTime() - s_time;
		if (run_time / 1000.0 > a.s_info.f_runtime)
		{
			break;
		}
		//stnu res;
		Candidate res = Candidate(whole_stnu, v_value_dis_var);

		SubSolution tmp_sub_solution;
		tmp_sub_solution.v_value_dis_var = v_value_dis_var;
		tmp_sub_solution.v_conflicts.clear();

		// debug: print the stnu into a dot file
		char buffer[100];

		// print the assignment of the current stnu
		if (a.s_info.id_debug >= DEBUG_C)
		{
			printf(ANSI_COLOR_GREEN); // change the color of output
			printf("Solve STNU with assignments: ");
			for (int k = 0; k < v_value_dis_var.size(); k++)
				printf(" V%d=%d	", k, v_value_dis_var[k].second);
			printf(ANSI_COLOR_RESET"\n");
		}

		// get conflict-resolution constraints of the current stnu
		//CDRU2(res, tmp_sub_solution); //debug
		Controllability current_checking = Controllability();

		//if(first)
		tmp_sub_solution.v_conflicts = current_checking.DCConflict(whole_stnu,
				res);
		if (tmp_sub_solution.v_conflicts.size() == 0)
			tmp_sub_solution.b_solved = true;

		n_leaves++;
		if (a.s_info.id_debug >= DEBUG_C)
		{
			printf(ANSI_COLOR_YELLOW); // change the color of output
			printf("solution to combine:\n");
			tmp_sub_solution.print();
			printf(ANSI_COLOR_RESET);
		}

		if (tmp_sub_solution.b_solvable)
		{
			if (tmp_sub_solution.b_solved)
			{
				//res.print_dot("sol.dot");
				run_time = getSystemTime() - s_time;
				generateResult(1);
				return 1; // the current assignment is DC
			}
			//this->simpli
			tmp_sub_solution.simplifyConflicts(whole_stnu);
		}

		std::vector<pair<int, int> > next_value_dis_var;
		next_value_dis_var = getNextAssignment(v_value_dis_var);

		bool fl_stop = false;
		if (tmp_sub_solution.b_solvable)
			for (int i_order = a.n_vars - 1; i_order >= 0; i_order--)
			{

				int i = whole_stnu.order_dv[i_order];
				if (next_value_dis_var.size() > 0 && fl_stop
						&& next_value_dis_var[i].second
								== v_value_dis_var[i].second)
					break;
				if (next_value_dis_var.size() > 0
						&& next_value_dis_var[i].second
								!= v_value_dis_var[i].second)
					fl_stop = true;
				//std::cerr<<i<<std::endl;
				if (tmp_sub_solution.b_solvable == false)
					break;
				//combine2STNU(sol[i_order], tmp_sub_solution, i);
				combineSubSol(sol[i_order], tmp_sub_solution, i);
				n_combines++;
				// keep a record of the layer
				if (mlayer < a.n_vars - i_order)
					mlayer = a.n_vars - i_order + 1;
				// end keep a record of the layer 20160421
				tmp_sub_solution = sol[i_order];
				//tmp_sub_solution.print();
				if (whole_stnu.s_info.id_debug >= DEBUG_C)
				{
					printf(ANSI_COLOR_YELLOW); // change the color of output
					printf("solution in order %d:\n", i_order);
					printf(ANSI_COLOR_RESET);
					sol[i_order].print();
				}
				if (sol[i_order].b_solved == true)
				{
					run_time = getSystemTime() - s_time;
					generateResult(a.n_vars - i_order + 1);
					return a.n_vars - i_order + 1;
					//return true; // the combined result is DC
				}

				if (i_order != a.n_vars - 1)
				{

					sol[i_order + 1] = SubSolution();
				}
			}
		run_time = getSystemTime() - s_time;
		if (run_time / 1000.0 > a.s_info.f_runtime)
		{
			break;
		}
		if (next_value_dis_var.size() == 0)
			break;
		else
			v_value_dis_var = next_value_dis_var;
	} while (fl == 0);
	run_time = getSystemTime() - s_time;
	generateResult(-1);
	return 0;
}

Candidate Cds::CdsCCTPU(const stnu &a)
{
	//use envelope to resolve conflicts
	whole_stnu = a;

	Candidate new_candidate = Candidate(a);
	long long s_time = getSystemTime();

	q_candidate.push(new_candidate);

	while (q_candidate.size() != 0)
	{
		std::vector<Conflict> v_current_conflicts;
		Candidate current_candidate = q_candidate.top();
		int tmp_resolve = current_candidate.v_resolved_conflicts.size();
		whole_stnu.setCtg(current_candidate.value_preference);
		q_candidate.pop();
		run_time = getSystemTime() - s_time;
		if (run_time / 1000.0 > a.s_info.f_runtime)
		{
			break;
		}
		//v_current_conflicts.clear();
		if (a.s_info.id_debug >= DEBUG_T)
		{
			printf("\nDequeue Candidate: size(%d)\n", q_candidate.size());
			if (a.s_info.id_debug == DEBUG_T)
				current_candidate.print(a);
			else
				current_candidate.print(a);
		}

		v_current_conflicts = resolveConflicts(current_candidate);
		if (current_candidate.b_checked == true
				&& v_current_conflicts.size() != 0)
		{
			//continue;
			// if the aim is to check dynamic controllability of CCTPU
			// without making choices dynamically, we don't expand on
			// conflicts.
			if (a.s_info.id_obj == O_DC_CHECK && a.s_info.id_variable == V_SC)
			{
				if (a.s_info.id_debug == DEBUG_C)
				{
					current_candidate.print(a);
				}
				continue;
			}
			if(whole_stnu.s_info.id_obj != O_DC_CHECK)
				this->expandOnConflict(current_candidate, v_current_conflicts);
		}
		else
		{
			if (isComplete(current_candidate) == false)
			{
				this->expandOnVariables(current_candidate);
			}
			else
			{

				Controllability current_checking = Controllability();
				if (a.s_info.id_debug == DEBUG_T)
				{
					printf("cd checking\n", q_candidate.size());

				}

				//	if(whole_stnu.s_info.id_constrs == C_DC)
				{
				  v_current_conflicts = current_checking.DCConflict(whole_stnu,
						current_candidate);
				}
				//	else
				{
				  //v_current_conflicts = current_checking.SCConflict(whole_stnu,
				//						    current_candidate);
				}
				if (current_candidate.b_checked && v_current_conflicts.size())
				{
					printf("error(%s)!\n", this->whole_stnu.file_name.c_str());
					current_candidate.print(whole_stnu);

					printf("\n\nUnknown conflicts:\n");
					for (int i = 0; i < v_current_conflicts.size(); i++)
						v_current_conflicts[i].print(whole_stnu);
					exit(0);
				}
				current_candidate.b_checked = true;

				//v_current_conflicts = current_checking.DCConflict(current_candidate.stnu_);
				if (a.s_info.id_debug == DEBUG_T)
				{
					printf("dc checking completed (%d conflicts)\n",
							v_current_conflicts.size());

				}
				if (v_current_conflicts.size())
				{
					for (int i = 0; i < v_current_conflicts.size(); i++)
					{
						v_current_conflicts[i].id = v_conflicts.size();
						v_conflicts.push_back(v_current_conflicts[i]);
					}

					//check duplicated conflicts
					if (whole_stnu.s_info.id_debug == DEBUG_T)
					{
						for (int i = 0; i < v_conflicts.size(); i++)
						{
							v_conflicts[i].print(whole_stnu);
							for (int j = 0; j < i; j++)
							{
								if (v_conflicts[i].isSame(v_conflicts[j]))
								{
									printf("Same conflicts!\n");
									//	exit(0);
								}
							}
						}
					}

					if (a.s_info.id_obj == O_DC_CHECK
							&& a.s_info.id_variable == V_SC)
					{
						if (a.s_info.id_debug == DEBUG_C)
						{
							//current_candidate.print(a);
							v_current_conflicts[0].print(a);
						}
						continue;
					}
					q_candidate.push(current_candidate);

				}
				else
				{
					run_time = getSystemTime() - s_time;
					if (a.s_info.id_debug == DEBUG_C)
						printf("find a feasible solution!\n");
					current_candidate.b_solvable = true;
					generateResult(current_candidate);
					if (a.s_info.id_debug == DEBUG_T)
						printf("Relaxed Bounds: %d\n",
								current_candidate.v_relaxed_bounds.size());
					for (int i = 0;
							i < current_candidate.v_relaxed_bounds.size(); i++)
					{
						bool b = current_candidate.v_relaxed_bounds[i].first;
						if (b == false)
							continue;
						LINK tmp_edge;
						if (i < whole_stnu.n_ctg * 2)
						{
							tmp_edge = whole_stnu.ctg[i / 2];
						}
						else
							tmp_edge = whole_stnu.rqm[(i - 2 * whole_stnu.n_ctg)
									/ 2];

						double ori, rel;
						rel = current_candidate.v_relaxed_bounds[i].second;
						if (i % 2)
							ori = tmp_edge.ub;
						else
							ori = tmp_edge.lb;

					}
					//continue;
					return current_candidate;
				}
			}
		}
	}

	//fail to find a feasible solution
	if (a.s_info.id_debug == DEBUG_T)
		printf("fail to find a feasible solution!\n");

	new_candidate.value_preference = -1;
	new_candidate.b_solvable = false;
	run_time = getSystemTime() - s_time;

	generateResult(new_candidate);
	return new_candidate;
}

void Cds::CDRU1(const stnu & a)
{
	long long s_time = getSystemTime();

	candidate_value.clear();
	Candidate new_candidate(a);
	new_candidate.setPreference(a);

	//printf("%s\n",a.file_name.c_str());
	//printf("Initial Candidate: %.3lf\n", new_candidate.value_preference);
	q_candidate.push(new_candidate);
	candidate_value.insert(new_candidate.value_preference);
	double ub_candidate = g_inf;
	whole_stnu = a;

	while (q_candidate.size() != 0)
	{
		Candidate current_candidate = q_candidate.top();
		q_candidate.pop();
		//	if(current_candidate.stnu_.value_preference +1e-4>= ub_candidate)
		//		continue;

		//ub_candidate = current_candidate.stnu_.value_preference;
		if (whole_stnu.s_info.id_debug >= DEBUG_C)
		{
			printf("\nDequeue Candidate: size(%d)\n", q_candidate.size());

			current_candidate.print();
		}
		Controllability current_checking = Controllability();
		v_conflicts.clear();

		if(isComplete(current_candidate) == false)
		{
		  expandOnVariables(current_candidate);
		  continue;
		}

		//	SimStnu simple_stnu;
		v_conflicts = current_checking.DCConflict1(whole_stnu,
				current_candidate);
		//printf("\nDequeue Candidate: size(%d)\n", v_conflicts.size());
		if (v_conflicts.size() == 0)
		{
			run_time = getSystemTime() - s_time;
			//sol.v_conflicts = current_candidate.v_resolved_conflicts;
			//sol.v_relaxed_bounds = current_candidate.v_relaxed_bounds;

			//	sol.simplifyConflicts(current_candidate);
			//sol.b_solve = true;
			while (q_candidate.size())
				q_candidate.pop();
			if (whole_stnu.s_info.id_debug >= DEBUG_C)
				printf("feasible stnu\n");
			this->generateResult(current_candidate);
			return;
			//return current_candidate.v_resolved_conflicts;
		}
		else
		{
			if(whole_stnu.s_info.id_obj != O_DC_CHECK)
				expandOnConflict1(current_candidate, v_conflicts);
		}
		//exit(0);
		new_candidate = current_candidate;
	}

	//fail to find a feasible solution
	run_time = getSystemTime() - s_time;
	new_candidate.value_preference = -1;
	new_candidate.b_solvable =false;
	this->generateResult(new_candidate);
	//sol.b_solve = false;

}

void Cds::CDRU2(const Candidate & o_candidate, SubSolution & sol)
{
	long long s_time = getSystemTime();

	candidate_value.clear();
	Candidate new_candidate = o_candidate;

	//printf("Initial Candidate: %.3lf\n", t_stnu.get_value_preference());
	q_candidate.push(new_candidate);
	candidate_value.insert(new_candidate.value_preference);
	double ub_candidate = g_inf;
	bool first = true;

	while (q_candidate.size() != 0)
	{
		Candidate current_candidate = q_candidate.top();
		q_candidate.pop();
		//	if(current_candidate.stnu_.value_preference +1e-4>= ub_candidate)
		//		continue;

		//ub_candidate = current_candidate.stnu_.value_preference;
		if (whole_stnu.s_info.id_debug)
		{
			printf("\nDequeue Candidate: size(%d)\n", q_candidate.size());

			current_candidate.print();
		}
		Controllability current_checking = Controllability();
		v_conflicts.clear();

		//	SimStnu simple_stnu;
		//if(first)
		v_conflicts = current_checking.DCConflict(whole_stnu,
				current_candidate);
		//	else
		//v_conflicts=current_checking.DCConflict1(whole_stnu,current_candidate);
		//current_candidate.print();
		//printf(" ");
		//std::cerr<<first<<" "<<v_conflicts.size()<<endl;
		if (v_conflicts.size() == 0)
		{
			run_time = getSystemTime() - s_time;
			sol.v_conflicts = current_candidate.v_resolved_conflicts;
			sol.v_relaxed_bounds = current_candidate.v_relaxed_bounds;

			//	sol.simplifyConflicts(current_candidate);
			sol.b_solved = true;
			while (q_candidate.size())
				q_candidate.pop();
			if (whole_stnu.s_info.id_debug >= DEBUG_C)
				printf("feasible stnu\n");

			return;
			//return current_candidate.v_resolved_conflicts;
		}
		else
		{
			//cerr<<first<<" "<<v_conflicts.size()<<endl;

			if (first == false)
			{

				v_conflicts[0].print(whole_stnu);
				printf("%s neg\n", whole_stnu.file_name.c_str());
				exit(0);
			}
			expandOnConflict(current_candidate, v_conflicts);

		}
		//exit(0);
		new_candidate = current_candidate;
		first = false;
	}

	//fail to find a feasible solution
	run_time = getSystemTime() - s_time;
	sol.b_solved = false;

}

void Cds::expandOnVariables(Candidate & a)
{
	// expand current candidate by assigning next variable
	int var_id = a.v_unassigned_vars.front();

	for (int i = 0; i < whole_stnu.v_dis_var[var_id].n_value; i++)
	{
		Candidate new_candidate = a;
		new_candidate.v_unassigned_vars.pop();
		new_candidate.b_checked = false;
		new_candidate.v_assignments[var_id] = i;
		new_candidate.vec_ass.push_back(make_pair(var_id, i));
//				whole_stnu.v_dis_var[var_id].v_values[i];

		// update 20160601: replace stnu from candidate by assignments
		bool fl_new = true;

		if (fl_new == true)
		{
			// the previous idea is to add the max
			//new_candidate.reward_vars +=
			//		whole_stnu.v_dis_var[var_id].v_utilities[i]
			//				- whole_stnu.v_dis_var[var_id].m_utility;
			if (whole_stnu.s_info.id_obj == O_MAX_DELAY)
			{
				new_candidate.reward_vars +=
						whole_stnu.v_dis_var[var_id].v_utilities[i];
				new_candidate.setPreference(whole_stnu);

				//new_candidate.value_preference = a.value_preference;
			}
			else
			{
				if(whole_stnu.s_info.id_obj == O_RELAX_COST)
					new_candidate.reward_vars +=
						whole_stnu.v_dis_var[var_id].v_utilities[i];
				else
					new_candidate.reward_vars = 0;

				new_candidate.value_preference = new_candidate.reward_vars;
			}
			if (whole_stnu.s_info.id_debug >= DEBUG_C)
			{
				printf("Add New Candidate: size(%d)\n", q_candidate.size());
				a.print();
				new_candidate.print(whole_stnu);
			}
			q_candidate.push(new_candidate);
		}
	}
}

void Cds::expandOnConflict(const Candidate & candidate,
		const std::vector<Conflict> v_c)
{
	// expand current candidate according to the vector of candidate

	/*std::set<Conflict> s_conflict;
	s_conflict.clear();

	Candidate new_candidate = candidate;
	for (int i = 0; i < v_c.size(); i++)
	{
		if (s_conflict.find(v_c[i]) != s_conflict.end())
			continue;
		//printf("\nConflict set size: %d\n", s_conflict.size());

		if (v_c[i].v_conflict_edge_id.size())
		{
			//CreateCandidateOnConflict(candidate, v_c[i]);
			new_candidate.v_resolved_conflicts.push_back(v_c[i]);
			//new_candidate.v_resolved_conflicts[new_candidate.v_resolved_conflicts.size() - 1].b_var = false;

			s_conflict.insert(v_c[i]);
			if (whole_stnu.s_info.id_debug == DEBUG_T)
			{
				printf("add conflict\n");
				v_conflicts[i].print(whole_stnu);
			}
		}
	}

	//LP lp_solver;
	//conflict.print();
	//lp_solver.solve(new_candidate, whole_stnu);
	//new_candidate.b_solvable = true;

	//calculate the preference value and push the new candidate to the q
	if (new_candidate.b_solvable)
	{
		if (whole_stnu.s_info.id_obj == O_MAX_DELAY
				|| whole_stnu.s_info.id_obj == O_MAX_EARLINESS)
		{
			new_candidate.setPreference(whole_stnu);
			whole_stnu.setCtg(new_candidate.value_preference);
			if (candidate_value.find(new_candidate.value_preference)
					== candidate_value.end())
			{
				q_candidate.push(new_candidate);

			}
		}
		else if (whole_stnu.s_info.id_obj == O_DC_CHECK
				|| whole_stnu.s_info.id_obj%10 == O_RELAX_COST)
		{
			new_candidate.setPreference(whole_stnu);
			if (whole_stnu.s_info.id_debug == DEBUG_T)
				printf("New Candidate: %.3lf (%.3lf)\n",
						new_candidate.value_preference,
						new_candidate.reward_vars);
			q_candidate.push(new_candidate);
		}
		else
		{
			std::cerr << "error: the objective type (-o) has not been set!"
					<< std::endl;
		}
	}
	s_conflict.clear();*/
}

void Cds::expandOnConflict1(const Candidate & candidate,
		const std::vector<Conflict> v_c)
{
	// expand current candidate according to the vector of candidate

	std::set<Conflict> s_conflict;
	s_conflict.clear();

	Candidate new_candidate = candidate;
	for (int i = 0; i < v_c.size(); i++)
	{
		if (s_conflict.find(v_c[i]) != s_conflict.end())
			continue;
		//printf("\nConflict set size: %d\n", s_conflict.size());
		if (whole_stnu.s_info.id_debug == DEBUG_T)
		{
			printf("Create Candidate \n");
			v_conflicts[i].print(whole_stnu);
		}

		CreateCandidateOnConflict(candidate, v_c[i]);

		s_conflict.insert(v_c[i]);

	}

	s_conflict.clear();
}

void Cds::CreateCandidateOnConflict(const Candidate & candidate,
		const Conflict & conflict)
{
	// create a new candidate from candidate and conflict,
	// push the new candidate into candidate queue

/*	Candidate new_candidate;
	new_candidate = candidate;
	//new_candidate.v_resolved_conflicts;
	new_candidate.v_resolved_conflicts.push_back(conflict);

	//encode a LP model for all conflicts and the new conflict
	//solve the problem and add relaxation to the new candidate
	//LP lp_solver;
	if (whole_stnu.s_info.id_debug >= DEBUG_T)
	{
		printf("new conflict\n");
		conflict.print();
	}
	//lp_solver.solve(new_candidate, whole_stnu);

	//calculate the preference value and push the new candidate to the q
	if (new_candidate.b_solvable)
	{
		if (whole_stnu.s_info.id_obj == O_MAX_DELAY
				|| whole_stnu.s_info.id_obj == O_MAX_EARLINESS)
		{
			new_candidate.setPreference(whole_stnu);
			new_candidate.setMaxDelay(whole_stnu,
					new_candidate.value_preference);
			//double best_value = *(candidate_value.begin());
			//printf("%lf\n",new_candidate.value_preference);
			if (candidate_value.find(new_candidate.value_preference)
					== candidate_value.end())
			{
				q_candidate.push(new_candidate);
				if (whole_stnu.s_info.id_debug >= DEBUG_C)
					printf("New Candidate: %.3lf (%.3lf)\n",
							new_candidate.value_preference,
							new_candidate.reward_vars);
				candidate_value.insert(new_candidate.value_preference);

			}
			else
			{
				if (whole_stnu.s_info.id_debug >= DEBUG_T)
					printf("Repeated Candidate\n");
			}
		}
		else if (whole_stnu.s_info.id_obj == O_DC_CHECK
				|| whole_stnu.s_info.id_obj%10 == O_RELAX_COST)
		{
			new_candidate.setPreference(whole_stnu);
			if (whole_stnu.s_info.id_debug == DEBUG_T)
				printf("New Candidate: %.3lf (%.3lf)\n",
						new_candidate.value_preference,
						new_candidate.reward_vars);
			q_candidate.push(new_candidate);
		}
		else
		{
			std::cerr << "error: the objective type (-o) has not been set!"
					<< std::endl;
		}
	}
	if (whole_stnu.s_info.id_debug >= DEBUG_T)
		printf("Not feasible\n");
*/
}

std::vector<pair<int, int> > Cds::getNextAssignment(
		std::vector<pair<int, int> > v_as)
{
	// usage: return the next assignment
	// if v_as is the last assignment, return empty
	int id;
	int i_order;
	//id = v_as.size() - 1;
	if (v_as.size() <= 0)
		return v_as;
	i_order = whole_stnu.n_vars - 1;
	id = this->whole_stnu.order_dv[i_order];
	v_as[id].second++;

	while (i_order >= 0 && v_as[id].second >= whole_stnu.v_dis_var[id].n_value)
	{
		v_as[id].second = 0;
		i_order--;
		//id--;
		if (i_order < 0) // the current assignment is the last one.
			return std::vector<pair<int, int> >();
		id = this->whole_stnu.order_dv[i_order];
		v_as[id].second++;
	}
	return v_as;
}

//stnu Cds::getSTNU(const stnu &current_stnu, std::vector<int> v_as)
//{
//	// usage: return the stnu with assignment v_as
//	stnu res = current_stnu;
//	bool *b_ctg = new bool[res.n_ctg]; // mark the mute ctgs
//	bool *b_rqm = new bool[res.n_rqm]; // mark the mure rqms
//
//	res.ctg.clear();
//	res.rqm.clear();
//
//	for (int i = 0; i < res.n_ctg; i++)
//		b_ctg[i] = false;
//	for (int i = 0; i < res.n_rqm; i++)
//		b_rqm[i] = false;
//
//	// mute the links according to the value of v_as
//	for (int i_order = 0; i_order < v_as.size(); i_order++)
//	{
//		int i = res.order_dv[i_order];
//		for (int k = 0; k < res.v_dis_var[i].vv_id_mute_edge[v_as[i]].size();
//				k++)
//		{
//			int id_edge = res.v_dis_var[i].vv_id_mute_edge[v_as[i]][k];
//			if (id_edge < 0) // mute ctg
//			{
//				id_edge = -id_edge;
//				id_edge--;
//				b_ctg[id_edge] = true;
//			}
//			else // mute rqm
//			{
//				id_edge--;
//				b_rqm[id_edge] = true;
//			}
//		}
//	}
//
//	for (int i = 0; i < res.n_ctg; i++)
//	{
//		if (b_ctg[i] == false)
//			res.ctg.push_back(current_stnu.ctg[i]);
//	}
//
//	for (int i = 0; i < res.n_rqm; i++)
//	{
//		if (b_rqm[i] == false)
//			res.rqm.push_back(current_stnu.rqm[i]);
//	}
//
//	// get feasible STNU from the branching STNU
//
//	res.n_ctg = res.ctg.size();
//	res.n_rqm = res.rqm.size();
//	return res;
//}

bool Cds::combineSubSol(SubSolution &a, SubSolution b, int v_id)
{
	/* combine 2 subsolutions.
	 * for every possible T(v_id):
	 * 		a'=update(a,v_id)
	 * 		b'=update(b,v_id)
	 * 		if (a' or b' is not feasible)
	 * 			continue
	 * 		if (a' or b' solves the problem)
	 * 			return true
	 * 		push (a' or b') into new subsolution
	 *
	 * 	a=new subsolution
	 * */

	// if b is not solvable, keep a
	if (b.b_solvable == false)
		return false;
	b.v_value_dis_var[v_id].second = -1;
	// if b is solvable and covers all situations
	if (b.b_solved == true)
	{
		a = b;
		return false;
	}

	// if a covers all conditions
	if (a.b_solved)
		return true;

	vector<int> v_time;
	vector<vector<SimConflict> > com_simconf;
	SimConflict l_branch, u_branch;
	//a.v_union_conj_conflicts.clear();
	v_time = whole_stnu.v_dis_var[v_id].v_ctg_pri;
	if (whole_stnu.s_info.id_debug)
	{
		a.print();
		b.print();
	}
	SubSolution c = a;
	c.v_value_dis_var = b.v_value_dis_var;

	vector<vector<SimConflict> > preh_conf;
	for (int i = 0; i < v_time.size(); i++)
	{
		c.v_union_conj_conflicts.clear();
		CTG c_ctg = whole_stnu.ctg[v_time[i]];
		if (c_ctg.isConsist(b.v_value_dis_var) == false)
			continue;
		if (whole_stnu.s_info.id_debug)
			c_ctg.print();
		com_simconf.clear();

		upSimConfbyCtg(com_simconf, a.v_current_conflicts, v_time[i]);
		upSimConfbyCtg(com_simconf, b.v_union_conj_conflicts, v_time[i]);

		if (com_simconf.size() == 0)
			continue;
		intervals.clear();
		for (int j = 0; j < com_simconf.size(); j++)
		{
			T tmp_t;
			tmp_t.lb = -1e20;
			tmp_t.ub = 1e20;
			for (int k = 0; k < com_simconf[j].size(); k++)
			{
				if (com_simconf[j][k].conf_type == SIM_CONF_L)
				{
					tmp_t.lb = com_simconf[j][k].lb;
					if (l_branch.v_variable.size()
							< com_simconf[j][k].v_variable.size())
						l_branch = com_simconf[j][k];
				}

				else
				{
					tmp_t.ub = -com_simconf[j][k].lb;
					if (u_branch.v_variable.size()
							< com_simconf[j][k].v_variable.size())
						u_branch = com_simconf[j][k];
				}
				if (whole_stnu.s_info.id_debug)
					com_simconf[j][k].print();

			}
			intervals.push_back(tmp_t);
			//	printf("%lf %lf\n",tmp_t.lb, tmp_t.ub);
		}
		sort(intervals.begin(), intervals.end(), cmpl);

		l_branch.removeCtg(v_time[i] + 1); //pre
		u_branch.removeCtg(v_time[i] + 1); //pre

		double p1 = intervals[0].lb - c_ctg.lb;
		if (p1 < 0)
			p1 = 0;
		double cl, cu;
		cu = intervals[0].ub - p1;
		bool b_cover = false;
		if (cu > c_ctg.ub)
		{
			b_cover = true;
			//break;
		}
		//printf("int0 [%lf, %lf]\n",intervals[0].lb, intervals[0].ub);
		for (int j = 1; j < intervals.size(); j++)
		{
			//printf("int%d [%lf, %lf]\n",j,intervals[j].lb, intervals[j].ub);

			cl = intervals[j].lb - p1;
			if (cl > cu)
			{
				if (b_cover)
				{
					l_branch.lb = p1;
					u_branch.lb = -(cu + p1 - c_ctg.ub);
					if (l_branch.lb <= -u_branch.lb)
					{
						c.v_union_conj_conflicts.push_back(
								vector<SimConflict>());
						int tmp_id = c.v_union_conj_conflicts.size() - 1;
						if (l_branch.v_variable.size())
							c.v_union_conj_conflicts[tmp_id].push_back(
									l_branch);
						if (u_branch.v_variable.size())
							c.v_union_conj_conflicts[tmp_id].push_back(
									u_branch);
						//preh.push_back(make_pair(p1,cu+p1 -c_ctg.ub));
						b_cover = false;
					}
				}
				p1 = intervals[j].lb - c_ctg.lb;
				if (p1 < 0)
					p1 = 0;
			}
			cu = max(intervals[j].ub - p1, cu);
			if (cu > c_ctg.ub)
			{
				b_cover = true;
				//break;
			}
		}

		if (b_cover)
		{
			if (l_branch.lb <= -u_branch.lb)
			{
				l_branch.lb = p1;
				u_branch.lb = -(cu + p1 - c_ctg.ub);
				c.v_union_conj_conflicts.push_back(vector<SimConflict>());
				int tmp_id = c.v_union_conj_conflicts.size() - 1;
				if (l_branch.v_variable.size())
					c.v_union_conj_conflicts[tmp_id].push_back(l_branch);
				if (u_branch.v_variable.size())
					c.v_union_conj_conflicts[tmp_id].push_back(u_branch);
			}

			//preh.push_back(make_pair(p1,cu+p1-c_ctg.ub));
		}

		bool empty_case = false;
		if (c.v_union_conj_conflicts.size()
				&& c.v_union_conj_conflicts[0].size() == 0
				&& l_branch.lb <= -u_branch.lb)
			empty_case = true;

		if (empty_case
				|| (c.v_union_conj_conflicts.size() && checkPre(c, c_ctg.start)))
		{
			//c.print();
			a.v_union_conj_conflicts.clear();
			a.v_current_conflicts.clear();
			a.b_solved = true;
			return true;
		}

		preh_conf.insert(preh_conf.end(), c.v_union_conj_conflicts.begin(),
				c.v_union_conj_conflicts.end());
		//else

//		sort(intervals.begin(),intervals.end(),cmpl);
	}
	a.b_solved = false;
	a.b_solvable = true;
	a.v_current_conflicts.insert(a.v_current_conflicts.end(),
			b.v_union_conj_conflicts.begin(), b.v_union_conj_conflicts.end());
	a.v_union_conj_conflicts = preh_conf;
	a.v_value_dis_var = b.v_value_dis_var;
	return false;
}

bool Cds::combine2STNU(SubSolution& a, SubSolution b, int var_id)
{
	/* Combine 2 STNUs
	 * 1. find the latest common discrete variable
	 * 2. combine the stnus/constraints
	 * 3. replace the variable which is combined by its bounds
	 * */

	//return combineSubSol(a,b,var_id);
	// update the sub-solution by replacing variables which
	// are following the discrete variable to its bounds
	if (this->whole_stnu.s_info.id_debug)
	{
		printf(ANSI_COLOR_RED); // change the color of output
		printf("update solution_a:\n");
	}
	this->updateSimConflict(a, var_id);
	if (this->whole_stnu.s_info.id_debug)
	{
		printf("updated solution_a:\n");
		a.print();
		printf("\nupdate solution_b:\n");
		b.print();
	}
	this->updateSimConflict(b, var_id);
	if (this->whole_stnu.s_info.id_debug)
	{
		printf("updated solution_b\n");
		b.print();
		printf(ANSI_COLOR_RESET); // change the color of output
	}

	// if b is not solvable, keep a
	if (b.b_solvable == false)
		return false;

	// if b is solvable and covers all situations
	if (b.b_solved == true && b.v_union_conj_conflicts.size() == 0)
	{
		a = b;
	}

	// if a covers all conditions
	if (a.b_solved && a.v_union_conj_conflicts.size() == 0)
		return true;

	/*
	 if (a.b_solve == false)
	 {
	 //a = b;
	 a.b_solve = b.b_solve;
	 a.v_union_conj_conflicts = b.v_union_conj_conflicts;
	 return false;
	 }
	 */
	//removeRedundant(a, b);
	//a.v_union_conj_conflicts.insert(a.v_union_conj_conflicts.begin(),
	//		b.v_union_conj_conflicts.begin(), b.v_union_conj_conflicts.end());
	//if(tryCombine(a,b,var_id))
	//if(checkSubSolution(a))
	{
		// if a covers the original problem, return TRUE
		a.v_union_conj_conflicts.clear();
		a.b_solved = true;
	}
	return true;
}
// combine the last ctg

//
//bool Cds::tryCombine(SubSolution &a, SubSolution b,int var_id)
//{
//	T tmp_t;
//	intervals.clear();
//	a.v_value_dis_var = b.v_value_dis_var;
//	for (int i = 0; i < a.v_lb.size(); i++)
//	{
//	//	T tmp_t;
//		tmp_t.lb = a.v_lb[i];
//		tmp_t.ub = a.v_ub[i];
//		tmp_t.id = i;
//		//printf("%lf %lf\n",tmp_t.lb, tmp_t.ub);
//		intervals.push_back(tmp_t);
//	}
//
//	tmp_t.id = a.v_lb.size();
//	tmp_t.lb = 1e10;
//	tmp_t.ub = -1e10;
//
//	for (int i = 0; i < b.v_union_conj_conflicts.size(); i++)
//	{
//		for( int j = 0; j < b.v_union_conj_conflicts[i].size(); j++)
//		{
//			double tmp_lb = b.v_union_conj_conflicts[i][j].lb;
//			if(b.v_union_conj_conflicts[i][j].conf_type == SIM_CONF_U)
//			{
//				tmp_lb =-tmp_lb;
//				tmp_t.ub = tmp_t.ub>tmp_lb?tmp_t.ub:tmp_lb;
//			}
//			else if(b.v_union_conj_conflicts[i][j].conf_type == SIM_CONF_L)
//			{
//				tmp_t.lb = tmp_t.lb<tmp_lb?tmp_t.lb:tmp_lb;
//			}
//		}
//	}
//
//	if(tmp_t.lb == 1e10) tmp_t.lb = -1e10;
//	if(tmp_t.ub == -1e10) tmp_t.ub = 1e10;
//	a.v_lb.push_back(tmp_t.lb);
//	a.v_ub.push_back(tmp_t.ub);
//
//	intervals.push_back(tmp_t);
//	sort(intervals.begin(), intervals.end(), cmpl);
//	vector<int> v_time;
//	v_time=whole_stnu.v_dis_var[var_id].v_ctg_pri;
//	//getCombinePoints(vc,v_time);
//	bool b_cover = false;
//	vector<pair<double,double> > preh;
//	for(int i = 0; i < v_time.size(); i++)
//	{
//		CTG c_ctg = whole_stnu.ctg[v_time[i]];
//
//
//		double p1 = intervals[0].lb - c_ctg.lb;
//		if(p1<0)p1=0;
//		double cl, cu;
//		cu = intervals[0].ub - p1;
//		for (int i = 1; i < intervals.size(); i++)
//		{
//			cl = intervals[i].lb - p1;
//			if(cl > cu)
//			{
//				if(b_cover)
//				{
//					preh.push_back(make_pair(p1,cu+p1 -c_ctg.ub));
//					b_cover = false;
//				}
//				p1  = intervals[i].lb - c_ctg.lb;
//				if(p1 < 0) p1 = 0;
//			}
//			cu = max(intervals[i].ub - p1,cu);
//			if(cu > c_ctg.ub)
//			{
//				b_cover = true;
//				//break;
//			}
//		}
//
//		if(b_cover)
//		{
//			preh.push_back(make_pair(p1,cu+p1-c_ctg.ub));
//		}
//
//		if(checkPre(preh, a, c_ctg.end))
//			return true;
//	}
//
//	//next step:
//	// if [p1,p2]
//	return false;
//}

bool Cds::checkPre(SubSolution b, int en)
{
	Candidate new_cand = Candidate(whole_stnu, b.v_value_dis_var);
	Controllability contr;

	SimConflict simc = b.v_union_conj_conflicts[0][0];
	int st = b.v_union_conj_conflicts[0][0].start_node;
	//whole_stnu.add_rqm(st,en);
	//whole_stnu.n_rqm++;
	//for (it = v.begin(); it!= v.end(); it++)
	for (int i = 0; i < b.v_union_conj_conflicts.size(); i++)
	{
		RQM tmp_rqm(st, en);
		for (int j = 0; j < b.v_union_conj_conflicts[i].size(); j++)
		{
			if (b.v_union_conj_conflicts[i][j].conf_type == SIM_CONF_L)
				tmp_rqm.lb = b.v_union_conj_conflicts[i][j].lb;
			else
				tmp_rqm.ub = -b.v_union_conj_conflicts[i][j].lb;
			//if(tmp_rqm.start == 16  && tmp_rqm.end == 30 &&
			//		tmp_rqm.lb < 400)
			//tmp_rqm.print();

		}
		tmp_rqm.id = whole_stnu.n_rqm;
		whole_stnu.add_rqm(tmp_rqm);
		whole_stnu.n_rqm++;
		contr.InitDistanceGraph(whole_stnu, new_cand);
		//RQM tmp_rqm(st,en,it->first,it->second);
		if (whole_stnu.s_info.id_debug)
		{
			printf("added rqm ");
			tmp_rqm.print();
		}

		contr.edgeEdge(tmp_rqm);
		vector<Conflict> vc;
		vc.clear();
		contr.extractEnvelope(vc);
		if (vc.size() == 0)
		{
			return true;
		}
		whole_stnu.rqm.erase(whole_stnu.rqm.end() - 1);
		whole_stnu.n_rqm--;
	}
	return false;
}

void Cds::updateSimConflict(SubSolution &a, int var_id)
{
	// update subsolution a by replacing variables of links after var[var_id]
	// by their actual bounds
	//printf("%d %d\n",var_id, a.v_value_dis_var.size());
	if (a.v_value_dis_var.size() == 0)
		return;

	//printf("%d %d\n",var_id, a.v_value_dis_var.size());
	//if(a.v_value_dis_var.size() )
	//a.v_value_dis_var.erase(a.v_value_dis_var.begin()+var_id);

	for (int i = a.v_union_conj_conflicts.size() - 1; i >= 0; --i)
	{
		bool b_solvable = true;
		for (int j = a.v_union_conj_conflicts[i].size() - 1; j >= 0; --j)
		{
			if (!updateSingleSimConflict(a.v_union_conj_conflicts[i][j],
					var_id))
			{
				// this conflict is not solvable, the conjunction is infeasible
				a.v_union_conj_conflicts.erase(
						a.v_union_conj_conflicts.begin() + i);
				b_solvable = false;
				break;
			}
			else // solvable
			{
				if (a.v_union_conj_conflicts[i][j].v_variable.size() == 0)
					a.v_union_conj_conflicts[i].erase(
							a.v_union_conj_conflicts[i].begin() + j);
			}
		}

		if (b_solvable && a.v_union_conj_conflicts[i].size() == 0)
		{
			// the current conjunction covers the full set
			a.b_solved = true;
			a.v_union_conj_conflicts.clear();
			return;
			//a.v_union_conj_conflicts.erase(a.v_union_conj_conflicts.begin()+i);
		}

		// if a.v_union[i] and bounds is empty remove a.v_union[i]
	}
	a.v_value_dis_var[var_id].second = -1;
	if (a.v_union_conj_conflicts.size() == 0)
		a.b_solved = false;
}

bool Cds::updateSingleSimConflict(SimConflict &sim_conf, int var_id)
{
	//update a single conflict by replacing variables of links after var[var_id]
	// by their actual bounds
	// if the conflict is not solvable return false
	// otherwise return true
	double tmp_lb = 0.0;
	for (int k = sim_conf.v_variable.size() - 1; k >= 0; k--)
	{
		bool fl_follow = false; // if the edge follows the dv
		SimEdge tmp_edge;

		tmp_edge = sim_conf.v_variable[k];
		LINK tmp_link;
		pair<LINK, bool> res;
		res = whole_stnu.idToEdge(tmp_edge.edge_id);
		tmp_link = res.first;
		tmp_edge.b_ub = res.second;
		sim_conf.lb = sim_conf.olb;

		for (int l = 0; fl_follow == false && l < tmp_edge.dis_vars.size(); ++l)
		{
			if (tmp_edge.dis_vars[l] == var_id)
				fl_follow = true;
		}

		if (fl_follow == false)
		{

			if (tmp_link.b_control == false)
			{ //ctg
				if (tmp_edge.b_ub)
					sim_conf.lb += tmp_link.ub;
				else
					sim_conf.lb -= tmp_link.lb;
			}
			else //rqm
			{
				if (tmp_edge.b_ub)
					sim_conf.lb -= tmp_link.ub;
				else
					sim_conf.lb += tmp_link.lb;
			}
			//sim_conf.v_variable.erase( sim_conf.v_variable.begin() + k);
		}
		else
		{
			//tmp_lb is the relaxation value of the variables
			//if(sim_conf.conf_type != SIM_CONF_D)
			{
				if (tmp_link.b_control == false)			//ctg
				{
					if (tmp_edge.b_ub)
						tmp_lb -= tmp_link.lb;
					else
						tmp_lb += tmp_link.ub;
				}
				else // rqm
				{
					if (tmp_edge.b_ub)
						tmp_lb += tmp_link.lb;
					else
						tmp_lb -= tmp_link.ub;
				}
			}
//else return false;
			//SIM_CONF_D: sum(u-l) => relaxation is 0.
		}
	}

	if (tmp_lb < sim_conf.lb)
		return false;
	return true;
}

//bool Cds::isInfeasible(const SimConflict & a)
//{
//	// usage: if the conflict is not solvable with current bounds return true
//	double lf = 0;
//	for (int k = 0; k < a.v_variable.size(); k++)
//	{
//		SimEdge tmp_edge;
//		tmp_edge = a.v_variable[k];
//		if (tmp_edge.b_ub)
//			lf += whole_stnu.ctg[tmp_edge.ctg_id].lb;
//		else
//			lf -= whole_stnu.ctg[tmp_edge.ctg_id].ub;
//	}
//	if (whole_stnu.s_info.id_debug == DEBUG_T)
//	{
//		printf("lb = %lf\n",lf);
//		a.print();
//		//std::cerr << lf << " " << a.v_union_conj_conflicts[i][j].lb
//		//		<< std::endl;
//	}
//	if (lf > -a.lb)
//		return true;
//	return false;
//}

//std::vector<SimConflict> Cds::getUnion(std::vector<SimConflict> a,
//		std::vector<SimConflict> b, int v_id)
//{
//	/* return the union of a1 or a2 ... an or b1 or b2 ... bn
//	 * 1. If the negation of the union is empty, return empty set
//	 * */
//	std::vector<SimConflict> res = a;
//	res.insert(res.end(), b.begin(), b.end());
//
//	for (int i = 0; i < res.size(); i++)
//	{
//		for (int j = i + 1; j < res.size(); j++)
//		{
//			bool s_v = false;
//			for (int k = 0; k < res[i].v_variable.size(); k++)
//			{
//				s_v = false;
//				for (int l = 0; l < res[j].v_variable.size(); l++)
//				//	if(res[j].v_variable[l] == res[i].v_variable[k])
//				{
//					s_v = true;
//					break;
//				}
//				if (s_v == false)
//					break;
//			}
//			if (s_v == true)
//				return std::vector<SimConflict>();
//		}
//	}
//
//	return res;
//
//}
bool Cds::upSimConfbyCtg(vector<vector<SimConflict> >&res,
		vector<vector<SimConflict> > v_union_conj_conflicts, int cid)
{
	/* split the sim-conflict of subsolution A into two parts by ctg
	 * record the new sim-conflicts in res
	 * */

	//printf("%d %d\n",var_id, a.v_value_dis_var.size());
	if (v_union_conj_conflicts.size() == 0)
	{
		return false;
	}

	//printf("%d %d\n",var_id, a.v_value_dis_var.size());
	//if(a.v_value_dis_var.size() )
	//a.v_value_dis_var.erase(a.v_value_dis_var.begin()+var_id);
	//a.v_lb.clear();
	//a.v_ub.clear();
	SimConflict tmp_simconf;
	for (int i = v_union_conj_conflicts.size() - 1; i >= 0; --i)
	{
		bool b_solvable = true;
		SimConflict l_branch; //max(l), max(len(v))
		SimConflict u_branch; //min(u), max(len(v))

		for (int j = v_union_conj_conflicts[i].size() - 1; j >= 0; --j)
		{
			tmp_simconf = v_union_conj_conflicts[i][j];
			if (!upOneSimConfbyCtg(tmp_simconf, cid))
			{
				// this conflict is not solvable, the conjunction is infeasible
				v_union_conj_conflicts.erase(
						v_union_conj_conflicts.begin() + i);
				b_solvable = false;
				break;
			}
			else // solvable
			{
				if (tmp_simconf.v_variable.size() == 0)
				{
					v_union_conj_conflicts[i].erase(
							v_union_conj_conflicts[i].begin() + j);
					continue;
				}
				else
				{
					if (v_union_conj_conflicts[i][j].conf_type == SIM_CONF_L)
					{
						if (l_branch.v_variable.size() == 0) //l_branch is empty
							l_branch = tmp_simconf;
						//l_branch = a.v_union_conj_conflicts[i][j],
						//l_branch.lb = tmp_simconf.lb;
						else
						{
							if (l_branch.lb < tmp_simconf.lb)
								l_branch.lb = tmp_simconf.lb;
							if (l_branch.v_variable.size()
									< tmp_simconf.v_variable.size())
								l_branch.v_variable = tmp_simconf.v_variable;

						}
					}
					else
					{
						if (u_branch.v_variable.size() == 0)
						{
							u_branch = tmp_simconf;	//a.v_union_conj_conflicts[i][j];
							//u_branch.lb = tmp_simconf.lb;//lb = -ub
						}
						else
						{
							if (u_branch.lb < tmp_simconf.lb)
								u_branch.lb = tmp_simconf.lb;
							if (u_branch.v_variable.size()
									< tmp_simconf.v_variable.size())
								u_branch.v_variable = tmp_simconf.v_variable;
						}

					}
				}
			}
		}
		if (!b_solvable)
		{
			continue;
		}

		if (b_solvable && v_union_conj_conflicts[i].size() == 0)
		{
			// the current conjunction covers the full set
			//b_solve = true;
			v_union_conj_conflicts.clear();
			res.clear();
			return true;
			//a.v_union_conj_conflicts.erase(a.v_union_conj_conflicts.begin()+i);
		}

		if (l_branch.v_variable.size() || u_branch.v_variable.size())
		{
			res.push_back(vector<SimConflict>());
			int tmp_id = res.size() - 1;
			if (l_branch.v_variable.size())
				res[tmp_id].push_back(l_branch);
			if (u_branch.v_variable.size())
				res[tmp_id].push_back(u_branch);
		}

		// if a.v_union[i] and bounds is empty remove a.v_union[i]
	}
	//a.v_value_dis_var[var_id].second = -1;
	if (v_union_conj_conflicts.size() == 0)
		return false;
	return false;
}

bool Cds::upOneSimConfbyCtg(SimConflict &sim_conf, int cid)
{
	double tmp_lb = 0.0;
	//whole_stnu.ctg[cid].print();
	for (int k = sim_conf.v_variable.size() - 1; k >= 0; k--)
	{
		bool fl_follow = false; // if the edge follows the dv
		SimEdge tmp_edge;

		tmp_edge = sim_conf.v_variable[k];
		LINK tmp_link;
		pair<LINK, bool> res;
		res = whole_stnu.idToEdge(tmp_edge.edge_id);
		tmp_link = res.first;
		tmp_edge.b_ub = res.second;
		//sim_conf.lb = sim_conf.olb;

		int ctgs = whole_stnu.ctg[cid].start;
		bool b_self = false;
		//for(int l = 0; fl_follow == false && l < tmp_edge.dis_vars.size(); ++l)
		if (tmp_link.id == cid && tmp_link.b_control == false) //tmp_link is ctg[cid]
			b_self = true;
		else if (whole_stnu.i_precede[tmp_link.start][ctgs] < 0
				|| whole_stnu.i_precede[tmp_link.end][ctgs] < 0)
			fl_follow = true;
		//tmp_link.print();
		//printf("%d %d ",(int)fl_follow,int(res.second));

		if (fl_follow)
		{
			if (tmp_link.b_control == false)
			{		//ctg
				if (tmp_edge.b_ub)
					sim_conf.lb += tmp_link.ub;
				else
					sim_conf.lb -= tmp_link.lb;
			}
			else		//rqm
			{
				if (tmp_edge.b_ub)
					sim_conf.lb -= tmp_link.ub;
				else
					sim_conf.lb += tmp_link.lb;
			}
			sim_conf.v_variable.erase(sim_conf.v_variable.begin() + k);
			//printf("%lf\n",sim_conf.lb);
		}
		else
		{

			if (tmp_link.b_control == false)		//ctg
			{
				if (tmp_edge.b_ub)
					tmp_lb -= tmp_link.lb;
				else
					tmp_lb += tmp_link.ub;
			}
			else // rqm
			{
				if (tmp_edge.b_ub)
					tmp_lb += tmp_link.lb;
				else
					tmp_lb -= tmp_link.ub;
			}
		}

		//if(b_self)
		//	sim_conf.v_variable.erase( sim_conf.v_variable.begin() + k);

	}

	if (tmp_lb < sim_conf.lb)
		return false;
	return true;
	return true;
}
stnu Cds::getSTNU(const stnu &current_stnu, std::vector<pair<int, int> > v_as)
{
	// usage: return the stnu with assignment v_as
	stnu res = current_stnu;
	bool *b_ctg = new bool[res.n_ctg]; // mark the mute ctgs
	bool *b_rqm = new bool[res.n_rqm]; // mark the mure rqms

	res.ctg.clear();
	res.rqm.clear();

	for (int i = 0; i < res.n_ctg; i++)
		b_ctg[i] = false;
	for (int i = 0; i < res.n_rqm; i++)
		b_rqm[i] = false;

	// mute the links according to the value of v_as
	for (int i_order = 0; i_order < v_as.size(); i_order++)
	{
		int i = res.order_dv[i_order];
		for (int k = 0;
				k < res.v_dis_var[i].vv_id_mute_edge[v_as[i].second].size();
				k++)
		{
			int id_edge = res.v_dis_var[i].vv_id_mute_edge[v_as[i].second][k];
			if (res.s_info.id_debug == DEBUG_T)
				printf("%d ", id_edge);
			if (id_edge < 0) // mute ctg
			{
				id_edge = -id_edge;
				id_edge--;
				b_ctg[id_edge] = true;
			}
			else // mute rqm
			{
				id_edge--;
				b_rqm[id_edge] = true;
			}
		}
	}
	if (res.s_info.id_debug == DEBUG_T)
		printf("\n");

	for (int i = 0; i < res.n_ctg; i++)
	{
		if (b_ctg[i] == false)
			res.ctg.push_back(current_stnu.ctg[i]);
	}

	for (int i = 0; i < res.n_rqm; i++)
	{
		if (b_rqm[i] == false)
			res.rqm.push_back(current_stnu.rqm[i]);
	}

	// get feasible STNU from the branching STNU

	res.n_ctg = res.ctg.size();
	res.n_rqm = res.rqm.size();
	for (int i = 0; i < res.n_ctg; i++)
	{
		res.ctg[i].id_o = res.ctg[i].id;
		res.ctg[i].id = i;
	}
	for (int i = 0; i < res.n_rqm; i++)
	{
		res.rqm[i].id_o = res.rqm[i].id;
		res.rqm[i].id = i;
	}
	res.setIdLr();
	res.setIdUr();
	return res;
}

std::vector<Conflict> Cds::resolveConflicts(Candidate & a)
{
	std::vector<Conflict> res;
	bool * b_resolved = new bool[v_conflicts.size() + 1];
	for (int i = 0; i < v_conflicts.size(); i++)
		b_resolved[i] = false;
	for (int i = 0; i < a.v_resolved_conflicts.size(); i++)
	{
		b_resolved[a.v_resolved_conflicts[i].id] = true;
	}

	for (int i = 0; i < v_conflicts.size(); i++)
	{
		if (b_resolved[v_conflicts[i].id])
			continue;
		if (whole_stnu.s_info.id_debug == DEBUG_T)
			v_conflicts[i].print(whole_stnu);
		int v_res = resolved(v_conflicts[i], a);
		if (v_res == 0)
		{
			if (whole_stnu.s_info.id_debug == DEBUG_T)
				printf("add conf\n");
			res.push_back(v_conflicts[i]);
			//return res;
		}

	}

	return res;
}

int Cds::resolved(const Conflict & c, Candidate &a)
{
	/* usage: return true, if conflict (c) is not in candidate(a)
	 * 1. check assignments: if c and a have the same variable
	 * 	  but different assignment, return true.
	 * 2. check bounds: if the sum bounds of a according to the
	 * 	  record of edge_id in c is not negative, return true.
	 * 3. otherwise: return false
	 * */

	for (int i = 0; i < c.v_assignments.size(); i++)
	{
		int v_id = c.v_assignments[i].first;
		int v_value = c.v_assignments[i].second; // wrong
		int c_value;

		if (a.v_assignments.find(v_id) == a.v_assignments.end())
			c_value = -1;
		else
			c_value = a.v_assignments[v_id];

		if (whole_stnu.s_info.id_debug == DEBUG_T)
			printf("v(%d)  %d %d\n", v_id, v_value, c_value);

		if (v_value != c_value)
		{
			if (c_value != -1)
				return 2;
			else
				return -2;
		}
	}

	double sum = 0;
	for (int i = 0; i < c.v_conflict_edge_id.size(); i++)
	{
		int edge_id = c.v_conflict_edge_id[i];
		pair<LINK, bool> res;
		res = whole_stnu.idToEdge(edge_id);
		LINK tmp_edge;
		bool b_ub;
		tmp_edge = res.first;
		b_ub = res.second;

		// use default as a rqm (value = ub | -lb)
		double value;
		if (whole_stnu.s_info.id_debug == DEBUG_T)
			printf("%d->%d(%d) %d", tmp_edge.start, tmp_edge.end,
					tmp_edge.b_control, b_ub);
		if (b_ub)
		{
			if (a.v_relaxed_bounds[tmp_edge.id_ur].first)
			{
				value = a.v_relaxed_bounds[tmp_edge.id_ur].second;
			}
			else
				value = tmp_edge.ub;
		}
		else
		{
			if (a.v_relaxed_bounds[tmp_edge.id_lr].first)
			{
				value = -a.v_relaxed_bounds[tmp_edge.id_lr].second;
			}
			else
				value = -tmp_edge.lb;
		}

		if (tmp_edge.b_control == false)
			value = -value;
		sum += value;
		if (whole_stnu.s_info.id_debug == DEBUG_T)
			printf(" +%lf\n", value);

	}
	if (whole_stnu.s_info.id_debug == DEBUG_T)
		printf("=%lf\n", sum);

	if (sum >= -1e-2)
		return 1;
	return 0;
}

bool Cds::isComplete(const Candidate & a)
{
	if (a.v_assignments.size() < whole_stnu.n_vars)
		return false;
	return true;
}

bool Cds::checkSubSolution(SubSolution &a)
{
	// if the sub-solution covers the original problem, return TRUE
	// otherwise, return false;

	// method: enumerate the combinations of one constraint from each disjunction
	//		   if there is one combination which negation is feasible, the
	//         sub-solution does not cover the original problem

	// TODO: the negation of combination which is infeasible can be removed
	//       in order to take this advantage, the sub-solution should keep a
	//		 record of the expanding combinations instead or the union of
	//		 conjunctions

	// initial combination
	vector<SimConflict> v_comb_conflicts;
	vector<int> v_comb_id;
	int comb_size = 0;
	comb_size = a.v_union_conj_conflicts.size();
	for (int i = 0; i < a.v_union_conj_conflicts.size(); i++)
	{
		v_comb_id.push_back(0);
		v_comb_conflicts.push_back(a.v_union_conj_conflicts[i][0]);
	}

	bool b_end = false;
	while (!b_end)
	{

		bool b_check;
		b_check = checkCombConflict(v_comb_conflicts);
		if (b_check) // the negation of this combination is feasible
		{
			return false;
		}
		v_comb_id[comb_size - 1]++;
		int la = comb_size - 1;
		while (la >= 0)
		{
			if (v_comb_id[la] >= a.v_union_conj_conflicts[la].size())
			{
				v_comb_id[la] = 0;
				v_comb_conflicts[la] = a.v_union_conj_conflicts[la][0];
				la--;
				if (la < 0)
					break;
				v_comb_id[la]++;
			}
			else
				break;
		}

		if (la < 0)
		{
			b_end = true;
		}
		else
		{
			v_comb_conflicts[la] = a.v_union_conj_conflicts[la][v_comb_id[la]];
		}
	}

	// if the union conflicts covers the original set
	a.b_solved = true;
	return true;
}

bool Cds::checkCombConflict(vector<SimConflict> &a)
{
	// check the combination of conflicts, if the result is infeasible
	// return false

	//LP lp_solver;
	//return lp_solver.solve(a, whole_stnu);

	return true;
}

void Cds::removeRedundant(SubSolution &old, SubSolution &fresh)
{
	// remove redundant constraints from old sub-solution according to the
	// fresh sub-solution to be combined

	// remove redundant constraints from fresh sub-solution according to the
	// old sub-solutions

	bool **b_red_old, **b_red_fresh; // mark the redundant constraints
	b_red_old = new bool *[old.v_union_conj_conflicts.size()];
	b_red_fresh = new bool *[fresh.v_union_conj_conflicts.size()];

	for (int i = old.v_union_conj_conflicts.size() - 1; i >= 0; --i)
	{
		b_red_old[i] = new bool[old.v_union_conj_conflicts[i].size() + 1];
		b_red_old[i][0] = true;
		for (int j = old.v_union_conj_conflicts[i].size() - 1; j >= 0; --j)
		{
			if (isRedundant(old.v_union_conj_conflicts[i][j], fresh))
				b_red_old[i][j + 1] = true;
			else
				b_red_old[i][j + 1] = false, b_red_old[i][0] = false;
		}
	}

	for (int i = fresh.v_union_conj_conflicts.size() - 1; i >= 0; --i)
	{
		b_red_fresh[i] = new bool[fresh.v_union_conj_conflicts.size() + 1];
		b_red_fresh[i][0] = true;
		for (int j = fresh.v_union_conj_conflicts[i].size() - 1; j >= 0; --j)
		{
			if (isRedundant(fresh.v_union_conj_conflicts[i][j], old))
				b_red_fresh[i][j + 1] = true;
			else
				b_red_fresh[i][j + 1] = false, b_red_fresh[i][0] = false;
		}
	}

	// remove the redundant constraints according to the mark;
	for (int i = old.v_union_conj_conflicts.size() - 1; i >= 0; --i)
	{
		if (b_red_old[i][0])
		{
			old.v_union_conj_conflicts.erase(
					old.v_union_conj_conflicts.begin() + i);
			continue;
		}
		for (int j = old.v_union_conj_conflicts[i].size() - 1; j >= 0; --j)
		{
			if (b_red_old[i][j + 1])
				old.v_union_conj_conflicts[i].erase(
						old.v_union_conj_conflicts[i].begin() + j);
		}
	}
	for (int i = fresh.v_union_conj_conflicts.size() - 1; i >= 0; --i)
	{
		if (b_red_fresh[i][0])
		{
			fresh.v_union_conj_conflicts.erase(
					fresh.v_union_conj_conflicts.begin() + i);
			continue;
		}
		for (int j = fresh.v_union_conj_conflicts[i].size() - 1; j >= 0; --j)
		{
			if (b_red_fresh[i][j + 1])
				fresh.v_union_conj_conflicts[i].erase(
						fresh.v_union_conj_conflicts[i].begin() + j);
		}
	}
}

bool Cds::isRedundant(const SimConflict &a, const SubSolution &b)
{
	// check if conflict a is redundant in sub-solution b
	// method: if a conflict(linear constraint) does not change
	// 		   the solution space of a conjunction of constraints,
	//		   it is redundant for the conjunction of constraints.
	//   	   If a conflict is redundant for a branch (a conjunction)
	//         it is redundant for the sub-solution.

	for (int i = 0; i < b.v_union_conj_conflicts.size(); ++i)
	{
		bool fl_redundant = true;
		for (int j = 0; j < b.v_union_conj_conflicts[i].size(); ++j)
		{
			// if a and b don't have any same edge-variables
			// the a U b = full set
			if (a.compare_vars(b.v_union_conj_conflicts[i][j]))
				continue;

			if (a.compare(b.v_union_conj_conflicts[i][j]) >= 0)
			{
				// if compare return 0/1, a is not redundant
				fl_redundant = false;
				break;
			}
		}

		if (fl_redundant)
		{
			if (whole_stnu.s_info.id_debug >= DEBUG_T)
			{
				printf("Remove a conflict:");
				a.print();
			}
			return true;
		}
	}

	return false;
}

void Cds::generateResult(const Candidate &a)
{
  if(whole_stnu.s_info.id_debug)
  {
    //print the dc stnu
    Controllability tmp_check;
    tmp_check.InitDistanceGraph(whole_stnu, a).print_cstnu("result.cst");
  }
	result = whole_stnu.file_name;
	char buffer[200];
	if (whole_stnu.s_info.id_obj%10 == O_RELAX_COST)
	{
		sprintf(buffer, "\tr = %lf", a.reward_vars);
		result += string(buffer);
		sprintf(buffer, "\tu = %lf\tconf = %d\ttime = %lf", a.value_preference,
				a.v_resolved_conflicts.size(), this->run_time / 1000.0);
		result += string(buffer);
	}
	else if (whole_stnu.s_info.id_obj == O_DC_CHECK)
	{
		if (whole_stnu.n_vars)
			sprintf(buffer, " %d",
					whole_stnu.v_dis_var[whole_stnu.n_vars - 1].n_value);
		else
			sprintf(buffer, " 0");
		result += string(buffer);
		if (a.b_solvable)
			sprintf(buffer, " DC	time = %lf", this->run_time / 1000.0);
		else
			sprintf(buffer, " NDC	time = %lf", this->run_time / 1000.0);
		result += string(buffer);
	}
	else if (whole_stnu.s_info.id_obj == O_MAX_DELAY)
	{
		sprintf(buffer, "\tu = %lf\ttime = %lf", a.value_preference,
				this->run_time / 1000.0);
		result += string(buffer);

	}
	else
	{
		printf("UNKOWN OBJ!\n");
		exit(0);
	}
	//printf("%s\n",result.c_str());
}

void Cds::generateResult(const int &x)
{
	result = whole_stnu.file_name;
	char buffer[200];
	if (whole_stnu.s_info.id_obj%10 == O_RELAX_COST)
	{
		printf("Relaxation Problem haven't been completed!\n");
		exit(0);
	}
	else if (whole_stnu.s_info.id_obj == O_DC_CHECK)
	{
		if (x)
			sprintf(buffer, "\tresult=%d\tleaves=%d\tcombines=%d\ttime = %lf",
					x, n_leaves, n_combines, this->run_time / 1000.0);
		//else sprintf(buffer, " %d	time = %lf",x,this->run_time/1000.0);
		result += string(buffer);
	}
	else
	{
		printf("UNKOWN OBJ!\n");
		exit(0);
	}
}

Cds::~Cds()
{
	// TODO Auto-generated destructor stub
}

void SubSolution::print()
{
	for (int i = 0; i < this->v_value_dis_var.size(); i++)
	{
		printf("V%d = %d	", i, v_value_dis_var[i].second);
	}
	printf("b_solve = %d\n", this->b_solved);

	for (int i = 0; i < v_union_conj_conflicts.size(); i++)
	{
		if (i)
			printf("or\n");
		for (int j = 0; j < v_union_conj_conflicts[i].size(); j++)
		{
			printf("Union Conlict %d %d:\n", i + 1, j + 1);
			v_union_conj_conflicts[i][j].print();
		}
	}

	for (int i = 0; i < v_current_conflicts.size(); i++)
	{
		if (i)
			printf("or\n");
		for (int j = 0; j < v_current_conflicts[i].size(); j++)
		{
			printf("c Conlict %d %d:\n", i + 1, j + 1);
			v_current_conflicts[i][j].print();
		}
	}
}

void SubSolution::simplifyConflicts(const stnu & current_stnu)
{
	//transfer v_conflicts to v_union_conj_conflicts

	v_union_conj_conflicts.push_back(vector<SimConflict>());

	//divide type III
	int size_v = v_conflicts.size();
	for (int i = 0; i < size_v; i++)
	{
		Conflict new_sconf = v_conflicts[i];
		Conflict back_conf = v_conflicts[i];
		new_sconf.v_conflict_edge_id.clear();
		bool fl = false;
		//new_sconf.v_conflict_edge_id.push_back(v_conflicts[i].v_conflict_edge_id[0]);
		for (int j = v_conflicts[i].v_conflict_edge_id.size() - 1; j > 0; j--)
		{
			int tmp_idi = v_conflicts[i].v_conflict_edge_id[j];
			int tmp_idj = v_conflicts[i].v_conflict_edge_id[j - 1];
			if (fl)
			{
				new_sconf.v_conflict_edge_id.push_back(
						v_conflicts[i].v_conflict_edge_id[j]);
				v_conflicts[i].v_conflict_edge_id.erase(
						v_conflicts[i].v_conflict_edge_id.begin() + j);
			}
			if (tmp_idi + tmp_idj == 0)
			{
				fl = !fl;
			}
		}
		if (fl)
		{
			new_sconf.v_conflict_edge_id.push_back(
					v_conflicts[i].v_conflict_edge_id[0]);
			v_conflicts[i].v_conflict_edge_id.erase(
					v_conflicts[i].v_conflict_edge_id.begin());
		}
		if (new_sconf.v_conflict_edge_id.size())
		{
			v_conflicts.push_back(new_sconf);
		}
		else
			v_conflicts[i] = back_conf;
	}

	for (int i = 0; i < v_conflicts.size(); i++)
	{
		//std::vector<SimConflict> sim_conflict;
		//SimConflict tmp_sim_conflict_u;
		//SimConflict tmp_sim_conflict_l;
		SimConflict tmp_sim_conflict;
		tmp_sim_conflict.lb = 0;
		tmp_sim_conflict.olb = 0;
		//tmp_sim_conflict_l = tmp_sim_conflict_u = tmp_sim_conflict;
		tmp_sim_conflict.conf_type = SIM_CONF_O;
		tmp_sim_conflict.start_node = -1;
		//tmp_sim_conflict_l.conf_type = SIM_CONF_L;
		//tmp_sim_conflict_u.conf_type = SIM_CONF_U;

		for (int j = 0; j < v_conflicts[i].v_conflict_edge_id.size(); j++)
		{
			LINK tmp_link;
			SimEdge sim_edge;
			pair<LINK, bool> res;

			res = current_stnu.idToEdge(v_conflicts[i].v_conflict_edge_id[j]);
			tmp_link = res.first;

			//tmp_sim_conflict.sim_edge.ctg_id = tmp_ctg.id;
			sim_edge.edge_id = v_conflicts[i].v_conflict_edge_id[j];
			sim_edge.dis_vars = tmp_link.v_pri_vars;
			//if(sim_edge.dis_vars.size() && sim_edge.dis_vars[0] == 0)
			{
				//tmp_link.print();
				//printf("%d %d\n",tmp_link.s)
			}
			sim_edge.b_ub = res.second;
			if ((res.second && tmp_link.b_control) || //rqm_ub
					(!res.second && !tmp_link.b_control)) //ctg_lb
			{
				sim_edge.start = tmp_link.start_name;
				sim_edge.end = tmp_link.end_name;
				if (tmp_link.b_control == false)
					tmp_sim_conflict.conf_type |= SIM_CONF_L;
				//tmp_sim_conflict_l.v_variable.push_back(sim_edge);
				//tmp_sim_conflict.lb -= tmp_ctg.ub;
			}
			else //rqm_lb ctg_ub
			{
				sim_edge.start = tmp_link.start_name;
				sim_edge.end = tmp_link.end_name;
				if (tmp_link.b_control == false)
					tmp_sim_conflict.conf_type |= SIM_CONF_U;
				//tmp_sim_conflict_u.v_variable.push_back(sim_edge);
				//tmp_sim_conflict.lb += tmp_ctg.lb;
			}
			current_stnu.priNode(tmp_sim_conflict.start_node, tmp_link.start);
			current_stnu.priNode(tmp_sim_conflict.start_node, tmp_link.end);

			tmp_sim_conflict.v_variable.push_back(sim_edge);

		}
		//tmp_sim_conflict.lb = -tmp_sim_conflict.lb;
		if (current_stnu.s_info.id_debug == DEBUG_T)
		{
			if (tmp_sim_conflict.conf_type != SIM_CONF_D)
				tmp_sim_conflict.print();

		}
		//if (tmp_sim_conflict.conf_type != SIM_CONF_D && tmp_sim_conflict.v_variable.size())
		if (tmp_sim_conflict.v_variable.size())
		{
			if (current_stnu.s_info.extra
					>= 10&& tmp_sim_conflict.conf_type == SIM_CONF_D)
				continue;
			if (current_stnu.s_info.id_debug == DEBUG_T)
				tmp_sim_conflict.print();
			//sort(tmp_sim_conflict.v_variable.begin(),
			//		tmp_sim_conflict.v_variable.end());
			if (tmp_sim_conflict.conf_type != SIM_CONF_D)
			{
				if (!isCovered(tmp_sim_conflict, v_union_conj_conflicts[0]))
					v_union_conj_conflicts[0].push_back(tmp_sim_conflict);
			}
			else
			{
				v_union_conj_conflicts[0].push_back(tmp_sim_conflict);
			}
		}
	}
	//this->v_conflicts.clear();
}

bool SubSolution::isCovered(const SimConflict &a, const vector<SimConflict> &b)
{
	if (b.size() == 0)
		return false;
	return false;
	bool fl_redundant = true;
	for (int i = 0; i < b.size(); i++)
	{

//		if(a.compare_vars(b[i]))
//			continue;
//
//		if(a.compare(b[i]) >= 0)
//		{
//			// if compare return 0/1, a is not redundant
//			fl_redundant = false;
//			break;
//		}
	}

	return fl_redundant;
}
} /* namespace CDS */

long long CDS::getSystemTime()
{
	struct timeb t;
	ftime(&t);
	return 1000 * t.time + t.millitm;
}
