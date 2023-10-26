/*
 * testgenerator.cc
 *
 *  Created on: 16/05/2016
 *      Author: jing
 */

#include "../include/testgenerator.h"

namespace CSTNU
{
TestGenerator::TestGenerator()
{
	// TODO Auto-generated constructor stub

}
long long getSystemTime()
{
	struct timeb t;
	ftime(&t);
	return 1000 * t.time + t.millitm;
}

string TestGenerator::createTest(CDS::stnu & original_stnu)
{
	/* steps:
	 * 0. change lb and ub of missions
	 * 1. expand on variables
	 * 2. check dc for each branched stnu
	 * 2.1 if not dc, continue 2
	 * 2.2 if dc, go to 3
	 * 3. for each ctg with labels, lb = eps and ub = ub + lb
	 * 4. check dc
	 * 4.1 if dc: change ub = inf, continue 2
	 * 4.2 if not dc: go to 5
	 * 5. check relaxations: lb? ub? both?
	 * 5.1 lb: lb = new_lb - eps
	 * 5.2 ub: ub = new_ub + eps
	 * 5.3 both: check dlb dub (choose the larger one)
	 * */
	n_stnu = n_ori_ndc = n_mod_dc = n_lb = n_ub = 0;
	expandOnVariables(original_stnu);
	string res;
	char buffer[200];
	sprintf(buffer, "%d %d %d %d %d",n_stnu, n_ori_ndc, n_mod_dc, n_lb,n_ub);
	res = string(buffer);
	return res;

}

void TestGenerator::expandOnVariables(CDS::stnu & current_stnu)
{
	// create the first candidate
	std::vector<pair<int,int> > v_as;
	for(int i = 0; i < current_stnu.n_vars; i++)
		v_as.push_back(make_pair(i,0));
	CDS::Cds t_cds = CDS::Cds();
	long long s_time = getSystemTime();
	//t_cds.setWholeSTNU()
	do{
		long long run_time = getSystemTime() - s_time;
		if(run_time/1000.0 > current_stnu.s_info.f_runtime)
		{
			break;
		}
		if(current_stnu.s_info.id_debug == DEBUG_T){
			char cst_name[100];
			string cst_name_str;
			for(int i = 0; i < current_stnu.n_vars; i++)
			{
				printf("v%d=%d ",v_as[i].first, v_as[i].second);
				sprintf(cst_name, "%d_%d_",v_as[i].first, v_as[i].second);
				cst_name_str += string(cst_name);
			}
			sprintf(cst_name,".cst");
			cst_name_str += string(cst_name);
			//current_stnu.print_cstnu(cst_name_str);
			printf("\n");
		}

		CDS::stnu new_stnu = t_cds.getSTNU(current_stnu, v_as);
		CDS::Controllability dc_check;
		n_stnu++;
		if (dc_check.DCConflict(new_stnu).size() == 0 )
		{
			if(current_stnu.s_info.id_debug == DEBUG_T)
				printf("original stnu is DC!\n");
			CDS::stnu exctg_stnu = new_stnu;
			for(int i = 0; i < new_stnu.n_rqm; i++)
			{
				exctg_stnu.rqm[i].lb_cost_r = 1e11;
				exctg_stnu.rqm[i].ub_cost_r = 1e11;
			}
			for(int i = 0; i < new_stnu.n_ctg; i++)
			{
				exctg_stnu.ctg[i].lb_cost_r = 1e11;
				exctg_stnu.ctg[i].ub_cost_r = 1e11;
			//	if(exctg_stnu.ctg[i].v_labels.size())
				{
				//	if(exctg_stnu.ctg[i].v_labels[0].first - 1 == exctg_stnu.order_dv[0])
				//		continue;
					exctg_stnu.ctg[i].lb_lb = 2;
					exctg_stnu.ctg[i].ub_lb = exctg_stnu.ctg[i].lb;
					exctg_stnu.ctg[i].lb = 2;
					//exctg_stnu.ctg[i].lb_o = 0.1;
					exctg_stnu.ctg[i].lb_ub = exctg_stnu.ctg[i].ub;
					exctg_stnu.ctg[i].ub += 1000;
					exctg_stnu.ctg[i].ub_ub = exctg_stnu.ctg[i].ub;
				//	exctg_stnu.ctg[i].setCostB();
				//	exctg_stnu.ctg[i].ub_o = exctg_stnu.ctg[i].ub;
					exctg_stnu.ctg[i].lb_cost_r = 1;
					exctg_stnu.ctg[i].ub_cost_r = 1;
				}
			}

			if(dc_check.DCConflict(exctg_stnu).size() == 0)
			{
				n_mod_dc++;
				if(current_stnu.s_info.id_debug == DEBUG_T)
					printf("modified stnu is DC!\n");
				for(int i = 0; i < new_stnu.n_ctg; i++)
				{
					if(new_stnu.ctg[i].v_labels.size())
					{
						int ctg_id;
						ctg_id = new_stnu.ctg[i].id_o;
						current_stnu.ctg[ctg_id].ub += 1000;
					//	new_stnu.ctg[i].ub += 1000;
					}
				}
			}
			else
			{
				bool c_ub = 0;
				bool c_lb = 0;
				if(current_stnu.s_info.id_debug == DEBUG_T)
					printf("modified stnu is not DC!\n");
				//printf("modified stnu is not DC!\n");
				CDS::SubSolution sol;
				exctg_stnu.s_info.id_obj = O_RELAX_COST;

				t_cds.setWholeSTNU(exctg_stnu);
			//	exctg_stnu.print_dot("part_stnu.dot");
				CDS::Candidate candi = CDS::Candidate(exctg_stnu);
			//	exctg_stnu.print_dot("part_stnu.dot");
				t_cds.CDRU2(candi, sol);
				//exctg_stnu.print_dot("part_stnu1.dot");

				for(int i = 0; i < new_stnu.n_ctg; i++)
				{

					if(exctg_stnu.ctg[i].v_labels.size() == 0)
						continue;
					if(exctg_stnu.ctg[i].v_labels[0].first - 1 == exctg_stnu.order_dv[0])
						continue;
					double d_lb = 0, d_ub = 0;
					const double eps = 1e-4;
					
					double tmp_lb, tmp_ub;				
					int tmp_id = exctg_stnu.ctg[i].id_lr;
							//exctg_stnu.edgeToId(exctg_stnu.ctg[i], 0);
					if(sol.v_relaxed_bounds[tmp_id].first)
						tmp_lb = sol.v_relaxed_bounds[tmp_id].second;
					else tmp_lb = exctg_stnu.ctg[i].lb;

					d_lb = tmp_lb - exctg_stnu.ctg[i].lb;
					
					tmp_id = exctg_stnu.ctg[i].id_ur;
							//exctg_stnu.edgeToId(exctg_stnu.ctg[i], 1);
					if(sol.v_relaxed_bounds[tmp_id].first)
						tmp_ub = sol.v_relaxed_bounds[tmp_id].second;
					else tmp_ub = exctg_stnu.ctg[i].ub;
					d_ub = exctg_stnu.ctg[i].ub - tmp_ub;
				//	d_lb = new_stnu.ctg[i].lb - sol.v_ctg[i].lb;
				//	d_ub = sol.v_ctg[i].ub - new_stnu.ctg[i].ub;
					int id_ctg = new_stnu.ctg[i].id_o;
					if(d_lb<=0&&d_ub<=0)
						continue;

					if((d_lb < d_ub || c_ub) && c_lb == false)
				//	if(d_ub > 0)
					{
						c_ub = true;

						if(current_stnu.s_info.id_debug == DEBUG_T)
						{
							printf("chub %d->%d: %lf %lf(%d)\n",current_stnu.ctg[id_ctg].start,
									current_stnu.ctg[id_ctg].end, current_stnu.ctg[id_ctg].ub,
									tmp_ub+1,sol.v_relaxed_bounds[tmp_id].first);
						}
						current_stnu.ctg[id_ctg].ub	= tmp_ub +1;

					}
					else if(d_lb > 0 || c_lb)
					{
						c_lb = true;
						if(current_stnu.s_info.id_debug == DEBUG_T)
						{
							printf("chulb %d->%d: %lf %lf\n",current_stnu.ctg[id_ctg].start,
									current_stnu.ctg[id_ctg].end, current_stnu.ctg[id_ctg].lb,
									tmp_lb-1);
						}
						current_stnu.ctg[id_ctg].lb = tmp_lb -1;
						if(current_stnu.ctg[id_ctg].lb < 0)
							current_stnu.ctg[id_ctg].lb = 0.01;
					}

				}
				n_ub+=c_ub;
				n_lb+=c_lb;
			}

		}
		else
		{
			n_ori_ndc++;
		}
		t_cds.setWholeSTNU(current_stnu);
		v_as = t_cds.getNextAssignment(v_as);
	}while(v_as.size());

	//current_stnu.print_cstnu("output.cst");

}

TestGenerator::~TestGenerator()
{
	// TODO Auto-generated destructor stub
}

} /* namespace CSTNU */
