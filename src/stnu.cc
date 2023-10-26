/*
 * stnu.cpp
 *
 *  Created on: 18/07/2014
 *      Author: jing
 */

#include "stnu.h"
#include <stdio.h>
#include "tcn.h"

namespace CDS
{

stnu::stnu()
{
  // TODO Auto-generated constructor stub
  n_vars = 0;
}

stnu::stnu(std::string problem_file, S_PARAM temp_info)
{
//	this->readSTNU(problem_file);
  //this->initSTNU(temp_info);
}

//void stnu::readSTNU(const std::string & problem_file)
//{
//	//read the input
//	this->file_name = problem_file;
//	reader = new Reader(this);
//	reader->read(problem_file);
//}

void stnu::initSTNU(S_PARAM &s_info)
{
  //initiate
  this->s_info = s_info;
  //ctg.clear();
  //rqm.clear();
  n_rqm = 0;
  n_ctg = 0;
  //utility = 0;
  n_ctg = this->ctg.size();
  n_rqm = this->rqm.size();
  //b_ctg = new bool*[n_nodes];
  //s_map = new bool*[n_nodes];
  for (int i = 0; i < n_nodes; i++)
  {
    //b_ctg[i] = new bool[n_nodes];
    //	s_map[i] = new bool[n_nodes];
    //	for (int j = 0; j < n_nodes; j++)
    //	b_ctg[i][j] = false, s_map[i][j] = false;
  }
  //according to the input
  setIdLr();
  setIdUr();
//std::cerr<<this->n_nodes<<" "<<this->n_ctg<<" "<<this->n_rqm<<std::endl;
//exit(0);
  init(n_nodes);

  if (s_info.id_ctg != C_CHANCE)
  {
    for (int i = 0; i < n_ctg; ++i)
    {
      int st = ctg[i].start;
      int en = ctg[i].end;
      double lb = ctg[i].lb;
      double ub;
      if (s_info.id_obj == O_MAX_DELAY) ctg[i].ub = ctg[i].lb + 10000;
      ub = ctg[i].ub;
      double DEFAULT = 1e6;
      if (lb < -DEFAULT) lb = -DEFAULT;
      if (ub > DEFAULT) ub = DEFAULT;
      if (s_info.id_obj%10 != O_RELAX_COST)
      {
        this->set_min(st, en, lb);
        this->set_max(st, en, ub);
      }
      else
      {
        this->set_min(st, en, lb);
        this->set_max(st, en, ub);
      }
    }
  }

  for (int i = 0; i < n_rqm; ++i)
  {
    int st = rqm[i].start;
    int en = rqm[i].end;
    double lb = rqm[i].lb;
    double ub = rqm[i].ub;
    double DEFAULT = 1e6;
    if (lb < -DEFAULT) lb = -DEFAULT;
    if (ub > DEFAULT) ub = rqm[i].ub =DEFAULT;
    if (s_info.id_obj%10 != O_RELAX_COST)
    {
      this->set_min(st, en, lb);
      this->set_max(st, en, ub);
    }
    else
    {
      if (rqm[i].lb_cost_r < 1e10) lb = 0;
      if (rqm[i].ub_cost_r < 1e10) ub = DEFAULT;
      this->set_min(st, en, lb);
      this->set_max(st, en, ub);
    }
  }

  //add labels to the links

  for (int i = 0; i < n_vars; i++)
  {
    v_dis_var[i].m_utility = 0;
    for (int j = 0; j < v_dis_var[i].n_value; j++)
    {
      v_dis_var[i].m_utility =
          v_dis_var[i].m_utility > v_dis_var[i].v_utilities[j] ?
              v_dis_var[i].m_utility : v_dis_var[i].v_utilities[j];

    }
    //	this->utility += v_dis_var[i].m_utility;
    //printf("%lf\n",utility);
  }

  //value_preference = utility;
  this->TransGraph();

  // add the order of ctgs and discrete variables
  bool **pre_dv = new bool*[n_vars];
  for (int i = 0; i < n_vars; i++)
  {
    pre_dv[i] = new bool[n_vars];
    for (int j = 0; j < n_vars; j++)
      pre_dv[i][j] = false;
  }

  for (int ic = 0; ic < n_ctg; ic++)
  {
    CTG tc = ctg[ic];
    for (int i = 0; i < n_vars; i++)
    {
      bool b_ad = bEdgePriVar(ctg[ic], i);

      if (b_ad)
      {
        v_dis_var[i].v_ctg_pri.push_back(ic);
        ctg[ic].v_pri_vars.push_back(i);
        if (ctg[ic].v_labels.size())
        {
          std::vector<std::pair<int, int> >::iterator it;
          for (it = ctg[ic].v_labels.begin(); it != ctg[ic].v_labels.end();
              it++)
          {
            int id_pre_var;
            id_pre_var = it->first;
            //id_pre_var is prior to the ith variable
            pre_dv[id_pre_var][i] = true;
          }
        }
      }
    }
  }

  for (int ir = 0; ir < n_rqm; ir++)
  {
    RQM tr = rqm[ir];
    for (int i = 0; i < n_vars; i++)
    {
      bool b_ad = bEdgePriVar(rqm[ir], i);

      if (b_ad)
      {
        rqm[ir].v_pri_vars.push_back(i);
      }
//				if (ctg[ic].v_labels.size())
//				{
//					std::vector<std::pair<int, int> >::iterator it;
//					for (it = ctg[ic].v_labels.begin();
//							it != ctg[ic].v_labels.end(); it++)
//					{
//						int id_pre_var;
//						id_pre_var = it->first ;
//						//id_pre_var is prior to the ith variable
//						pre_dv[id_pre_var][i] = true;
//					}
//				}
//			}
    }
  }

  // topological sort the discrete variables
  // report error when a cycle exists in pre_dv
  bool *b_used = new bool[n_vars];
  this->order_dv = new int[n_vars];
  for (int i = 0; i < n_vars; i++)
    b_used[i] = false, order_dv[i] = -1;

  for (int i = 0; i < n_vars; i++)
  {
     for(int j = 0; j < v_dis_var[i].v_ctg_pri.size(); j++)
     {
       bool fl = true;
       int ctg1 = v_dis_var[i].v_ctg_pri[j];

       for(int k = 0; k < v_dis_var[i].v_ctg_pri.size(); k++)
       {
    	   if(j == k) continue;
    	   int ctg2 = v_dis_var[i].v_ctg_pri[k];

    	   if(i_precede[ctg[ctg1].end][ctg[ctg2].start] >=0
    			   && i_precede[ctg[ctg1].start][ctg[ctg2].start]>0)
    	   {
    		   fl = false;
    		   break;
    	   }
       }
       if(fl) v_dis_var[i].v_dt.push_back(ctg1);
     }
  }
  for (int i = 0; i < n_vars; i++)
  {

    for (int j = 0; j < n_vars; j++)
    {
      if (b_used[j]) continue;
      int k;
      for (k = 0; k < n_vars; k++)
      {
        if (b_used[k]) continue;
        if (pre_dv[k][j]) break;
      }
      //no discrete variable is prior to var_j
      if (k == n_vars)
      {
        b_used[j] = true;
        order_dv[i] = j;
        break;
      }
    }

    if (order_dv[i] == -1)
    {
      std::cerr << "error: a cycle exists in the precede order of "
          "discrete varisbles!" << std::endl;
    }
  }

  // add preference for DC checking of CCTPU
  if (this->s_info.id_obj == O_DC_CHECK && this->s_info.id_variable == V_DC)
  {
    resetCost();
  }
  return;
  //check the consistency of the network
  if (s_info.id_ctg != C_CHANCE)
  {
    if (this->consistent() == false)
    {
      printf("Inconsistent STNU!\n");
      exit(0);
    }
  }
}

bool stnu::bEdgePriVar(const LINK &a, const int &v_id)
{
  for (int j = 0; j < v_dis_var[v_id].n_value; j++)
  {
    for (int k = 0; k < v_dis_var[v_id].vv_id_mute_edge[j].size(); k++)
    {
      int id_edge = v_dis_var[v_id].vv_id_mute_edge[j][k];
      int start, end;
      if (id_edge < 0)
      {
        id_edge = -id_edge - 1;
        start = ctg[id_edge].start;
        end = ctg[id_edge].end;
        //ctg[id_edge].v_labels.push_back(std::make_pair(i+1, j+1));
      }
      else
      {
        id_edge--;
        start = rqm[id_edge].start;
        end = rqm[id_edge].end;
//      rqm[id_edge].v_labels.push_back(std::make_pair(i+1, j+1));
      }

      if (i_precede[a.end][start] == -1 || i_precede[a.end][end] == -1)
      {
    	//  char buffer[200];
    	//  sprintf(buffer, "p[%d][%d]=%d p[%d][%d]=%d",a.end,start,i_precede[a.end][start],
    	//		  a.end,end,i_precede[a.end][end]);
    	//  std::cerr<<buffer<<std::endl;
        return false;
      }
    }
  }
  return true;
}

void stnu::resetCost()
{
  // reset the cost for dc checking for CCTPU
  // the earlier ctgs have less cost
  for (int i = 0; i < n_ctg; i++)
  {
    ctg[i].lb_cost_r = 0.1;
    ctg[i].ub_cost_r = 0.1;

    for (int j = 0; j < n_nodes; j++)
    {
      if (i_precede[j][ctg[i].start] == 1)
      {
        ctg[i].lb_cost_r *= 2;
        ctg[i].ub_cost_r *= 2;
      }
    }
  }

  for (int i = 0; i < n_vars; i++)
  {
    int var_id = order_dv[i];

    for (int j = 0; j < v_dis_var[var_id].v_values.size(); j++)
    {
      for (int k = 0; k < v_dis_var[var_id].vv_id_mute_edge[j].size(); k++)
      {
        int id_edge = v_dis_var[var_id].vv_id_mute_edge[j][k];
        if (id_edge < 0) // mute ctg
        {
          id_edge = -id_edge;
          id_edge--;
          ctg[id_edge].lb_cost_r = (i + 1) * 10;
          ctg[id_edge].ub_cost_r = (i + 1) * 10;
        }
      }
    }
  }

}

void stnu::sortCtgs()
{
  vector<bool> v;
  for (int i = 0; i < ctg.size(); i++)
    v.push_back(false), v_order_ctg.push_back(-1);

  for (int k = 0; k < ctg.size(); k++)
  {
    for (int i = 0; i < ctg.size(); i++)
    {
      if (v[i]) continue;
      bool fl = false;
      for (int j = 0; j < ctg.size(); j++)
      {
        if (i == j || v[j]) continue;
        if (i_precede[ctg[j].end][ctg[i].start] == 1)
        {
          fl = true;
          break;
        }
      }

      if (fl == false)
      {
        v_order_ctg[i] = k;
	//  printf("c%d o%d ",i,k);
        v[i] = true;
        break;
      }
    }
  }
//  printf("\n");
}

void stnu::TransGraph() //Floyd
{
  if (i_precede != NULL)
  {
    //for(int i = 0; i < num_nodes_;i++)
    //	delete i_precede[i];
    //	delete i_precede;
  }

  i_precede = new int*[this->n_nodes];

  for (int i = 0; i < n_nodes; i++)
  {
    i_precede[i] = new int[n_nodes];

    for (int j = 0; j < n_nodes; j++)
    {
      i_precede[i][j] = 0;
    }
  }

  for (int i = 0; i < ctg.size(); i++)
  {
    int st = ctg[i].start;
    int en = ctg[i].end;
    i_precede[st][en] = 1;
    i_precede[en][st] = -1;
  }

  for (int i = 0; i < rqm.size(); i++)
    if (rqm[i].lb >= 0)
    {
      int st = rqm[i].start;
      int en = rqm[i].end;
      i_precede[st][en] = 1;
      i_precede[en][st] = -1;
    }

  for (int k = 0; k < n_nodes; k++)
  {
    for (int i = 0; i < n_nodes; i++)
    {
      if (i == k || i_precede[i][k] == 0) continue;
      for (int j = 0; j < n_nodes; j++)
      {
        if (i == j || j == k) continue;
        if (i_precede[j][k] == 0) continue;
        if (i_precede[i][k] == i_precede[k][j])
        {
          i_precede[i][j] = i_precede[i][k];
          i_precede[j][i] = -i_precede[i][j];
        }
      }
    }
  }

  sortCtgs();
  for (int i = 0; i < n_nodes; i++)
  {
    cnt_pri_node.push_back(0);
    for (int j = 0; j < n_nodes; j++)
    {
      if (i_precede[j][i] == 1)
      {
        cnt_pri_node[i]++;
      }
    }
    if(cnt_pri_node[i] == 0)
      start_node = i;
  }
//
//	for(int i = 0; i < n_nodes; i++)
//	{
//		for(int j = 0; j < n_nodes; j++)
//			printf("%2d ",i_precede[i][j]);
//		printf("\n");
//	}
}

void stnu::priNode(int &pri, int current) const
{
  if (pri == -1)
  {
    pri = current;
    return;
  }

  if (i_precede[pri][current] == 1) return;
  if (i_precede[current][pri] == 1)
  {
    pri = current;
    return;
  }

  if (cnt_pri_node[pri] > cnt_pri_node[current])
  {
    pri = current;
  }

}

void stnu::setCtgUb(double x)
{
  for (int i = 0; i < n_ctg; i++)
  {
    this->ctg[i].ub = this->ctg[i].lb + x;
  }
}

void stnu::setCtgLb(double x)
{
  for (int i = 0; i < n_ctg; i++)
  {
    this->ctg[i].lb = this->ctg[i].ub - x;
    if (this->ctg[i].lb < 0) this->ctg[i].lb = 0;
  }
}

void stnu::setCtg(double x)
{
  if (this->s_info.id_obj == O_MAX_DELAY)
    setCtgUb(x);
  else if (this->s_info.id_obj == O_MAX_EARLINESS) setCtgLb(x);
}

void stnu::setIdLr()
{
  for (int i = 0; i < n_ctg; i++)
  {
    ctg[i].id_lr = 2 * i;
  }
  for (int i = 0; i < n_rqm; i++)
  {
    rqm[i].id_lr = 2 * i + 2 * n_ctg;
  }
}

void stnu::setIdUr()
{
  for (int i = 0; i < n_ctg; i++)
  {
    ctg[i].id_ur = 2 * i + 1;
  }
  for (int i = 0; i < n_rqm; i++)
  {
    rqm[i].id_ur = 2 * i + 2 * n_ctg + 1;
  }
}

void stnu::add_ctg(CTG c)
{
  c.setOriginal();
  c.b_control = false;
  ctg.push_back(c);
  ctg[ctg.size() - 1].id = ctg.size() - 1;
}

void stnu::add_ctg(int st, int en, double p[])
{
  CTG c = CTG(st, en, p);
  this->add_ctg(c);
}

void stnu::add_ctg(int st, int en, double l, double u, double lc, double uc)
{
  CTG c = CTG(st, en, l, u, lc, uc);
  this->add_ctg(c);
}

void stnu::add_ctg(int st, int en)
{
  CTG c = CTG(st, en);
  this->add_ctg(c);
}

void stnu::add_rqm(RQM r)
{
  r.setOriginal();
  r.b_control = true;
  r.id = rqm.size();
  r.id_lr = 2 * r.id + n_ctg * 2;
  r.id_ur = r.id_lr + 1;
  rqm.push_back(r);
  //rqm[rqm.size() - 1].id = rqm.size() - 1;

}

void stnu::add_rqm(int st, int en)
{
  RQM r = RQM(st, en);
  this->add_rqm(r);
}

void stnu::add_rqm(int st, int en, double l, double u, double lc, double uc)
{
  RQM r = RQM(st, en, l, u, lc, uc);
  this->add_rqm(r);
}

void stnu::print()
{
  for (int i = 0; i < n_ctg; i++)
  {
    ctg[i].print();
  }
  for (int i = 0; i < n_rqm; i++)
    rqm[i].print();
}

void stnu::print_dot(std::string file_name)
{
  FILE *fdot;
  fdot = fopen(file_name.c_str(), "w");
  fprintf(fdot, "digraph G {\n nodesep = .45; \n "
      "size = 30;\nlabel=\"%s\";\n", this->file_name.c_str());

  for (int i = 0; i < this->n_ctg; i++)
  {
    CTG tmp_edge;
    tmp_edge = ctg[i];
    if (tmp_edge.b_mute) continue;
    if (tmp_edge.lb > tmp_edge.ub)
    {
      tmp_edge.print();
      exit(0);
    }

    if (tmp_edge.start_name.length() == 0)
      fprintf(fdot, "\"%d\"->\"%d\"[style=dotted label = \"[%.2lf,%.2lf] ",
          tmp_edge.start, tmp_edge.end, tmp_edge.lb, tmp_edge.ub);
    else fprintf(fdot, "\"%s\"->\"%s\"[style=dotted label = \"[%.2lf,%.2lf] ",
        tmp_edge.start_name.c_str(), tmp_edge.end_name.c_str(), tmp_edge.lb,
        tmp_edge.ub);
    for (int j = 0; j < tmp_edge.v_labels.size(); j++)
    {
      fprintf(fdot, "%d=%d ", tmp_edge.v_labels[j].first,
          tmp_edge.v_labels[j].second);
    }
    fprintf(fdot, "\nr %.2lf,%.2lf ", tmp_edge.lb_cost_r, tmp_edge.ub_cost_r);

    fprintf(fdot, "\"];\n");
  }

  for (int i = 0; i < this->n_rqm; i++)
  {
    RQM tmp_edge;
    tmp_edge = this->rqm[i];
    if (tmp_edge.b_mute) continue;
    if (tmp_edge.start_name.length() == 0)
      fprintf(fdot, "\"%d\"->\"%d\"[label = \"[%.2lf,%.2lf] ", tmp_edge.start,
          tmp_edge.end, tmp_edge.lb, tmp_edge.ub);
    else fprintf(fdot, "\"%s\"->\"%s\"[label = \"[%.2lf,%.2lf] %d ",
        tmp_edge.start_name.c_str(), tmp_edge.end_name.c_str(), tmp_edge.lb,
        tmp_edge.ub, 0);
    for (int j = 0; j < tmp_edge.v_labels.size(); j++)
    {
      fprintf(fdot, "%d=%d ", tmp_edge.v_labels[j].first,
          tmp_edge.v_labels[j].second);
    }
    fprintf(fdot, "\nr %.2lf,%.2lf ", tmp_edge.lb_cost_r, tmp_edge.ub_cost_r);
    fprintf(fdot, "\"];\n");

  }
  fprintf(fdot, "}\n");
  fclose(fdot);
}

void stnu::print_cstnu(std::string file_name)
{
  FILE * fp;
  fp = fopen(file_name.c_str(), "w");

  fprintf(fp, "%d	%d	%d %lf\n", n_nodes, n_ctg + n_rqm, n_vars, utility);

  for (int i = 0; i < n_vars; i++)
  {
    fprintf(fp, "%d	", v_dis_var[i].n_value);
    for (int j = 0; j < v_dis_var[i].n_value; j++)
    {
      fprintf(fp, "%d	", v_dis_var[i].v_values[j]);
    }
    for (int j = 0; j < v_dis_var[i].n_value; j++)
    {
      fprintf(fp, "[%lf]	", v_dis_var[i].v_utilities[j]);
    }
    fprintf(fp, "\n");
  }

  for (int i = 0; i < n_rqm; i++)
  {
    RQM tmp_edge = rqm[i];
    vector<int> labels;
    //labels = get_labels(tmp_edge);
    fprintf(fp, "%d	%d	R	%lf	%lf", tmp_edge.start, tmp_edge.end, tmp_edge.lb,
        tmp_edge.ub);
    fprintf(fp, "	%lf	%lf	%d\n", tmp_edge.lb_cost_r, tmp_edge.ub_cost_r,
        tmp_edge.v_labels.size());

    if (tmp_edge.v_labels.size())
    {
      for (int i = 0; i < tmp_edge.v_labels.size(); i = i + 1)
      {
        fprintf(fp, "%d	[%d]	", tmp_edge.v_labels[i].first,
            tmp_edge.v_labels[i].second);
        //fprintf(fp, "%d	[%d]	", labels[i], labels[i+1]);
      }
      fprintf(fp, "\n");
    }
  }

  for (int i = 0; i < n_ctg; i++)
  {
    CTG tmp_edge = ctg[i];
    vector<int> labels;
    //labels = get_labels(tmp_edge);
    fprintf(fp, "%d	%d	C	%lf	%lf", tmp_edge.start, tmp_edge.end, tmp_edge.lb,
        tmp_edge.ub);
    fprintf(fp, "	%lf	%lf	%d\n", tmp_edge.lb_cost_r, tmp_edge.ub_cost_r,
        tmp_edge.v_labels.size());

    if (tmp_edge.v_labels.size())
    {
      for (int i = 0; i < tmp_edge.v_labels.size(); i = i + 1)
      {
        fprintf(fp, "%d	[%d]	", tmp_edge.v_labels[i].first,
            tmp_edge.v_labels[i].second);
        //fprintf(fp, "%d	[%d]	", labels[i], labels[i+1]);
      }
      fprintf(fp, "\n");
    }
  }
  fclose(fp);

}

std::vector<int> stnu::get_labels(LINK tmp_edge)
{
  // return the labels of link tmp_edge
  // the labels are the assignments which activate tmp_edge
  // *tmp_ege.labels is a record of assignments which mute tmp_edge
  vector<int> res;
  res.clear();

  bool *b_var_in_label = new bool[n_vars];
  bool **b_value_label = new bool*[n_vars];
  int b_size = 0;
  for (int i = 0; i < n_vars; i++)
  {
    b_value_label[i] = new bool[100];
    memset(b_value_label[i], false, sizeof(b_value_label[i]));
  }

  memset(b_var_in_label, false, sizeof(b_var_in_label));
  //memset(b_value_label, false, sizeof(b_value_label)) ;

  for (int j = 0; j < tmp_edge.v_labels.size(); j++)
  {
    b_var_in_label[tmp_edge.v_labels[j].first] = true;
    b_value_label[tmp_edge.v_labels[j].first][tmp_edge.v_labels[j].second] =
        true;
    //fprintf(fp, "%d	[%d]	", tmp_edge.v_labels[j].first, tmp_edge.v_labels[j].second);
  }
  for (int j = 0; j < n_vars; j++)
  {
    if (b_var_in_label[j] == true)
    {
      for (int k = 0; k < v_dis_var[j].n_value; k++)
      {
        if (b_value_label[j][k] == false)
        {
          res.push_back(j);
          res.push_back(k);
        }
      }
    }
  }

  return res;
}

void stnu::print_grid(std::string file_name)
{
  FILE *fgrid;
  fgrid = fopen(file_name.c_str(), "w");
  fprintf(fgrid, "p dc %d	1 %d\n", this->n_nodes, this->n_ctg);
  for (int i = 0; i < this->n_ctg; i++)
  {
    CTG tmp_ctg;
    tmp_ctg = this->ctg[i];
    fprintf(fgrid, "a	%d	%d	%lf	%lf\n", tmp_ctg.start + 1, tmp_ctg.end + 1,
        tmp_ctg.lb, tmp_ctg.ub);
  }

  for (int i = 0; i < this->n_rqm; i++)
  {
    RQM tmp_rqm;
    tmp_rqm = this->rqm[i];
    fprintf(fgrid, "a	%d	%d	%lf	%lf\n", tmp_rqm.start + 1, tmp_rqm.end + 1,
        tmp_rqm.lb, tmp_rqm.ub);
  }

  for (int i = 0; i < this->n_ctg; i++)
  {
    CTG tmp_ctg;
    tmp_ctg = this->ctg[i];
    fprintf(fgrid, "g	%d	%d\n", tmp_ctg.start + 1, tmp_ctg.end + 1);
  }

  fclose(fgrid);
}

bool stnu::removeCTG(int x)
{
  for (int i = 0; i < n_ctg; i++)
  {
    if (ctg[i].id == x)
    {
      ctg.erase(ctg.begin() + i);
      n_ctg--;
      return true;
    }
  }
  return false;
}

bool stnu::removeRQM(int id)
{
  for (int i = 0; i < n_rqm; i++)
  {
    if (rqm[i].id == id)
    {
      rqm.erase(rqm.begin() + i);
      n_rqm--;
      return true;
    }
  }

  return false;
}

void stnu::relaxMission(const int &i_plan)
{
  // change key upper bounds for AUV file
  for (int i = 0; i < n_rqm; i++)
  {
    if (this->rqm[i].b_mission == true)
    {
      //printf("%d %d\n",rqm[i].start, rqm[i].end);
      changeUb(rqm[i], i, i_plan);
    }
  }
}

// change key upper bounds for AUV file
void stnu::changeUb(RQM & r, const int id, const int &i_plan)
{
  double *dis = new double[this->n_nodes];
  double *disl = new double[this->n_nodes];
  bool *v = new bool[n_nodes];
  bool *vl = new bool[n_nodes];
  for (int i = 0; i < n_nodes; i++)
  {
    dis[i] = g_inf;
    disl[i] = -1;
    v[i] = false;
    vl[i] = false;
  }

  dis[r.start] = 0;
  disl[r.start] = 0;
  bool fl = true;

  do
  {
    fl = false;
    int next = -1;
    int nextl = -1;
    for (int i = 0; i < n_nodes; i++)
    {
      if (v[i] == false && dis[i] < g_inf && (next == -1 || dis[next] > dis[i]))
      {
        next = i;
        fl = true;
      }
      if (vl[i] == false && disl[i] >= 0
          && (nextl == -1 || dis[nextl] < disl[i]))
      {
        nextl = i;
        fl = true;
      }

    }
    //printf("%d %d\n",next, nextl);
    //if(next == r.end)
    //	break;

    if (next >= 0) v[next] = true;
    if (nextl >= 0) vl[nextl] = true;
    if (nextl < 0 && next < 0) break;
    for (int i = 0; i < n_rqm; i++)
    {
      if (next != -1)
        if (i != id && rqm[i].start == next && v[rqm[i].end] == false)
        {
          if (dis[rqm[i].end] > dis[rqm[i].start] + rqm[i].lb)
            dis[rqm[i].end] = dis[rqm[i].start] + rqm[i].lb;
        }

      if (nextl != -1)
        if (i != id && rqm[i].start == nextl && vl[rqm[i].end] == false)
        {
          if (disl[rqm[i].end] < disl[rqm[i].start] + rqm[i].ub)
            disl[rqm[i].end] = disl[rqm[i].start] + rqm[i].ub;
        }
    }

    for (int i = 0; i < n_ctg; i++)
    {
      if (next != -1)
        if (ctg[i].start == next && v[ctg[i].end] == false)
        {
          if (dis[ctg[i].end] > dis[ctg[i].start] + ctg[i].ub)
            dis[ctg[i].end] = dis[ctg[i].start] + ctg[i].ub;
        }

      if (nextl != -1)
        if (ctg[i].start == nextl && vl[ctg[i].end] == false)
        {
          if (disl[ctg[i].end] < disl[ctg[i].start] + ctg[i].lb)
            disl[ctg[i].end] = disl[ctg[i].start] + ctg[i].lb;
        }

    }

//		for(int i = 0; i < n_nodes;i++)
//		{
//			printf("%lf ",dis[i]);
//		}
//		printf("\n");
//		for(int i = 0; i < n_nodes; i++)
//			printf("%lf ",disl[i]);
//		printf("\n");

  } while (fl);
  double ku[] = { 1, 1.2, 1.5 };
  //int i_u = i_plan%3
  r.ub = dis[r.end] * ku[i_plan % 3];
  r.lb = (disl[r.end]) / ku[i_plan / 3];
  //r.ub_o = r.ub;
  //r.lb_o = r.lb;
//	printf("%lf %lf\n", r.lb, r.ub);

  //exit(0);
  //delete dis;
  //delete v;
}
// end change key upper bounds for AUV file 20160421
int stnu::edgeToOid(const LINK &l, bool b_ub) const
{
  int id = l.id_o + 1;	//ctg: 1..n_ctg
  if (l.b_control == true) //rqm: n_ctg+1 .. n_ctg+n_rqm
  {
    id += n_ctg;
  }

  if (b_ub) id = -id;
  //printf("id=%d lid=%d ctg=%d rqm=%d\n",id, l.id, n_ctg,n_rqm);
  return id;
}

int stnu::edgeToId(const LINK &l, bool b_ub) const
{
  int id = l.id + 1; //ctg: 1..n_ctg
  if (l.b_control == true) //rqm: n_ctg+1 .. n_ctg+n_rqm
  {
    id += n_ctg;
  }

  if (b_ub) id = -id;
//	if(l.start == 7 && l.end == 10)
//	printf("id = %d oid=%d %d %d\n",id, l.id,l.start, l.end);
  return id;
}

std::pair<LINK, bool> stnu::idToEdge(int id) const
{
  bool b_ub = false;
  if (id < 0) id = -id, b_ub = true;
  id--;
  LINK res(-1, -1);
  if (id >= n_ctg)
  {
    if (id - n_ctg > n_rqm)
    {
      printf("error: id_to_edge got invalid parameters!\n");
      exit(0);
    }
    if (id - n_ctg != n_rqm) res = rqm[id - n_ctg];
  }
  else
  {
    return std::make_pair(ctg[id], b_ub);
    //res = ctg[id];
  }

  return std::make_pair(res, b_ub);
}

void stnu::exp4mip()
{
  /* expand the networks to deal with mip model with dynamic options
   * For c in DisVars:
   * 	For n in Nodes:
   * 		If n is not related to A(c) and not before D(c):
   * 			make |D(c)|-1 repetitions
   *  For e in E:
   *  	If e is not related to A(c) and not before D(c):
   *  		make |D(c)|-1 repetitions and attached to the new nodes
   * */

  // init labels
  for (int i = 0; i < n_nodes; i++)
  {
    node_label.push_back(vector<int>());
    for (int j = 0; j < n_vars; j++)
      node_label[i].push_back(0);
    father_node.push_back(i);
    node_rep.push_back(vector<int>());
    node_rep[i].push_back(i);
  }

  for (int i = 0; i < n_vars; i++)
  {
	  int dt_node = ctg[v_dis_var[i].v_dt[0]].end;
	  for(int j = 0 ;j < n_nodes; j++)
	  {
		  if(i_precede[j][dt_node] >= 0) //j is before dt
			  node_label[j][i] = -1;
	  }
  }

  for (int i = 0; i < n_ctg; i++)
  {
    for (int j = 0; j < ctg[i].v_labels.size(); j++)
    {
      int v_id = ctg[i].v_labels[j].first;
      int value = ctg[i].v_labels[j].second;
      node_label[ctg[i].start][v_id] |= 1 << value;
      node_label[ctg[i].end][v_id] |= 1 << value;
    }
  }

  for (int i = 0; i < n_rqm; i++)
  {
    for (int j = 0; j < rqm[i].v_labels.size(); j++)
    {
      int v_id = rqm[i].v_labels[j].first;
      int value = rqm[i].v_labels[j].second;
      node_label[rqm[i].start][v_id] |= 1 << value;
      node_label[rqm[i].end][v_id] |= 1 << value;
    }
  }
  if(s_info.id_debug == DEBUG_T && n_vars)
  {
	  printf("n_nodes=%d\n",n_nodes);
	  for(int j = 0 ; j < n_nodes; j++)
		  printf("ol[%d]=%d ",j,node_label[j][0]);
	  printf("\n");
  }
  // main process
  for (int i = 0; i < n_vars; i++)
  {
	int a_n_nodes=n_nodes;
    for (int j = 0; j < a_n_nodes; j++)
    {
      if (node_label[j][i] == 0) //|| v_dis_var[i]._dt == j)
      {
        for (int k = 1; k < v_dis_var[i].n_value; k++)
        {
          add_node(j, i, k);
        }
        node_label[j][i] = 1;
      }
      else if(node_label[j][i] > 0)
      {
    	  int first_value = -1;
          for (int k = 0; k < v_dis_var[i].n_value; k++)
          {
        	  if(first_value>=0 && (node_label[j][i]&(1<<k)))
        	  {
        		  add_node(j, i, k);

        	  }
        	  if(first_value<0 && (node_label[j][i]&(1<<k)))
        		  first_value = k;
          }
          node_label[j][i] = (1<<first_value);
      }
    }
   // n_nodes+=a_n_nodes;
//    printf("n_nodes=%d\n",n_nodes);
 //   for(int j = 0 ; j < n_nodes; j++)
//    	printf("l[%d]=%d ",j,node_label[j][i]);
 //   printf("\n");
  }

  for(int i = 0; i < n_vars; i++)
  {
    for(int j = 0; j < v_dis_var[i].v_dt.size(); j++)
    {
    	int ori_dt = ctg[v_dis_var[i].v_dt[j]].end;
    	int la = node_rep[ori_dt].size();
    	for(int l = 0; l < la; l++)
    	{
    		int tmp_dt = node_rep[ori_dt][l];
    		node_label[tmp_dt][i] = 0;
    		for(int k = 1; k < v_dis_var[i].n_value; k++)
    			add_node(tmp_dt, i, k);
    		node_label[tmp_dt][i] = 1;
    	}
    }
  }

  for (int i = 0; i < n_ctg; i++)
  {
    ctg_rep.push_back(vector<int> ());
    for (int j = 0; j < node_rep[ctg[i].start].size(); j++)
    {
      for (int k = 0; k < node_rep[ctg[i].end].size(); k++)
      {
	//	if(!(j|k)) continue;
        int st = node_rep[ctg[i].start][j];
        int en = node_rep[ctg[i].end][k];
        if (check_labels(node_label[st], node_label[en]))
        {

          add_rep_link(ctg[i], st, en);
          ctg_rep[i].push_back(ctg.size()-1-n_ctg);
        }
      }
    }
  }

  for (int i = 0; i < n_rqm; i++)
  {
    for (int j = 0; j < node_rep[rqm[i].start].size(); j++)
    {
      for (int k = 0; k < node_rep[rqm[i].end].size(); k++)
      {
	//	if(!(j|k))continue;
        int st = node_rep[rqm[i].start][j];
        int en = node_rep[rqm[i].end][k];
        if (check_labels(node_label[st], node_label[en]))
        {
          add_rep_link(rqm[i], st, en);
        }
      }
    }
  }
  ctg.erase(ctg.begin(), ctg.begin()+n_ctg);
  rqm.erase(rqm.begin(), rqm.begin()+n_rqm);
  n_ctg = ctg.size();
  n_rqm = rqm.size();
  for(int i = 0; i < n_ctg; i++)
    ctg[i].id = i;
  for(int i = 0; i < n_rqm; i++)
    rqm[i].id = i;
  // printf("n_ctg = %d n_rqm = %d\n",n_ctg, n_rqm);
  
}

void stnu::add_node(int ori_id, int v_id, int value)
{
  father_node.push_back(father_node[ori_id]);
  node_label.push_back(vector<int>());
  int tmp = father_node[ori_id];
  node_rep[tmp].push_back(n_nodes);

//  node_rep[tmp].push_back(n_nodes);
  for (int i = 0; i < n_vars; i++)
  {
    node_label[n_nodes].push_back(node_label[ori_id][i]);
  }

  node_label[n_nodes][v_id] = 1 << value;

//  printf("add node %d[%d]", n_nodes, father_node[n_nodes]);
//  for(int i = 0; i < n_vars; i++)
//   printf(" %d=%d",i,node_label[n_nodes][i]);
//  printf("\n");

  ++n_nodes;
}

void stnu::add_rep_link(LINK a, int st, int en)
{
  a.start = st;
  a.end = en;
  a.start_name = a.end_name = "";
//a.print();
  vector<int> ori_labels;
  for(int i = 0; i < n_vars; i++)
    ori_labels.push_back(0);
  for(int i = 0; i < a.v_labels.size(); i++)
    ori_labels[a.v_labels[i].first] = 1;
  for(int i = 0; i < n_vars; i++)
  {
    int v_exp = node_label[st][i] & node_label[en][i];
 //   printf("%d %d %d\n",node_label[st][i],node_label[en][i],v_exp);
    if(v_exp == 0 || ori_labels[i] || v_exp < 0) continue;
    int value = -1;
    int now_value = 0;
    while(v_exp >= 1)
    {
      if((1&v_exp) && value >=0 )
           break;
      if(1&v_exp)
	  value = now_value;
      now_value++;
      v_exp=v_exp>>1;
    }
    if(v_exp < 1)
     a.v_labels.push_back(make_pair(i, value));
  }
 // a.print();
  if(a.b_control)
    add_rqm(a);
  else add_ctg(a);
}

bool stnu::check_labels(vector<int> a, vector<int> b)
{
  if (a.size() != b.size())
  {
    printf("error! the labels of two nodes are in different length!\n");
    printf("l(a) = %d l(b) = %d\n", a.size(), b.size());
    exit(0);
  }
  
  for (int i = 0; i < a.size(); i++)
  {
    if (a[i] == 0 || b[i] == 0) // don't have labels on c_i
    continue;
    if ((a[i] & b[i]) == 0) // different
      return false;
  }
  return true;
}

stnu::~stnu()
{
  // TODO Auto-generated destructor stub
}

} /* namespace CDS */
