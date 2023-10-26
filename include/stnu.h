/*
 * stnu.h
 *
 *  Created on: 18/07/2014
 *      Author: jing
 */

#ifndef STNU_H_
#define STNU_H_

#include "numeric_type.h"
#include "tcn.h"
#include "head.h"

const double g_inf = std::numeric_limits<double>::infinity();

//define the objective functions
#define O_MAX_DELAY 00
#define O_TOT_DELAY 10

#define O_MAX_EARLINESS 20

#define O_MAKESPAN 1
#define O_RELAX_COST 2
#define O_RELAX_COST_NR 12 //relax cost without rewards
#define O_MIN_FLEX 3
#define O_MIN_RISK 4
#define O_DC_CHECK 5
#define O_EP_FLEX 7

#define O_MAX_FLEX 60
#define O_MAX_TOT_FLEX 61
#define O_MAX_MIN_FLEX 62

//define the types of contingent links
#define C_LU 0
#define C_CHANCE 1

//solvers
#define S_CDS 0
#define S_SNOPT_M 1
#define S_SNOPT_S 2
#define S_DC_CHECK 3
#define S_BASELINE 5
#define S_BISEARCH 6
#define S_CDS_E 10
//initializations
#define I_LOOSE 0
#define I_TIGHT 1
#define I_RAND 2

//constraints
#define C_DC 0
#define C_SC 1
#define C_CONSIST 2
#define C_SCDC 10

//debug
#define DEBUG_F 0
#define DEBUG_C 1 //only show the combining process
#define DEBUG_T 2
#define DEBUG_K -1

//variable
#define V_DC 0
#define V_SC 1

namespace CDS
{

class S_PARAM
{
public:
  int id_obj;
  int id_solver;
  int id_ctg;
  int id_constrs;
  int id_init;
  int id_debug;
  double extra;
  int id_variable;
  double f_runtime;
  S_PARAM()
  {
    id_obj = id_solver = id_ctg = id_constrs = id_init = extra = 0;
    id_debug = 0;
    f_runtime = 1000000000;
    id_variable = 0;
  }
  ;
};

//define the structure of links (controllable & uncontrollable)
typedef class LINK
{
public:
  std::string start_name, end_name;
  int start, end;
  double lb_cost_r, ub_cost_r;
  double lb, ub;
  double lb_lb, ub_lb; //lower and upper bounds of lb
  double lb_ub, ub_ub; //lower and upper bounds of ub
  //double lb_o, ub_o;
  //int i_lr,i_ur;
//	int i_lt,i_ut;
  int id; //edge_id in current stnu
  int id_o; //edge_id in original cctpu
  bool b_control; //ctg-false
  bool b_mute;
  int id_lr, id_ur;

  // in the AVU file, the constraint with name "Missionx-xxx" is a
  // key upper bound
  bool b_mission;

  std::vector<std::pair<int, int> > v_labels; //from 0
  std::vector<int> v_pri_vars; //20160201 the list of variables following the ctg

  int get_start()
  {
    return start;
  }
  int get_end()
  {
    return end;
  }

  //double param[5];
  double (*G)(double, double *);
  double (*F)(double, double *);
  double (*dFdp)(double, double *);

  LINK()
  {
  }
  ;
  void setOriginal()
  {
    b_mute = false;
    lb_lb = lb;
    ub_lb = ub;
    lb_ub = lb;
    ub_ub = ub;
    //lb_o = lb;
    //	ub_o = ub;
  }

  //basic construct
  LINK(int st, int en)
  {
    start = st, end = en;
    lb_cost_r = 1e20;
    ub_cost_r = 1e20;
    lb = -1e20;
    ub = 1e20;
    v_labels.clear();
    b_control = true;
    this->start_name = "";
    this->end_name = "";
    setOriginal();
  }
  ;

  //construct with parameters
  LINK(int st, int en, double p[])
  {
    start = st, end = en;
    //for(int i = 0; i < 5; i++)
    //	param[i] = p[i];
    setOriginal();
  }

  //construct with bounds
  LINK(int st, int en, double l, double u)
  {
    start = st;
    end = en;
    if (l < -1e10) l = -1e10;
    lb = l;
    if (u > 1e10) u = 1e10;
    ub = u;
    lb_cost_r = 1e20;
    ub_cost_r = 1e20;
    this->start_name = "";
    this->end_name = "";
    setOriginal();
  }
  ;

  //construct with bounds and cost ratios
  LINK(int st, int en, double l, double u, double lc, double uc)
  {
    start_name = "s";
    end_name = "e";
    start = st;
    end = en;
    if (l < -1e10) l = -1e10;
    lb = l;
    if (u > 1e10) u = 1e10;
    ub = u;
    lb_cost_r = lc;
    ub_cost_r = uc;
    setOriginal();
  }
  ;

  void print()
  {
    printf("%d->%d [%lf, %lf]", this->start, this->end, lb, ub);
    for(int i = 0; i < this->v_labels.size(); i++)
    	printf(" v(%d)=%d",v_labels[i].first, v_labels[i].second);
    printf(" %c",b_control?'r':'c');
    printf("\n");
  }

  bool isConsist(vector<pair<int, int> > v) const
  {
    for (int i = 0; i < v_labels.size(); i++)
    {
      for (int j = 0; j < v.size(); j++)
      {
        if (v_labels[i].first == v[j].first
            && v_labels[i].second != v[j].second) return false;
      }
    }
    return true;
  }

  //	LINK operator= (LINK a)
//	{
//		this->start = a.start;
//		this->end = a.end;
//		this->
//		return a;
//	};
} CTG, RQM;

class DiscreteVariable
{
public:
  int id;	//node id of the assignment time-point
  std::string name; //node name
  std::string ID;
  int n_value;	//the number of values
  double m_utility; // the maximum utility
  std::vector<int> v_values;
  std::vector<double> v_utilities;
  std::vector<std::vector<int> > vv_id_mute_edge; // the vector of mute edge
  std::vector<int> v_ctg_pri; // the list of ctg prior to the variable
  vector<int> v_dt; // decision node, used for mip_baseline solving dynamic option

  DiscreteVariable()
  {
    v_values.clear();
    v_utilities.clear();
    vv_id_mute_edge.clear();
    n_value = 0;
  }
};

class stnu: public hsps::STN
{

public:
  //double value_preference;

  std::string file_name;
  S_PARAM s_info; //input parameters
  int n_nodes;	//the number of nodes
  int n_ctg;		//the number of contingent links
  int n_rqm;		//the number of requirement links
  int n_vars;		//the number of discrete variables
  int start_node; // the start node

  double utility;	//for cost perference
  double fldt, value_mddc, value_mdsc; //fldt, value of mddc and mdsc;

  int **i_precede; //1-precede, -1-follow, 0-unknown
  int *order_dv; // the order of discrete variables
  std::vector<vector<int> > node_pre_var; // the id of nodes prior to var_i
  std::vector<int> v_order_ctg; //the id of topological order of ctgs
  std::vector<CTG> ctg; //the vector of contingent links (start_node, end_node)
  std::vector<RQM> rqm;	//the vector of rqm
  std::vector<DiscreteVariable> v_dis_var; // the vector of discrete variables
  std::vector<int> cnt_pri_node;
  std::vector<int> father_node; //used for expanded nodes
  vector<vector<int> > node_label;
  vector<vector<int> > node_rep; // recording repetitions
  vector<vector<int> > ctg_rep; // repretitions of ctg

  stnu();
  stnu(std::string problem_file, S_PARAM s_info);
  void initSTNU(S_PARAM &s_info);
  void resetCost();

  void TransGraph(); //get i_precede;
  void sortCtgs(); //get v_order_ctg

  void add_ctg(CTG a);
  void add_ctg(int st, int en);
  void add_ctg(int st, int en, double p[]);
  void add_ctg(int st, int en, double l, double u, double lc, double uc);
  void add_rqm(RQM c);
  void add_rqm(int st, int en);
  void add_rqm(int st, int en, double l, double u, double lc, double uc);

  void setCtgUb(double x);
  void setCtgLb(double x);
  void setCtg(double x);
  void setIdLr();
  void setIdUr();
  bool bEdgePriVar(const LINK &a, const int &v_id);

  // change key upper bounds for AUV file
  void changeUb(RQM & r, const int id, const int & pl);
  void relaxMission(const int& id);
  // end change key upper bounds for AUV file 20160421

  bool removeCTG(int id);
  bool removeRQM(int id);
  void print();

  void print_dot(std::string file_name);
  void print_cstnu(std::string file_name);
  std::vector<int> get_labels(LINK tmp_edge);

  void print_grid(std::string file_name);

  int edgeToId(const LINK & l, bool b_lb) const;
  int edgeToOid(const LINK &l, bool b_lb) const;
  std::pair<LINK, bool> idToEdge(int id) const;

  void priNode(int & pri, int current) const;

  virtual ~stnu();

  void exp4mip();
  void add_node(int ori_id, int v_id, int value);
  void add_rep_link(LINK a, int st, int en);
  bool check_labels(vector<int> a, vector<int> b);
};

} /* namespace CDS */
#endif /* STNU_H_ */
