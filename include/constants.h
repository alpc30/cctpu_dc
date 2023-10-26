#ifndef CONSTANS_H_
#define CONSTANS_H_

#include <utility>
#include <set>
#include <vector>
#include <map>
#include "updatequeue.h"
#include <climits>
#include <cfloat>
#include <string>
#include <map>
#include <vector>
#include <cassert>
//#include "search_delta.h"

namespace TRANSPLAN
{

  //GOAL == dummy gaol transition
  //INIT == dummy init transition
  //EFFT == effet transitions
  //PRV  == prevail transitions
  enum   EFFT_TYPE {GOAL, INIT, EFFT, PRV};
  
  //SV_T == state variable resource transition
  //U_T  == unary resource transiton
  //B_T == borrow resource // needes to be decomposed to U_T
  //NUM_T == numeric resource transition //needs to be decomposed
  //DECOM_T == decomposed numeric resource transition //special type of U_T, where it has a parent transition
  //SEQ_NUM_T == sequential numeric resoure transition
  //PROD_T == produce transitions
  //CONS_T == consume transitions
  enum   TRANS_TYPE {SV_T, U_T, B_T, NUM_T, DECOM_T, SEQ_NUM_T};

  //Resource types: SV-State variable, 
  //U-Unary, SD-Sequential Discrete, SR-Sequential Reservior, B-Parallel Borrow, D-Discrete, R-Reservior
  enum   RES_TYPE {SV_RES, U_RES, SD_RES, SR_RES, B_RES, D_RES, R_RES};

  //requirement types: A: Assignments, P: Production, C:Consumption, NA:not applicable, INV: invariant prevail
  enum   REQ_TYPE {NA, A, P, C, INV};

  enum   CONSTR {exactly,
		 min_change,
		 min_achieve,
		 max_change,
		 max_achieve,
		 change_before,
		 change_after,
		 achieve_before,
		 achieve_after,
		 before,
		 after,
		 start_before,
		 start_after,
		 end_before,
                 end_after,
		 non_goal};

  enum   CSP_MOD_EVENT {IGNORED, SUCCESS, FAILURE};

  enum   PROP_STATUS {FIX, FAIL};

  enum   CSP_SOL {SOLVED, NO_SOLUTION};

  enum   PROP_TYPE {SUPP, TEMP, CHANGE, PREC, BOUND, INCL_EXCL, LINEAR, BOUND_INFER};

  enum   NUM_OP {EQ, LQ, GQ, LE, GR}; //==, <=, >=, <, >

  //this is to identify which type of interval is it, T-TIME, C-COST, and L-Level
  enum   INTERVAL_VAR_TYPE {TIME, COST, LEVEL};

  typedef std::set<int>                 IntSet;
  typedef std::vector<int>              IntVector;
  typedef std::vector<double>			DblVector;
  typedef std::vector<IntSet>           IntSetVector;
  typedef std::vector<IntVector>        IntMatrix;
  typedef std::vector<DblVector>		DblMatrix;
  typedef std::pair<int, int>           IntPair;
  typedef std::set<IntPair>             IntPairSet;
  typedef std::vector<IntPair>          IntPairVector;
  typedef std::pair<std::string, int>   NameIndexPair;
  typedef std::map<int, int>            IntIntMap;
  typedef std::map<int, IntSet>         IntIntSetMap;

  typedef std::pair<CONSTR, std::pair<int,int> > BinaryConstraint;
  typedef std::pair<CONSTR, int>                 UnaryConstraint;

  typedef UpdateQueue<IntPair> IntPairStack ;
  typedef UpdateQueue<int>     IntStack;

  
  class CSPIntVar;
  class CSPIntervalVar;
  class CSPBoolVar;
  class CSPMaxCounterVar;
  class CSPAuxIntVar;
  class CSPAuxSetVar;
  class CSPAuxSetMonoIncVar;
  class CSPAuxSetMonoDecVar;
  class CSPAuxPairSet;

  typedef std::set<CSPIntVar>              CSPIntSet;
  typedef std::set<CSPIntervalVar>         CSPIntervalSet;
  typedef std::set<CSPBoolVar>             CSPBooleanSet;
  typedef std::set<CSPMaxCounterVar>       CSPMaxCounterSet;
  typedef std::set<CSPAuxSetVar>           CSPAuxSetSet;
  typedef std::set<CSPAuxIntVar>           CSPAuxIntSet;
  typedef std::set<CSPAuxSetMonoIncVar>    CSPAuxSetMonoIncSet;
  typedef std::set<CSPAuxSetMonoDecVar>    CSPAuxSetMonoDecSet;
  typedef std::set<CSPAuxPairSet>          CSPAuxPairSetSet;
  
  typedef std::vector<CSPIntVar>           CSPIntVector;
  typedef std::vector<CSPIntervalVar>      CSPIntervalVector;
  typedef std::vector<CSPBoolVar>          CSPBooleanVector;
  typedef std::vector<CSPMaxCounterVar>    CSPMaxCounterVector;
  typedef std::vector<CSPAuxIntVar>        CSPAuxIntVector;
  typedef std::vector<CSPAuxSetVar>        CSPAuxSetVector;
  typedef std::vector<CSPAuxSetMonoIncVar> CSPAuxSetMonoIncVector;
  typedef std::vector<CSPAuxSetMonoDecVar> CSPAuxSetMonoDecVector;
  typedef std::vector<CSPAuxPairSet>       CSPAuxPairSetVector;
  


  //special CSP domain values
  const int DUMMY_VAL        = -1;
  const int NOT_IN_PLAN      = -2;

  //csp assignment status
  const int CSP_NOT_ASSIGNED = -5; //if a varaible is not yet assigned
  const int CSP_ASSIGNED     = -6; //if a varaible is already assigned

  //values retrun by branching
  const int CSP_NO_NEXT_VAR  = -12; // all assigned
  const int CSP_NO_NEXT_VAL  = -13; // all vlaue exhausted==fail


  //dependencies between transitions, for heuristics
  const int TRANS1_INDEPENDENT_TRANS2 = 0;
  const int TRANS1_PRECEDE_TRANS2     = 1;
  const int TRANS1_SUCCEED_TRANS2     = 2;

  //dependencies between state variables, for heuristics
  const int VARS_INDEPENDENT          = 0;
  const int VARS_INTER_DEPENDENT      = 1;
  const int VAR1_SUPPORTS_VAR2        = 2;
  const int VAR1_DEPENDS_ON_VAR2      = 3;

  //Default planning horizon
  const int DEFAULT_RELEASE_DATE      = 0;
  const int DEFAULT_DEADLINE          = INT_MAX;
  const int DEFAULT_MAX_COST          = INT_MAX;
  const int DEFAULT_MAX_LEVEL         = INT_MAX;

  //parent for temporal trees, est-tree and let-tree
  const int NULL_PARENT               = -1;


  const int EMPTY_HEAP                = -100;

}
#endif
