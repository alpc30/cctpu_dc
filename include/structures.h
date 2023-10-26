#ifndef STRUCTURES_H_
#define STRUCTURES_H_

#include <vector>
#include <map>
#include <set>
#include <stack>
#include <string>
#include <deque>
#include <iostream>
#include "constants.h"


namespace TRANSPLAN
{

#define DEFAULT_BOUND IntPair(0, INT_MAX)

  struct Constraint
  {
    CONSTR constraint;
    int left;
    int right;
    bool unary;

    Constraint(CONSTR constr, int l, int r)
    :constraint(constr)
    ,left(l)
    ,right(r)
    {
      //if unary
      if ( r == -1 )
	unary = true;
      else
	unary = false;
    }

    //copy constructor
    Constraint(const Constraint& c)
    :constraint(c.constraint)
    ,left(c.left)
    ,right(c.right)
    ,unary(c.unary)
    {
      //
    }

    //= operator
    Constraint& operator=(const Constraint& c)
    {
      this->constraint = c.constraint;
      this->left       = c.left;
      this->right      = c.right;
      this->unary      = c.unary;

      return *this;
    }
      

    void print()
    {
      std::cout <<left << ":";
      
      if (constraint == TRANSPLAN::exactly)
	std::cout << "exactly";
      else if (constraint == TRANSPLAN::min_change)
	std::cout << "min_change";
      else if (constraint == TRANSPLAN::min_achieve)
	std::cout << "min_achieve";
      else if (constraint == TRANSPLAN::max_change)
	std::cout << "max_change";
      else if (constraint == TRANSPLAN::max_achieve)
	std::cout << "max_achieve";
      else if (constraint == TRANSPLAN::change_before)
	std::cout << "changes_before";
      else if (constraint == TRANSPLAN::change_after)
	std::cout << "changes_after";
      else if (constraint == TRANSPLAN::achieve_before)
	std::cout << "achieves_before";
      else if (constraint == TRANSPLAN::achieve_after)
	std::cout << "achieves_after";
      else if (constraint == TRANSPLAN::before)
	std::cout << "before";
      else if (constraint == TRANSPLAN::after)
	std::cout << "after";
      else if (constraint == TRANSPLAN::start_before)
	std::cout << "starts_before";
      else if (constraint == TRANSPLAN::start_after)
	std::cout << "starts_after";
      else if (constraint == TRANSPLAN::end_before)
	std::cout << "ends_before";
      else if (constraint == TRANSPLAN::end_after)
	std::cout << "ends_after";

      if (!unary)
	std::cout <<":" << right << std::endl;
    }

  };

  class Matrix
  {
  public:
    int dimention;
    std::vector<std::vector<int> > matrix;

  Matrix(int n_state, int def_val = 0)
    :dimention(n_state)
    {
      for (int i = 0; i < dimention; ++i)
      {
	matrix.push_back(std::vector<int>());
	for (int j = 0; j < dimention; ++j)
	{
	  matrix[i].push_back(def_val);
	}
      }
    }

    //copy
    Matrix( const Matrix& m)
    :dimention(m.dimention)
    {
      matrix.clear();
      for (int i = 0; i < dimention; ++i)
      {
	matrix.push_back(std::vector<int>());
	for (int j = 0; j < dimention; ++j)
	{
	  matrix[i].push_back(m.matrix[i][j]);
	}
      }
    }

    //operator =
    Matrix& operator=(const Matrix& m)
    {
      assert(this->dimention == m.dimention);

      for (int i = 0; i < dimention; ++i)
      {
	for (int j = 0; j < dimention; ++j)
	{
	  this->matrix[i][j] = m.matrix[i][j];
	}
      }

      return *this;
    }

    //operator []
    std::vector<int>& operator[](int col)
    {
      //std::cout << "col = " << col << std::endl;
      assert(col < dimention);
      return matrix[col];
    }

    void assign(int row, int col, int val)
    {
      assert(col < dimention);
      assert(row < dimention);

      this->matrix[row][col] = val;
    }

  };

  /*!
    Store Search Statistics
  */
  struct SearchStatistics
  {
    //Some statistics
    int n_backtrack;
    int n_fails;
    int n_unit_propagation_call;
    int cpu_time, start_time, end_time;
    int n_nodes;
    int n_derived_prec;
    int make_span;

    double init_fldt_b;
    double init_fldt_a;
    double init_flex_b;
    double init_flex_a;
    int init_total_slack;
    double plan_fldt;
    double plan_flex;
    double plan_total_slack;
    

    SearchStatistics()
    {
      n_backtrack = 0;
      n_fails = 0;
      n_derived_prec = 0;
      n_unit_propagation_call = 0;
      n_nodes = 0;
      cpu_time = 0;
      make_span = 0;

      plan_fldt = 0;
      plan_flex = 0;
      plan_total_slack = 0;

    }

    void print()
    {
      std::cout << "Backtrack = " << n_backtrack << std::endl;
      std::cout << "Expanded Nodes = " << n_nodes << std::endl;
      std::cout << "Derived Precedence Rel = " << n_derived_prec << std::endl;
    }

    void print_in_one_line(std::string& pref)
    {
      
      std::cout << pref << " " << cpu_time << " " << n_backtrack << " " << n_nodes  << " " << make_span 
		<< " " << plan_flex  << " " << init_flex_b <<" " << init_flex_a 
		<< " " << plan_fldt<< " " << init_fldt_b << " " << init_fldt_a 
		<< " " << plan_total_slack << " " << init_total_slack
		<< std::endl;
    }
  };

  struct Triplet
  {
    int lhs;
    int rhs;
    bool type;

    Triplet()
    {
      lhs = -1;
      rhs = -1;
      type = false;
    }

    Triplet(const Triplet& t)
    :lhs(t.lhs)
    ,rhs(t.rhs)
    ,type(t.type)
    {
      //
    }

    Triplet& operator=(const Triplet& t)
    {
      this->lhs = t.lhs;
      this->rhs = t.rhs;
      this->type = t.type;

      return *this;
    }
  };
}
#endif
