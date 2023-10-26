/*
 * distance_graph.h
 *
 *  Created on: 16/09/2016
 *      Author: jing
 */

#ifndef DISTANCE_GRAPH_H_
#define DISTANCE_GRAPH_H_

#define TYPE_LOWER 1
#define TYPE_UPPER -1
#define TYPE_ORDINARY 0
#include "stnu.h"
#include "candidate.h"
#include "conflict.h"

namespace CDS
{
	class DistanceEdge
	{
	public:

		std::string start_name, end_name;
		int start_; //start node id
		int end_;  //end node id
		//int id_; // edge id
		double lb_;	//distance
		double ub_;	//upper bound of the distance (for ctg)
		int label_;	//low = -id-1, upper = id+1
		bool b_block_;	//1-can be blocked;

		std::vector<int> v_support_; // 0~i:ctg_l, -1~-i-1:ctg_u, b+2*i:r_u, b+2*i+1:r_l

		std::vector<std::pair<int, int> > v_label; //first-id, second-value
		//std::vector<std::string> v_support_names; // keep a record of reductions

		DistanceEdge()
		{
			start_ = 0;
			end_ = 0;
			lb_ = -g_inf;
			ub_ = g_inf;
			label_ = 0;
			b_block_ = false;
			v_support_.clear();
		}
		DistanceEdge(const LINK& a)
		{
			start_ = a.start;
			end_ = a.end;
			start_name = a.start_name;
			end_name = a.end_name;
			lb_ = a.lb;
			ub_ = a.ub;
			b_block_ = false;
			v_support_.clear();
			label_ = 0;
			v_label = a.v_labels;

		}
		DistanceEdge(const LINK& a, bool re) :
				DistanceEdge(a)
		{
			if (re) // reverse link
			{
				start_ = a.end;
				end_ = a.start;
			}

			if (a.b_control) // rqm
			{
				lb_ = re ? (-a.lb) : a.ub;
				label_ = 0;
			}
			else //ctg
			{
				label_ = a.end + 1;
				label_ = re ? (label_) : -label_;
				lb_ = re ? (-a.ub) : lb_;
			}
		}
		bool operator<(const DistanceEdge a) const
		{
			return lb_ > a.lb_;
			if (start_ != a.start_)
				return start_ < a.start_;
			if (end_ != a.end_)
				return end_ < a.end_;
			return label_ < a.label_;
		}

		bool operator!=(const DistanceEdge a) const
		{
			if (*this < a)
				return true;
			if (a < *this)
				return true;
			return false;
		}

		DistanceEdge operator+(const DistanceEdge a) const;
		void print() const;

	};

	class DistanceGraph //distance graph with lower/upper case labels
	{
		int num_nodes_;

	public:
		bool b_updated;
		vector<DistanceEdge> v_edge_;	//all edges
		vector<int> v_ordinary_id;	//the index of ordinary edges
		vector<int> v_lower_id;	// the index of lower case edges
		vector<int> v_upper_id;	// the index of upper case edges
		vector<vector<int> > v_end_id;// the index of edges whose end node is i
		vector<vector<int> > v_start_id;// the index of edges whose start node is i
		map<pair<int, int>, int> mp_index_edge_;	// the map to ordinary edges
		map<pair<int, int>, int> mp_start_upper;
		map<int, int> mp_start_lower;
		void clear()
		{
			v_edge_.clear();
			v_ordinary_id.clear();
			v_lower_id.clear();
			v_upper_id.clear();
			v_end_id.clear();
			v_start_id.clear();
			mp_index_edge_.clear();
			mp_start_upper.clear();
			mp_start_lower.clear();
		}

		int add_edge(const DistanceEdge &e);
		void addLowerCaseEdge(const DistanceEdge &e);
		void addUpperCaseEdge(const DistanceEdge &e);
		void addOrdinaryEdge(const DistanceEdge &e);

		DistanceEdge getEdge(const int &id);
		DistanceEdge getLower(const int &id);
		DistanceEdge getUpper(const int &id);
		DistanceEdge getOrdinary(const int &id);

		DistanceEdge getStart(const int &st, const int &id);
		DistanceEdge getEnd(const int &st, const int &id);

		int get_num_nodes() const
		{
			return num_nodes_;
		}
		;
		void set_num_nodes(int num)
		{
			num_nodes_ = num;
		}
		;

		//void TransGraph(stnu * p_stnu);//get i_precede;
		void PrintDot(std::string output_file);

		// the bound between ids of lower edges and ordinary edges in support vector
		int low_ord_edge_bound;

		int edge_to_id(DistanceEdge edge);

	};

}

#endif /* DISTANCE_GRAPH_H_ */
