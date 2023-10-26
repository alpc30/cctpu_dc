/*
 * controllability.cc
 *
 *  Created on: 07/09/2015
 *      Author: jing
 */

#include "controllability.h"

namespace CDS
{

Controllability::Controllability()
{
	// TODO Auto-generated constructor stub
}

stnu Controllability::InitDistanceGraph(const stnu & a, const Candidate &b)
{
	stnu_ = a;

	//mute the edges according to assignments
	for (int id = 0; id < b.vec_ass.size(); id++)
	{
		int i = b.vec_ass[id].first;
		int val = b.vec_ass[id].second;
		if (a.s_info.id_debug == DEBUG_T)
			printf("mute v%d=%d\n", i, val);
		int tmpr;
		if (val == -1)
			val = 0, tmpr = 2;
		else
			tmpr = 1;
		while (tmpr--)
		{
			for (int k = 0; k < stnu_.v_dis_var[i].vv_id_mute_edge[val].size();
					k++)
			{
				int id_edge = stnu_.v_dis_var[i].vv_id_mute_edge[val][k];
				if (id_edge < 0) // mute ctg
				{
					id_edge = -id_edge;
					id_edge--;
					stnu_.ctg[id_edge].b_mute = true;
					if (a.s_info.id_debug == DEBUG_T)
						stnu_.ctg[id_edge].print();
				}
				else // mute rqm
				{
					id_edge--;
					stnu_.rqm[id_edge].b_mute = true;
					if (a.s_info.id_debug == DEBUG_T)
						stnu_.rqm[id_edge].print();
				}
			}
			val++;
		}
	}
	if (a.s_info.id_debug == DEBUG_T)
		stnu_.print_dot("stnu.dot");
	//use the relaxed bounds of candidate b
	for (int i = 0; i < stnu_.n_ctg; i++)
	{
		int id = stnu_.ctg[i].id_lr;
		if (b.v_relaxed_bounds[id].first)
		{
			stnu_.ctg[i].lb = b.v_relaxed_bounds[id].second;
		}
		id = stnu_.ctg[i].id_ur;
		if (b.v_relaxed_bounds[id].first)
		{
			stnu_.ctg[i].ub = b.v_relaxed_bounds[id].second;
		}
	}

	for (int i = 0; i < stnu_.n_rqm; i++)
	{
		int id = stnu_.rqm[i].id_lr;
		if (b.v_relaxed_bounds[id].first)
		{
			stnu_.rqm[i].lb = b.v_relaxed_bounds[id].second;
		}
		id = stnu_.rqm[i].id_ur;
		if (b.v_relaxed_bounds[id].first)
		{
			stnu_.rqm[i].ub = b.v_relaxed_bounds[id].second;
		}
	}

	//buid distance graph
	InitDistanceGraph(stnu_);
	return stnu_;

}

void Controllability::InitDistanceGraph(const stnu &a)
{
	//Initialization of the distance graph
	stnu_ = a;
	distance_graph_.clear();
	distance_graph_.set_num_nodes(stnu_.n_nodes);
	distance_graph_.low_ord_edge_bound = stnu_.n_ctg + 1;
	for (int i = 0; i < distance_graph_.get_num_nodes(); i++)
	{
		std::vector<int> tmp_v_id;
		distance_graph_.v_end_id.push_back(tmp_v_id);
		distance_graph_.v_start_id.push_back(tmp_v_id);
	}

	for (int i = 0; i < stnu_.n_ctg; i++)
	{
		if (stnu_.ctg[i].b_mute) // the ctg is muted
			continue;

		// add A->B [b:l]
		DistanceEdge tmp_edge(stnu_.ctg[i], false);
		tmp_edge.v_support_.push_back(stnu_.edgeToId(stnu_.ctg[i], false));
		distance_graph_.addLowerCaseEdge(tmp_edge);

		// add A<-B [B:-u]
		tmp_edge = DistanceEdge(stnu_.ctg[i], true);
		tmp_edge.v_support_.push_back(stnu_.edgeToId(stnu_.ctg[i], true));
		distance_graph_.addUpperCaseEdge(tmp_edge);
	}

	for (int i = 0; i < stnu_.n_rqm; i++)
	{
		if (stnu_.rqm[i].b_mute)
			continue;
		//add A->B [u]
		DistanceEdge tmp_edge(stnu_.rqm[i], false);
		tmp_edge.v_support_.push_back(stnu_.edgeToId(stnu_.rqm[i], true));
		if (tmp_edge.lb_ > -g_inf && tmp_edge.lb_ < g_inf)
			distance_graph_.addOrdinaryEdge(tmp_edge);

		//add A<-B [-l]
		tmp_edge = DistanceEdge(stnu_.rqm[i], true);
		tmp_edge.v_support_.push_back(stnu_.edgeToId(stnu_.rqm[i], false));
		if (tmp_edge.lb_ > -g_inf && tmp_edge.lb_ < g_inf)
			distance_graph_.addOrdinaryEdge(tmp_edge);
	}

	//generate source nodes
	int cnt_nodes = 0;
	vector<int> v_nodes_;
	for(int i = 0; i < stnu_.n_nodes; i++)
	  v_nodes_.push_back(0);
	while(cnt_nodes < stnu_.n_nodes)
	{
	  int now = 0;
	  for(now = 0; now < stnu_.n_nodes; now++)
	  {
	    if(v_nodes_[now] == 0) // find the first unmarked node
	      break;
	  }

	  v_src_nodes.push_back(now); // add now as a source node
	  queue<int> q_nodes;
	  q_nodes.push(now);
	  while(q_nodes.size()){
	    int crt = q_nodes.front();
	    q_nodes.pop();
	    if(v_nodes_[crt]) continue;
	    v_nodes_[crt] = 1; cnt_nodes++;
	    for(int i = 0; i < distance_graph_.v_start_id[crt].size(); i++)
	    {
	      int next = distance_graph_.v_edge_[distance_graph_.v_start_id[crt][i]].end_;
	      if(v_nodes_[next])continue;
	      q_nodes.push(next);
	    }
	  }
	}
	

	//distance_graph_.PrintDot("dis2.dot");
}

void Controllability::edgeEdge(LINK a)
{
	DistanceEdge tmp_edge(a, false); //ub
	if (a.b_control) //rqm
	{
		//if(tmp_edg.ub)

		tmp_edge.v_support_.push_back(stnu_.edgeToId(a, true));

		if (distance_graph_.mp_index_edge_.find(
		{ tmp_edge.start_, tmp_edge.end_ })
				== distance_graph_.mp_index_edge_.end())
		{
			distance_graph_.addOrdinaryEdge(tmp_edge);
		}
		else if (tmp_edge.lb_ > -g_inf && tmp_edge.lb_ < g_inf)
		{
			int id_new_edge;
			id_new_edge = distance_graph_.mp_index_edge_[
			{ tmp_edge.start_, tmp_edge.end_ }] - 1;
			distance_graph_.v_edge_[id_new_edge] = tmp_edge;
		}

		tmp_edge = DistanceEdge(a, true); //lb
		tmp_edge.v_support_.push_back(stnu_.edgeToId(a, false));
		if (distance_graph_.mp_index_edge_.find(
		{ tmp_edge.start_, tmp_edge.end_ })
				== distance_graph_.mp_index_edge_.end())
		{
			distance_graph_.addOrdinaryEdge(tmp_edge);
		}
		else if (tmp_edge.lb_ > -g_inf && tmp_edge.lb_ < g_inf)
		{
			int id_new_edge;

			id_new_edge = distance_graph_.mp_index_edge_[
			{ tmp_edge.start_, tmp_edge.end_ }] - 1;
			distance_graph_.v_edge_[id_new_edge] = tmp_edge;
		}

	}

}

bool Controllability::CheckingSC()
{
  // store ctg_id with index of start and end nodes
  vector<int> v_end_lower, v_start_upper;
  for(int i = 0; i < distance_graph_.get_num_nodes(); i++)
  {
    v_end_lower.push_back(0);
    v_start_upper.push_back(0);
  }

  for(int i = 0; i < distance_graph_.v_lower_id.size(); i++)
  {
    int end = distance_graph_.getLower(i).end_;
    v_end_lower[end] = i+1;
  }

  for(int i = 0; i < distance_graph_.v_upper_id.size(); i++)
  {
    int start = distance_graph_.getUpper(i).start_;
    v_start_upper[start] = i + 1;
  }
	//checking strong controllability for the current stnu
	//stnu_.reader->PrintDot("2.dot");


  // std::cerr<<"Before CheckingSC, num_edge="<<distance_graph_.v_ordinary_id.size()<<std::endl;

  int id_edge = 0;

	while(id_edge < distance_graph_.v_ordinary_id.size())
	{
	  DistanceEdge tmp_rqm_edge = distance_graph_.getOrdinary(id_edge);
	  // tmp_rqm_edge.print();
	  //1. A --b:ctg_l--> B --rqm_(-l)--> C
	  //add A --(ctg_l-l)--> C
	  if(v_end_lower[tmp_rqm_edge.start_] !=0)
	  {
	    int id_lower = v_end_lower[(tmp_rqm_edge.start_)]-1;
	    DistanceEdge tmp_ctg_edge = distance_graph_.getLower(id_lower);


	    DistanceEdge tmp_add_edge;
	    tmp_add_edge = tmp_ctg_edge + tmp_rqm_edge;
//	    tmp_add_edge.start_ = tmp_ctg_edge.start_;
//	    tmp_add_edge.end_ = tmp_rqm_edge.end_;
//	    tmp_add_edge.lb_ = tmp_ctg_edge.lb_ + tmp_rqm_edge.lb_;
//	    tmp_add_edge.support = tmp_ctg_edge.support;
//	    tmp_
	    tmp_add_edge.label_ = 0;
	    if(isTighterOrdinaryEdge(tmp_add_edge))
	      addOrdinaryEdge(tmp_add_edge);

	   


	    // if(distance_graph_.mp_index_edge_[{tmp_ctg_edge.start_,tmp_rqm_edge.end_}]==0)
	    //{
	    //  distance_graph_.mp_index_edge_[{tmp_ctg_edge.start_,
//		    tmp_rqm_edge.end_}] = distance_graph_.v_ordinary_id.size()+1;
//	      distance_graph_.v_ordinary_id.push_back(tmp_add_edge);
//	    }
//	    else
//	    {
//	      int tmp_id_edge = distance_graph_.mp_index_edge_[{tmp_ctg_edge.start_,
//								tmp_rqm_edge.end_}]-1;

//	      if(distance_graph_.v_ordinary_id[tmp_id_edge].lb_ > tmp_add_edge.lb_)
//		distance_graph_.v_ordinary_id[tmp_id_edge].lb_ = tmp_add_edge.lb_;
//	    }
	  }

		//2. A <--B:(-ctg_u)-- B <--rqm_u-- C
		// add C --(-ctg_u+rqm_u)--> A

	  if(v_start_upper[tmp_rqm_edge.end_] != 0)
	  {
	    int id_ctg = v_start_upper[tmp_rqm_edge.end_]-1;
	    DistanceEdge tmp_ctg_edge = distance_graph_.getUpper(id_ctg);
	    DistanceEdge tmp_add_edge;
	    tmp_add_edge =  tmp_rqm_edge + tmp_ctg_edge;
	    tmp_add_edge.label_ = 0;
	    if(isTighterOrdinaryEdge(tmp_add_edge))
	      addOrdinaryEdge(tmp_add_edge);
//	    tmp_add_edge.start_ = tmp_rqm_edge.start_;
//	    tmp_add_edge.end_ = tmp_ctg_edge.end_;
//	    tmp_add_edge.lb_ = tmp_ctg_edge.lb_ + tmp_rqm_edge.lb_;

			// if(distance_graph_.mp_index_edge_[{tmp_ctg_edge.start_,tmp_rqm_edge.end_}]==0)
			// {
			// 	distance_graph_.mp_index_edge_[{tmp_ctg_edge.start_,
			// 		tmp_rqm_edge.end_}] = distance_graph_.v_ordinary_id.size()+1;
			// 	distance_graph_.v_ordinary_id.push_back(tmp_add_edge);
			// }
			// else
			// {
			// 	int tmp_id_edge = distance_graph_.mp_index_edge_[{tmp_ctg_edge.start_,
			// 				tmp_rqm_edge.end_}]-1;
			// 	if(distance_graph_.v_ordinary_id[tmp_id_edge].lb_ > tmp_add_edge.lb_)
			// 		distance_graph_.v_ordinary_id[tmp_id_edge].lb_ = tmp_add_edge.lb_;
			// }

	  }
	  ++id_edge;
	}

//	std::vector<Conflict> v_neg_cycle;
//	v_neg_cycle =
	v_conf = BellmanFord();
	
	//distance_graph_.PrintDot("dis2.dot");

	if (v_conf.size())
	{
    
//	  ng_path = v_neg_cycle[0];
		return false;
	}
	else
       
		return true;
}

bool Controllability::CheckingDC()
{
	//checking dynamic controllability for the current stnu
	// O(n^4)

	std::vector<Conflict> v_neg_cycle;
	v_neg_cycle.clear();

	int K = distance_graph_.get_num_nodes();
	for (int i = 0; i < K; i++)
	{
		std::vector<Conflict> tmp_v_neg_cycle;
		tmp_v_neg_cycle = BellmanFord();
		v_neg_cycle.insert(v_neg_cycle.begin(), tmp_v_neg_cycle.begin(),
				tmp_v_neg_cycle.end());

		if (v_neg_cycle.size())
		{
			return false;
		}

		for (int j = 0; j < distance_graph_.v_lower_id.size(); ++j)
		{
			v_neg_cycle = Propagation(distance_graph_.getLower(j));
		}
	}
	return true;
}

int Controllability::isTighterOrdinaryEdge(DistanceEdge new_edge)
{
	// return: 	true - the new edge changed the distance graph
	//			false - the new edge did not change the graph

	if (containNegLoop(new_edge))
		return 0;
	if (new_edge.start_ == new_edge.end_)
		return 0;
	int id_new_edge;

	id_new_edge = distance_graph_.mp_index_edge_[
	{ new_edge.start_, new_edge.end_ }];

	if (id_new_edge != 0)
	{
		id_new_edge--;
		DistanceEdge tmp_new_edge;
		tmp_new_edge = distance_graph_.v_edge_[id_new_edge];
		if (tmp_new_edge.lb_ > new_edge.lb_ + 1e-3)
		{
			return 2;
		}
	}
	else
	{
		return 1;
	}
	return 0;
}

int Controllability::isTighterUpperEdge(DistanceEdge new_edge)
{
	// return: 	1 - the new edge changed the distance graph (add edge)
	//			2 - the new edge changed the distance graph (change lb)
	//			0- the new edge did not change the graph

	if (containNegLoop(new_edge))
		return 0;

	bool fl_update = false;
	int id_update = -1;
	if (distance_graph_.mp_start_upper.find(
			make_pair(new_edge.start_, new_edge.end_))
			!= distance_graph_.mp_start_upper.end())
	{
		int id_update = distance_graph_.mp_start_upper[make_pair(
				new_edge.start_, new_edge.end_)] - 1;
		if (distance_graph_.v_edge_[id_update].lb_ > new_edge.lb_)
		{
			//distance_graph_.v_edge_[id_update].lb_ = new_edge.lb_;
			return 2;
		}
		return 0;
	}
	//distance_graph_.addUpperCaseEdge(new_edge);
	return 1;

//	for(int k = 0; k < distance_graph_.v_upper_id.size() && !fl_update; k++)
//	{
//		DistanceEdge tmp_upper_edge;
//		tmp_upper_edge = distance_graph_.getUpper(k);
//		if(tmp_upper_edge.start_ == new_edge.start_ &&
//				tmp_upper_edge.end_ == new_edge.end_)
//		{
//			id_update = k;
//			if(tmp_upper_edge.lb_ <= new_edge.lb_ + 1e-6)
//			{
//				return true;
//				fl_update = true;
//			}
//		}
//	}
	if (fl_update == false)
	{
		return true;
	}
	else if (id_update >= 0
			&& distance_graph_.getUpper(id_update).lb_ > new_edge.lb_ + 1e-4)
	{
		return true;
	}
	return false;
}

bool Controllability::addOrdinaryEdge(DistanceEdge new_edge)
{
	// add a new ordinary edge to the distance graph

	int id_new_edge;
//	new_edge.print();

	id_new_edge = distance_graph_.mp_index_edge_[
	{ new_edge.start_, new_edge.end_ }];

	if (id_new_edge != 0)
	{
		id_new_edge--;
		DistanceEdge tmp_new_edge;
		tmp_new_edge = distance_graph_.v_edge_[id_new_edge];
		if (tmp_new_edge.lb_ > new_edge.lb_ + 1e-3)
		{
			distance_graph_.b_updated = true;
			distance_graph_.v_edge_[id_new_edge] = new_edge;
			return true;
		}
	}
	else
	{
		distance_graph_.b_updated = true;
		distance_graph_.addOrdinaryEdge(new_edge);
		return true;
	}
	return false;
}

bool Controllability::addUpperEdge(DistanceEdge new_edge)
{

	bool fl_update = false;
	int id_update = -1;
	if (distance_graph_.mp_start_upper.find(
			make_pair(new_edge.start_, new_edge.end_))
			!= distance_graph_.mp_start_upper.end())
	{
		int id_update = distance_graph_.mp_start_upper[make_pair(
				new_edge.start_, new_edge.end_)] - 1;
		if (distance_graph_.v_edge_[id_update].lb_ > new_edge.lb_)
		{
			distance_graph_.v_edge_[id_update].lb_ = new_edge.lb_;
			distance_graph_.b_updated = true;
			return true;
		}
		return false;
	}
	distance_graph_.b_updated = true;
	distance_graph_.addUpperCaseEdge(new_edge);
	return true;

	if (fl_update == false)
	{
		distance_graph_.b_updated = true;
		//new_edge.print();
		distance_graph_.addUpperCaseEdge(new_edge);
		return true;
	}
	else if (id_update >= 0
			&& distance_graph_.v_edge_[id_update].lb_ > new_edge.lb_ + 1e-4)
	{
		distance_graph_.b_updated = true;
		///	distance_graph_.v_upper_edge_[id_update].print();
		///	new_edge.print();
		distance_graph_.v_edge_[id_update].lb_ = new_edge.lb_;
		return true;
	}

	return false;
}

bool Controllability::isLabelConsistent(const DistanceEdge &e1,
		const DistanceEdge &e2)
{
	for (int i = 0; i < e1.v_label.size(); i++)
	{
		for (int j = 0; j < e2.v_label.size(); j++)
		{
			if (e1.v_label[i].first == e2.v_label[j].first
					&& e1.v_label[i].second != e2.v_label[j].second)
				return false;
		}
	}

	return true;
}

bool Controllability::containNegLoop(DistanceEdge edge)
{
	/*If the support vector contains negative loop return true*/

	return false;
	std::set<int> tmp_s;
	tmp_s.clear();
	for (int i = 0; i < edge.v_support_.size(); i++)
	{
		if (edge.v_support_[i] >= distance_graph_.low_ord_edge_bound)
			continue;
		if (tmp_s.find(edge.v_support_[i]) != tmp_s.end())
			return true;
		tmp_s.insert(edge.v_support_[i]);
	}
	tmp_s.clear();
	return false;
}

std::vector<Conflict> Controllability::Propagation(DistanceEdge low_edge)
{
	return PropagationDijkstra(low_edge);
}

std::vector<Conflict> Controllability::PropagationDijkstra(
		DistanceEdge low_edge)
{
	//Propagation through a lower-case edge
	/* 0. Initialize
	 * 1. Find potential starts of moat edge (connect to low_edge and not from the same ctg)
	 * 2. push potential start edges into vector
	 * 3. propagate the edges in the vector (as BFS) to find moat edges
	 * 4. do lower-case and cross-case reductions with low_edge and moat edges
	 * */

	//0. Initialization
	//std::vector<DistanceEdge> v_edge; // vector of edges to propagate
	std::priority_queue<DistanceEdge> q_edge;
	std::vector<DistanceEdge> v_moat_edge; // vector of moat edges
	std::vector<DistanceEdge> v_dis;
	std::vector<Conflict> res; //return result
	std::vector<std::set<int> > precede_set;

	vector<bool> b_visit;
	//bool *b_visit = new bool[distance_graph_.get_num_nodes()+1];
	for (int i = 0; i < distance_graph_.get_num_nodes() + 1; i++)
	{
		b_visit.push_back(false);
		//b_visit[i] = false;
		DistanceEdge tmp_edge;
		tmp_edge.start_ = low_edge.end_;
		tmp_edge.end_ = i + 1;
		tmp_edge.lb_ = g_inf;
		v_dis.push_back(tmp_edge);
		std::set<int> tmp_set;
		tmp_set.clear();
		precede_set.push_back(tmp_set);

	}

	//1. Find potential starts of moat edge
	for (int i = 0; i < distance_graph_.v_start_id[low_edge.end_].size(); i++)
	{
		DistanceEdge tmp_edge;
		tmp_edge = distance_graph_.getStart(low_edge.end_, i);
		if (tmp_edge.label_ + low_edge.label_ != 0 && tmp_edge.label_ >= 0) //not from the same ctg
		{
			q_edge.push(tmp_edge);
			v_dis[tmp_edge.end_] = tmp_edge;
		}
	}

	//3. Propagation
	while (q_edge.size())
	{
		DistanceEdge tmp_edge;

		tmp_edge = q_edge.top();
		q_edge.pop();
		printf("Expand: %d ", q_edge.size());
		tmp_edge.print();

		if (tmp_edge.start_ == tmp_edge.end_)
			continue;

		int id_node;
		id_node = tmp_edge.end_;

		//if tmp_edge < 0, it is a moat edge;
		if (v_dis[tmp_edge.end_].lb_ <= tmp_edge.lb_ + 1e-5 && b_visit[id_node])
			continue;
		if (tmp_edge.label_ == 0)
		{
			b_visit[id_node] = true;
			v_dis[tmp_edge.end_].lb_ = tmp_edge.lb_;
		}
		if (tmp_edge.lb_ < -1e-5 && tmp_edge.label_ + low_edge.label_ != 0)
		{
			v_moat_edge.push_back(tmp_edge);
			continue;
		}
		if (tmp_edge.label_ > 0) //upper case
		{
			//label removal
			for (int i = 0; i < distance_graph_.v_lower_id.size(); i++)
			{
				DistanceEdge next_edge;
				next_edge = distance_graph_.getLower(i);

				if (tmp_edge.end_ == next_edge.start_
						&& tmp_edge.label_ + next_edge.label_ == 0
						&& tmp_edge.lb_ >= -next_edge.lb_)
				{
					tmp_edge.label_ = 0;
					if (isTighterOrdinaryEdge(tmp_edge))
					{
						q_edge.push(tmp_edge);
					}
				}
			}
		}
		else // ordinary
		{
			//upper-case
			for (int i = 0;
					i < distance_graph_.v_start_id[tmp_edge.end_].size(); i++)
			{
				DistanceEdge next_edge;
				next_edge = distance_graph_.getStart(tmp_edge.end_, i);
				if (next_edge.label_ > 0)
				{
					DistanceEdge new_edge;
					new_edge = tmp_edge + next_edge;

					if (isTighterUpperEdge(new_edge))
					{
						//to avoid negative loop
						if (precede_set[tmp_edge.end_].find(next_edge.end_)
								!= precede_set[tmp_edge.end_].end())
							continue;
						precede_set[new_edge.end_] = precede_set[tmp_edge.end_];
						precede_set[new_edge.end_].insert(tmp_edge.end_);

						q_edge.push(new_edge);
					}
				}
			}

			//no-case
			for (int i = 0;
					i < distance_graph_.v_start_id[tmp_edge.end_].size(); i++)
			{

				DistanceEdge next_edge;
				next_edge = distance_graph_.getStart(tmp_edge.end_, i);
				if (next_edge.label_ == 0)
				{
					DistanceEdge new_edge;
					new_edge = tmp_edge + next_edge;

					if (isTighterOrdinaryEdge(new_edge))
					{
						//to avoid negative loop
						if (precede_set[tmp_edge.end_].find(next_edge.end_)
								!= precede_set[tmp_edge.end_].end())
							continue;
						precede_set[new_edge.end_] = precede_set[tmp_edge.end_];
						precede_set[new_edge.end_].insert(tmp_edge.end_);

						q_edge.push(new_edge);
					}
				}
			}
		}
	}
//	printf("finish\n");

	for (int i = 0; i < v_moat_edge.size(); i++)
	{
		DistanceEdge tmp_edge;
		tmp_edge = v_moat_edge[i];
		tmp_edge.b_block_ = true;

		DistanceEdge new_edge;
		new_edge = low_edge + tmp_edge;

		bool fl_update = false;

		if (tmp_edge.label_ == 0)
		{
			fl_update = addOrdinaryEdge(new_edge);
		}
		else if (tmp_edge.label_ + low_edge.label_ != 0)
		{
			fl_update = addUpperEdge(new_edge);
		}

		if (fl_update)
		{
//			printf("Lower-case Edge: ");
//			low_edge.print();
//			printf("Moat Edge: ");
//			tmp_edge.print();
			printf("New Edge(%d): ", distance_graph_.b_updated);
			new_edge.print();
			if (tmp_edge.v_support_.size())
				res.push_back(getBlockedRelaxation(tmp_edge));
		}
	}
	return res;
}

std::vector<Conflict> Controllability::BellmanFord()
{
//#ifdef DEBUG
//	stnu_.reader->PrintDot("1.dot");
//	std::cerr<<"Before Bellmanford, num_edge="<<distance_graph_.v_ordinary_edge_.size()<<std::endl;
//	for(int i = 0; i < distance_graph_.v_ordinary_edge_.size(); ++i)
//	{
//		DistanceEdge tmp_edge = distance_graph_.v_ordinary_edge_[i];
//		if( tmp_edge.lb_ != g_inf)
//		std::cerr<<tmp_edge.start_+1<<"--("<<tmp_edge.lb_<<")-->"<<tmp_edge.end_+1<<std::endl;
//	}
//#endif
	std::vector<double> tmp_dis;
	std::vector<std::vector<DistanceEdge> > v_v_path;
	std::vector<int> predecessor;
	std::vector<Conflict> res;
//	std::set<Conflict> set_conflicts;
	res.clear();
	v_v_path.clear();
	tmp_dis.clear();
	predecessor.clear();

	for (int i = 0; i < distance_graph_.get_num_nodes(); ++i)
	{

		tmp_dis.push_back(g_inf);
		predecessor.push_back(-1);

		std::vector<DistanceEdge> tmp_v_e;
		tmp_v_e.clear();
		v_v_path.push_back(tmp_v_e);
	}

	//set the source node
	for(int i = 0; i < v_src_nodes.size(); i++)
	  tmp_dis[v_src_nodes[i]] = 0;
	//the graph may be disconnected 

	for (int i = 1; i < distance_graph_.get_num_nodes(); ++i)
	{
		for (int j = 0; j < distance_graph_.v_ordinary_id.size(); ++j)
		{
			DistanceEdge tmp_edge = distance_graph_.getOrdinary(j);
			if (tmp_dis[tmp_edge.start_] < g_inf && tmp_edge.lb_ < g_inf)
			{
				if (tmp_dis[tmp_edge.end_]
						> tmp_dis[tmp_edge.start_] + tmp_edge.lb_ + 1e-3)
				{

					/*printf("tighten %d(%lf) to %d(%lf)+(%lf) = %lf\n",
					 tmp_edge.end_,tmp_dis[tmp_edge.end_],
					 tmp_edge.start_,tmp_dis[tmp_edge.start_],
					 tmp_edge.lb_,tmp_edge.lb_+tmp_dis[tmp_edge.start_]);*/

					tmp_dis[tmp_edge.end_] = tmp_dis[tmp_edge.start_]
							+ tmp_edge.lb_;
					v_v_path[tmp_edge.end_] = v_v_path[tmp_edge.start_];
					v_v_path[tmp_edge.end_].push_back(tmp_edge);
					predecessor[tmp_edge.end_] = tmp_edge.start_;
				}
			}
		}

		//upper edge
		for (int j = 0; j < distance_graph_.v_upper_id.size(); ++j)
		{
			DistanceEdge tmp_edge = distance_graph_.getUpper(j);
			if (tmp_dis[tmp_edge.start_] < g_inf && tmp_edge.lb_ < g_inf)
			{
				if (tmp_dis[tmp_edge.end_]
						> tmp_dis[tmp_edge.start_] + tmp_edge.lb_ + 1e-3)
				{
					tmp_dis[tmp_edge.end_] = tmp_dis[tmp_edge.start_]
							+ tmp_edge.lb_;
					predecessor[tmp_edge.end_] = tmp_edge.start_;
					v_v_path[tmp_edge.end_] = v_v_path[tmp_edge.start_];
					v_v_path[tmp_edge.end_].push_back(tmp_edge);
				}
			}
		}
	}

	for (int j = 0; j < distance_graph_.v_ordinary_id.size(); ++j)
	{
		DistanceEdge tmp_edge = distance_graph_.getOrdinary(j);
		
		if (tmp_dis[tmp_edge.start_] < g_inf)
		{
			if (tmp_dis[tmp_edge.end_]
					> tmp_dis[tmp_edge.start_] + tmp_edge.lb_ + 1e-3)
			{
			  
				Conflict tmp_conflict = getNegCycle(v_v_path[tmp_edge.start_]);

				//	if (set_conflicts.find(tmp_conflict) == set_conflicts.end())
				{
					if (stnu_.s_info.id_debug == DEBUG_T)
						printf("conf = %lf\n",
								tmp_dis[tmp_edge.end_]
										- tmp_dis[tmp_edge.start_]);

					res.push_back(tmp_conflict);
					//	set_conflicts.insert(tmp_conflict);
					return res;
				}
			}
		}
	}

	for (int j = 0; j < distance_graph_.v_upper_id.size(); ++j)
	{
		DistanceEdge tmp_edge = distance_graph_.getUpper(j);
		if (tmp_dis[tmp_edge.start_] < g_inf)
		{
			if (tmp_dis[tmp_edge.end_]
					> tmp_dis[tmp_edge.start_] + tmp_edge.lb_ + 1e-3)
			{
				Conflict tmp_conflict = getNegCycle(v_v_path[tmp_edge.start_]);
//				tmp_conflict.print();
				///	v_v_path[tmp_edge.start_].print();
				//	if (set_conflicts.find(tmp_conflict) == set_conflicts.end())
				{
					res.push_back(tmp_conflict);
					if (stnu_.s_info.id_debug == DEBUG_T)
						printf("conf = %lf\n",
								tmp_dis[tmp_edge.end_]
										- tmp_dis[tmp_edge.start_]);
					//	set_conflicts.insert(tmp_conflict);
					return res;
				}
			}
		}
	}

	return res;
}

std::vector<Conflict> Controllability::DCConflict(stnu &n)
{
	std::vector<Conflict> conflicts;
	conflicts.clear();

	//init distance graph
	InitDistanceGraph(n);

	//dc checking return conflict
	conflicts = this->CheckingDCNegCycle();
	return conflicts;
}

std::vector<Conflict> Controllability::DCConflict1(const stnu &n,
		const Candidate &b)
{
	std::vector<Conflict> conflicts;
	conflicts.clear();

	InitDistanceGraph(n, b);
	if (n.s_info.id_debug == DEBUG_T)
		printf("Initialization Completed!\n");
	if(n.s_info.id_constrs == C_DC)
	{
	  CubicCheckingDC();

        //reverse the ng_path 
	vector<DistanceEdge> rev_ng_path;
	while(ng_path.size())
	{
          rev_ng_path.push_back(ng_path.back());
          ng_path.pop_back();
	}

	if (rev_ng_path.size())
	{
          conflicts.push_back(getNegCycle(rev_ng_path));
			//conflicts.push_back(getNegCycle(ng_path));
			//conflicts[0].print();
	} //	else this->extractEnvelope(conflicts);
	}
	else
	{
	  CheckingSC();
	  conflicts = v_conf;
	}
	

	return conflicts;
}

std::vector<Conflict> Controllability::DCConflict(const stnu &n,
		const Candidate &b)
{

	std::vector<Conflict> conflicts;
	conflicts.clear();

	InitDistanceGraph(n, b);
	if (n.s_info.id_debug == DEBUG_T)
		printf("Initialization Completed!\n");

	//use new structure
	extractEnvelope(conflicts);
//	CubicCheckingDC();
//	vector<DistanceEdge> rev_ng_path;
//	while(ng_path.size())
//	{
//		rev_ng_path.push_back(ng_path.back());
//		ng_path.pop_back();
//	}
//	conflicts.push_back(getNegCycle(rev_ng_path));
//
	//remove repeated conflicts
	for (int i = conflicts.size() - 1; i >= 0; i--)
	{
		for (int j = i - 1; j >= 0; j--)
		{
			if (conflicts[i].isSame(conflicts[j]))
			{
				conflicts.erase(conflicts.begin() + i);
				break;
			}
		}
	}

	return conflicts;

}

std::vector<Conflict> Controllability::CheckingDCNegCycle()
{
	//checking DC and return negative cycles

	std::vector<Conflict> v_neg_cycle;
	std::vector<Conflict> v_neg_cycle_prop;
	v_neg_cycle.clear();
	v_neg_cycle_prop.clear();
	std::vector<Conflict> tmp_v_neg_cycle;
	int K = distance_graph_.get_num_nodes();

	for (int i = 0; i < K; i++)
	{
//		printf("K = %d, edges = %d %d\n",i, this->distance_graph_.v_ordinary_edge_.size()
//				, this->distance_graph_.v_upper_edge_.size());
////
//		printf("Bellman:\n");
		tmp_v_neg_cycle = BellmanFord();

		if (tmp_v_neg_cycle.size())
		{
			v_neg_cycle = tmp_v_neg_cycle;
			//	v_neg_cycle.insert(v_neg_cycle.begin(),
			//			tmp_v_neg_cycle.begin(), tmp_v_neg_cycle.end());
			//	v_neg_cycle.insert(v_neg_cycle.begin(),
			//		v_neg_cycle_prop.begin(), v_neg_cycle_prop.end());
			break;
		}
		v_neg_cycle.clear();
		//	v_neg_cycle_prop.clear();
		distance_graph_.b_updated = false;
		for (int j = 0; j < distance_graph_.v_lower_id.size(); ++j)

		{
			//	printf("Propagation:\n");
			tmp_v_neg_cycle = Propagation(distance_graph_.getLower(j));
			//	if(tmp_v_neg_cycle.size())
			//	v_neg_cycle_prop.insert(v_neg_cycle_prop.begin(),
			//				tmp_v_neg_cycle.begin(), tmp_v_neg_cycle.end());
		}

		if (distance_graph_.b_updated == false)
			break;
	}

	//v_neg_cycle = BellmanFord();
	return v_neg_cycle;
}
std::vector<Conflict> Controllability::CheckingDCNegCycle1()
{
	//checking DC and return negative cycle

	std::vector<Conflict> v_neg_cycle;
	std::vector<Conflict> v_neg_cycle_prop;
	v_neg_cycle.clear();
	v_neg_cycle_prop.clear();
	std::vector<Conflict> tmp_v_neg_cycle;
	int K = distance_graph_.get_num_nodes();

	for (int i = 0; i < K; i++)
	{
		//	printf("K = %d, edges = %d %d\n",i, this->distance_graph_.v_ordinary_edge_.size()
		//		, this->distance_graph_.v_upper_edge_.size());
////
		printf("Bellman: %d\n", i);
		tmp_v_neg_cycle = BellmanFord();

		if (tmp_v_neg_cycle.size())
		{
			//v_neg_cycle = tmp_v_neg_cycle;
			v_neg_cycle.insert(v_neg_cycle.begin(), tmp_v_neg_cycle.begin(),
					tmp_v_neg_cycle.end());
			v_neg_cycle.insert(v_neg_cycle.begin(), v_neg_cycle_prop.begin(),
					v_neg_cycle_prop.end());
			break;
		}
		v_neg_cycle.clear();
		v_neg_cycle_prop.clear();
		distance_graph_.b_updated = false;
		for (int j = 0; j < distance_graph_.v_lower_id.size(); ++j)
		{
			printf("Propagation: %d\n", j);
			tmp_v_neg_cycle = Propagation(distance_graph_.getLower(j));
			if (tmp_v_neg_cycle.size())
				v_neg_cycle_prop.insert(v_neg_cycle_prop.begin(),
						tmp_v_neg_cycle.begin(), tmp_v_neg_cycle.end());
		}

		if (distance_graph_.b_updated == false)
			break;
	}
	return v_neg_cycle;
}

Conflict Controllability::getNegCycle(std::vector<DistanceEdge> v_neg_edge_o)
{
	Conflict tmp_conflict;
//	tmp_conflict.lb = 0;
	tmp_conflict.v_conflict_edge_id.clear();
	vector<DistanceEdge> v_neg_edge;
	vector<bool> visit;
	vector<double> sum;

	//bool* visit = new bool[this->stnu_.n_nodes];
	if(v_neg_edge_o.size() > 1) // check the direction
	{
	  bool fl_r = false;
	  for(int i = 1; i < v_neg_edge_o.size(); i++)
	  {
	    if(v_neg_edge_o[i-1].start_ != v_neg_edge_o[i].end_)
	    {
	      fl_r = true;
	      break;
	    }
	  }
	  
	  if(fl_r)
	  {
	    int len = v_neg_edge_o.size();
	    for(int i = 0; i < len/2; i++)
	    {
	      DistanceEdge tmp_edge = v_neg_edge_o[i];
	      v_neg_edge_o[i] = v_neg_edge_o[len - i - 1];
	      v_neg_edge_o[len-i-1] = tmp_edge;
	    }
	  }
	    
	}

	for (int i = 0; i < stnu_.n_nodes; i++)
	{
		visit.push_back(false);
		sum.push_back(0);
	}
	for (int i = 0; i < v_neg_edge_o.size(); i++)
	{
	  //	printf("fa: ");
	  //	v_neg_edge_o[i].print();
		for (int j = 0; j < v_neg_edge_o[i].v_support_.size(); j++)
		{
			pair<LINK, bool> tmp;
			tmp = stnu_.idToEdge(v_neg_edge_o[i].v_support_[j]);
			//if()
			if (tmp.first.b_control)		//rqm
				tmp.second = !tmp.second;
			DistanceEdge x(tmp.first, tmp.second);
			x.v_support_.push_back(v_neg_edge_o[i].v_support_[j]);
			//	x.print();
			v_neg_edge.push_back(x);
		}
	}
	for (int i = 0; i < v_neg_edge.size(); i++)
	{

		DistanceEdge tmp_edge;
		tmp_edge = v_neg_edge[i];
		//	tmp_edge.print();
		//continue;
		if (tmp_edge.lb_ < 0)
			visit[tmp_edge.end_] = 1;

		if (visit[tmp_edge.start_] == true
				&& sum[tmp_edge.end_] + tmp_edge.lb_ - sum[tmp_edge.start_] < 0)
		{

			int id_node_cycle;
			id_node_cycle = tmp_edge.start_;

			int t_i = 0;
			while (tmp_edge.end_ != id_node_cycle && t_i < i)
			{
				t_i++;
				if (tmp_edge.b_block_)
					tmp_conflict.v_alter_blocked_moats.push_back(
							tmp_edge.v_support_);
				tmp_conflict.v_conflict_edge_id.insert(
						tmp_conflict.v_conflict_edge_id.begin(),
						tmp_edge.v_support_.begin(), tmp_edge.v_support_.end());
				//tmp_conflict.v_conflict_edge.push_back(tmp_edge);
				//tmp_conflict.lb += tmp_edge.lb_;
				tmp_edge = v_neg_edge[i - t_i];
			}
			//tmp_conflict.v_conflict_edge.push_back(tmp_edge);
			tmp_conflict.v_conflict_edge_id.insert(
					tmp_conflict.v_conflict_edge_id.begin(),
					tmp_edge.v_support_.begin(), tmp_edge.v_support_.end());
			if (tmp_edge.b_block_)
				tmp_conflict.v_alter_blocked_moats.push_back(
						tmp_edge.v_support_);

			//	tmp_conflict.lb += tmp_edge.lb_;

			break;

		}
		sum[tmp_edge.start_] = sum[tmp_edge.end_] + tmp_edge.lb_;
		//visit[tmp_edge.end_] = 1;
	}
//	tmp_conflict.print();

	set<pair<int, int> > uniq;
	for (int i = 0; i < tmp_conflict.v_conflict_edge_id.size(); i++)
	{
		int id = tmp_conflict.v_conflict_edge_id[i];
		LINK tmp_edge = this->stnu_.idToEdge(id).first;
		if (tmp_edge.v_labels.size())
		{
			for (int j = 0; j < tmp_edge.v_labels.size(); j++)
			{
				if (uniq.find(tmp_edge.v_labels[j]) == uniq.end())
				{
					uniq.insert(tmp_edge.v_labels[j]);
					tmp_conflict.v_assignments.push_back(tmp_edge.v_labels[j]);
				}
			}
		}
	}
	return tmp_conflict;
}

Conflict Controllability::getBlockedRelaxation(DistanceEdge tmp_edge)
{
	Conflict tmp_conflict;

//	res.v_conflict_edge.clear();

	// replace the distance_graph_edge into edge_id of original networks 20160531
	//res.v_conflict_edge.push_back(tmp_edge);
	tmp_conflict.v_conflict_edge_id = tmp_edge.v_support_;
	// end change in 20160531
	set<pair<int, int> > uniq;
	for (int i = 0; i < tmp_conflict.v_conflict_edge_id.size(); i++)
	{
		int id = tmp_conflict.v_conflict_edge_id[i];
		LINK tmp_edge = this->stnu_.idToEdge(id).first;
		if (tmp_edge.v_labels.size())
		{
			for (int j = 0; j < tmp_edge.v_labels.size(); j++)
			{
				if (uniq.find(tmp_edge.v_labels[j]) == uniq.end())
				{
					uniq.insert(tmp_edge.v_labels[j]);
					tmp_conflict.v_assignments.push_back(tmp_edge.v_labels[j]);
				}
			}
		}
	}

	//res.print();
	//res.v_conflict_edge[0].print();
	//res.lb = tmp_edge.lb_;
	//res.ub = tmp_edge.ub_;
	return tmp_conflict;
}

bool Controllability::CubicCheckingDC()
{
	backprop_depth.clear();
	for (int i = 0; i < distance_graph_.v_edge_.size(); i++)
		backprop_depth.push_back(0); //

	for (int i = 0; i < distance_graph_.get_num_nodes(); i++)
	{
		v_q_neg_edge.push_back(vector<depth_edge>());
	}

	for (int i = 0; i < distance_graph_.v_upper_id.size(); i++)
	{
		// enumerate upper-case edge (they are negative)
		if (DCbackpropCheck(distance_graph_.getUpper(i)) == false)
		{
			return false;
		}
	}

	for (int i = 0; i < distance_graph_.v_ordinary_id.size(); i++)
	{
		if (distance_graph_.getOrdinary(i).lb_ < -1e-5)
		{
			if (DCbackpropCheck(distance_graph_.getOrdinary(i)) == false)
			{

				return false;
			}
		}
	}
	if (stnu_.s_info.id_debug)
		printf("cubic: feasible\n");

	return true;
}

bool Controllability::DCbackpropCheck(DistanceEdge e)
{
	int back_id = get_id_backprop(e);
	//printf("%d ",backprop_depth.size());
	if (stnu_.s_info.id_debug >= DEBUG_C)
	{
		printf("back prop: %d->%d [%d]%lf d[%d]=%d\n", e.start_, e.end_,
				e.label_, e.lb_, back_id, backprop_depth[back_id]);
	}

	if (backprop_depth[back_id] < 0) // ancestor call
	{
		if (stnu_.s_info.id_debug >= DEBUG_C)
			printf("anc: end back prop %d %d d[%d]=%d\n", e.start_, e.end_,
					back_id, backprop_depth[back_id]);

		return false;
	}

	if (backprop_depth[back_id] == 1) //terminated call
	{
		if (stnu_.s_info.id_debug >= DEBUG_C)
			printf("end back prop %d %d d[%d]=%d\n", e.start_, e.end_, back_id,
					backprop_depth[back_id]);
		return true;
	}

	backprop_depth[back_id] = -1;
	std::priority_queue<DistanceEdge> q_dis;
	std::vector<DistanceEdge> v_dis;
	std::vector<bool> v;
	for (int i = 0; i < distance_graph_.get_num_nodes(); i++)
	{
		DistanceEdge tmp_edge;
		tmp_edge.start_ = -1;
		tmp_edge.lb_ = 1e10; //+inf
		v_dis.push_back(tmp_edge);
		v.push_back(false);
	}

	v_dis[e.start_] = e;
	//v[e.end_] = true;
	q_dis.push(e);

	while (q_dis.size())
	{
		DistanceEdge now = q_dis.top();
		q_dis.pop();
		int id = now.start_;
		// now.print();
		if (v[id])
			continue;
		v[id] = true;

		if (now.lb_ >= -1e-5) //add new edge
		{
			if (now.label_)
				now.label_ = 0;
			if (this->isTighterOrdinaryEdge(now))
			{
				addOrdinaryEdge(now);

				if (stnu_.s_info.id_debug >= DEBUG_C)
				{
					//tmp_edge.print();
					//	dis[id].print();
					printf("new edge: ");
					now.print();
				}
			}

			continue;
		}
		for (int i = 0; i < distance_graph_.v_end_id[id].size(); i++)
		{
			DistanceEdge tmp_edge;
			tmp_edge = distance_graph_.getEnd(id, i);
			if (stnu_.s_info.id_debug == DEBUG_T)
			{
				printf("propa: ");
				v_dis[id].print();
				printf("through: ");
				tmp_edge.print();
			}

			if (tmp_edge.lb_ < -1e-5) //negative link
			{
				if (DCbackpropCheck(tmp_edge) == false)
				{
					ng_path.push_back(v_dis[id]);
					//v_dis[id].print();
					return false;
				}
				continue;
			}

			if (tmp_edge.label_ + v_dis[id].label_ == 0 && tmp_edge.label_) //same label different cases
				continue;
			DistanceEdge new_edge = tmp_edge + v_dis[id];
			DistanceEdge rep_edge;
			rep_edge = v_dis[tmp_edge.start_];
			//	printf("newprop: ");
			//	new_edge.print();
			//	rep_edge.print();
			if (rep_edge.lb_ > new_edge.lb_)
			{

				v_dis[tmp_edge.start_] = new_edge;
				q_dis.push(new_edge);
			}
		}

	}

	backprop_depth[back_id] = 1;
	if (stnu_.s_info.id_debug >= DEBUG_C)
		printf("end back prop %d %d d[%d]=%d\n", e.start_, e.end_, back_id,
				backprop_depth[back_id]);

	return true;
}

void Controllability::printResult(FILE *fp)
{
	for (int i = 0; i < v_result_.size(); i++)
	{
		fprintf(fp, "%.2lf	", v_result_[i]);
	}
	fprintf(fp, "\n");
}

long long Controllability::get_run_time()
{
	return run_time_;
}

void Controllability::extractEnvelope(vector<Conflict>& v_neg_cycle)
{
	// this function tries to enumerate negative cycles within the graph
//	std::vector<Conflict> v_neg_cycle;
	//backprop_depth = new int[distance_graph_.get_num_nodes()];
	backprop_depth.clear();
	for (int i = 0; i < distance_graph_.v_edge_.size(); i++)
		backprop_depth.push_back(0); //
	//backprop_depth.resize(distance_graph_.v_upper_id.size()
	//		+distance_graph_.v_ordinary_id.size(), 0);
	being_back.clear();
	for (int i = 0; i < distance_graph_.get_num_nodes(); i++)
	{
		this->v_q_neg_edge.push_back(vector<depth_edge>());
		being_back.push_back(false);
	}

	//vector<DistanceEdge> v_inc_neg; //restore incompleted negative cycles
	//v_inc_neg[3*i] is the source edge, 3*i+1 is the parent neg link, 3*i+2 is the path

	//being_back.resize(distance_graph_.get_num_nodes(),false);

	for (int i = 0; i < distance_graph_.v_upper_id.size(); i++)
	{
		// enumerate upper-case edge (they are negative)
		DCbackprop(distance_graph_.getUpper(i), v_neg_cycle
		//, v_inc_neg
				, 0);
		for (int j = 0; j < backprop_depth.size(); j++)
			if (backprop_depth[j] < -1)
				backprop_depth[j] = 0;
	}

	for (int i = 0; i < distance_graph_.v_ordinary_id.size(); i++)
	{
		if (distance_graph_.getOrdinary(i).lb_ < 0)
		{
			DCbackprop(distance_graph_.getOrdinary(i), v_neg_cycle,
			//v_inc_neg,
					0);

		}
	}
	if (stnu_.s_info.id_debug == DEBUG_T)
	{
		printf("print neg cycle %d\n", v_neg_cycle.size());
		for (int i = 0; i < v_neg_cycle.size(); i++)
		{
			v_neg_cycle[i].print(this->stnu_);
		}
	}
	if (stnu_.s_info.id_debug == DEBUG_T)
		printf("End of extracting conflicts.\n\n");
}

int Controllability::DCbackprop(DistanceEdge e, vector<Conflict> &neg_cycle,
//vector<DistanceEdge >& inc_neg,
		const int &depth)
{
	int back_id = get_id_backprop(e);
	if (stnu_.s_info.id_debug >= DEBUG_T)
	{
		//printf("%d[%d] %d\n",back_id,backprop_depth[back_id], backprop_depth.size());
		//printf("%d[%d] %d",back_id,backprop_depth[back_id],backprop_depth.size());
		printf("back prop: %d->%d [%d]%lf d[%d]=%d\n", e.start_, e.end_,
				e.label_, e.lb_, back_id, depth);
	}
	//cerr<<back_id<<"["<<backprop_depth[back_id]<<"]"<<endl;
	// propagation backwards
	bool b_back = true;
	if (backprop_depth[back_id] < 0)
	{
		//extract conflict
		if (stnu_.s_info.id_debug == DEBUG_T)
			printf("end back prop %d %d d[%d]=%d\n", e.start_, e.end_, back_id,
					backprop_depth[back_id]);
		b_back = false;
		return backprop_depth[back_id];
	}

	if (backprop_depth[back_id] == 1)
	{
		if (stnu_.s_info.id_debug == DEBUG_T)
			printf("end back prop %d %d d[%d]=%d\n", e.start_, e.end_, back_id,
					backprop_depth[back_id]);
		return 1;
	}

	//if(backprop_depth[back_id] == 0)
	backprop_depth[back_id] = depth - 1;
	//being_back[e.start_] = true;
	being_back[e.end_] = true;
	std::vector<DistanceEdge> dis;
	std::vector<bool> v;
	for (int i = 0; i < distance_graph_.get_num_nodes(); i++)
	{
		DistanceEdge tmp_edge;
		tmp_edge.start_ = -1;
		dis.push_back(tmp_edge);
		v.push_back(false);
	}
	//dis[e.start_] = 0;
	dis[e.start_] = e;
	v[e.end_] = true;
	//printf("inc_neg.size = %d\n",inc_neg.size());
	if (dfs(e.start_, dis, v, e, neg_cycle/*, inc_neg*/))
	{
		backprop_depth[back_id] = 1;
		being_back[e.end_] = false;
	}

	//printf("inc_neg.size = %d\n",inc_neg.size());
	//completeNegCycle(e,neg_cycle,inc_neg,e);

	if (stnu_.s_info.id_debug == DEBUG_T)
		printf("end back prop %d %d d[%d]=%d\n", e.start_, e.end_, back_id,
				backprop_depth[back_id]);

	if (backprop_depth[back_id] == 1)
		return 1;
	else
		return backprop_depth[back_id];
}

bool Controllability::dfs(int id, vector<DistanceEdge> &dis, vector<bool>& v,
		const DistanceEdge &source, vector<Conflict> &neg_cycle
		//, vector<DistanceEdge > &b
		)
{
	//enumerate negative cycle
	//bool b_return = false;
	bool return_value = true;
	if (v[id])
	{

		if (dis[id].lb_ < -1e-4)
		{
			vector<DistanceEdge> v_edge;
			v_edge.push_back(dis[id]);
			neg_cycle.push_back(this->getNegCycle(v_edge));
			if (stnu_.s_info.id_debug == DEBUG_T)
			{
				printf("ng v %d\n", id);
				dis[id].print();
			}
		}
		//negative cyle
		return return_value;
		//	b_return = true;
	}

	if (being_back[id])
	{

		//	if(dis[id].lb_ < -1e-4)
		{
			DistanceEdge tmp_edge;
			tmp_edge.start_ = tmp_edge.end_ = id;
			//b.push_back(tmp_edge); // push the end of the cycle
			//b.push_back(source); //push the previous negative edge
			//b.push_back(dis[id]); // push the path
			//dis[id].print();
			if (stnu_.s_info.id_debug == DEBUG_T)
			{
				printf("ng v back %d\n", id);
				printf("endpoint ");
				tmp_edge.print();
				printf("pre ");
				source.print();
				printf("path ");
				dis[id].print();
			}
			//cerr<<id<<" "<<being_back[id]<<" "<<v_q_neg_edge[id].size()<<endl;
			//printf("%d %\n",id,being_back[id]);
			for (int i = 0; i < this->v_q_neg_edge[id].size(); i++)
			{
				//if(v_q_neg_edge[id][i].edge.start_ == dis[id].end_) // ng
				{
					//printf("qc: ");
					//(v_q_neg_edge[id][i].edge+dis[id]).print();
					vector<DistanceEdge> v;
					//v.push_back();
					v.push_back(v_q_neg_edge[id][i].edge + dis[id]);

					Conflict tmp_conf;
					tmp_conf = this->getNegCycle(v);
					if (tmp_conf.v_conflict_edge_id.size())
					{
						//	tmp_conf.print();
						neg_cycle.push_back(this->getNegCycle(v));
						continue;
					}
				}
				//v_q_neg_edge[id].push_back(depth_edge(0,v_q_neg_edge[id][i].edge+dis[id]));
				//b.push_back(tmp_edge); // push the end of the cycle
				//b.push_back(source); //push the previous negative edge
				//b.push_back(v_q_neg_edge[id][i].edge+dis[id]); // push the path
				//v_q_neg_edge[id][i].edge.print();
				//dis[id].print();
				//printf("q: %d ",id);
				//	DistanceEdge new_edge = v_q_neg_edge[id][i].edge+dis[id];
				//v_q_neg_edge[new_edge.end_].push_back(depth_edge(1,new_edge));
				//v_q_neg_edge[id][i].edge.print();
				//	dis[id].print();
				//v_q_neg_edge[new_edge.end_][v_q_neg_edge[new_edge.end_].size()-1].edge.print();

				//printf("%d %d\n",new_edge.end_,id);
				//	printf("\t new%d: ",new_edge.end_);
				//new_edge.print();
			}

		}

		//return;
	}
	//if(b_return)
	//	return return_value;

	if (stnu_.s_info.id_debug == DEBUG_T)
		printf("id=%d\n", id);
	v[id] = true;

	for (int i = 0; i < distance_graph_.v_end_id[id].size(); i++)
	{
		DistanceEdge tmp_edge;
		tmp_edge = distance_graph_.getEnd(id, i);
		if (tmp_edge.b_block_)
		{
			distance_graph_.v_edge_[distance_graph_.v_end_id[id][i]].b_block_ =
					false;
			continue;
		}
		if (stnu_.s_info.id_debug == DEBUG_T)
		{
			printf("propa: ");
			dis[id].print();
			printf("through: ");
			tmp_edge.print();
		}
		//printf("%d %d %lf\n",tmp_edge.start_, tmp_edge.end_,tmp_edge.lb_);
		if (tmp_edge.lb_ < 0) //negative link
		{
			int tmp_return_value = DCbackprop(tmp_edge, neg_cycle, // b ,
					backprop_depth[get_id_backprop(source)]);
			if (tmp_return_value <= 0)
			{
				return_value = false;
				being_back[tmp_edge.end_] = true;

				//b.push_back(tmp_edge); // push the end of the cycle
				//b.push_back(source); //push the previous negative edge
				//b.push_back(dis[id]); // push the path

				//dis[id].print();
//				if(stnu_.s_info.id_debug >= DEBUG_T){
//					printf("end ");
//					tmp_edge.print();
//					printf("pre ");
//					source.print();
//					printf("path ");
//					dis[id].print();
//				}

				//dis[id].print();
				//if(tmp_return_value == -1)
				{
					v_q_neg_edge[dis[id].end_].push_back(
							depth_edge(tmp_return_value, dis[id]));
					//	continue;
				}
				for (int i = 0; i < v_q_neg_edge[dis[id].start_].size(); i++)
				{
					if (v_q_neg_edge[dis[id].start_][i].edge.start_
							== v_q_neg_edge[dis[id].start_][i].edge.end_)
						continue;
					//if(v_q_neg_edge[dis[id].start_][i].edge.start_ == dis[id].end_)
					//if(v_q_neg_edge[id][i].edge.start_ == dis[id].end_) // ng
					{
						//printf("qc: ");
						//(v_q_neg_edge[dis[id].start_][i].edge+dis[id]).print();
						vector<DistanceEdge> v;
						//v.push_back();
						v.push_back(
								v_q_neg_edge[dis[id].start_][i].edge + dis[id]);

						Conflict tmp_conf;
						tmp_conf = this->getNegCycle(v);
						if (tmp_conf.v_conflict_edge_id.size())
						{
							neg_cycle.push_back(this->getNegCycle(v));
							continue;
						}
					}
					v_q_neg_edge[dis[id].end_].push_back(
							depth_edge(tmp_return_value,
									v_q_neg_edge[dis[id].start_][i].edge
											+ dis[id]));
					//	v_q_neg_edge[dis[id].start_][i].edge.print();
					//	dis[id].print();
					//	v_q_neg_edge[dis[id].end_][v_q_neg_edge[dis[id].end_].size()-1].edge.print();

				}

				if (stnu_.s_info.id_debug >= DEBUG_T)
				{

					for (int i = 0; i < v_q_neg_edge[dis[id].end_].size(); i++)
					{
						//printf("id source: %d ",dis[id].end_);
						//v_q_neg_edge[dis[id].end_][i].edge.print();
					}
					//printf("ng back %d\n",id);
				}

			}
			if (tmp_return_value)
			{
//			for (int j = 0; 3*j < b.size(); j++)
//			{
//
//				if(b[3*j+1] != tmp_edge) // the previous negative edge of inc_neg is the current edge
//					continue;
//				if(stnu_.s_info.id_debug == DEBUG_T){
//					b[3*j+2].print();
//					(dis[id]).print();
//				}
//				b[3*j+2] = b[3*j+2] + dis[id];
//				b[3*j+1] = source; // update the previous neg edge
//
//				if(stnu_.s_info.id_debug == DEBUG_T){
//					printf("added: ");
//					b[j*3+2].print();
//					printf("pre ");
//					source.print();
//				}
//
//			}
			}
			continue;
		}

		if (tmp_edge.label_ + source.label_ == 0 && tmp_edge.label_) //same label different cases
			continue;
		DistanceEdge new_edge = tmp_edge + dis[id];

		bool bloop = false;
		for (int it = 0; it < tmp_edge.v_support_.size(); it++)
		{

			for (int il = 0; il < dis[id].v_support_.size(); il++)
				if (dis[id].v_support_[il] == tmp_edge.v_support_[it])
				{
					bloop = true;
					break;
				}
			if (bloop)
				break;
		}
		if (bloop)
			continue;
		if (new_edge.lb_ >= 0)
		{
			if (new_edge.label_)
				new_edge.label_ = 0;
			if (this->isTighterOrdinaryEdge(new_edge))
			{
				//if(new_edge.lb_ > source.ub_)
				//	new_edge.label_ = 0;
				addOrdinaryEdge(new_edge);
				//backprop_depth.resize(distance_graph_.v_edge_.size(), 0);

				if (stnu_.s_info.id_debug == DEBUG_T)
				{
					tmp_edge.print();
					dis[id].print();
					printf("new edge: ");
					new_edge.print();
				}

			}
		}
		else
		{
			DistanceEdge rep_edge;
			rep_edge = dis[tmp_edge.start_];
			dis[tmp_edge.start_] = new_edge;
			if (v[tmp_edge.start_] == false || tmp_edge.start_ == source.end_)
			{
				if (dfs(tmp_edge.start_, dis, v, source, neg_cycle
				//, b
						) == false)
					return_value = false;
			}
			dis[tmp_edge.start_] = rep_edge;
		}
	}
	v[id] = false;
	return return_value;
}

void Controllability::completeNegCycle(DistanceEdge e, vector<Conflict>&c,
		vector<DistanceEdge>& inc_neg, DistanceEdge& dis)
{

	for (int i = inc_neg.size() / 3 - 1; i >= 0; i--)
	{
		//inc_neg[i*3+2].print();
		//printf("f1=%d f2=%d\n",int(inc_neg[3*i].start_ == inc_neg[3*i].end_),
		//		int(inc_neg[3*i] != e && inc_neg[i*3+2].start_ != inc_neg[i*3+2].end_));
		//inc_neg[i].second.push_back(dis);
		if (inc_neg[3 * i].start_ == inc_neg[3 * i].end_
				&& inc_neg[i * 3 + 2].start_ != inc_neg[i * 3 + 2].end_)//source is not the current link
		//point
		{
			if (inc_neg[3 * i].start_ != e.end_
					&& inc_neg[3 * i].end_ != e.start_)
				continue;
		}
		else if (inc_neg[3 * i] != e
				&& inc_neg[i * 3 + 2].start_ != inc_neg[i * 3 + 2].end_)//source is not the current link
		{
			continue;
		}
		vector<DistanceEdge> v;
		v.push_back(inc_neg[i * 3 + 2]);
		Conflict tmp_conf;
		tmp_conf = this->getNegCycle(v);
		if (tmp_conf.v_conflict_edge_id.size())
		{
			//printf("comp ng: ");
			//inc_neg[i*3+2].print();
			c.push_back(this->getNegCycle(v));
		}

		if (inc_neg[3 * i].start_ != inc_neg[3 * i].end_
				|| inc_neg[3 * i + 2].lb_ > 0)
			inc_neg.erase(inc_neg.begin() + 3 * i,
					inc_neg.begin() + 3 * (i + 1));

	}
	//exit(0);
}

int Controllability::get_id_backprop(const DistanceEdge &edge)
{
	int res;
	res = distance_graph_.edge_to_id(edge);

	return res;
}

Controllability::~Controllability()
{
	// TODO Auto-generated destructor stub
}

void Conflict::print() const
{
	printf("Printing Conflict\n");
	printf("Assignments (%d):", this->v_assignments.size());
	for (int i = 0; i < v_assignments.size(); i++)
		printf(" v%d=%d", v_assignments[i].first, v_assignments[i].second);
	if (v_assignments.size())
		printf("\n");
	printf("Number of Edges: %d\n", this->v_conflict_edge_id.size());
	for (int i = 0; i < this->v_conflict_edge_id.size(); i++)
	{
		printf("%d 	", v_conflict_edge_id[i]);
	}
	printf("Complete Conflict Printing\n");
}

void Conflict::print(const stnu& a) const
{
	//print();
	printf("Printing Conflict\n");
	printf(ANSI_COLOR_RED); // change the color of output

	printf("Assignments (%d):", this->v_assignments.size());

	for (int i = 0; i < v_assignments.size(); i++)
		printf(" v%d=%d", v_assignments[i].first, v_assignments[i].second);
	printf(ANSI_COLOR_RESET); // change the color of output

	if (v_assignments.size())
		printf("\n");
	printf("Number of Edges: %d\n", this->v_conflict_edge_id.size());
	for (int i = 0; i < this->v_conflict_edge_id.size(); i++)
	{
		printf("%d ", v_conflict_edge_id[i]);
		printf("%d ", a.idToEdge(v_conflict_edge_id[i]).second);
		LINK s = a.idToEdge(v_conflict_edge_id[i]).first;
		if (s.start >= 0)
			s.print();

		//printf("%d 	",v_conflict_edge_id[i]);
	}
	printf("Complete Conflict Printing\n");
}

void SimConflict::removeCtg(int cid)
{
	for (int i = 0; i < v_variable.size(); i++)
	{
		if (v_variable[i].edge_id == cid || v_variable[i].edge_id == -cid)
		{
			v_variable.erase(v_variable.begin() + i);
			return;
		}
	}
}

bool SimConflict::isNew(const vector<SimConflict> a) const
{
	//if the current constraint is implied by a, return false

	for (int i = 0; i < a.size(); i++)
	{

		if (a[i].conf_type != conf_type)
		{
			// if the current conflict is in a different type with a[i]
			// then a[i] cannot imply the current conflict
			continue;
		}

		for (int j = 0; j < v_variable.size(); j++)
		{

		}

	}
	return true;
}

} /* namespace CDS */
