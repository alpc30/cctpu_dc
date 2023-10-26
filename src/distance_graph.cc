/*
 * distance_graph.cc
 *
 *  Created on: 16/09/2016
 *      Author: jing
 */

#include "distance_graph.h"

namespace CDS
{

int DistanceGraph::add_edge(const DistanceEdge &e)
{
	if (v_end_id.size() <= e.end_)
	{
		e.print();
		printf("error: the v_end_id is smaller than the end node.\n");
	}
	if (v_start_id.size() <= e.start_)
	{
		e.print();
		printf("error: the v_start_id is smaller than the start node.\n");
	}

	int edge_id = v_edge_.size();
	v_edge_.push_back(e);
	v_end_id[e.end_].push_back(edge_id);
	v_start_id[e.start_].push_back(edge_id);
	if (e.label_ == 0)
		mp_index_edge_[make_pair(e.start_, e.end_)] = edge_id + 1;
	else if (e.label_ > 0) //upper
	{
		mp_start_upper[make_pair(e.start_, e.end_)] = edge_id + 1;
	}
	else
		mp_start_lower[e.start_] = edge_id + 1;

	return edge_id;
}

void DistanceGraph::addLowerCaseEdge(const DistanceEdge &e)
{
	v_lower_id.push_back(add_edge(e));
}

void DistanceGraph::addUpperCaseEdge(const DistanceEdge &e)
{
	v_upper_id.push_back(add_edge(e));
}
void DistanceGraph::addOrdinaryEdge(const DistanceEdge &e)
{
	v_ordinary_id.push_back(add_edge(e));
}

DistanceEdge DistanceGraph::getEdge(const int & id)
{
	if (id >= v_edge_.size())
	{
		printf("Error: the edge id is out of v_edge!\n");
		exit(0);
	}
	return v_edge_[id];
}

DistanceEdge DistanceGraph::getLower(const int & id)
{
	if (id >= v_lower_id.size())
	{
		printf("Error: the id is out of v_lower_id!\n");
		exit(0);
	}
	return getEdge(v_lower_id[id]);
}

DistanceEdge DistanceGraph::getUpper(const int & id)
{
	if (id >= v_upper_id.size())
	{
		printf("Error: the id is out of v_upper_id!\n");
		exit(0);
	}
	return getEdge(v_upper_id[id]);
}

DistanceEdge DistanceGraph::getOrdinary(const int & id)
{
	if (id >= v_ordinary_id.size())
	{
		printf("Error: the id is out of v_ordinary_id!\n");
		exit(0);
	}
	return getEdge(v_ordinary_id[id]);
}

DistanceEdge DistanceGraph::getStart(const int &st, const int &id)
{
	if (st >= v_start_id.size())
	{
		printf("Error: the start node is out of v_start_id!\n");
		exit(0);
	}

	if (id >= v_start_id[st].size())
	{
		printf("Error: the id is out of v_start_id[start]!\n");
		exit(0);
	}

	return getEdge(v_start_id[st][id]);
}
DistanceEdge DistanceGraph::getEnd(const int &en, const int &id)
{
	if (en >= v_end_id.size())
	{
		printf("Error: the end node is out of v_end_id!\n");
		exit(0);
	}

	if (id >= v_end_id[en].size())
	{
		printf("Error: the id is out of v_end_id[end]!\n");
		exit(0);
	}

	return getEdge(v_end_id[en][id]);
}

DistanceEdge DistanceEdge::operator +(const DistanceEdge tmp_edge) const
{
	DistanceEdge new_edge;
	new_edge.start_ = this->start_;
	new_edge.end_ = tmp_edge.end_;
	new_edge.lb_ = this->lb_ + tmp_edge.lb_;
	new_edge.v_support_ = this->v_support_;
	new_edge.v_support_.insert(new_edge.v_support_.begin(),
			tmp_edge.v_support_.begin(), tmp_edge.v_support_.end());

	new_edge.label_ = tmp_edge.label_;
	return new_edge;
}

void DistanceGraph::PrintDot(std::string output_file)
{
	FILE *fdot;
	fdot = fopen(output_file.c_str(), "w");
	fprintf(fdot, "digraph G {\n rankdir = LR;\n nodesep = .45; \n "
			"size = 30;\nlabel=\"result\";\n");

	for (int i = 0; i < v_edge_.size(); ++i)
	{
		DistanceEdge tmp_edge;
		tmp_edge = v_edge_[i];
		if (tmp_edge.lb_ != g_inf)
			fprintf(fdot, "\"%d\"->\"%d\"[label = \"[%.2lf] %d\"];\n",
					tmp_edge.start_, tmp_edge.end_, tmp_edge.lb_, 0);
	}

	fprintf(fdot, "}\n");
	fclose(fdot);
}

int DistanceGraph::edge_to_id(DistanceEdge edge)
{
	int res;
	if (edge.label_ == 0)
	{
	  if(mp_index_edge_.find(make_pair(edge.start_, edge.end_)) ==
	     mp_index_edge_.end())
	  {
	    printf("error: distance_graph_ doesn't have %d-%d\n",edge.start_,
		   edge.end_);
	  }
		res = mp_index_edge_[make_pair(edge.start_, edge.end_)] - 1;
	}
	else if (edge.label_ > 0) //upper
	{
	   if(mp_start_upper.find(make_pair(edge.start_, edge.end_)) ==
	     mp_start_upper.end())
	  {
	    printf("error: distance_graph_ doesn't have %d-%d\n",edge.start_,
		   edge.end_);
	  }

		res = mp_start_upper[make_pair(edge.start_, edge.end_)] - 1;
	}
	else
	{
	  if(mp_start_lower.find(edge.start_) == mp_start_lower.end())
	  {
	    printf("error: distance_graph_\n");
	  }
		res = mp_start_lower[edge.start_] - 1;
	}

//	if(edge.label_ == 0)//ordinary
//	{
//		res = mp_index_edge_[make_pair(edge.start_, edge.end_)] - 1;
//	}
//	else if(edge.label_  < 0) // lower
//	{
//		res = mp_end_lower[edge.end_] - 1;
//	}
//	else //upper
//	{
//		res = mp_start_upper[edge.start_] - 1;
//	}
	return res;

}

void DistanceEdge::print() const
{
	if (this->start_name.length() == 0 || true)
		printf("%d -> %d: [%lf, %lf] lab = %d\n", start_, end_, lb_, ub_,
				label_);
	else
		printf("%s -> %s: [%lf, %lf] lab = %d\n", start_name.c_str(),
				end_name.c_str(), lb_, ub_, label_);
	if (this->v_support_.size())
	{
		printf("\t support links (%d): ", v_support_.size());
		for (int i = 0; i < this->v_support_.size(); i++)
		{
			//int id = this->v_support_[i];
			printf("%d ", this->v_support_[i]);
		}
		printf("\n");

//		for(int i = this->v_support_names.size()-1; i >=0; i--)
//			printf("%s ",this->v_support_names[i].c_str());
//		printf("\n");
	}

	if (this->v_label.size())
	{
		printf("\tLabels (%d):", this->v_label.size());
		for (int i = 0; i < this->v_label.size(); i++)
		{
			printf(" V%d=%d ", v_label[i].first, v_label[i].second);
		}
		printf("\n");

	}
}
}

