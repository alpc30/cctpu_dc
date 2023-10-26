#include "reader.h"
#include "structures.h"
#include "utils.h"
#include <fstream>
#include <string>
#include <vector>
#include <map>
#define CT_LINK 0
#define CT_VAR 1
#define CTGL 1
#define RQML 0

namespace CDS
{

const double INF = g_inf;

void Reader::read(std::string problem_file, stnu & st)
{
	std::ifstream reader(problem_file.data());
	if (!reader.good())
	{
		printf("No such file!\n");
		exit(0);
	}
	reader.close();
	std::vector<std::string> tokens;
	TRANSPLAN::Utils::tokenize(problem_file, tokens, std::string("."));
	t = &st;
	//printf("%s\n",tokens[0].c_str());
	//if (st.s_info.id_debug)
	//	printf("%s\n", tokens[tokens.size() - 1].c_str());
	//judge the type of the input file
	if (tokens[tokens.size() - 1].compare("pos") == 0)
	{
		this->readpos(problem_file);
	}
	else if (tokens[tokens.size() - 1].compare("cctp") == 0)
	{
		this->readcctp(problem_file);
	}/*
	 else if(tokens[tokens.size()-1].compare("out")==0)
	 {
	 this->readout(problem_file);
	 }*/
	else if (tokens[tokens.size() - 1].compare("gr") == 0)
	{
		this->readgrid(problem_file);
	}/*
	 else if(tokens[tokens.size()-1].compare("PPSTP")==0)
	 {
	 this->readppstp(problem_file);
	 }
	 */
	else if (tokens[tokens.size() - 1].compare("stnu") == 0)
		this->readevacplan(problem_file);
	else if (tokens[tokens.size() - 1].compare("cstnu") == 0)
		this->readcstnu(problem_file);
	else if (tokens[tokens.size() - 1].compare("cst") == 0)
		readcst(problem_file);
	else
	{
		printf("Wrong Input File!\n");
		exit(0);
	}

}

void Reader::readcst(std::string problem_file)
{
	std::ifstream reader(problem_file.data());
	std::string line;
	t->file_name = problem_file;
	t->n_nodes = 0;
	int n_links = 0;
	while (std::getline(reader, line))
	{
		std::vector<std::string> tokens;
		//std::cerr<<"line"<<line<<std::endl;
		TRANSPLAN::Utils::tokenize(line, tokens, std::string(" []\t"));
		if (tokens.size() == 0)
			continue;
		int id_t = 0;
		if (t->n_nodes == 0) // the amount has not been read
		{
			t->n_nodes = atoi(tokens[id_t++].c_str());
			n_links = atoi(tokens[id_t++].c_str());
			t->n_vars = atoi(tokens[id_t++].c_str());
			t->utility = atof(tokens[id_t].c_str());
			continue;
		}

		if (t->n_vars > t->v_dis_var.size()) // the vars have not been read
		{
			DiscreteVariable tmp_var;
			tmp_var.n_value = atoi(tokens[id_t++].c_str());
			for (int i = 0; i < tmp_var.n_value; i++)
			{
				tmp_var.v_values.push_back(atoi(tokens[id_t++].c_str()));
				tmp_var.vv_id_mute_edge.push_back(std::vector<int>());
			}
			for (int i = 0; i < tmp_var.n_value; i++)
			{
				tmp_var.v_utilities.push_back(atof(tokens[id_t++].c_str()));
			}

			tmp_var.id = t->v_dis_var.size();
			t->v_dis_var.push_back(tmp_var);
			continue;
		}

		LINK tmp_edge;
		tmp_edge.start_name = tokens[id_t];
		tmp_edge.start = atoi(tokens[id_t++].c_str());
		tmp_edge.end_name = tokens[id_t];
		tmp_edge.end = atoi(tokens[id_t++].c_str());
		bool b_controllable = false;
		if (tokens[id_t][0] == 'r' || tokens[id_t][0] == 'R')
			b_controllable = true;
		id_t++;
		tmp_edge.lb = atof(tokens[id_t++].c_str());
		tmp_edge.ub = atof(tokens[id_t++].c_str());
		tmp_edge.lb_cost_r = atof(tokens[id_t++].c_str());
		tmp_edge.ub_cost_r = atof(tokens[id_t++].c_str());

		int n_labels = atoi(tokens[id_t].c_str());
		if (n_labels)
		{
			std::getline(reader, line);
			tokens.clear();
			TRANSPLAN::Utils::tokenize(line, tokens, std::string(" []\t"));
			id_t = 0;
			while (id_t < tokens.size())
			{
				int id_var = atoi(tokens[id_t++].c_str());
				int id_value = atoi(tokens[id_t++].c_str());
				tmp_edge.v_labels.push_back(make_pair(id_var, id_value));

				for (int i = 0; i < t->v_dis_var[id_var].v_values.size(); i++)
				{
					if (t->v_dis_var[id_var].v_values[i] != id_value)
					{
						int id_link = 0;
						if (b_controllable != CTGL)
							id_link = -t->ctg.size() - 1;
						else
							id_link = t->rqm.size() + 1;
						t->v_dis_var[id_var].vv_id_mute_edge[i].push_back(
								id_link);
					}
				}
			}

		}
		if (b_controllable)
			t->add_rqm(tmp_edge);
		//t->rqm.push_back(tmp_edge);
		else
			t->add_ctg(tmp_edge);
		//t->ctg.push_back(tmp_edge);

	}
	t->init(t->n_nodes);
	t->n_ctg = t->ctg.size();
	t->n_rqm = t->rqm.size();

}

void Reader::readcstnu(std::string problem_file)
{
	std::ifstream reader(problem_file.data());
	std::string line;

	t->n_ctg = 0;
	t->n_rqm = 0;
	t->ctg.clear();
	t->rqm.clear();
	t->file_name = problem_file;
	std::map<std::string, int> nid;
	std::map<std::string, int>::iterator it;
	nid.clear();

	int laid = 0;
	while (std::getline(reader, line))
	{
		std::vector<std::string> tokens;
		// std::cerr<<"line"<<line<<std::endl;
		TRANSPLAN::Utils::tokenize(line, tokens, std::string(" \t"));
		if (tokens.size() == 0)
			continue;

		int content_type = atoi(tokens[0].c_str());
		int pos_token = 1;

		if (content_type == CT_LINK)
		{
			std::string s_st, s_en;
			s_st = tokens[pos_token++];
			s_en = tokens[pos_token++];
			if (nid.find(s_st) == nid.end())
			{
				nid[s_st] = laid++;
			}
			if (nid.find(s_en) == nid.end())
			{
				nid[s_en] = laid++;
			}

			int st = nid[s_st];
			int en = nid[s_en];
			printf("%d %d\n", st, en);
			double lb = atof(tokens[pos_token++].c_str());
			double ub = atof(tokens[pos_token++].c_str());

			int type = getType(tokens[pos_token++]);
			if (type)
			{
				CTG tmp_ctg = LINK(st, en);
				tmp_ctg.lb = lb;
				tmp_ctg.ub = ub;
				//debug 20151207 add name of the node
				tmp_ctg.start_name = s_st;
				tmp_ctg.end_name = s_en;
				tmp_ctg.setOriginal();
				int n_pri_var = atoi(tokens[pos_token++].c_str());
				for (int i = 0; i < n_pri_var; i++)
				{
					tmp_ctg.v_pri_vars.push_back(
							atoi(tokens[pos_token++].c_str()) - 1);
				}
				//end debug 20151207
				t->add_ctg(tmp_ctg);
				//t->ctg.push_back(tmp_ctg);
			}
			else
			{
				RQM tmp_rqm = LINK(st, en);
				tmp_rqm.lb = lb;
				tmp_rqm.ub = ub;
				//debug 20151207 add name of the nodes to the link
				tmp_rqm.start_name = s_st;
				tmp_rqm.end_name = s_en;
				//end of debug 20151207
				t->rqm.push_back(tmp_rqm);
			}
		}
		else if (content_type == CT_VAR)
		{
			// DisVar tmp_var;
			// tmp_var.id = t->v_dis_vars.size();

			//	 for(int i = 0; i < atoi(tokens[pos_token].c_str()); i++)
//			 {
//				 tmp_var.v_dc.push_back(i);
//			 }
//			 pos_token++;
			// t->n_vars = atoi(tokens[pos_token++].c_str());
			// for(int i = 0; i < t->n_vars; i++)
			//{
			DiscreteVariable tmp_dis_var;
			tmp_dis_var.name = tokens[pos_token++];
			tmp_dis_var.id = nid[tmp_dis_var.name];
			tmp_dis_var.n_value = atoi(tokens[pos_token++].c_str());
			for (int i = 0; i < tmp_dis_var.n_value; i++)
			{
				int n_tmp = atoi(tokens[pos_token++].c_str());
				tmp_dis_var.vv_id_mute_edge.push_back(std::vector<int>());
				for (int j = 0; j < n_tmp; j++)
					tmp_dis_var.vv_id_mute_edge[i].push_back(
							atoi(tokens[pos_token++].c_str()));
			}

			while (tmp_dis_var.v_utilities.size() < tmp_dis_var.n_value)
				tmp_dis_var.v_utilities.push_back(0);
			int tmp = 0;
			while (tmp_dis_var.v_values.size() < tmp_dis_var.n_value)
				tmp_dis_var.v_values.push_back(tmp++);
			tmp_dis_var.id = t->v_dis_var.size();
			t->v_dis_var.push_back(tmp_dis_var);
			// }
		}
	}

//	FILE *fname = fopen("node_id_80.out", "w");
//	for(it = nid.begin(); it != nid.end(); it++)
//	{
//		fprintf(fname, "%s %d\n", (it->first).c_str(), it->second);
//	}
//	fclose(fname);

	t->n_ctg = t->ctg.size();
	t->n_rqm = t->rqm.size();
	t->n_nodes = nid.size();
	t->n_vars = t->v_dis_var.size();
	//printf("stat %d	%d	%d\n",t->n_nodes, t->n_rqm, t->n_ctg);
	//exit(0);
	t->init(t->n_nodes);

	for (int i = 0; i < t->n_ctg; i++)
	{
		int st = t->ctg[i].start;
		int en = t->ctg[i].end;
		double lb = t->ctg[i].lb;
		double ub = t->ctg[i].ub;
		printf("c %d %d %.2lf %.2lf\n", st, en, lb, ub);
		//	t->set_min(st, en, lb);
		//	t->set_max(st, en, ub);
	}

	for (int i = 0; i < t->n_rqm; i++)
	{
		int st = t->rqm[i].start;
		int en = t->rqm[i].end;
		double lb = t->rqm[i].lb;
		double ub = t->rqm[i].ub;

		printf("r %d %d %.2lf %.2lf\n", st, en, lb, ub);
		//t->set_min(st, en, lb);
		//t->set_max(st, en, ub);
	}
}

void Reader::readpos(std::string problem_file)
{

	std::ifstream reader(problem_file.data());
	std::string line;
	std::map<std::pair<int, int>, int> mp_nodes_rqm;
	mp_nodes_rqm.clear();
	getline(reader, line);

	std::vector<std::string> tokens;
	TRANSPLAN::Utils::tokenize(line, tokens, std::string(" \t"));

	int n_act = atoi(tokens[0].c_str());
	t->fldt = atof(tokens[1].c_str());
	t->n_ctg = n_act - 2; //2 dummy activities (source and sink)
	t->ctg.clear();
	t->rqm.clear();
	t->n_nodes = t->n_ctg * 2 + 2;
	t->file_name = problem_file;
	//printf("%s\n",problem_file.c_str());

	if (tokens.size() > 2)
	{
		t->value_mdsc = atof(tokens[2].c_str());
		t->value_mddc = atof(tokens[3].c_str());
	}
	else
	{
		if (t->s_info.id_obj == O_MAX_DELAY)
			t->value_mddc = t->value_mdsc = g_inf;
		else
			t->value_mddc = t->value_mdsc = 0;
	}

	for (int i = 0; i < t->n_ctg; i++)
	{
		CDS::CTG tmp_ctg;
		tmp_ctg.start = 2 * i + 1;
		tmp_ctg.end = 2 * i + 2;
		t->add_ctg(tmp_ctg);
		//t->ctg.push_back(tmp_ctg);
	}
	std::vector<std::vector<int> > v_succ;
	std::vector<std::vector<double> > v_time_lag;

	for (int i = 0; i < n_act; ++i)
	{
		std::getline(reader, line);
		std::vector<std::string> tokens;
		TRANSPLAN::Utils::tokenize(line, tokens, std::string(" \t"));
		std::vector<int> succ;
		succ.clear();

		int start = 1;
		int dur = atoi(tokens[start++].c_str());
		if (i != 0 && i != n_act - 1)
		{
			t->ctg[i - 1].lb = dur;
			t->ctg[i - 1].ub = dur;

			if (t->s_info.id_obj == O_MAX_DELAY)
			{
				t->ctg[i - 1].ub = 1e6;
			}
		}

		int est = atoi(tokens[start++].c_str());

		int n_succ = atoi(tokens[start++].c_str());
		for (int s = 0; s < n_succ; ++s)
		{
			int suc = atoi(tokens[start++].c_str());
			succ.push_back(suc);
			//Transplan::activities[task].successors.push_back(suc);
		}

		v_succ.push_back(succ);
		std::vector<double> t_lag;
		t_lag.clear();
		for (int s = 0; s < n_succ; ++s)
		{
			std::vector<std::string> time_lag;

			TRANSPLAN::Utils::tokenize(tokens[start++], time_lag,
					std::string(" []"));
			assert(time_lag.size() >= 1);
			int lag = atoi(time_lag[0].c_str());
			t_lag.push_back(lag);
		}
		v_time_lag.push_back(t_lag);
	}

	for (int i = 0; i < n_act; i++)
	{
		CTG tmp_ctg;
		if (i == 0)
		{
			tmp_ctg.start = tmp_ctg.end = tmp_ctg.lb = 0;
		}
		else if (i == n_act - 1)
		{
			tmp_ctg.start = tmp_ctg.end = t->n_nodes - 1;
			tmp_ctg.lb = 0;
		}
		else
			tmp_ctg = t->ctg[i - 1];
		for (int j = 0; j < v_succ[i].size(); j++)
		{
			CTG tmp_succ_ctg;
			if (v_succ[i][j] == 0)
			{
				tmp_succ_ctg.start = tmp_succ_ctg.end = tmp_succ_ctg.lb = 0;
			}
			else if (v_succ[i][j] == n_act - 1)
			{
				tmp_succ_ctg.start = tmp_succ_ctg.end = t->n_nodes - 1;
				tmp_succ_ctg.lb = 0;
			}
			else
				tmp_succ_ctg = t->ctg[v_succ[i][j] - 1];

			CDS::RQM tmp_rqm;
			if (v_time_lag[i][j] >= 0)
			{
				tmp_rqm.start = tmp_ctg.end;
				tmp_rqm.end = tmp_succ_ctg.start;
				tmp_rqm.lb = v_time_lag[i][j] - tmp_ctg.lb;
				tmp_rqm.ub = INF;
			}
			else
			{
				tmp_rqm.start = tmp_succ_ctg.end;
				tmp_rqm.end = tmp_ctg.start;
				tmp_rqm.ub = -v_time_lag[i][j] - tmp_succ_ctg.lb;
				tmp_rqm.lb = -INF;
			}

			if (mp_nodes_rqm[
			{ tmp_rqm.start, tmp_rqm.end }] == 0)
			{
				//t->rqm.push_back(tmp_rqm);
				t->add_rqm(tmp_rqm);
				mp_nodes_rqm[
				{ tmp_rqm.start, tmp_rqm.end }] = t->rqm.size();
			}
			else
			{
				int tmp_id_rqm = mp_nodes_rqm[
				{ tmp_rqm.start, tmp_rqm.end }] - 1;
				if (t->rqm[tmp_id_rqm].lb < tmp_rqm.lb)
					t->rqm[tmp_id_rqm].lb = tmp_rqm.lb;
				if (t->rqm[tmp_id_rqm].ub > tmp_rqm.ub)
					t->rqm[tmp_id_rqm].ub = tmp_rqm.ub;
			}
		}
	}

}

void Reader::readevacplan(std::string problem_file)
{
	std::ifstream reader(problem_file.data());
	std::string line;

	t->n_ctg = 0;
	t->n_rqm = 0;
	t->ctg.clear();
	t->rqm.clear();
	t->file_name = problem_file;
	std::map<std::string, int> nid;
	std::map<std::string, int>::iterator it;
	nid.clear();

	int laid = 0;
	while (std::getline(reader, line))
	{
		std::vector<std::string> tokens;
		// std::cerr<<"line"<<line<<std::endl;
		TRANSPLAN::Utils::tokenize(line, tokens, std::string(" \t"));
		std::string s_st, s_en;
		s_st = tokens[0];
		s_en = tokens[1];
		if (nid.find(s_st) == nid.end())
		{
			nid[s_st] = laid++;
		}
		if (nid.find(s_en) == nid.end())
		{
			nid[s_en] = laid++;
		}

		int st = nid[s_st];
		int en = nid[s_en];
		printf("%d %d\n", st, en);
		double lb = atof(tokens[2].c_str());
		double ub = atof(tokens[3].c_str());

		int type = atoi(tokens[4].c_str());

		// t->set_min(st, en, lb);
		// t->set_max(st, en, ub);

		if (type)
		{
			CTG tmp_ctg = LINK(st, en);
			tmp_ctg.lb = lb;
			tmp_ctg.ub = ub;
			t->ctg.push_back(tmp_ctg);
		}
		else
		{
			RQM tmp_rqm = LINK(st, en);
			tmp_rqm.lb = lb;
			tmp_rqm.ub = ub;
			t->rqm.push_back(tmp_rqm);

		}
	}

//	FILE *fname = fopen("node_id_80.out", "w");
//	for(it = nid.begin(); it != nid.end(); it++)
//	{
//		fprintf(fname, "%s %d\n", (it->first).c_str(), it->second);
//	}
//	fclose(fname);

	t->n_ctg = t->ctg.size();
	t->n_rqm = t->rqm.size();
	t->n_nodes = nid.size();
	//printf("stat %d	%d	%d\n",t->n_nodes, t->n_rqm, t->n_ctg);
	//exit(0);
	t->init(t->n_nodes);

	for (int i = 0; i < t->n_ctg; i++)
	{
		int st = t->ctg[i].start;
		int en = t->ctg[i].end;
		double lb = t->ctg[i].lb;
		double ub = t->ctg[i].ub;
		printf("c %d %d %.2lf %.2lf\n", st, en, lb, ub);
		//	t->set_min(st, en, lb);
		//	t->set_max(st, en, ub);
	}

	for (int i = 0; i < t->n_rqm; i++)
	{
		int st = t->rqm[i].start;
		int en = t->rqm[i].end;
		double lb = t->rqm[i].lb;
		double ub = t->rqm[i].ub;

		printf("r %d %d %.2lf %.2lf\n", st, en, lb, ub);
		//t->set_min(st, en, lb);
		//t->set_max(st, en, ub);
	}
}
/*
 void Reader::readppstp(std::string problem_file)
 {
 char dotfile[]="2.dot";
 printf("%s\n",dotfile);
 FILE *fdot =fopen(dotfile,"w");
 fprintf(fdot,"digraph G {\n rankdir = LR;\n nodesep = .45; \n size = 30;\nlabel=\"CCTP\";\n");
 std::map<std::string,int> name;
 name.clear();
 CcpstpProb* prob = ccpstpParser(problem_file.c_str());
 //printf("%s %s\n",prob->start_name,prob->objective_fn);
 stnu *t = this->getOwner();
 int n_nodes = 1;
 int *ind = new int[2000];
 for(int i = 0; i < 2000; i++)ind[i] = 0;
 if(t)
 {
 t->n_nodes = prob->n_nodes;
 t->init(t->n_nodes);
 t->rqm.clear();
 CcpstpContArc * t_arc;
 for(t_arc = prob->arc_list; t_arc; t_arc = t_arc->next)
 {
 int st,en;

 if(name.find((std::string)t_arc->start)!=name.end())
 {

 }
 else
 {
 st = name[(std::string)t_arc->start] = n_nodes++;
 }

 if(name.find((std::string)t_arc->end)!=name.end())
 {
 en = name[(std::string)t_arc->end];
 }
 else
 {
 name[(std::string)t_arc->end] = n_nodes++;
 }
 st = name[(std::string)t_arc->start]-1;
 en = name[(std::string)t_arc->end]-1;
 if(t_arc->low < -1e6)t_arc->low = -1e6;
 if(t_arc->upp > 1e6) t_arc->upp = 1e6;
 if(st == 0)
 printf("%s %s %d %d [%.2lf,%.2lf]\n",t_arc->start,t_arc->end,st,en,
 t_arc->low,t_arc->upp);
 int i;
 for(i = 0;i < t->rqm.size(); i++)
 if(t->rqm[i].start == st && t->rqm[i].end == en)
 break;
 if(i == t->rqm.size())
 t->add_rqm(st,en,t_arc->low,t_arc->upp,1e20,1e20);
 t->set_min(st,en,t_arc->low);
 t->set_max(st,en,t_arc->upp);
 fprintf(fdot,"\"%d\"->\"%d\"[ label = \"[%lf,%lf] %d\"];\n",st,en,t_arc->low,t_arc->upp,0);
 ind[en]++;
 }

 CcpstpNode * t_node;
 for(t_node = prob->node_list; t_node; t_node = t_node->next)
 {
 if(strlen(t_node->parent) == 0)continue;
 int st,en;
 if(name.find((std::string)t_node->parent) != name.end())
 {

 }
 else
 {
 st = name[(std::string)t_node->parent] = n_nodes++;
 }
 if(name.find((std::string)t_node->name) != name.end())
 {

 }
 else
 {
 en = name[(std::string)t_node->name] = n_nodes++;
 }

 st = name[(std::string)t_node->parent]-1;
 en = name[(std::string)t_node->name]-1;
 printf("%s %s %d %d\n",t_node->parent,t_node->name,st,en);

 if(st < t->n_nodes && en < t->n_nodes)
 {t->add_ctg(st,en,t_node->param);
 double lb = t_node->param[2]-10*sqrt(t_node->param[1]);
 double ub = t_node->param[3]+10*sqrt(t_node->param[1]);
 t->set_min(st,en,lb);
 t->set_max(st,en,ub);}

 if(strcmp(prob->objective_fn,t_node->name) == 0)
 strcpy(prob->objective_fn,t_node->parent);
 fprintf(fdot,"\"%d\"->\"%d\"[ label = \"[%.2lf,%.2lf] %d\"];\n",st,en,t_node->param[2],t_node->param[3],1);
 ind[en]++;
 }

 //t->n_nodes = n_nodes -1;
 t->end = name[(std::string)prob->objective_fn]-1;
 t->start = name[(std::string)prob->start_name]-1;
 for(int i = 0; i < t->n_nodes; i++)
 if(i != t->start && ind[i] == 0)
 {
 //	t->add_rqm((t->start), i, 0.0, 1e10,1e20,1e20);
 //	fprintf(fdot,"\"%d\"->\"%d\"[ label = \"[%.2lf,%.2lf] %d\"];\n",(t->start),i,0.0,1e10,2);

 }
 delete[] ind;
 t->n_ctg = t->ctg.size();
 t->n_rqm = t->rqm.size();
 t->p_chance = prob->chance_con;

 printf("%d %d %d\n",t->n_nodes,t->n_ctg,t->n_rqm);


 for(int i = 0; i < t->n_rqm; i++)
 {
 int st = t->rqm[i].start;
 int en = t->rqm[i].end;


 printf("rqm %d %d %.3lf %.3lf\n",st, en, t->min_distance_d(st,en), t->max_distance_d(st,en));
 }

 for(int i = 0; i < t->n_ctg; i++)
 {
 int  st = t->ctg[i].start;
 int en = t->ctg[i].end;
 printf("ctg %d %d ",st,en);

 printf("[%lf,%lf]",t->min_distance_d(st,en),t->max_distance_d(st,en));
 for(int j = 0; j < 5; j++)
 printf("%.3lf ",t->ctg[i].param[j]);
 printf("\n");
 }
 }

 /*cmp sc

 FILE *scf;
 scf = fopen("test.out","r");
 char s[100];
 while(fscanf(scf," %s",s)!=EOF)
 {
 if(strcmp(s,"Event") == 0)
 {
 fscanf(scf," %s",s);
 int id = name[(std::string)s]-1;
 char tmp[100];
 for(int i = 0; i < 3; i++)
 fscanf(scf," %s", tmp);
 double st;
 fscanf(scf, " %lf", &st);
 fprintf(fdot,"\"%d\" [ label = \"%d %lf\"];\n",id,id,st);

 }
 else if(strcmp(s,"Duration") == 0)
 {

 }
 else
 {

 }
 }

 /*cmp sc*/
/*
 fprintf(fdot,"\"start=%d\"->\"end=%d\";\n",t->start,t->end);
 fprintf(fdot,"}\n");
 fclose(fdot);
 printf("%s %s\n",prob->start_name,prob->objective_fn);
 printf("%d %d %d\n",t->n_nodes,t->ctg.size(),t->rqm.size());
 //  exit(0);
 }
 /*
 void Reader::readout(std::string problem_file)
 {
 std::ifstream reader(problem_file.data());
 std::string line;
 getline(reader,line);
 std::vector<std::string> tokens;
 Utils::tokenize(line, tokens, std::string(" \t"));
 Transplan::n_activities = atoi(tokens[0].c_str());
 int nline = atoi(tokens[1].c_str());
 for (int i = 0; i < Transplan::n_activities; ++i)
 {
 Activity act(i,0);
 Transplan::activities.push_back(act);
 }

 for(int i = 0; i < nline; ++i)
 {
 std::getline(reader, line);
 std::vector<std::string> tokens;
 Utils::tokenize(line, tokens, std::string(" :,[]"));
 int start = 0;
 int _from = atoi(tokens[start++].c_str());
 int _to = atoi(tokens[start++].c_str());
 int _type = atoi(tokens[start++].c_str());
 Transplan::activities[_from].successors.push_back(_to);
 Transplan::activities[_to].successors.push_back(_from);

 //  Utils::tokenize(tokens[start++], time_lag, std::string(" ,[]"));
 //	assert(time_lag.size()>=1);
 int minlag = atoi(tokens[start++].c_str());
 int maxlag = atoi(tokens[start++].c_str());
 Transplan::activities[_from].time_lags.push_back(minlag);
 Transplan::activities[_to].time_lags.push_back(-maxlag);
 std::cerr<<_from<<_to<<minlag<<maxlag<<nline<<i<<std::endl;
 }
 }

 */
void Reader::readcctp(std::string problem_file)
{

	//stnu *t = this->getOwner();

	std::ifstream reader(problem_file.data());
	std::string line;
	std::string v_name[10];
	std::string v_value;
	std::vector<std::string> tokens;
	std::vector<std::string> nodes;
	DiscreteVariable tmp_dis_var;

	std::map<std::string, int> mp_name_var;
	std::map<std::string, int> mp_id_var;
	//map from the value of the discrete variable to integer
	std::vector<std::map<std::string, int> > v_mp_value_int;
	std::map<std::string, int> mp_value_int;

	int fl = 0; //1-constraint, 2-discrete vars
	int st, en, ty;
	double lb, ub;
	int g_v_id, g_value;
	t->n_ctg = 0;
	double lbr, ubr;
	bool b_urelax = false, b_lrelax = false;
	bool b_mission; // change the key upper bound which cause ng
	int cntl = 0, cntcl = 0;
	//char dotfile[] = "2.dot";
	vector<pair<int, int> > labels;
	t->utility = 0;
	t->file_name = problem_file;
//	printf("%s\n", dotfile);
	//FILE *fdot = fopen(dotfile, "w");
	//fprintf(fdot,
	//		"digraph G {\n nodesep = .45; \n size = 30;\nlabel=\"CCTP\";\n");
	while (getline(reader, line))
	{
		tokens.clear();
		TRANSPLAN::Utils::tokenize(line, tokens, std::string(" <>\t"));

		if (tokens.size() == 1)
		{
			if (tokens[0] == std::string("CONSTRAINT"))
			{
				fl = 1;
				cntl++;
				g_v_id = -1;
				g_value = -1;
				b_mission = false;
				labels.clear();
				continue;
			}
			if (tokens[0] == std::string("DECISION-VARIABLE"))
			{
				fl = 2;
				cntl++;
				tmp_dis_var = DiscreteVariable();
				mp_value_int.clear();
				continue;
			}
			
			if (tokens[0] == std::string("UNCERTAIN-DURATION-CONSTRAINT"))//extra constraints
			{
			  fl = 3;
			  cntl++;
			  
			}
			if (tokens[0] == std::string("/CONSTRAINT"))
			{
				fl = 0;
				if (b_lrelax == false)
					lbr = INF;
				if (b_urelax == false)
					ubr = INF;
				if (ty == CTGL)
				{

					if (b_lrelax == false)
						lbr = 1e10 + 1;
					if (b_urelax == false)
						ubr = 1e10 + 1;
					CTG tmp_ctg = CTG(st, en, lb, ub, lbr, ubr);
					char buffer[10];
					sprintf(buffer, "%d", st);
					tmp_ctg.start_name = std::string(buffer);
					sprintf(buffer, "%d", en);
					tmp_ctg.end_name = std::string(buffer);
					tmp_ctg.v_labels = labels;
					t->add_ctg(tmp_ctg);
//					t->add_ctg(st, en, lb, ub, lbr, ubr);
				}
				else
				{
					RQM tmp_rqm = RQM(st, en, lb, ub, lbr, ubr);
					char buffer[10];
					sprintf(buffer, "%d", st);
					tmp_rqm.start_name = std::string(buffer);
					sprintf(buffer, "%d", en);
					tmp_rqm.end_name = std::string(buffer);
					tmp_rqm.b_mission = b_mission;
					tmp_rqm.v_labels = labels;
					b_mission = false;
					t->add_rqm(tmp_rqm);
					//				t->add_rqm(st, en, lb, ub, lbr, ubr);
				}

				//	std::cout << st << "->" << en << " [" << lb << ", " << ub
				//			<< "] " << ty << " " << lbr << " " << ubr << std::endl;
				//	fprintf(fdot,
				//			"\"%d\"->\"%d\"[ label = \"[%.2lf,%.2lf] %d\"];\n", st,
				///			en, lb, ub, ty);
				lbr = ubr = 0;
				b_urelax = false;
				b_lrelax = false;

				if (g_v_id >= 0)
				{
					for (int i = 0; i < t->v_dis_var[g_v_id].v_values.size();
							i++)
					{
						if (t->v_dis_var[g_v_id].v_values[i] != g_value)
						{
							int id_link = 0;
							if (ty == CTGL)
								id_link = -t->ctg.size();
							else
								id_link = t->rqm.size();
							t->v_dis_var[g_v_id].vv_id_mute_edge[i].push_back(
									id_link);
						}
					}
				}
			}

			if (tokens[0] == std::string("/DECISION-VARIABLE"))
			{
				fl = 0;
				tmp_dis_var.n_value = tmp_dis_var.v_values.size();
				if (tmp_dis_var.n_value > 1) //if more than one choices exist
				{
					t->v_dis_var.push_back(tmp_dis_var);
					v_mp_value_int.push_back(mp_value_int);
				}
				else
				{
					if (tmp_dis_var.n_value == 1)
						t->utility += tmp_dis_var.v_utilities[0];
					mp_id_var[tmp_dis_var.ID] = 0;
					mp_name_var[tmp_dis_var.name] = 0;
				}
			}


		}

		if (fl == 1) //constraint
		{
			if (tokens[0] == std::string("START"))
			{
				int i = 0;
				for (i = 0; i < nodes.size(); i++)
					if (nodes[i] == tokens[1])
						break;
				if (i < nodes.size())
					st = i;
				else
					nodes.push_back(tokens[1]), st = i;
			}
			if (tokens[0] == std::string("END"))
			{
				int i = 0;
				for (i = 0; i < nodes.size(); ++i)
					if (nodes[i] == tokens[1])
						break;
				if (i < nodes.size())
					en = i;
				else
				{
					nodes.push_back(tokens[1]), en = i;
				}
			}
			else if (tokens[0] == std::string("LOWERBOUND"))
			{
				lb = atof(tokens[1].c_str());
			}
			else if (tokens[0] == std::string("UPPERBOUND"))
			{
				ub = atof(tokens[1].c_str());
			}
			else if (tokens[0] == std::string("TYPE"))
			{
				if (tokens[1].find("Uncontrollable") != tokens[1].npos)
					ty = CTGL;
				else
					ty = RQML;
			}
			else if (tokens[0].find("LB-RELAX-COST-RATIO") != tokens[1].npos)
			{
				lbr = atof(tokens[1].c_str());
			}
			else if (tokens[0].find("UB-RELAX-COST-RATIO") != tokens[1].npos)
			{
				ubr = atof(tokens[1].c_str());
			}
			else if (tokens[0].find("UBRELAXABLE") != tokens[1].npos)
			{
				b_urelax = true;
			}
			else if (tokens[0].find("LBRELAXABLE") != tokens[1].npos)
			{
				b_lrelax = true;
			}
			else if (tokens[0] == std::string("GUARD-VARIABLE"))
			{
//				if(mp_name_var.find(tokens[1]) != mp_name_var.end())
//					g_v_id = mp_name_var[tokens[1]] - 1;
//				else
				if (mp_id_var.find(tokens[1]) != mp_id_var.end())
					g_v_id = mp_id_var[tokens[1]] - 1;
				else
					g_v_id = -1;
			}
			else if (tokens[0] == std::string("GUARD-VALUE"))
			{
				//std::vector<std::string> elements;
				//Utils::tokenize(tokens[1], elements, std::string("-"));
				//the format of value-name is xx-yyyyyyy, and xx is a unique int.
				if (g_v_id != -1)
				{
					g_value = v_mp_value_int[g_v_id][tokens[1]] - 1;
					labels.push_back(make_pair(g_v_id, g_value));
				}
			}
			else if (tokens[0] == std::string("CONTROLLABLE"))
			{
			  if (tokens[1].find("false") != tokens[1].npos)
					ty = CTGL;

			}

			if (tokens[0] == std::string("NAME"))
			{
				if (tokens[1].find("Mission") != -1)
				{
					b_mission = true;
				}
			}

		}

		if (fl == 2) //discrete variable
		{
			if (tokens[0] == std::string("ID"))
			{
				if (tmp_dis_var.ID.length())
					continue;
				tmp_dis_var.ID = tokens[1];
				mp_id_var[tokens[1]] = t->v_dis_var.size() + 1;
			}
			if (tokens[0] == std::string("DECISION-NAME"))
			{
				tmp_dis_var.name = tokens[1];
				mp_name_var[tokens[1]] = t->v_dis_var.size() + 1;
			}
			if (tokens[0] == std::string("VALUE-NAME"))
			{
				//std::vector<std::string> elements;
				//Utils::tokenize(tokens[1], elements, std::string("-"));
				//the format of value-name is xx-yyyyyyy, and xx is a unique int.

				if (mp_value_int.find(tokens[1]) == mp_value_int.end())
					mp_value_int[tokens[1]] = mp_value_int.size();
				int value_id;
				value_id = mp_value_int[tokens[1]] - 1;
				tmp_dis_var.v_values.push_back(value_id);
				tmp_dis_var.vv_id_mute_edge.push_back(std::vector<int>());
			}

			if (tokens[0] == std::string("VALUE-UTILITY"))
			{
				tmp_dis_var.v_utilities.push_back(atof(tokens[1].c_str()));
			}

		}
		//if (tokens[0] == std::string("VALUE-UTILITY")) {
		//t->utility += atof(tokens[1].c_str());

	}
	//fprintf(fdot, "}\n");
//	fclose(fdot);
	//printf("%d %d %d %d\n", nodes.size(), cntl, t->ctg.size(), t->rqm.size());
//  exit(0);
	t->n_ctg = t->ctg.size();
	t->n_rqm = t->rqm.size();
	t->n_vars = t->v_dis_var.size();
	t->n_nodes = nodes.size();

	//Transplan::n_activities = nodes.size();
	if (t->s_info.id_debug == DEBUG_T)
	{
		FILE *fp;
		fp = fopen("nodes.txt", "w");
		for (int i = 0; i < nodes.size(); i++)
		{
			fprintf(fp, "%d	%s\n", i, nodes[i].c_str());
		}
		fclose(fp);
	}

	if (t != NULL)
	{
		t->init(nodes.size());
		t->n_ctg = t->ctg.size();
		t->n_rqm = t->rqm.size();
		t->n_nodes = nodes.size();

		/*  for(int i = 0; i < t->n_ctg; i++)
		 {
		 int st = t->ctg[i].start;
		 int en = t->ctg[i].end;
		 t->set_min(st,en,t->ctg[i].lb);
		 t->set_max(st,en,t->ctg[i].ub);
		 }

		 for(int i = 0; i < t->n_rqm; i++)
		 {
		 int st = t->rqm[i].start;
		 int en = t->rqm[i].end;

		 if(t->id_obj != O_RELAXCOST)
		 {
		 t->set_min(st,en,t->rqm[i].lb);
		 t->set_max(st,en,t->rqm[i].ub);
		 }
		 else
		 {
		 t->set_min(st,en,0);
		 t->set_max(st,en,1e6);
		 }
		 }
		 t->compute_minimal();
		 */
	}

}
/*

 void Reader::readpos(std::string problem_file)
 {
 //read from file
 std::ifstream reader(problem_file.data());
 std::string line;
 getline(reader,line);

 stnu* t=this->getOwner();
 std::vector<std::string> tokens;
 TRANSPLAN::Utils::tokenize(line, tokens, std::string(" \t"));
 Transplan::n_activities = atoi(tokens[0].c_str());

 Transplan::FLDT_b = atof(tokens[1].c_str());

 Transplan::activities.clear();
 for (int i = 0; i < Transplan::n_activities; ++i)
 {
 Activity act(i,0);
 Transplan::activities.push_back(act);
 Transplan::activities[i].successors.clear();
 Transplan::activities[i].time_lags.clear();
 }

 //read project precedences
 for (int task = 0 ; task < Transplan::n_activities; ++task)
 {
 std::getline(reader, line);
 std::vector<std::string> tokens;
 Utils::tokenize(line, tokens, std::string(" \t"));

 int start = 1;
 int dur = atoi(tokens[start++].c_str());
 Transplan::activities[task].duration = dur;
 int est = atoi(tokens[start++].c_str());

 int n_succ = atoi(tokens[start++].c_str());
 for (int s = 0; s < n_succ; ++s)
 {
 int suc = atoi(tokens[start++].c_str());
 Transplan::activities[task].successors.push_back(suc);
 }

 for (int t = 0; t < n_succ; ++t)
 {
 std::vector<std::string> time_lag;
 Utils::tokenize(tokens[start++], time_lag, std::string(" []"));
 assert(time_lag.size()>=1);
 int lag = atoi(time_lag[0].c_str());
 Transplan::activities[task].time_lags.push_back(lag);
 }
 }

 Transplan::n_nodes = (Transplan::n_activities-1)*2;

 //transfer the input to the STNU
 if(t!=NULL)
 {
 t->init(Transplan::n_nodes);
 t->n_ctg = Transplan::n_activities - 2;
 t->n_nodes = Transplan::n_nodes;

 for(int i = 1;i < Transplan::n_activities-1; i++)
 t->add_ctg(i*2-1,2*i);
 Transplan::min_time_lags.clear();
 Transplan:: max_time_lags.clear();

 for (int i = 0; i < Transplan::n_activities; ++i)
 {
 Transplan::min_time_lags.push_back(DblVector());
 Transplan::max_time_lags.push_back(DblVector());

 for (int j = 0; j < Transplan::n_activities; ++j)
 {
 Transplan::min_time_lags[i].push_back(-1*Transplan::DEFAULT);
 Transplan::max_time_lags[i].push_back(Transplan::DEFAULT);
 }
 Transplan:: activities[i].ub = Transplan::DEFAULT;
 }

 for (int i = 0; i < Transplan::n_activities; ++i)
 {
 int act_index = i;
 //time lag posting
 IntVector succ = Transplan::activities[act_index].successors;
 DblVector lags = Transplan::activities[act_index].time_lags;

 for (int j = 0; j < succ.size(); ++j)
 {
 int succ_index = succ[j];
 int time_lag = lags[j];
 RQM t_rqm(0,0);

 double lb = Transplan::min_time_lags[act_index][succ_index];
 double ub = Transplan::max_time_lags[succ_index][act_index];
 if (time_lag >= 0 && lb <= time_lag)
 {
 Transplan::min_time_lags[act_index][succ_index] = time_lag;
 Transplan::activities[act_index].min_time_fw_const[succ_index] = time_lag;
 Transplan::activities[succ_index].min_time_bw_const[act_index] = time_lag;
 t_rqm.start = act_index*2;
 t_rqm.end = succ_index*2-1;
 if(t_rqm.end<0)t_rqm.end++;
 if(t_rqm.start > t->n_nodes) t_rqm.start--;
 }
 else if(time_lag < 0 && ub >= - time_lag)
 {
 Transplan::max_time_lags[succ_index][act_index] = -1*time_lag;
 Transplan::activities[succ_index].max_time_fw_const[act_index] = -1*time_lag;
 Transplan::activities[act_index].max_time_bw_const[succ_index] = -1*time_lag;
 t_rqm.start = succ_index*2;
 t_rqm.end = act_index*2-1;
 if(t_rqm.end<0)t_rqm.end++;
 if(t_rqm.start > t->n_nodes) t_rqm.start--;
 }
 int k;

 for(k = 0; k < (t->rqm).size(); k++)
 if((t->rqm[k]).start == t_rqm.start
 && (t->rqm[k]).end == t_rqm.end)
 {
 break;
 }
 if(k == (t->rqm).size())
 {
 (t->rqm).push_back(t_rqm);
 printf("rqm %d %d\n",t_rqm.start,t_rqm.end);
 }
 }
 }
 t->n_rqm = (t->rqm).size();

 for(int i = 0; i < Transplan::n_nodes; ++i)
 {
 if (i%2 == 0 && i!=0)
 {
 t->set_min(i-1,i,Transplan::activities[i/2].duration);
 }
 }
 t->b_ctg = new bool*[Transplan::n_nodes];
 t->s_map = new bool*[Transplan::n_nodes];

 for(int i = 0; i < Transplan::n_nodes; i++)
 {
 t->b_ctg[i] = new bool[Transplan::n_nodes];
 t->s_map[i] = new bool[Transplan::n_nodes];
 for(int j = 0;j < Transplan::n_nodes; j++)
 {	t->b_ctg[i][j] = false,
 t->s_map[i][j] = false;

 if(i%2 == 0 && j%2 == 1 && j+1 != i)//end2start
 {

 if(Transplan::min_time_lags[i/2][(j+1)/2] != -Transplan::DEFAULT)
 t->set_min(i,(j), Transplan::min_time_lags[i/2][(j+1)/2]
 - Transplan::activities[(i)/2].duration);


 if(Transplan::max_time_lags[i/2][(j+1)/2] != Transplan::DEFAULT)
 t->set_max(i,(j),Transplan::max_time_lags[i/2][(j+1)/2]
 - Transplan::activities[(i)/2].duration);
 }
 }
 }

 printf("%d %d\n",t->ctg.size(),t->rqm.size());

 t->compute_minimal();
 for(int i = 0; i < t->n_rqm; i++)
 {
 int st = t->rqm[i].start;
 int en = t->rqm[i].end;

 t->rqm[i].lb = t->min_distance_d(st,en);
 t->rqm[i].ub = t->max_distance_d(st,en);
 }

 for(int i = 0; i < t->n_ctg; i++)
 {
 int st = t->ctg[i].start;
 int en = t->ctg[i].end;

 t->ctg[i].lb = t->min_distance_d(st,en);
 t->ctg[i].ub = 1e6;;
 }
 }

 }
 */
void Reader::readgrid(std::string problem_file)
{
	std::ifstream reader(problem_file.data());
	std::string line;

	while (std::getline(reader, line))
	{
		std::vector<std::string> tokens;
		std::cerr << line << std::endl;
		TRANSPLAN::Utils::tokenize(line, tokens, std::string(" \t"));
		if (tokens.size() < 1)
			continue;
		if (tokens[0] == "p")
		{
			t->n_nodes/* = Transplan::n_nodes */= atoi(tokens[2].c_str());
			t->n_rqm = atoi(tokens[3].c_str());
			t->n_ctg/* = Transplan::n_activities*/= atoi(tokens[4].c_str());
			t->init(t->n_nodes);
		}
		else if (tokens[0] == "a")
		{
			int st, en;
			double lb, ub;
			st = atoi(tokens[1].c_str()) - 1;
			en = atoi(tokens[2].c_str()) - 1;
			lb = atof(tokens[3].c_str());
			ub = atof(tokens[4].c_str());
			t->set_min(st, en, lb);
			t->set_max(st, en, ub);
			//	 t->compute_minimal();
			if (ub > 0)
				t->add_rqm(st, en, lb, ub, 1, 1);
			else
				t->add_rqm(en, st, -ub, -lb, 1, 1);
		}
		else if (tokens[0] == "g")
		{
			int st, en, lb, ub;
			st = atoi(tokens[1].c_str()) - 1;
			en = atoi(tokens[2].c_str()) - 1;
			t->add_ctg(st, en, t->min_distance_d(st, en),
					t->max_distance_d(st, en), 1e20, 1e20);
			for (std::vector<LINK>::iterator it = t->rqm.begin();
					it != t->rqm.end(); it++)
			{
				if ((it->start) == st && (it->end) == en)
				{
					t->rqm.erase(it);
					break;
				}
				if ((it->start) == en && (it->end) == st)
				{
					t->rqm.erase(it);
					break;
				}
			}
			t->n_rqm = t->rqm.size();
		}
	}

}

void Reader::PrintDot(std::string problem_file)
{
	FILE *fdot;
	fdot = fopen(problem_file.c_str(), "w");
	fprintf(fdot, "digraph G {\n rankdir = LR;\n nodesep = .45; \n "
			"size = 30;\nlabel=\"result\";\n");

	for (int i = 0; i < t->n_ctg; i++)
	{
		CTG tmp_edge;
		tmp_edge = t->ctg[i];

		fprintf(fdot,
				"\"%d\"->\"%d\"[style=dotted label = \"[%.2lf,%.2lf] %d\"];\n",
				tmp_edge.start, tmp_edge.end, tmp_edge.lb, tmp_edge.ub, 1);
	}

	for (int i = 0; i < t->n_rqm; i++)
	{
		RQM tmp_edge;
		tmp_edge = t->rqm[i];
		fprintf(fdot, "\"%d\"->\"%d\"[label = \"[%.2lf,%.2lf] %d\"];\n",
				tmp_edge.start, tmp_edge.end, tmp_edge.lb, tmp_edge.ub, 0);

	}
	fprintf(fdot, "}\n");
	fclose(fdot);
}

int Reader::getType(std::string x)
{
	if (x[0] == 'r' || x[0] == 'R')
		return 0;
	if (x[0] == 'c' || x[0] == 'C')
		return 1;
	int res = atoi(x.c_str());
	//  if(res == 'r')return 0;
	// if(res == 'c')return 1;
	return res;
}

std::string Reader::getName(std::string input_file)
{
	std::vector<std::string> tokens;
	TRANSPLAN::Utils::tokenize(input_file, tokens, std::string("./"));

	return tokens[tokens.size() - 2];

}
}
