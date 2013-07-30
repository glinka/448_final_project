/*
 * Dynamic_Voter.cpp
 *
 *  Created on: Apr 18, 2011
 *      Author: Bill
 */

#include "Dynamic_Voter.h"
#include "global_variables.h"
//#include <engine.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <queue>
#include <math.h>
#include <algorithm>

using namespace std;

//Generate a regular random graph of the given degree
Dynamic_Voter::Dynamic_Voter(long int number_of_nodes, int degree, int number_of_opinions)
: max_edge_ID(0), population(number_of_nodes), sites(number_of_opinions)
{
	bool need_new_neigh,need_new_edge;
	long int i,j,k,i1,i2,p,q,r,s,u,v,e1;
	long int number_of_edges=number_of_nodes*degree/2;
	vector<Node>::iterator person_it;
	vector<Edge>::iterator new_edge;
	vector<long int> nodequeue(number_of_nodes,0); // all the nodes of degree less than d
	edges.reserve(number_of_edges);
	for(i=0;i<number_of_opinions;i++)
		sites[i].reserve(number_of_nodes);
	edge_boundary.reserve(number_of_edges);
	edge_boundary01.reserve(number_of_edges);

	for(person_it=population.begin(),j=0; person_it!=population.end(); person_it++,j++)
	{
		person_it->ID=j;
		person_it->state=0;
		person_it->myself=person_it;
	}

	for (i=0; i<number_of_nodes; i++)
		nodequeue[i] = i;
	if ((degree % 2) != 0)
		if ((number_of_nodes % 2) != 0)
		{
			cout<<"if degree is odd then number_of_nodes must be even"<<endl;
			return;
		}
	if (number_of_nodes <= degree)
	{
		cout<<"number_of_nodes should be greater than degree"<<endl;
		return;
	}

	while (nodequeue.empty()==false)
	{
		u=-1;
		v=-1;
		if (nodequeue.size()==1)
		{
			// node r has less than the required degree,
			// find two adjacent nodes p and q non-adjacent to r.
			r=nodequeue[0];
			need_new_edge=true;
			while (need_new_edge)
			{
				e1 = random_number.integer(edges.size());
				new_edge = edges[e1].myself;
				need_new_edge = (new_edge->person1->ID==r || new_edge->person2->ID==r);
				if (need_new_edge==false)
				{
					p = new_edge->person1->ID;
					q = new_edge->person2->ID;
					need_new_edge = is_neighbor(p, r) || is_neighbor(q, r);
				}
			}
			delete_edge(new_edge);
			add_edge(p, r);
			add_edge(q, r);
			if ((int)population[r].edge_list.size()>=degree)
			{
				nodequeue.pop_back();
			}
		}

		// find two non-adjacent nodes each has less than required degree
		if ((int)nodequeue.size()>degree)
		{
			need_new_neigh=true;
			while (need_new_neigh)
			{
				i1=random_number.integer(nodequeue.size());
				i2 = i1;
				while (i2 == i1)
				{
					i2=random_number.integer(nodequeue.size());
				}
				need_new_neigh=is_neighbor(nodequeue[i1],nodequeue[i2]);
			}
			add_edge(nodequeue[i1],nodequeue[i2]); // add this edge
			if ((int)population[nodequeue[i1]].edge_list.size()>=degree && (int)population[nodequeue[i2]].edge_list.size()>=degree)
			{
				if (i1==(int)nodequeue.size()-1)
				{
					nodequeue.pop_back();
					nodequeue[i2]=nodequeue.back();
					nodequeue.pop_back();
				}
				else
				{
					nodequeue[i2]=nodequeue.back();
					nodequeue.pop_back();
					nodequeue[i1]=nodequeue.back();
					nodequeue.pop_back();
				}
			}
			else if((int)population[nodequeue[i1]].edge_list.size()>=degree)
			{
				nodequeue[i1]=nodequeue.back();
				nodequeue.pop_back();
			}
			else if ((int)population[nodequeue[i2]].edge_list.size()>=degree)
			{
				nodequeue[i2]=nodequeue.back();
				nodequeue.pop_back();
			}
		}
		else
		{
			for (i=0; i<(int)nodequeue.size(); i++)  // randomly permute node queue
			{
				j = (int)(i + random_number.real() * (nodequeue.size() - i));
				k = nodequeue[i];
				nodequeue[i] = nodequeue[j];
				nodequeue[j] = k;
			}
			for (i1=0; i1<(int)nodequeue.size()-1; i1++)
			{
				for (i2=i1+1; i2<(int)nodequeue.size(); i2++)
				{
					if (!is_neighbor(nodequeue[i1],nodequeue[i2]))
					{
						add_edge(nodequeue[i1],nodequeue[i2]);  // add this edge
						if ((int)population[nodequeue[i1]].edge_list.size()>=degree && (int)population[nodequeue[i2]].edge_list.size()>=degree)
						{
							if (i1==(int)nodequeue.size()-1)
							{
								nodequeue.pop_back();
								nodequeue[i2]=nodequeue.back();
								nodequeue.pop_back();
							}
							else
							{
								nodequeue[i2]=nodequeue.back();
								nodequeue.pop_back();
								nodequeue[i1]=nodequeue.back();
								nodequeue.pop_back();
							}
						}
						else if((int)population[nodequeue[i1]].edge_list.size()>=degree)
						{
							nodequeue[i1]=nodequeue.back();
							nodequeue.pop_back();
						}
						else if ((int)population[nodequeue[i2]].edge_list.size()>=degree)
						{
							nodequeue[i2]=nodequeue.back();
							nodequeue.pop_back();
						}
						v=-1;
						break;
					}
					else
					{
						u=i1;
						v=i2;
						r = nodequeue[i1];
						s = nodequeue[i2];
					}
				}
				if (v==-1)
					break;
			}
		}
		if (v!=-1)
		{
			// nodes r and s of less than required degree, find two
			// adjacent nodes p & q such that (p,r) & (q,s) are not edges.
			need_new_edge=true;
			while (need_new_edge)
			{
				e1 = random_number.integer(edges.size());
				new_edge = edges[e1].myself;
				need_new_edge = (new_edge->person1->ID==r || new_edge->person2->ID==r || new_edge->person1->ID==s || new_edge->person2->ID==s);
				if (need_new_edge==false)
				{
					p = new_edge->person1->ID;
					q = new_edge->person2->ID;
					need_new_edge = is_neighbor(p, r) || is_neighbor(q, s);
					if (need_new_edge==true)
					{
						q = new_edge->person1->ID;
						p = new_edge->person2->ID;
						need_new_edge = is_neighbor(p, r) || is_neighbor(q, s);
					}
				}
			}
			delete_edge(new_edge);
			add_edge(p, r);
			if((int)population[r].edge_list.size()>=degree)
			{
				nodequeue[u]=nodequeue.back();
				nodequeue.pop_back();
			}
			add_edge(q, s);
			if((int)population[s].edge_list.size()>=degree)
			{
				nodequeue[v]=nodequeue.back();
				nodequeue.pop_back();
			}
		}
	}
	return;
}

//Generate a ER random graph with given number of nodes and edges
Dynamic_Voter::Dynamic_Voter(long int number_of_nodes, long int number_of_edges, int number_of_opinions)
: max_edge_ID(0), population(number_of_nodes),sites(number_of_opinions)
{
	bool need_new_neigh;
	long int j,i1,i2;
	int i;
	vector<Node>::iterator person_it;
	degree=2.0*number_of_edges/number_of_nodes;
	edges.reserve(number_of_edges);
	for(i=0;i<number_of_opinions;i++)
		sites[i].reserve(number_of_nodes);
	edge_boundary.reserve(number_of_edges);
	edge_boundary01.reserve(number_of_edges);
//	igraph = gsl_matrix_calloc(number_of_nodes, number_of_nodes);
//	igraph_vector_t edgelist;
//	igraph_vector_init(&edgelist,2*number_of_edges);

	for(person_it=population.begin(),j=0; person_it!=population.end(); person_it++,j++)
	{
		person_it->ID=j;
		person_it->state=0;
		person_it->myself=person_it;
	}

	for(j=0; j<number_of_edges; j++)
	{
		// choose two nodes
		need_new_neigh=true;
		while (need_new_neigh)
		{
			i1=random_number.integer(number_of_nodes);
			i2 = i1;
			while (i2 == i1)
			{
				i2=random_number.integer(number_of_nodes);
			}
			need_new_neigh=is_neighbor(i1,i2);
		}
		add_edge(i1, i2);  // add this edge
		//		VECTOR(edgelist)[2*j]=i1;
		//		VECTOR(edgelist)[2*j+1]=i2;
//		gsl_matrix_set(igraph,i1,i2,1);
//		gsl_matrix_set(igraph,i2,i1,1);
	}
	//	igraph_create(&igraph, &edgelist, number_of_nodes, IGRAPH_UNDIRECTED);
	//	igraph_vector_destroy(&edgelist);

	return;
}


Dynamic_Voter::~Dynamic_Voter()
{
	// TODO Auto-generated destructor stub
}

bool Dynamic_Voter::is_neighbor(long int i1, long int i2)
{
	list<vector<Edge>::iterator>::iterator neigh_edge_it;

	for (neigh_edge_it=population[i1].edge_list.begin(); neigh_edge_it!=population[i1].edge_list.end(); neigh_edge_it++)
	{
		if ((*neigh_edge_it)->person1->ID == i2)
			return true;
		if ((*neigh_edge_it)->person2->ID == i2)
			return true;
	}
	return false;
}

vector<Edge>::iterator Dynamic_Voter::add_edge(long int i1, long int i2)
{
	Edge new_edge;
	vector<Edge>::iterator edge_it;

	new_edge.ID = max_edge_ID++;
	edge_it = edges.insert(edges.end(), new_edge);
	edge_it->person1 = population[i1].myself;
	edge_it->person2 = population[i2].myself;
	edge_it->link1 = population[i1].edge_list.insert(population[i1].edge_list.begin(), edge_it);
	edge_it->link2 = population[i2].edge_list.insert(population[i2].edge_list.begin(), edge_it);
	edge_it->myself=edge_it;
	if (population[i1].state != population[i2].state)
	{
		edge_it->state = 1;
		edge_it->boundary_place = edge_boundary.insert(edge_boundary.end(),	edge_it);
		if (population[i1].state==0 && population[i2].state==1 || population[i1].state==1 && population[i2].state==0)
		{
			edge_it->state01=1;
			edge_it->boundary01_place = edge_boundary01.insert(edge_boundary01.end(),edge_it);
		}
		else
		{
			edge_it->state01=0;
			edge_it->boundary01_place = NULL_EDGE_BOUNDARY_ITERATOR;
		}
	}
	else
	{
		edge_it->state=0;
		edge_it->boundary_place = NULL_EDGE_BOUNDARY_ITERATOR;
		edge_it->state01=0;
		edge_it->boundary01_place = NULL_EDGE_BOUNDARY_ITERATOR;
	}

	return edge_it;
}

// swap the edge with the last edge in the edgelist and delete it
void Dynamic_Voter::delete_edge(vector<Edge>::iterator edge_it)
{
	if (edge_it->state == 1)
	{
		swap_delete(edge_it, edge_boundary);
		if (edge_it->state01 == 1)
			swap_delete01(edge_it, edge_boundary01);
	}
	edge_it->person1->edge_list.erase(edge_it->link1);
	edge_it->person2->edge_list.erase(edge_it->link2);
	if (edge_it->ID != (edges.back()).ID && edges.size()!=1)
	{
		*edge_it=edges.back();
		*(edge_it->link1) = edge_it;
		*(edge_it->link2) = edge_it;
		edge_it->myself = edge_it;
		if (edge_it->state==1)
		{
			*(edge_it->boundary_place) = edge_it;
			if (edge_it->state01==1)
				*(edge_it->boundary01_place) = edge_it;
		}
	}
	edges.pop_back();

	return;
}

// delete two edges at the same time
void Dynamic_Voter::delete_edge(vector<Edge>::iterator edge1_it, vector<Edge>::iterator edge2_it)
{
	if (edge1_it->ID == edges.back().ID)
	{
		delete_edge(edge1_it);
		delete_edge(edge2_it);
	}
	else
	{
		delete_edge(edge2_it);
		delete_edge(edge1_it);
	}
}

// swap the node pointer with the last node pointer in the list and remove it
bool Dynamic_Voter::swap_delete(vector<Node>::iterator person_it, vector<vector<Node>::iterator> & site)
{
	if(person_it->sites_place == NULL_SITES_ITERATOR)
		return false;
	if (site.size()==1 || person_it == site.back())
		;
	else
	{
		*(person_it->sites_place) = site.back();
		(*site.back()).sites_place = person_it->sites_place;
	}
	site.pop_back();
	person_it->sites_place = NULL_SITES_ITERATOR;
	return true;
}

// swap the edge pointer with the last edge pointer in the list and remove it
bool Dynamic_Voter::swap_delete(vector<Edge>::iterator edge_it, vector<vector<Edge>::iterator> & site)
{
	if (edge_it->boundary_place == NULL_EDGE_BOUNDARY_ITERATOR)
		return false;
	if (site.size()==1 || edge_it == site.back())
		;
	else
	{
		*(edge_it->boundary_place) = site.back();
		(*site.back()).boundary_place = edge_it->boundary_place;
	}
	site.pop_back();
	edge_it->boundary_place=NULL_EDGE_BOUNDARY_ITERATOR;
	return true;
}

// swap the edge pointer with the last edge pointer in the list and remove it, restricted to 0-1 edges
bool Dynamic_Voter::swap_delete01(vector<Edge>::iterator edge_it, vector<vector<Edge>::iterator> & site)
{
	if (edge_it->boundary01_place == NULL_EDGE_BOUNDARY_ITERATOR)
		return false;
	if (site.size()==1 || edge_it == site.back())
		;
	else
	{
		*(edge_it->boundary01_place) = site.back();
		(*site.back()).boundary01_place = edge_it->boundary01_place;
	}
	site.pop_back();
	edge_it->boundary01_place=NULL_EDGE_BOUNDARY_ITERATOR;
	return true;
}

void Dynamic_Voter::assign_states(vector<float> &initial_density)
{
	vector<Node>::iterator person_it;
	vector<Edge>::iterator edge_it;
	double u;
	int i;
	vector<float> cmf;
	cmf.reserve(initial_density.size()+1);
	cmf.push_back(0);

	for (i=0;i<(int)initial_density.size();i++)
	{
		cmf.push_back(initial_density[i]+cmf[i]);
	}

	for (person_it = population.begin(); person_it != population.end(); person_it++)
	{
		u=random_number.real();
		for (i=0;i<(int)cmf.size()-1;i++)
			if (u > cmf[i] && u<=cmf[i+1])
			{
				person_it->state = i;
				person_it->sites_place = (sites[i]).insert((sites[i]).end(),person_it);
				break;
			}
	}

	for (edge_it = edges.begin(); edge_it != edges.end(); edge_it++)
	{
		if (edge_it->person1->state != edge_it->person2->state)
		{
			edge_it->state=1;
			edge_it->boundary_place=edge_boundary.insert(edge_boundary.end(),edge_it);
			if (edge_it->person1->state==0 && edge_it->person2->state==1 || edge_it->person1->state==1 && edge_it->person2->state==0)
			{
				edge_it->state01=1;
				edge_it->boundary01_place=edge_boundary01.insert(edge_boundary01.end(),edge_it);
			}
		}
	}

	return;
}

double Dynamic_Voter::simulate(float alpha, float beta, int mode)
{
	long int e1,e2;
	long int N101=0,N010=0,N100=0,N110=0;
	double step=0, real_time=0;
	int counter=0,isrewire;
	int dt,i,j,k,l;
	vector<Edge>::iterator edge_it,edge_it2;
	vector<Node>::iterator person1_it, person2_it, person3_it, person4_it;
	vector<long int>::iterator comp_it;
	vector<long int> degreedist,degreedist1;
	stringstream filename;
	bool need_new_edge = true;
	dt=(int)population.size();

	// create file to save the frequencies of opinions throughout the process
	filename<<"mode"<<mode<<"n"<<population.size()<<"k"<<(edges.size())*2/(population.size())<<"b"<< setiosflags(ios::fixed) << setprecision(2)<<beta<<"u";
	for (j=0; j<(int)sites.size()-1; j++)
		filename << (float)sites[j].size()/population.size()<<"_";
	filename<< (float)sites[j].size()/population.size()<<".process";
	ifstream testfile(filename.str().c_str());
	i=1;
	while (testfile.is_open())
	{
		i++;
		filename.str("");
		filename<<"mode"<<mode<<"n"<<population.size()<<"k"<<(edges.size())*2/(population.size())<<"b"<<beta<<"u";
		for (j=0; j<(int)sites.size()-1; j++)
			filename <<(float)sites[j].size()/population.size()<<"_";
		filename<< (float)sites[j].size()/population.size()<<"("<< i <<").process";
		testfile.close();
		testfile.open(filename.str().c_str());
	}
	ofstream pFile_process(filename.str().c_str());

//	create file to save the component sizes throughout the process
//	filename.str("");
//	filename<<"comp_size_mode"<<mode<<"n"<<population.size()<<"k"<<(edges.size())*2/(population.size())<<"a"<<alpha<<"u";
//	for (j=0; j<(int)sites.size()-1; j++)
//		filename << setiosflags(ios::fixed) << setprecision(1)<<(float)sites[j].size()/population.size()<<"_";
//	filename<<(float)sites[j].size()/population.size();
//	if (i>1)
//		filename<<"("<< i <<")";
//	ofstream pFile_comp(filename.str().c_str());

//	create file to save the component sizes of 1-nodes throughout the process
//	filename.str("");
//	filename<<"comp_size_1_mode"<<mode<<"n"<<population.size()<<"k"<<(edges.size())*2/(population.size())<<"a"<<alpha<<"u";
//	for (j=0; j<(int)sites.size()-1; j++)
//		filename << setiosflags(ios::fixed) << setprecision(1)<<(float)sites[j].size()/population.size()<<"_";
//	filename<<(float)sites[j].size()/population.size();
//	if (i>1)
//		filename<<"("<< i <<")";
//	ofstream pFile_comp_1(filename.str().c_str());

//	create file to save the degree distributions throughout the process
//	filename.str("");
//	filename<<"mode"<<mode<<"n"<<population.size()<<"k"<<(edges.size())*2/(population.size())<<"a"<< setiosflags(ios::fixed) << setprecision(2)<<alpha<<"u";
//	for (j=0; j<(int)sites.size()-1; j++)
//		filename << (float)sites[j].size()/population.size()<<"_";
//	filename<< (float)sites[j].size()/population.size();
//	if (i>1)
//		filename<<"("<< i <<")";
//	filename<<".ddist";
//	ofstream pFile_degree(filename.str().c_str());

//	create file to save the degree distributions of 1-nodes throughout the process
//	filename.str("");
//	filename<<"degree_dist_1_mode"<<mode<<"n"<<population.size()<<"k"<<(edges.size())*2/(population.size())<<"a"<<alpha<<"u";
//	for (j=0; j<(int)sites.size()-1; j++)
//		filename << setiosflags(ios::fixed) << setprecision(1)<<(float)sites[j].size()/population.size()<<"_";
//	filename<<(float)sites[j].size()/population.size();
//	if (i>1)
//		filename<<"("<< i <<")";
//	ofstream pFile_degree1(filename.str().c_str());

//	fill in the initial data
	pFile_process<<"step "<<"time ";
	for (j=0; j< (int)sites.size(); j++)
		pFile_process<<"N"<<j<<" ";
	pFile_process<<"N_neq ";
	for (i=0;i<(int)sites.size();i++)
		for (j=i;j<(int)sites.size();j++)
			pFile_process<<"N"<<i<<j<<" ";
	for (i=0;i<(int)sites.size();i++)
		for (j=0;j<(int)sites.size();j++)
			for (k=i;k<(int)sites.size();k++)
				pFile_process<<"N"<<i<<j<<k<<" ";
	pFile_process<<"E_d Std_d ";
	for (i=0;i<(int)sites.size();i++)
		pFile_process<<"E_d"<<i<<" Std_d"<<i<<" ";
//	pFile_process<<"Node_Between ";
//	for (j=0; j< (int)sites.size(); j++)
//		pFile_process<<"Node_Between"<<j<<" ";
//	pFile_process<<"Node_Close ";
//	for (j=0; j< (int)sites.size(); j++)
//		pFile_process<<"Node_Close"<<j<<" ";
	pFile_process<<"Node_Eig ";
	for (j=0; j< (int)sites.size(); j++)
		pFile_process<<"Node_Eig"<<j<<" ";
//	pFile_process<<"Edge_Between ";
//	for (i=0;i<(int)sites.size();i++)
//		for (j=i;j<(int)sites.size();j++)
//			pFile_process<<"Edge_Between"<<i<<j<<" ";
	for (i=0;i<(int)sites.size();i++)
		for (j=0;j<=11;j++){
			k=i;
			for(l=j;l<=11;l++)
				pFile_process<<"E"<<i<<"_"<<j<<"_"<<k<<"_"<<l<<" ";
			for (k=i+1;k<(int)sites.size();k++)
				for (l=0;l<=11;l++)
					pFile_process<<"E"<<i<<"_"<<j<<"_"<<k<<"_"<<l<<" ";
		}
	pFile_process<<endl;
	pFile_process<<step<<" "<<real_time<<" ";
	print_statistics(pFile_process);

//	component();
//	for(comp_it=comp_size.begin()+1;comp_it!=comp_size.end();comp_it++){
//		pFile_comp<<*comp_it<<" ";
//	}
//	pFile_comp<<endl;
//	for(comp_it=comp_size_1.begin()+1;comp_it!=comp_size_1.end();comp_it++){
//		pFile_comp_1<<*comp_it<<" ";
//	}
//	pFile_comp_1<<endl;

//	degree_dist(degreedist,degreedist1);
//	for(comp_it=degreedist.begin();comp_it!=degreedist.end();comp_it++){
//		pFile_degree<<*comp_it<<" ";
//	}
//	pFile_degree<<endl;
//	for(comp_it=degreedist1.begin();comp_it!=degreedist1.end();comp_it++){
//		pFile_degree1<<*comp_it<<" ";
//	}
//	pFile_degree1<<endl;

//	printEdgelist();

//  visualization using matlab
//	vector<long int> sort_comp, sort_comp_1;
//	sort_comp=comp_size;
//	sort_comp_1=comp_size_1;
//	sort(sort_comp.begin(),sort_comp.end());
//	sort(sort_comp_1.begin(),sort_comp_1.end());
//	Engine *ep =engOpen(NULL);
//	stringstream matlab_cmd;
//	engEvalString(ep, "figure;hold on;");
//	matlab_cmd<<"plot("<<step<<","<<sort_comp[sort_comp.size()-1]<<",'b.');"
//			"plot("<<step<<","<<sort_comp_1[sort_comp_1.size()-1]<<",'r.');";
//	engEvalString(ep, matlab_cmd.str().c_str());

//	stringstream matlab_cmd;
//	Engine *ep =engOpen(NULL);
//	engEvalString(ep, "figure;");
//	matlab_cmd << "h(1)=subplot(4,1,1);"
//			"title('alpha = " << alpha << ", initial fraction = " <<setprecision(2)<<(float)sites[1].size()/population.size()<<"');"
//			"ylabel('number of 0-1 edges');"
//			"hold on;"
//			"plot(" << step << "," << (float)edge_boundary.size()/edges.size()	<< ",'.');"
//			"h(2)=subplot(4,1,2);"
//			"ylabel('fraction of 0s');"
//			"hold on;"
//			"plot(" << step << "," << (float)sites[0].size()/population.size() << ",'.');"
//			"h(3)=subplot(4,1,3);"
//			"ylabel('fraction of 1s');"
//			"hold on;"
//			"plot(" << step << "," << (float)sites[1].size()/population.size() << ",'.');"
//			"h(4)=subplot(4,1,4);"
//			"ylabel('fraction of 2s');"
//			"xlabel('step');"
//			"hold on;"
//			"plot(" << step << "," << (float)sites[2].size()/population.size() << ",'.');";
//	engEvalString(ep, matlab_cmd.str().c_str());


	if (mode==1) // rewire to random
	{
		while (edge_boundary.empty() == false)
		{
			real_time+=(float)edges.size()/edge_boundary.size();
			e1 = random_number.integer(edge_boundary.size());
			edge_it=edge_boundary[e1];
			if (random_number.real()<alpha)
			{
				rand_rewire(edge_it);
			}
			else
			{
				adopt_state(edge_it);
			}
			step++;

			if ((long int)step%dt == 0){
				pFile_process<<step<<" "<<real_time<<" ";
				print_statistics(pFile_process);
//				component();
//				for(comp_it=comp_size.begin()+1;comp_it!=comp_size.end();comp_it++)
//				{
//					pFile_comp<<*comp_it<<" ";
//				}
//				pFile_comp<<endl;
//				for(comp_it=comp_size_1.begin()+1;comp_it!=comp_size_1.end();comp_it++)
//				{
//					pFile_comp_1<<*comp_it<<" ";
//				}
//				pFile_comp_1<<endl;
//				degree_dist(degreedist,degreedist1);
//				for(comp_it=degreedist.begin();comp_it!=degreedist.end();comp_it++)
//				{
//					pFile_degree<<*comp_it<<" ";
//				}
//				pFile_degree<<endl;
//				for(comp_it=degreedist1.begin();comp_it!=degreedist1.end();comp_it++)
//				{
//					pFile_degree1<<*comp_it<<" ";
//				}
//				pFile_degree1<<endl;
//				printEdgelist();
			}

//			visulization using matlab
//			if ((long int)step%dt == 0)
//			{
//				sort_comp=comp_size;
//				sort_comp_1=comp_size_1;
//				sort(sort_comp.begin(),sort_comp.end());
//				sort(sort_comp_1.begin(),sort_comp_1.end());
//				matlab_cmd.str("");
//				matlab_cmd<<"plot("<<step<<","<<sort_comp[sort_comp.size()-1]<<",'b.');"
//						"plot("<<step<<","<<sort_comp_1[sort_comp_1.size()-1]<<",'r.');";
//				engEvalString(ep, matlab_cmd.str().c_str());
//				matlab_cmd.str("");
//				matlab_cmd << "axes(h(1)); plot(" << step << "," << (float)edge_boundary.size()/edges.size()<< ",'.');"
//						"axes(h(2)); plot(" << step << "," << (float)one_sites.size()/population.size() << ",'.');";
//				engEvalString(ep, matlab_cmd.str().c_str());
//			}
		}
	}
	else if (mode==2) // rewire to same
	{
		while (edge_boundary.empty() == false)
		{
			e1 = random_number.integer(edge_boundary.size());
			edge_it=edge_boundary[e1];
			if (random_number.real()<alpha)
			{
				if (pref_rewire(edge_it)==false)
					counter++;
			}
			else
			{
				adopt_state(edge_it);
			}
			step++;

			if ((long int)step%dt == 0)
			{
				pFile_process<<step<<" "<<real_time<<" ";
				print_statistics(pFile_process);
//				degree_dist(degreedist,degreedist1);
//				for(comp_it=degreedist.begin();comp_it!=degreedist.end();comp_it++)
//				{
//					pFile_degree<<*comp_it<<" ";
//				}
//				pFile_degree<<endl;
//				component();
//				for(comp_it=comp_size.begin()+1;comp_it!=comp_size.end();comp_it++)
//				{
//					pFile_comp<<*comp_it<<" ";
//				}
//				pFile_comp<<endl;
//				for(comp_it=comp_size_1.begin()+1;comp_it!=comp_size_1.end();comp_it++)
//				{
//					pFile_comp_1<<*comp_it<<" ";
//				}
//				pFile_comp_1<<endl;
			}

//			if (step%dt == 0)
//			{
//				matlab_cmd.str("");
//				matlab_cmd << "axes(h(1)); plot(" << step << "," << (float)edge_boundary.size()/edges.size()<< ",'.');"
//						"axes(h(2)); plot(" << step << "," << (float)one_sites.size()/population.size() << ",'.');";
//				engEvalString(ep, matlab_cmd.str().c_str());
//			}
		}
	}
	else if(mode==3) // rewire to random and preserve degrees
	{
		while (edge_boundary.empty() == false)
		{
			e1 = random_number.integer(edge_boundary.size());
			edge_it=edge_boundary[e1];
			if (random_number.real()<alpha)
			{
				rand_rewire_DP(edge_it);
			}
			else
			{
				adopt_state(edge_it);
			}
			step++;

			//			if (step%dt == 0)
			//			{
			//				matlab_cmd.str("");
			//				matlab_cmd << "axes(h(1)); plot(" << step << "," << (float)edge_boundary.size()/edges.size()<< ",'.');"
			//						"axes(h(2)); plot(" << step << "," << (float)one_sites.size()/population.size() << ",'.');";
			//				engEvalString(ep, matlab_cmd.str().c_str());
			//			}
			//			outfile << step <<" "<< counter <<" "<<(float)one_sites.size()/population.size()<<" "<<edge_boundary.size()<<endl;
		}
	}
	else if(mode==4) // rewire to random and preserve degrees
	{
		while (edge_boundary.empty() == false)
		{
			e1 = random_number.integer(edge_boundary.size());
			edge_it=edge_boundary[e1];
			if (random_number.real()<alpha)
			{
				if (pref_rewire_DP(edge_it)==false)
					counter++;
			}
			else
			{
				adopt_state(edge_it);
			}
			step++;

//			if (step%dt == 0)
//			{
//				matlab_cmd.str("");
//				matlab_cmd << "axes(h(1)); plot(" << step << "," << (float)edge_boundary.size()/edges.size()<< ",'.');"
//						"axes(h(2)); plot(" << step << "," << (float)one_sites.size()/population.size() << ",'.');";
//				engEvalString(ep, matlab_cmd.str().c_str());
//			}
//			outfile << step <<" "<< counter <<" "<<(float)one_sites.size()/population.size()<<" "<<edge_boundary.size()<<endl;
		}
	}
	else if(mode==5)
	{
		while (edge_boundary.empty() == false)
		{
			real_time+=1.0/edge_boundary.size();
			if (edge_boundary.size()==1)
			{
				edge_it=edge_boundary[0];
				if (random_number.real()<alpha)
				{
					pref_rewire(edge_it);
				}
				else
				{
					adopt_state(edge_it);
				}
				step++;
			}
			else
			{
				e1 = random_number.integer(edge_boundary.size());
				e2=e1;
				counter=0;
				need_new_edge = true;
				edge_it = edge_boundary[e1];
				person1_it = edge_it->person1;
				person2_it = edge_it->person2;
				while (need_new_edge && counter<10)
				{
					e2 = random_number.integer(edge_boundary.size());
					need_new_edge = (e1==e2);
					if (need_new_edge==false)
					{
						edge_it2 = edge_boundary[e2];
						if (edge_it2->person1->state == person1_it->state)
						{
							person3_it = edge_it2->person1;
							person4_it = edge_it2->person2;
						}
						else
						{
							person3_it = edge_it2->person2;
							person4_it = edge_it2->person1;
						}
						need_new_edge = is_neighbor(person1_it->ID, person3_it->ID) || is_neighbor(person2_it->ID, person4_it->ID);
					}
					counter++;
				}
				if (need_new_edge == false)
				{
					if (random_number.real()<alpha)
					{
						delete_edge(edge_it,edge_it2);
						add_edge(person1_it->ID, person3_it->ID);
						add_edge(person2_it->ID, person4_it->ID);
					}
					else
					{
						adopt_state(edge_it);
						adopt_state(edge_it2);
					}
					step++;
				}
				else
				{
					if (random_number.real()<alpha)
					{
						pref_rewire(edge_it);
					}
					else
					{
						adopt_state(edge_it);
					}
					step++;
				}
			}

			if ((long int)step%dt == 0){
				pFile_process<<step<<" "<<real_time<<" ";
				print_statistics(pFile_process);
			}
		}
	}
	else if (mode==6)
	{
		while (edge_boundary.empty() == false && step<10000000)
		{
			real_time+=1.0/edge_boundary.size();
			e1 = random_number.integer(edge_boundary.size());
			edge_it=edge_boundary[e1];
			if (random_number.real()<alpha)
			{
				isrewire=1;
				rand_rewire(edge_it);
			}
			else
			{
				isrewire=0;
				swap_states(edge_it);
			}
			step++;

			if ((long int)step%dt == 0){
				pFile_process<<step<<" "<<real_time<<" ";
				print_statistics(pFile_process);
			}
		}
	}
	else if (mode==7)
	{
		while (edge_boundary.empty() == false)
		{
			real_time+=1.0/edge_boundary.size();
			e1 = random_number.integer(edge_boundary.size());
			edge_it=edge_boundary[e1];
			if (random_number.real()<alpha){
				rand_rewire(edge_it);
			}
			else{
				adopt_state2(edge_it);
			}
			step++;

			if ((long int)step%dt == 0){
				pFile_process<<step<<" "<<real_time<<" ";
				print_statistics(pFile_process);
			}
		}
	}

	pFile_process<<step<<" "<<real_time<<" ";
	print_statistics(pFile_process);
//	component();
//	for(comp_it=comp_size.begin()+1;comp_it!=comp_size.end();comp_it++)
//	{
//		pFile_comp<<*comp_it<<" ";
//	}
//	pFile_comp<<endl;
//	for(comp_it=comp_size_1.begin()+1;comp_it!=comp_size_1.end();comp_it++)
//	{
//		pFile_comp_1<<*comp_it<<" ";
//	}
//	pFile_comp_1<<endl;
//	degree_dist(degreedist,degreedist1);
//	for(comp_it=degreedist.begin();comp_it!=degreedist.end();comp_it++)
//	{
//		pFile_degree<<*comp_it<<" ";
//	}
//	pFile_degree<<endl;
//	for(comp_it=degreedist1.begin();comp_it!=degreedist1.end();comp_it++)
//	{
//		pFile_degree1<<*comp_it<<" ";
//	}
//	pFile_degree1<<endl;

	pFile_process.close();
//	pFile_comp.close();
//	pFile_comp_1.close();
//	pFile_degree.close();
//	pFile_degree1.close();
//	print();
//	matlab_cmd.str("");
//	matlab_cmd << "axes(h(1)); plot(" << step << "," << (float)edge_boundary.size()/edges.size()<< ",'.');"
//		"axes(h(2)); plot(" << step << "," << (float)sites[0].size()/population.size() << ",'.');"
//		"axes(h(3)); plot(" << step << "," << (float)sites[1].size()/population.size() << ",'.');"
//		"axes(h(4)); plot(" << step << "," << (float)sites[2].size()/population.size() << ",'.');";
//	engEvalString(ep, matlab_cmd.str().c_str());
//	engClose(ep);
	return step;
}

void Dynamic_Voter::adopt_state(vector<Edge>::iterator edge_it)
{
	vector<Node>::iterator person1_it, person2_it, person3_it;
	list<vector<Edge>::iterator>::iterator neigh_edge_it;

	person1_it = edge_it->person1;
	person2_it = edge_it->person2;
	if (random_number.real() < 0.5) // swap roles of the two people
	{
		person1_it = edge_it->person2;
		person2_it = edge_it->person1;
	}
	// person 1 adopts the state of person 2.

	swap_delete(person1_it, sites[person1_it->state]);
	person1_it->sites_place = sites[person2_it->state].insert(sites[person2_it->state].end(),person1_it);
	person1_it->state = person2_it->state;

	// need to go through all of person 1's neighbors:
	// each edge that was concordant is now discordant and vice-versa
	for (neigh_edge_it = person1_it->edge_list.begin(); neigh_edge_it!=person1_it->edge_list.end(); neigh_edge_it++)
	{
		edge_it = *neigh_edge_it;
		if (edge_it->state==0) // "0" means that this edge was concordant
		{
			edge_it->state=1;
			edge_it->boundary_place=edge_boundary.insert(edge_boundary.end(),edge_it);
			if (edge_it->person1->state==0 && edge_it->person2->state==1 || edge_it->person1->state==1 && edge_it->person2->state==0)
			{
				edge_it->state01=1;
				edge_it->boundary01_place=edge_boundary01.insert(edge_boundary01.end(),edge_it);
			}
		}
		else
		{ // edge was discordant, it may be concordant now
			if (edge_it->person1 == person1_it)
				person3_it=edge_it->person2;
			else
				person3_it=edge_it->person1;
			if (person1_it->state == person3_it->state)
			{
				edge_it->state=0;
				swap_delete(edge_it, edge_boundary);
			}
			if (edge_it->state01==0)
			{
				if (person1_it->state==0 && person3_it->state==1 || person1_it->state==1 && person3_it->state==0)
				{
					edge_it->state01=1;
					edge_it->boundary01_place=edge_boundary01.insert(edge_boundary01.end(),edge_it);
				}
			}
			else
			{
				edge_it->state01=0;
				swap_delete01(edge_it, edge_boundary01);
			}
		}
	}

	return;
}

void Dynamic_Voter::adopt_state2(vector<Edge>::iterator edge_it)
{
	vector<Node>::iterator person1_it, person2_it, person3_it;
	list<vector<Edge>::iterator>::iterator neigh_edge_it;

	if (edge_it->person1->state!=2 && edge_it->person2->state!=2)
	{
		if (edge_it->person1->state==lastadopt)
		{
			person1_it = edge_it->person1;
			person2_it = edge_it->person2;
			lastadopt=edge_it->person2->state;
		}
		else
		{
			person1_it = edge_it->person2;
			person2_it = edge_it->person1;
			lastadopt=edge_it->person1->state;
		}
	}
	else
	{
		person1_it = edge_it->person1;
		person2_it = edge_it->person2;
		if (random_number.real() < 0.5) // swap roles of the two people
		{
			person1_it = edge_it->person2;
			person2_it = edge_it->person1;
		}
	}
	// person 1 adopts the state of person 2.

	swap_delete(person1_it, sites[person1_it->state]);
	person1_it->sites_place = sites[person2_it->state].insert(sites[person2_it->state].end(),person1_it);
	person1_it->state = person2_it->state;

	// need to go through all of person 1's neighbors:
	// each edge that was concordant is now discordant and vice-versa
	for (neigh_edge_it = person1_it->edge_list.begin(); neigh_edge_it!=person1_it->edge_list.end(); neigh_edge_it++)
	{
		edge_it = *neigh_edge_it;
		if (edge_it->state==0) // "0" means that this edge was concordant
		{
			edge_it->state=1;
			edge_it->boundary_place=edge_boundary.insert(edge_boundary.end(),edge_it);
			if (edge_it->person1->state==0 && edge_it->person2->state==1 || edge_it->person1->state==1 && edge_it->person2->state==0)
			{
				edge_it->state01=1;
				edge_it->boundary01_place=edge_boundary01.insert(edge_boundary01.end(),edge_it);
			}
		}
		else
		{ // edge was discordant, it may be concordant now
			if (edge_it->person1 == person1_it)
				person3_it=edge_it->person2;
			else
				person3_it=edge_it->person1;
			if (person1_it->state == person3_it->state)
			{
				edge_it->state=0;
				swap_delete(edge_it, edge_boundary);
			}
			if (edge_it->state01==0)
			{
				if (person1_it->state==0 && person3_it->state==1 || person1_it->state==1 && person3_it->state==0)
				{
					edge_it->state01=1;
					edge_it->boundary01_place=edge_boundary01.insert(edge_boundary01.end(),edge_it);
				}
			}
			else
			{
				edge_it->state01=0;
				swap_delete01(edge_it, edge_boundary01);
			}
		}
	}

	return;
}

void Dynamic_Voter::swap_states(vector<Edge>::iterator edge_it)
{
	vector<Node>::iterator person1_it, person2_it, person3_it;
	list<vector<Edge>::iterator>::iterator neigh_edge_it;
	int tempstate;

	person1_it = edge_it->person1;
	person2_it = edge_it->person2;

	//swap states of person1 and person2
	tempstate=person1_it->state;
	swap_delete(person1_it, sites[person1_it->state]);
	person1_it->sites_place = sites[person2_it->state].insert(sites[person2_it->state].end(),person1_it);
	person1_it->state = person2_it->state;
	swap_delete(person2_it, sites[person2_it->state]);
	person2_it->sites_place = sites[tempstate].insert(sites[tempstate].end(),person2_it);
	person2_it->state = tempstate;

	// need to go through all of person 1's neighbors:
	// each edge that was concordant is now discordant and maybe vice-versa
	for (neigh_edge_it = person1_it->edge_list.begin(); neigh_edge_it!=person1_it->edge_list.end(); neigh_edge_it++)
	{
		edge_it = *neigh_edge_it;
		if (edge_it->state==0) // "0" means that this edge was concordant
		{
			edge_it->state=1;
			edge_it->boundary_place=edge_boundary.insert(edge_boundary.end(),edge_it);
			if (edge_it->person1->state==0 && edge_it->person2->state==1 || edge_it->person1->state==1 && edge_it->person2->state==0)
			{
				edge_it->state01=1;
				edge_it->boundary01_place=edge_boundary01.insert(edge_boundary01.end(),edge_it);
			}
		}
		else
		{ // edge was discordant, it may be concordant now
			if (edge_it->person1 == person1_it)
				person3_it=edge_it->person2;
			else
				person3_it=edge_it->person1;
			if (person1_it->state == person3_it->state)
			{
				edge_it->state=0;
				swap_delete(edge_it, edge_boundary);
			}
			if (edge_it->state01==0)
			{
				if (person1_it->state==0 && person3_it->state==1 || person1_it->state==1 && person3_it->state==0)
				{
					edge_it->state01=1;
					edge_it->boundary01_place=edge_boundary01.insert(edge_boundary01.end(),edge_it);
				}
			}
			else
			{
				edge_it->state01=0;
				swap_delete01(edge_it, edge_boundary01);
			}
		}
	}
	// need to go through all of person 2's neighbors:
	// each edge that was concordant is now discordant and maybe vice-versa
	for (neigh_edge_it = person2_it->edge_list.begin(); neigh_edge_it!=person2_it->edge_list.end(); neigh_edge_it++)
	{
		edge_it = *neigh_edge_it;
		if (edge_it->state==0) // "0" means that this edge was concordant
		{
			edge_it->state=1;
			edge_it->boundary_place=edge_boundary.insert(edge_boundary.end(),edge_it);
			if (edge_it->person1->state==0 && edge_it->person2->state==1 || edge_it->person1->state==1 && edge_it->person2->state==0)
			{
				edge_it->state01=1;
				edge_it->boundary01_place=edge_boundary01.insert(edge_boundary01.end(),edge_it);
			}
		}
		else
		{ // edge was discordant, it may be concordant now
			if (edge_it->person1 == person2_it)
				person3_it=edge_it->person2;
			else
				person3_it=edge_it->person1;
			if (person2_it->state == person3_it->state)
			{
				edge_it->state=0;
				swap_delete(edge_it, edge_boundary);
			}
			if (edge_it->state01==0)
			{
				if (person2_it->state==0 && person3_it->state==1 || person2_it->state==1 && person3_it->state==0)
				{
					edge_it->state01=1;
					edge_it->boundary01_place=edge_boundary01.insert(edge_boundary01.end(),edge_it);
				}
			}
			else
			{
				edge_it->state01=0;
				swap_delete01(edge_it, edge_boundary01);
			}
		}
	}

	return;
}

bool Dynamic_Voter::rand_rewire(vector<Edge>::iterator edge_it)
{
	vector<Node>::iterator person1_it, person2_it;
	long int i1, i2;
	bool need_new_neigh=true;

	person1_it = edge_it->person1;
	person2_it = edge_it->person2;
	if (random_number.real() < 0.5) // swap roles of the two people, but need to note this!
	{
		person1_it = edge_it->person2;
		person2_it = edge_it->person1;
	}
	i1 = person1_it->ID;
	while (need_new_neigh)
	{
		i2 = random_number.integer(population.size());
		need_new_neigh = (i1==i2) || is_neighbor(i1, i2);
	}
	delete_edge(edge_it);
	add_edge(i1, i2);
	// update igraph
	//	igraph_integer_t eid;
	//	igraph_get_eid(&igraph,&eid,person1_it->ID,person2_it->ID,IGRAPH_UNDIRECTED,0);
	//	igraph_delete_edges(&igraph,igraph_ess_1(eid));
	//	igraph_add_edge(&igraph,i1,i2);
//	gsl_matrix_set(igraph,person1_it->ID,person2_it->ID,0);
//	gsl_matrix_set(igraph,person2_it->ID,person1_it->ID,0);
//	gsl_matrix_set(igraph,i1,i2,1);
//	gsl_matrix_set(igraph,i2,i1,1);
	return true;
}

bool Dynamic_Voter::pref_rewire(vector<Edge>::iterator edge_it)
{
	vector<Node>::iterator person1_it, person2_it;
	long int i1, i2;
	int counter=0;
	bool need_new_neigh=true;

	person1_it = edge_it->person1;
	person2_it = edge_it->person2;
	if (random_number.real() < 0.5) // swap roles of the two people, but need to note this!
	{
		person1_it = edge_it->person2;
		person2_it = edge_it->person1;
	}

	// remove edge between person 1 and person 2
	// then form list between person 1 and person x
	// x chosen at random from the same group as person 1

	i1 = person1_it->ID;
	while (need_new_neigh && counter<10)
	{
		i2 = random_number.integer(sites[person1_it->state].size());
		i2 = (sites[person1_it->state][i2])->ID;
		need_new_neigh = (i1==i2) || is_neighbor(i1, i2);
		counter++;
	}

	if (need_new_neigh == false)
	{
		delete_edge(edge_it);
		add_edge(i1, i2);
		return true;
	}
	else
		return false;
}

bool Dynamic_Voter::pref_rewire_DP(vector<Edge>::iterator edge_it)
{
	vector<Node>::iterator person1_it, person2_it, person3_it, person4_it;
	vector<Edge>::iterator new_edge;
	int counter=0;
	long int e1;
	bool need_new_edge = true;

	if (edge_boundary.size()<=1)
		return false;

	person1_it = edge_it->person1;
	person2_it = edge_it->person2;
	// choose another discordant edge that orients the same way as this edge, then rewire p1 to p3, p2 to p4.
	while (need_new_edge && counter< 10)
	{
		e1 = random_number.integer(edge_boundary.size());
		need_new_edge = (edge_it==edge_boundary[e1]);
		if (need_new_edge==false)
		{
			new_edge = edge_boundary[e1];
			if (new_edge->person1->state == person1_it->state)
			{
				person3_it = new_edge->person1;
				person4_it = new_edge->person2;
			}
			else
			{
				person3_it = new_edge->person2;
				person4_it = new_edge->person1;
			}
			need_new_edge = is_neighbor(person1_it->ID, person3_it->ID) || is_neighbor(person2_it->ID, person4_it->ID);
		}
		counter++;
	}
	if (need_new_edge == false)
	{
		delete_edge(edge_it,new_edge); //need to delete the two edges at same time, because deleting one edge changes the edge list.
		add_edge(person1_it->ID, person3_it->ID);
		add_edge(person2_it->ID, person4_it->ID);
		return true;
	}
	return false;
}

bool Dynamic_Voter::rand_rewire_DP(vector<Edge>::iterator edge_it)
{
	vector<Node>::iterator person1_it, person2_it, person3_it, person4_it;
	vector<Edge>::iterator new_edge;
	bool need_new_edge = true;
	long int e1;

	person1_it = edge_it->person1;
	person2_it = edge_it->person2;
	if (random_number.real() < 0.5) // swap roles of the two people, but need to note this!
	{
		person1_it = edge_it->person2;
		person2_it = edge_it->person1;
	}
	while (need_new_edge)
	{
		e1 = random_number.integer(edges.size());
		new_edge = edges[e1].myself;
		need_new_edge = (edge_it==new_edge);
		if (need_new_edge==false)
		{
			person3_it = new_edge->person1;
			person4_it = new_edge->person2;
			if (random_number.real() < 0.5) // swap roles of the two people, but need to note this!
			{
				person3_it = new_edge->person2;
				person4_it = new_edge->person1;
			}
			need_new_edge = is_neighbor(person1_it->ID, person3_it->ID) || is_neighbor(person2_it->ID, person4_it->ID);
		}
	}
	delete_edge(edge_it,new_edge); //need to delete the two edges at same time, because deleting one edge changes the edge list.
	add_edge(person1_it->ID, person3_it->ID);
	add_edge(person2_it->ID, person4_it->ID);
	return true;
}

void Dynamic_Voter::component(vector<int> &comp, vector<int> &comp_size)
{
	queue<int> que;
	vector<Node>::iterator person_it;
	vector<Edge>::iterator edge_it;
	list<vector<Edge>::iterator>::iterator neigh_edge_it;
	long int num_comp=0, n1, n2;

	comp_size.push_back(0);
	for (person_it = population.begin(); person_it != population.end(); person_it++)
	{
		if(comp[person_it->ID]==0)
		{
			que.push(person_it->ID);
			while(que.empty()==false)
			{
				n1=que.front();
				que.pop();
				if(comp[n1]==0)
				{
					comp[n1]=++num_comp;
					comp_size.push_back(1);
				}
				for (neigh_edge_it = population[n1].edge_list.begin(); neigh_edge_it!=population[n1].edge_list.end(); neigh_edge_it++)
				{
					edge_it = *neigh_edge_it;
					n2=edge_it->person1->ID;
					if(n2==n1)
						n2=edge_it->person2->ID;
					if(comp[n2]==0)
					{
						comp[n2]=comp[n1];
						que.push(n2);
						(comp_size[num_comp])++;
					}
				}
			}
		}
	}
}

void Dynamic_Voter::component()
{
	queue<long int> que;
	vector<Node>::iterator person_it;
	vector<Edge>::iterator edge_it;
	list<vector<Edge>::iterator>::iterator neigh_edge_it;
	long int num_comp=0, n1, n2;
	comp.assign(population.size(),0);
	comp_size.clear();
	comp_size_1.clear();
	comp_size.reserve((int)(population.size()/log((float)population.size())));
	comp_size_1.reserve((int)(population.size()/log((float)population.size())));
	comp_size.push_back(0);
	comp_size_1.push_back(0);

	for (person_it = population.begin(); person_it != population.end(); person_it++)
	{
		if(comp[person_it->ID]==0)
		{
			que.push(person_it->ID);
			while(que.empty()==false)
			{
				n1=que.front();
				que.pop();
				if(comp[n1]==0)
				{
					comp[n1]=++num_comp;
					comp_size.push_back(1);
					if (population[n1].state==1)
						comp_size_1.push_back(1);
					else
						comp_size_1.push_back(0);
				}
				for (neigh_edge_it = population[n1].edge_list.begin(); neigh_edge_it!=population[n1].edge_list.end(); neigh_edge_it++)
				{
					edge_it = *neigh_edge_it;
					n2=edge_it->person1->ID;
					if(n2==n1)
						n2=edge_it->person2->ID;
					if(comp[n2]==0)
					{
						comp[n2]=comp[n1];
						que.push(n2);
						(comp_size[num_comp])++;
						if (population[n2].state==1)
							(comp_size_1[num_comp])++;
					}
				}
			}
		}
	}
}

void Dynamic_Voter::printEdgelist()
{
	vector<Edge>::iterator edge_it;
	vector<Node>::iterator person_it;
	ofstream edgelist;
	ofstream status;

	edgelist.open("edgelist1",ios_base::app);
	for (edge_it = edges.begin(); edge_it != edges.end(); edge_it++)
	{
		edgelist<<edge_it->person1->ID<<" "<<edge_it->person2->ID<<endl;
	}
	edgelist.close();
	status.open("status1",ios_base::app);
	for (person_it = population.begin(); person_it!=population.end(); person_it++)
	{
		status << person_it->state <<endl;
	}
	status.close();
}

void Dynamic_Voter::print_statistics(ofstream &pFile_process)
{
	int i,j,k,l;
	vector<Node>::iterator person_it,person_it1,person_it2;
	vector<Edge>::iterator edge_it;
	list<vector<Edge>::iterator>::iterator neigh_edge_it1, neigh_edge_it2;
	long int maxd=0,first_moment_d_whole=0,second_moment_d_whole=0;
	long int *d=new long int [sites.size()];
	long int *first_moment_d=new long int [sites.size()];
	long int *second_moment_d=new long int [sites.size()];
	for (i=0;i<(int)sites.size();i++)
	{
		first_moment_d[i]=0;
		second_moment_d[i]=0;
	}
	long int **N_edges=new long int* [sites.size()];
	for (i=0;i<(int)sites.size();i++)
	{
		N_edges[i]=new long int [sites.size()];
		for (j=0;j<(int)sites.size();j++)
			N_edges[i][j]=0;
	}
	long int ***N_triples=new long int** [sites.size()];
	for (i=0;i<(int)sites.size();i++)
	{
		N_triples[i]=new long int* [sites.size()];
		for (j=0;j<(int)sites.size();j++)
		{
			N_triples[i][j]=new long int [sites.size()];
			for (k=0;k<(int)sites.size();k++)
				N_triples[i][j][k]=0;
		}
	}
//	double node_between_whole=0,*node_between=new double[sites.size()];
//	double node_close_whole=0,*node_close=new double[sites.size()];
	double node_eig_whole=0,*node_eig=new double[sites.size()];
	for (i=0;i<(int)sites.size();i++)
	{
//		node_between[i]=0;
//		node_close[i]=0;
		node_eig[i]=0;
	}
//	double edge_between_whole=0,**edge_between=new double*[sites.size()];
//	for (i=0;i<(int)sites.size();i++)
//	{
//		edge_between[i]=new double[sites.size()];
//		for (j=0;j<(int)sites.size();j++)
//			edge_between[i][j]=0;
//	}
	long int ****e=new long int***[sites.size()];
	for (i=0;i<(int)sites.size();i++)
	{
		e[i]=new long int**[12];
		for (j=0;j<=11;j++)
		{
			e[i][j]=new long int*[sites.size()];
			for (k=0;k<(int)sites.size();k++)
			{
				e[i][j][k]=new long int[12];
				for (l=0;l<=11;l++)
					e[i][j][k][l]=0;
			}
		}
	}

	//	igraph_vector_t iresult;
	//	igraph_arpack_options_t ioptions;
//	gsl_vector* iresult;

//	Calculate N_ij, N_ijk, average degree and second moment of the degree distribution.
//	Note: N11, N00, N111, N000, N101, N010 are double counted
	for (person_it = population.begin(); person_it != population.end(); person_it++)
	{
		for (i=0;i<(int)sites.size();i++)
			d[i]=0;
		if (maxd<(int)person_it->edge_list.size())
			maxd=person_it->edge_list.size();
		for (neigh_edge_it1 = person_it->edge_list.begin(); neigh_edge_it1!=person_it->edge_list.end(); neigh_edge_it1++)
		{
			if ((*neigh_edge_it1)->person1==person_it)
				person_it1=(*neigh_edge_it1)->person2;
			else
				person_it1=(*neigh_edge_it1)->person1;
			d[person_it1->state]++;
		}
		for (i=person_it->state;i<(int)sites.size();i++)
			N_edges[person_it->state][i]+=d[i];
		for (i=0;i<(int)sites.size();i++)
			N_triples[i][person_it->state][i]+=d[i]*(d[i]-1);
		for (i=0;i<(int)sites.size();i++)
			for (j=i+1;j<(int)sites.size();j++)
				N_triples[i][person_it->state][j]+=(d[i])*(d[j]);
		first_moment_d_whole+=person_it->edge_list.size();
		second_moment_d_whole+=(person_it->edge_list.size())*(person_it->edge_list.size());
		first_moment_d[person_it->state]+=person_it->edge_list.size();
		second_moment_d[person_it->state]+=(person_it->edge_list.size())*(person_it->edge_list.size());
	}

//	Calculate the number of edges connecting nodes with state i and degree j to nodes nodes with state k and degree l
	for (edge_it = edges.begin(),i=0; edge_it != edges.end(); edge_it++,i++)
	{
		person_it1=edge_it->person1;
		person_it2=edge_it->person2;
		e[person_it1->state][min((int)person_it1->edge_list.size(),11)][person_it2->state][min((int)person_it2->edge_list.size(),11)]++;
		e[person_it2->state][min((int)person_it2->edge_list.size(),11)][person_it1->state][min((int)person_it1->edge_list.size(),11)]++;
	}

//	Initialize igraph vector
//	igraph_integer_t i1,i2;
//	igraph_vector_init(&iresult, 0);
//	igraph_arpack_options_init(&ioptions);

//	Calculate edge betweenness
//	igraph_edge_betweenness(&igraph, &iresult, IGRAPH_UNDIRECTED,0);
//	gsl_matrix* EBC;
//	EBC = bct::edge_betweenness_bin(igraph, &iresult);
//	for (edge_it = edges.begin(),i=0; edge_it != edges.end(); edge_it++,i++)
//	{
//		person_it1=edge_it->person1;
//		person_it2=edge_it->person2;
//		if (person_it1->state>person_it2->state)
//		{
//			person_it1=edge_it->person2;
//			person_it2=edge_it->person1;
//		}
//		edge_between_whole+=VECTOR(iresult)[i];
//		igraph_edge(&igraph,i, &i1, &i2);
//		if (i1>i2)
//			edge_between[population[i2].state][population[i1].state]+=VECTOR(iresult)[i];
//		else
//			edge_between[population[i1].state][population[i2].state]+=VECTOR(iresult)[i];
//		edge_between_whole+=gsl_matrix_get(EBC,person_it1->ID,person_it2->ID);
//		edge_between[person_it1->state][person_it2->state]+=gsl_matrix_get(EBC,person_it1->ID,person_it2->ID);
//	}
//	gsl_matrix_free(EBC);

//	Calculate node betweenness
//	igraph_betweenness(&igraph, &iresult, igraph_vss_all(), IGRAPH_UNDIRECTED,0,1);
//	for (person_it = population.begin(),i=0; person_it != population.end(); person_it++,i++)
//	{
//		node_between_whole+=VECTOR(iresult)[i];
//		node_between[person_it->state]+=VECTOR(iresult)[i];
//		node_between_whole+=gsl_vector_get(iresult,i);
//		node_between[person_it->state]+=gsl_vector_get(iresult,i);
//	}
//	gsl_vector_free(iresult);

//	for (person_it = population.begin(),i=0; person_it != population.end(); person_it++,i++)
//	{
//		node_close_whole+=VECTOR(iresult)[i];
//		node_close[person_it->state]+=VECTOR(iresult)[i];
//	}

//	Calculate eigenvector centrality
//	igraph_eigenvector_centrality(&igraph, &iresult, 0, IGRAPH_UNDIRECTED, 0, 0,&ioptions);
//	iresult=bct::eigenvector_centrality(igraph);
	for (person_it = population.begin(),i=0; person_it != population.end(); person_it++,i++)
	{
	  //		node_eig_whole+=VECTOR(iresult)[i];
	  //	node_eig[person_it->state]+=VECTOR(iresult)[i];
//		node_eig_whole+=gsl_vector_get(iresult,i);
//		node_eig[person_it->state]+=gsl_vector_get(iresult,i);
	}
//	gsl_vector_free(iresult);

//	Output results
	for (j=0; j<(int)sites.size(); j++)
		pFile_process<<(float)sites[j].size()/population.size()<<" ";
	pFile_process<<(float)edge_boundary.size()/edges.size()<<" ";
	for (i=0; i<(int)sites.size(); i++)
		for (j=i; j<(int)sites.size(); j++)
			pFile_process<<(float)N_edges[i][j]/edges.size()/2<<" ";
	for (i=0; i<(int)sites.size(); i++)
		for (j=0; j<(int)sites.size(); j++)
			for (k=i; k<(int)sites.size(); k++)
				pFile_process<<(float)N_triples[i][j][k]/edges.size()/(edges.size()-1)*2<<" ";
	pFile_process<<(float)first_moment_d_whole/population.size()/degree<<" ";
	pFile_process<<sqrt(second_moment_d_whole/(float)population.size()-(first_moment_d_whole/(float)population.size())*(first_moment_d_whole/(float)population.size()))/degree<<" ";
	for (i=0; i<(int)sites.size(); i++)
	{
		if (sites[i].size()!=0)
		{
			pFile_process<<(float)first_moment_d[i]/sites[i].size()/degree<<" ";
			pFile_process<<sqrt(second_moment_d[i]/(float)sites[i].size()-(first_moment_d[i]/(float)sites[i].size())*(first_moment_d[i]/(float)sites[i].size()))/degree<<" ";
		}
		else
			pFile_process<<"0 0 ";
	}
//	pFile_process<<node_between_whole/population.size()<<" ";
//	for (i=0; i<(int)sites.size(); i++)
//		pFile_process<<node_between[i]/sites[i].size()<<" ";
//	pFile_process<<node_close_whole<<" ";
//	for (i=0; i<(int)sites.size(); i++)
//		pFile_process<<node_close[i]<<" ";
	pFile_process<<node_eig_whole/population.size()<<" ";
	for (i=0; i<(int)sites.size(); i++)
		if (sites[i].size()!=0)
			pFile_process<<node_eig[i]/sites[i].size()<<" ";
		else
			pFile_process<<"0 ";
//	pFile_process<<edge_between_whole/edges.size()<<" ";
//	for (i=0; i<(int)sites.size(); i++)
//	{
//		edge_between[i][i]*=2;
//		for (j=i; j<(int)sites.size(); j++)
//			pFile_process<<edge_between[i][j]/N_edges[i][j]<<" ";
//	}
	for (i=0;i<(int)sites.size();i++)
		for (j=0;j<=11;j++)
		{
			k=i;
			for(l=j;l<=11;l++)
				pFile_process<<(float)e[i][j][k][l]/edges.size()/2<<" ";
			for (k=i+1;k<(int)sites.size();k++)
				for (l=0;l<=11;l++)
					pFile_process<<(float)e[i][j][k][l]/edges.size()/2<<" ";
		}
	pFile_process<<endl;

	delete[] d;
	delete[] first_moment_d;
	delete[] second_moment_d;
//	delete[] node_between;
//	delete[] node_close;
	delete[] node_eig;
	for (i=0;i<(int)sites.size();i++)
	{
		delete[] N_edges[i];
//		delete[] edge_between[i];
		for (j=0;j<(int)sites.size();j++)
		{
			delete[] N_triples[i][j];
		}
		delete[] N_triples[i];
	}
	delete[] N_edges;
	delete[] N_triples;
//	delete[] edge_between;
//	igraph_vector_destroy(&iresult);
	for (i=0;i<(int)sites.size();i++)
	{
		for (j=0;j<=11;j++)
		{
			for (k=0;k<(int)sites.size();k++)
				delete[] e[i][j][k];
			delete[] e[i][j];
		}
		delete[] e[i];
	}
	delete[] e;
}

void Dynamic_Voter::degree_dist(vector<long int> &d, vector<long int> &d1)
{
	vector<Node>::iterator person_it;
	d.clear();
	d1.clear();
	d.reserve(population.size());
	d1.reserve(population.size());

	for (person_it = population.begin(); person_it != population.end(); person_it++)
	{
		if(d.size()<=person_it->edge_list.size())
		{
			d.resize(person_it->edge_list.size()+1);
			d1.resize(person_it->edge_list.size()+1);
		}
		d[(person_it->edge_list.size())]++;
		if (person_it->state==1)
			d1[(person_it->edge_list.size())]++;
	}
}

