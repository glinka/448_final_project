/*
 * main.cpp
 *
 *  Created on: Apr 18, 2011
 *      Author: Bill
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <time.h>
#include <vector>
#include "Dynamic_Voter.h"
#include <algorithm>
#include <math.h>
#include <cstdlib>
#include<cstring>
using namespace std;

int main(int argc, char* argv[]){ // -n number_of_nodes -m number_of_edges -d degree -u number_of_opinions initial_densities -a alpha -r mode -b beta
	long int number_of_nodes, number_of_edges;
	double steps;
	float alpha,beta;
	int mode, realization, i,j, degree, number_of_opinions;
	vector<float> initial_density;
	clock_t t_start;
	stringstream filename;

//	Initialize parameters
	number_of_nodes=1000;
	number_of_edges=2000;
	degree=4;
	number_of_opinions=2;
	initial_density.push_back(1.0/2);
	initial_density.push_back(1.0/2);
	mode=1;
	alpha=0.3;
	beta=alpha/(1-alpha);
	if(argc>1){
		for (i = 1; i < argc; i++){
			if (!strcmp(argv[i],"-n")){
				i++;
				number_of_nodes=atol(argv[i]);
			}
			else if (!strcmp(argv[i],"-m")){
				i++;
				number_of_edges=atol(argv[i]);
				degree=number_of_edges*2/number_of_nodes;
			}
			else if (!strcmp(argv[i],"-d")){
				i++;
				degree=atoi(argv[i]);
				number_of_edges=number_of_nodes*degree/2;
			}
			else if (!strcmp(argv[i],"-u")){
				i++;
				number_of_opinions=atoi(argv[i]);
				for (int j=0; j<number_of_opinions; j++){
					i++;
					initial_density.push_back(atof(argv[i]));
				}
			}
			else if (!strcmp(argv[i],"-a")){
				i++;
				alpha=atof(argv[i]);
				beta=alpha/(1-alpha);
			}
			else if (!strcmp(argv[i],"-b")){
				i++;
				beta=atof(argv[i]);
				alpha=beta/(1+beta);
			}
			else if (!strcmp(argv[i],"-r")){
				i++;
				mode=atoi(argv[i]);
			}
			else{
				cout<<"unknown parameter: "<<argv[i]<<endl;
				exit(1);
			}
		}
	}

//	summary of the simulation
	filename<<"mode"<<mode<<"n"<<number_of_nodes<<"k"<<degree<<"b"<< setiosflags(ios::fixed) << setprecision(2)<<beta<<"u";
	for (j=0; j<number_of_opinions-1; j++)
		filename<< initial_density[j]<<"_";
	filename<< initial_density[j]<<".summary";

//	test if the file already exists
	ifstream testfile(filename.str().c_str());
	i=1;
	while (testfile.is_open()){
		i++;
		filename.str("");
		filename<<"mode"<<mode<<"n"<<number_of_nodes<<"k"<<degree<<"b"<<beta<<"u";
		for (j=0; j<number_of_opinions-1; j++)
			filename <<initial_density[j] <<"_";
		filename<<initial_density[j] <<"("<< i <<").summary";
		testfile.close();
		testfile.open(filename.str().c_str());
	}
	ofstream outfile(filename.str().c_str());
	outfile<<"realization "<<"alpha ";
	for (i=0; i<number_of_opinions; i++)
		outfile<<"u"<<i<<" ";
	outfile<< "steps " << "time " << "discordant_edges " << "number_of_components "<<"largest_comp_size"<< endl;

	for(realization=1;realization<=1;realization++){
		t_start=clock();
//		Dynamic_Voter *dv = new Dynamic_Voter(number_of_nodes,degree,number_of_opinions);
		Dynamic_Voter *dv = new Dynamic_Voter(number_of_nodes,number_of_edges,number_of_opinions);
		dv->assign_states(initial_density);
		steps = dv->simulate(alpha,beta,mode);
		dv->component();
//		outfile <<realization<<" "<< alpha <<" "<< min(dv->zero_sites.size(),dv->one_sites.size())/(float)number_of_nodes <<" "<< steps <<" "<< (float)(clock()-t_start)/CLOCKS_PER_SEC <<" "<< dv->edge_boundary.size()<<" "<<dv->comp_size.size()-1<<" "<<*max_element(dv->comp_size.begin(),dv->comp_size.end())/(float)number_of_nodes<< endl;
		outfile <<realization<<" "<< alpha <<" ";
		for (i=0; i<number_of_opinions; i++)
			outfile << dv->sites[i].size()/(float)number_of_nodes <<" ";
		outfile << steps <<" "<< (float)(clock()-t_start)/CLOCKS_PER_SEC <<" "<< dv->edge_boundary.size()<<" "<<dv->comp_size.size()-1<<" "<<*max_element(dv->comp_size.begin(),dv->comp_size.end())/(float)number_of_nodes<< endl;
		delete dv;
	}

	outfile.close();
	return 0;
}

