/*
 * Node.h
 *
 *  Created on: Apr 17, 2011
 *      Author: Bill
 */

#ifndef NODE_H_
#define NODE_H_

#include <list>
#include <vector>
class Edge;
#include "Edge.h"
using namespace std;

class Node
{
public:
	long int ID;
	int state;
	vector<Node>::iterator myself;
	list<vector<Edge>::iterator> edge_list; //iterators to the edges that join to/from this node
	vector<vector<Node>::iterator>::iterator sites_place; // iterator to where it is in the list of i nodes

public:
	Node();
	virtual ~Node();
};

#endif /* NODE_H_ */
