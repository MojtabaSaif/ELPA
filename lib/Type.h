#pragma once
#include<map>
#include<set>
#include<unordered_map>
#include<unordered_set>
#include<vector>
#include<deque>
typedef  int           CLUSTER_ID;
typedef  int           NODE_ID;
typedef  unsigned int  NODE_NAME;
typedef  double        WEIGHT;
typedef  double        MEMBERSHIPE;
typedef  unsigned int  LABEL_ID;

//struct LABEL_ID { unsigned int Id; int Mbr; };
using std::vector;
using std::map;
using std::set;
using std::unordered_map;
using std::unordered_set;
using std::string;
typedef map<NODE_NAME, NODE_ID> MAP_NAME_ID;

/*template<class myKEY, class myVALUE>
class myMAP {
public:
	std::map<myKEY, myVALUE> m;
	int FindByValue() { return -1; };
};*/
class OverlappingWeightedGraph {
public:
	int numEdges=0;
	int numNodes=0;
	vector<map<NODE_ID, WEIGHT>>   edges;
	MAP_NAME_ID                    nodeNames;
	vector<WEIGHT>                 inflNodes;
	vector<WEIGHT>                 capecity;
};

typedef   OverlappingWeightedGraph OWG;
typedef   set<NODE_ID> CLUSTER;
typedef   vector<CLUSTER>        CLUSTERING;
typedef  std::pair<LABEL_ID, MEMBERSHIPE> PACKET;
//typedef   vector<map<LABEL_ID,MEMBERSHIPE>>  MBR_CLUSTERING;
class NodeLabeles: public std::vector<PACKET>
{
public:
	vector<PACKET>::iterator Find(LABEL_ID id);
};
typedef   vector<NodeLabeles>  MBR_CLUSTERING;
