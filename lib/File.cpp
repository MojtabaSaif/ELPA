#include"File.h"

bool File::LoadClustering(string fName, const MAP_NAME_ID& nodeNames, CLUSTERING& cl)
{
	std::ifstream file;
	file.open(fName);
	if (!file)
	{
		std::cerr << "Can not open : " << fName << "\n";
		return false;
	}
	else
	{
		std::string line;
		cl.clear();
		int label = 0;
		while (std::getline(file, line))
		{
			std::istringstream str(line);
			NODE_NAME name;
			cl.resize(label + 1);
			while (str >> name)
			{
				NODE_ID nID;
				MAP_NAME_ID::const_iterator it = nodeNames.find(name);
				if (it != nodeNames.end())
					nID = it->second;
				else
					continue;

				cl[label].insert(nID);
			}
			label++;
		}
	}

	return true;
}
bool File::LoadFileEdgeList(string fName, OWG& graph)
{

	std::ifstream inFile;
	inFile.open(fName);
	if (!inFile)
	{
		std::cerr << "Can not open file :" << fName << std::endl;
		return false;
	}
	else
	{
		string line;
		map<NODE_NAME, set<NODE_NAME>> tmp;
		map<NODE_NAME, set<NODE_NAME>>::iterator it_tmp;
		int id = 0;
		while (std::getline(inFile, line))
		{
			if (line[0] != '#')
			{
				std::istringstream str(line);
				NODE_NAME n1, n2;
				NODE_ID   n1_id, n2_id;
				str >> n1 >> n2;
				map<NODE_NAME, NODE_ID>::iterator it = graph.nodeNames.find(n1);
				if (it == graph.nodeNames.end())
				{
					it = graph.nodeNames.insert(std::make_pair(n1, id)).first;
					id++;
				}
				n1_id = it->second;
				it = graph.nodeNames.find(n2);
				if (it == graph.nodeNames.end())
				{
					it = graph.nodeNames.insert(std::make_pair(n2, id)).first;
					id++;
				}
				n2_id = it->second;

				it_tmp = tmp.find(n1_id);
				if (it_tmp == tmp.end())
				{
					std::set<NODE_NAME> value;
					value.insert(n2_id);
					tmp.insert(std::make_pair(n1_id, value));
				}
				else
					it_tmp->second.insert(n2_id);
				////////////////////////////////////////
				it_tmp = tmp.find(n2_id);
				if (it_tmp == tmp.end())
				{
					std::set<NODE_NAME> value;
					value.insert(n1_id);
					tmp.insert(std::make_pair(n2_id, value));
				}
				else
					it_tmp->second.insert(n1_id);

			}
		}
		/// <summary>
		/// Init variable
		/// </summary>
		graph.numNodes = tmp.size();
		graph.numEdges = 0;
		graph.edges.resize(graph.numNodes);

		it_tmp = tmp.begin();

		for (; it_tmp != tmp.end(); it_tmp++)
		{
			for (NODE_ID nbr : it_tmp->second)
				graph.edges[it_tmp->first].insert(std::make_pair(nbr, -1));
			graph.numEdges += it_tmp->second.size();
		}
		graph.numEdges = graph.numEdges / 2;
		inFile.close();
	}
	return true;

}
bool File::SaveClustering(string fName, const MAP_NAME_ID& nodeNames, const CLUSTERING& cl)
{
	std::ofstream outFile{ fName,std::ios::out };
	if (!outFile) {
		std::cerr << "File could not open " << fName << std::endl;
		//exit(EXIT_FAILURE);
		return false;
	}
	// Create map NODE_ID to NODE_NAME
	std::map<NODE_ID, NODE_NAME> listNodeID;
	MAP_NAME_ID::const_iterator it_listNodeName = nodeNames.begin();
	for (; it_listNodeName != nodeNames.end(); it_listNodeName++)
		listNodeID.insert(std::make_pair(it_listNodeName->second, it_listNodeName->first));
	CLUSTERING::const_iterator it_clustering = cl.begin();
	for (; it_clustering != cl.end(); it_clustering++)
	{
		CLUSTER::iterator it_cluster = it_clustering->begin();
		for (; it_cluster != it_clustering->end(); it_cluster++)
		{
			std::map<NODE_ID, NODE_NAME>::iterator it_listNodeID = listNodeID.find(*it_cluster);
			if (it_listNodeID != listNodeID.end())
				outFile << it_listNodeID->second << "\t";
		}
		outFile << "\n";
	}

	return true;
}