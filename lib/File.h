#pragma once
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <map>
#include <set>
#include"Type.h"

namespace File
{
	using std::string;
	using std::map;
	using std::set;
	enum FileType
	{
		EDGE_LIST,
		CLUSTERING_FILE
	};
	
	bool LoadClustering(string fName,const MAP_NAME_ID& gra,CLUSTERING& cl);
	bool LoadFileEdgeList(string fName, OWG& graph);
	
	//bool SaveEdgeListFile(string fName, FileType type);
	bool SaveClustering(string fName, const MAP_NAME_ID& gra,const CLUSTERING& cl);
}
