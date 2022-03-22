#pragma once
#include <iostream>
#include<deque>
#include<algorithm>
#include<chrono>
#include<random>
#include"Type.h"
#include"..\NMI_Package\mutual.h"
namespace Utilities 
{
	/// <summary>
	/// Random function
	/// </summary>
	
	double RandomDouble(double s, double e);
	int    RandomInt(int s, int e);
	
	void PrintClustering(const CLUSTERING& cl);
	void PrintClustering(const CLUSTERING& cl,const MAP_NAME_ID& nodeNames);
	void MbrClusteringToClustering(const MBR_CLUSTERING& mbrCL,CLUSTERING& cl);
	void ClusteringToDeque(const CLUSTERING& cl, std::deque<std::deque<int>>& dCLU);
	double OverlappingModularity(const OWG& gra,const MBR_CLUSTERING& mbrCL);
	void UpdateClusteringWithNodeName(CLUSTERING& cl, const MAP_NAME_ID& nodeNames);
	void IntersectionTwoCluster(const CLUSTER& cl1, const CLUSTER& cl2, CLUSTER& interCL);
	double F1Score(const CLUSTERING& trueComms, const CLUSTERING& preComms);
	void PrintInfoClustering(const MBR_CLUSTERING& mbrCl);
	void CalCriteria(const OWG& gra,const MBR_CLUSTERING& cl, const CLUSTERING& gt_cl, double& nmi, double& q, double& f1);
}


