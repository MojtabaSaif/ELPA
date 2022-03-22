#pragma once
#include<set>
#include<iostream>
#include<deque>
#include"Type.h"
#include"Alg.h"
#include"Utilities.h"
using std::pair;
namespace HALPA {
	typedef pair<NODE_ID, PACKET> SEED;
	void Weighting(OWG& gra);
	//void propagation(OWG& gra, vector<SEED> seeds, MEMBERSHIP_CLUSTERING& mbrCL);
	void Run(const OWG& gra,vector<SEED> seeds, MBR_CLUSTERING& mbrCL,vector<MEMBERSHIPE>& maxMbr);
	void Run(const OWG& gra, MBR_CLUSTERING& mbrCL);
	void Normalizing(MBR_CLUSTERING& mbrCL);
	
}
namespace LPA {
	void Run(OWG& gra,CLUSTERING& cl);
	LABEL_ID GetMaxLabe(OWG& gra,NODE_ID n,const vector<LABEL_ID>& labels);
}

