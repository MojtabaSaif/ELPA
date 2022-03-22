#include"Alg.h"
void HALPA::Weighting(OWG& gra)
{
	for (int src = 0; src < gra.edges.size(); src++)
	{
		map<NODE_ID, WEIGHT>::iterator it_dst = gra.edges[src].begin();
		for (; it_dst != gra.edges[src].end(); it_dst++)
		{
			if (it_dst->second != -1) continue;

			// Find Common neighbor (src,dst)
			std::set<NODE_ID> cmnNbrs;
			for (pair<NODE_ID, WEIGHT> nbr_src : gra.edges[src])
				if (gra.edges[it_dst->first].find(nbr_src.first) != gra.edges[it_dst->first].end())
					cmnNbrs.insert(nbr_src.first);

			//
			int m4 = 0;
			std::set<NODE_ID>::iterator it1 = cmnNbrs.begin();
			for (; it1 != cmnNbrs.end(); it1++)
			{
				std::set<NODE_ID>::iterator it2 = it1;
				it2++;
				for (; it2 != cmnNbrs.end(); it2++)
					if (gra.edges[*it1].find(*it2) != gra.edges[*it1].end())
						m4++;
			}
			int w = cmnNbrs.size() + m4;
			it_dst->second = w;
			gra.edges[it_dst->first].find(src)->second = w;
		}
	}
	// Calculating influance of nodes
	gra.inflNodes.resize(gra.numNodes);
	gra.capecity.resize(gra.numNodes);
	for (int n = 0; n < gra.numNodes; n++)
	{
		WEIGHT influ = 0.0, capecity = -1;
		map<NODE_ID, WEIGHT>::iterator it = gra.edges[n].begin();
		for (; it != gra.edges[n].end(); it++)
		{
			influ += it->second;
			if (it->second > capecity)
				capecity = it->second;
		}
		gra.inflNodes[n] = influ;//(gra.edges[n].size() != 0) ? influ / static_cast<double>(gra.edges[n].size()) : 0.0;
		gra.capecity[n]  = capecity;
	}
}
void HALPA::Run(const OWG& gra, MBR_CLUSTERING& mbrCL)
{
	vector<MEMBERSHIPE> maxMbr(gra.numNodes,0.0);
	mbrCL.clear();
	mbrCL.resize(gra.numNodes);

	LABEL_ID label = 0;
	while (true)
	{
		vector<SEED> seed;
		seed.resize(1);
		seed[0].first = -1;
		seed[0].second.first = label;
		seed[0].second.second = -1;
		double maxCap = 0.0;
		for (int n = 0; n < gra.numNodes; n++) //Find node with max influance 
		{
			if ((gra.inflNodes[n] > seed[0].second.second) && (mbrCL[n].size() == 0))
			{
				seed[0].first = n;
				seed[0].second.second = gra.inflNodes[n];
				maxCap = gra.capecity[n];
			}
		}
		if ((seed[0].first == -1) || (seed[0].second.second == 0)) break;
		//std::cout << "Seed " << seed[0].first << "\n";
		seed[0].second.second = maxCap;
		Run(gra, seed, mbrCL,maxMbr);
		Normalizing(mbrCL);
		label++;
	}
}

void HALPA::Run(const OWG& gra, vector<SEED> seeds, MBR_CLUSTERING& mbrCL, vector<MEMBERSHIPE>& maxMbrNodes)
{
	//static int label = 0;
		
	
	for (int s = 0; s < seeds.size(); s++)
	{
		PACKET p;
		set<NODE_ID> nComm;
		MEMBERSHIPE maxMbr, minMbr;
		p = seeds[s].second;
		// The membershipe of seed is set max edige weight
		//p.second = -1;// set defult max=-1
		//for (pair<NODE_ID, WEIGHT> n : gra.edges[seeds[s].first])
		//	if (n.second > p.second)
		//		p.second = n.second;
		// 
		mbrCL[seeds[s].first].push_back(p);
		std::deque<NODE_ID> list;
		list.push_back(seeds[s].first);
		nComm.insert(seeds[s].first);
		maxMbr = p.second;
		minMbr = p.second;
		while (!list.empty())
		{
			NODE_ID nID = list.front();
			list.pop_front();
			p.second = mbrCL[nID].Find(p.first)->second;
			for (pair<NODE_ID, WEIGHT> nbr : gra.edges[nID])
			{
				MEMBERSHIPE tmpMBR = fmin(nbr.second, p.second);
				if (tmpMBR< (gra.capecity[nbr.first] / 3.0)) continue;
				//Penalty
				//if (tmpMBR < nbr.second) continue; 
				//if (tmpMBR < nbr.second) tmpMBR -= 1;
				//if (tmpMBR < 0.00) continue;
				NodeLabeles::iterator it = mbrCL[nbr.first].Find(p.first);
				if (it == mbrCL[nbr.first].end())
				{
					mbrCL[nbr.first].push_back(std::make_pair(p.first, tmpMBR));
					if (tmpMBR > maxMbr) maxMbr = tmpMBR;
					else if (tmpMBR < minMbr) minMbr = tmpMBR;
					nComm.insert(nbr.first);
					if (gra.capecity[nbr.first]!=0)
						list.push_back(nbr.first);
				}
				else
					if (tmpMBR > it->second)
					{
						it->second = tmpMBR;
						if (tmpMBR > maxMbr) maxMbr = tmpMBR;
						else if (tmpMBR < minMbr) minMbr = tmpMBR;
						list.push_back(nbr.first);
					}
			}

		}
		vector<pair<NODE_ID, MEMBERSHIPE>> tmpNodeExt;
		for (NODE_ID n : nComm)
		{
			NodeLabeles::iterator it = mbrCL[n].Find(p.first);
			it->second =(maxMbr==minMbr)?1.0: (it->second) / static_cast<double>(maxMbr);
			// update maxMBr here
			if (it->second > maxMbrNodes[n]) maxMbrNodes[n] = it->second;
		}
		
	}
}
void HALPA::Normalizing(MBR_CLUSTERING& mbrCL)
{
	
	vector<MEMBERSHIPE> maxNode(mbrCL.size(), -1);
	
	for (int n = 0; n < mbrCL.size(); n++)
	{
		NodeLabeles::iterator item = mbrCL[n].begin();
		for (; item != mbrCL[n].end(); item++)
		{
			if (item->second > maxNode[n])
				maxNode[n] = item->second;
		}
	}
	for (int n = 0; n < mbrCL.size(); n++)
	{
		NodeLabeles::iterator item = mbrCL[n].begin();
		for (; item != mbrCL[n].end(); )
		{
			MEMBERSHIPE mbr = (maxNode[n] > 0) ? pow(item->second / maxNode[n], 6.0) : 0.0;
			if (mbr < 0.5)
				item = mbrCL[n].erase(item);
			else
				item++;
		}
	}
}

void LPA::Run(OWG& gra, CLUSTERING& cl)
{
	vector<LABEL_ID> labels(gra.numNodes);
	for (int i = 0; i < gra.numNodes; i++)
		labels[i] = i;
	bool flgChange = true;
	int itr_count = 0;
	while (flgChange)
	{
		flgChange = false;
		for (int n = 0; n < gra.numNodes; n++)
		{
			LABEL_ID maxL = GetMaxLabe(gra, n, labels);
			if (labels[n] != maxL)
			{
				labels[n] = maxL;
				flgChange = true;
			}
		}
		itr_count++;
	}// End of LPA
	map<LABEL_ID, CLUSTER> clutering;
	for (int n = 0; n < labels.size(); n++)
	{
		map<LABEL_ID, CLUSTER>::iterator it = clutering.find(labels[n]);
		if (it != clutering.end())
			it->second.insert(n);
		else
		{
			CLUSTER tmp;
			tmp.insert(n);
			clutering.insert(std::make_pair(labels[n],tmp));
		}	
	}
	std::cout <<"Number of iteration\t" << itr_count << "\n";
	for (pair<LABEL_ID, CLUSTER> c : clutering)
		cl.push_back(c.second);
}

LABEL_ID LPA::GetMaxLabe(OWG& gra, NODE_ID n, const vector<LABEL_ID>& labels)
{
	map<NODE_ID, WEIGHT> nbr = gra.edges[n];
	map<LABEL_ID, int> tmpL;
	LABEL_ID maxL_ID=-1;
	int      maxL_count = -1;
	for (pair<NODE_ID, WEIGHT> n : nbr)
	{
		map<LABEL_ID, int>::iterator it = tmpL.find(labels[n.first]);
		if (it != tmpL.end())
			it->second++;
		else
			it = tmpL.insert(std::make_pair(labels[n.first], 1)).first;
		if (it->second > maxL_count)
		{
			maxL_ID = labels[n.first];
			maxL_count = it->second;
		}
		
	}
	vector<LABEL_ID> ids;
	for (pair<LABEL_ID, int> item : tmpL)
		if (item.second == maxL_count)
			ids.push_back(item.first);

	return ids[Utilities::RandomInt(0,ids.size()-1)];
}
