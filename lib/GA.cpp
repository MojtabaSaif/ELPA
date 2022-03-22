#include "GA.h"
using namespace  EA_Alg;
std::string OptionGA::ToString()
{
	std::string str;
	str = "<OptionGA>\n";
	str += "lenSol\t" + std::to_string(lenSol) + "\n";
	str += "pSeed\t" + std::to_string(pSeeds) + "\n";
	str += "pCrossOver\t" + std::to_string(pCrossOver) + "\n";
	str += "pMutation\t" + std::to_string(pMutation) + "\n";
	str += "popSize\t" + std::to_string(popSize) + "\n";
	//str += "numComm\t" + std::to_string(numComm) + "\n";
	str += "MaxGenCount\t" + std::to_string(MaxGenCount) + "\n";
	str += "TournamentSize\t" + std::to_string(TournamentSize) + "\n";
	str += "Max Thread\t" + std::to_string(maxThread) + "\n";
	str += "\n</OptionGA>";
	return str;
}
void  OptionGA::setFromFile(std::string fName) 
{
	
	std::ifstream inFile{ fName,std::ios::in };
	if (!inFile) {
		std::cerr << "File could not open " << fName << std::endl;
		return;
	}

	std::string key;
	std::string value;
	while (inFile >> key >> value) {

		std::transform(key.begin(), key.end(), key.begin(), ::tolower);
		std::transform(value.begin(), value.end(), value.begin(), ::tolower);
		if (key == "lensol")				lenSol = std::stoi(value);
		else if ((key == "pseeds"))			pSeeds = std::stod(value);
		else if ((key == "pcrossover"))		pCrossOver = std::stod(value);
		else if ((key == "pmutation"))		pMutation = std::stod(value);
		else if ((key == "popsize"))		popSize = std::stoi(value);
		//else if ((key == "numcomm"))		numComm = std::stoi(value);
		else if ((key == "maxgencount"))    MaxGenCount = std::stoi(value);
		else if ((key == "tournamentsize")) TournamentSize = std::stoi(value);
		else if ((key == "maxthread"))      maxThread = std::stoi(value);
		else if ((key == "showreprot"))		ShowReprot = (value == "true") ? true : false;
		else if ((key == "multi"))			Multi = (value == "true") ? true : false;
	}
}
void OptionGA::SetMaxThread(int v)
{
	maxThread = v;
}

void EA_Alg::GA::CreateChromosome(CHROMOSOME& chr, FITNESS& fit)
{
	
	vector<bool> flg_label(ga_gra.numNodes, false);
	
	MBR_CLUSTERING mbrCL;
	vector<MEMBERSHIPE> maxMbr;
	vector<pair<double, NODE_ID>> tmpCumP = cumP;
	mbrCL.resize(ga_gra.numNodes);
	maxMbr.resize(ga_gra.numNodes, 0);
	
	while (tmpCumP.size() >0)
	{
		NODE_ID selectN = -1;
		while ((selectN==-1)&&(tmpCumP.size()>0))
		{
			double p=RandomP();
			int it_tmpCumP = 0;
			while ((it_tmpCumP < tmpCumP.size()) && (p > tmpCumP[it_tmpCumP].first) )
				it_tmpCumP++;
			if (it_tmpCumP >= tmpCumP.size()) continue;
			//while ((id > 0) && flg_label[id])
			//	id--;
			NODE_ID id = tmpCumP[it_tmpCumP].second;
			if (!flg_label[id])
			{
				if ((mbrCL[id].size() == 0))
					selectN = id;
				tmpCumP.erase(tmpCumP.begin() + it_tmpCumP);
				flg_label[id] = true;
				
			}
			
		}
		if (selectN == -1) continue;
		HALPA::SEED seed;
		
		seed.first = selectN;
		seed.second.first = GenerateRandomLabel();
		seed.second.second = ga_gra.inflNodes[selectN];
		vector<HALPA::SEED> ss;
		ss.push_back(seed);
		chr.push_back(seed);
		
	    HALPA::Run(ga_gra, ss, mbrCL,maxMbr);	
		HALPA::Normalizing(mbrCL);
		
	}
	
	fit.value.push_back(1-Utilities::OverlappingModularity(ga_gra,mbrCL));
	if (opt.Multi)
	{
		double sIN, sOUT;
		SumEdges(ga_gra, mbrCL, sIN, sOUT);
		fit.value.push_back(sOUT);
	}
		
}

void EA_Alg::GA::CreateChromosome2(CHROMOSOME& chr, FITNESS& fit)
{
	MBR_CLUSTERING mbrCL(ga_gra.numNodes);
	vector<MEMBERSHIPE> maxMbr(ga_gra.numNodes, 0.0);
	//vector<NODE_ID> nodes(ga_gra.numNodes);
	vector<pair<double, NODE_ID>> tmpCumD = cumP;
	//for (int i = 0; i < nodes.size(); i++)
	//	nodes[i] = i;
	//nodes = gSeeds;
	while (tmpCumD.size()>0)
	{
		NODE_ID selectN = -1;
		
		while ((selectN == -1))
		{
			double p = RandomP();
			int it_CumP = 0;
			while ((it_CumP < tmpCumD.size()) && (p > tmpCumD[it_CumP].first))
				it_CumP++;
			if (it_CumP >= tmpCumD.size()) it_CumP=tmpCumD.size()-1;
			//if (maxMbr[it_CumP] > 0.5) continue;
			int i = 0;
			selectN = tmpCumD[it_CumP].second;
			while ((i < chr.size()))
			{
				if ((chr[i].first == selectN))
					selectN = -1;
				i++;
			}
		}
		
		HALPA::SEED seed;
		seed.first = selectN;
		seed.second.first = GenerateRandomLabel();
		seed.second.second = ga_gra.inflNodes[selectN];
		vector<HALPA::SEED> ss;
		ss.push_back(seed);
		chr.push_back(seed);

		HALPA::Run(ga_gra, ss, mbrCL, maxMbr);
		HALPA::Normalizing(mbrCL);

		/*vector<NODE_ID>::iterator it_nodes = nodes.begin();
		for (; it_nodes != nodes.end(); )
			if (maxMbr[*it_nodes] > 0.5)
				it_nodes = nodes.erase(it_nodes);
			else
				it_nodes++;*/
		vector<pair<double, NODE_ID>>::iterator it = tmpCumD.begin();
		for(;it!=tmpCumD.end();)
			if (mbrCL[it->second].size()!=0)
				it = tmpCumD.erase(it);
			else
				it++;
	}

	//Utilities::PrintInfoClustering(mbrCL);
	//std::cout << "\n";

	fit.value.push_back(1 - Utilities::OverlappingModularity(ga_gra, mbrCL));
	if (opt.Multi)
	{
		double sIN, sOUT;
		SumEdges(ga_gra, mbrCL, sIN, sOUT);
		fit.value.push_back(sOUT);
	}

}

void EA_Alg::GA::CreateChromosome3(CHROMOSOME& chr, FITNESS& fit)
{
	MBR_CLUSTERING mbrCL(ga_gra.numNodes);
	vector<MEMBERSHIPE> maxMbr(ga_gra.numNodes, 0.0);
	vector<NODE_ID> nodes(ga_gra.numNodes);
	//for (int i = 0; i < nodes.size(); i++)
	//	nodes[i] = i;
	nodes = gSeeds;
	while (nodes.size() > 0)
	{
		NODE_ID selectN = -1;
		
		while ((selectN == -1)&& (nodes.size()>0))
		{
			selectN =nodes[Utilities::RandomInt(0, nodes.size() - 1)];
			//selectN = RandomInt(0, ga_gra.numNodes - 1);
			
			int i = 0;
			while ((selectN!=-1)&&(i < chr.size()))
			{
				if ((chr[i].first == selectN))
				{
					selectN = -1;
					break;
				}
				i++;
			}
		}

		
		HALPA::SEED seed;
		seed.first = selectN;
		seed.second.first = GenerateRandomLabel();
		seed.second.second = ga_gra.capecity[selectN];
		vector<HALPA::SEED> ss;
		ss.push_back(seed);
		chr.push_back(seed);

		HALPA::Run(ga_gra, ss, mbrCL, maxMbr);
		/*int count = 0;
		for (int i = 0; i < mbrCL.size(); i++)
			if (mbrCL[i].Find(seed.second.first) != mbrCL[i].end())
				count++;
		std::cout <<"A: " << count << "\n\n\n";
		if (count == 1)
		{
			map<NODE_ID, WEIGHT>::const_iterator it = ga_gra.edges[selectN].begin();
			while ((it != ga_gra.edges[selectN].end()))
			{
				std::cout << ga_gra.capecity[selectN] << "\n";
				std::cout << it->second << "\n";
				std::cout << (ga_gra.capecity[it->first] / 2.0) << "\n";
				std::cout << "----------------\n";
				if ((ga_gra.capecity[selectN] >= it->second) && (it->second >= (ga_gra.capecity[it->first] / 2.0)))
					break;
				it++;
			}
			std::cout <<"D =\t" << ga_gra.edges[selectN].size() << "\n";
			std::cout << "Influ =\t" << ga_gra.inflNodes[selectN] << "\n";
			std::cout << "Capecity =\t" << ga_gra.capecity[selectN] << "\n";
		}*/
		HALPA::Normalizing(mbrCL);
		/*count = 0;
		for (int i = 0; i < mbrCL.size(); i++)
			if (mbrCL[i].Find(seed.second.first) != mbrCL[i].end())
				count++;
		std::cout <<"B: " << count << "\n\n\n";*/
		vector<NODE_ID>::iterator it_nodes = nodes.begin();
		for (; it_nodes != nodes.end(); )
			if (maxMbr[*it_nodes]>=0.7)
				it_nodes = nodes.erase(it_nodes);
			else
				it_nodes++;
	}
	Utilities::PrintInfoClustering(mbrCL);
	std::cout << "\n\n\n";
	
	fit.value.push_back(1 - Utilities::OverlappingModularity(ga_gra, mbrCL));
	if (opt.Multi)
	{
		double sIN, sOUT;
		SumEdges(ga_gra, mbrCL, sIN, sOUT);
		fit.value.push_back(sOUT / (sIN + sOUT));
	}
}

void EA_Alg::GA::CreateChromosome4(CHROMOSOME& chr, FITNESS& fit)
{
	MBR_CLUSTERING mbrCL(ga_gra.numNodes);
	vector<MEMBERSHIPE> maxMbr(ga_gra.numNodes, 0.0);
	vector<pair<double, NODE_ID>> tmpCumD = cumP;
	vector<HALPA::SEED> ss;
	int num = Utilities::RandomInt(2, 20);

	while (num > 0)
	{
		NODE_ID selectN = -1;

		while ((selectN == -1))
		{
			double p = RandomP();
			int it_CumP = 0;
			while ((it_CumP < tmpCumD.size()) && (p > tmpCumD[it_CumP].first))
				it_CumP++;
			if (it_CumP >= tmpCumD.size()) continue;
			//if (maxMbr[it_CumP] > 0.5) continue;
			int i = 0;
			selectN = tmpCumD[it_CumP].second;
			while ((i < chr.size()))
			{
				if ((chr[i].first == selectN))
				{
					selectN = -1;
					break;
				}
				i++;
			}
		}
		
		HALPA::SEED seed;
		seed.first = selectN;
		seed.second.first = GenerateRandomLabel();
		seed.second.second = ga_gra.inflNodes[selectN];
		ss.push_back(seed);
		chr.push_back(seed);
		num--;
	}

	HALPA::Run(ga_gra, ss, mbrCL, maxMbr);
	HALPA::Normalizing(mbrCL);

	Utilities::PrintInfoClustering(mbrCL);
	std::cout << "\n";

	fit.value.push_back(1 - Utilities::OverlappingModularity(ga_gra, mbrCL));
	if (opt.Multi)
	{
		double sIN, sOUT;
		SumEdges(ga_gra, mbrCL, sIN, sOUT);
		fit.value.push_back(sOUT);
	}


}

POP_ITR EA_Alg::GA::TournamentSelection(int size)
{
	POP_ITR index=pop.begin();
	std::advance(index, Utilities::RandomInt(0, opt.popSize - 1));
	for (int i = 0; i < size; i++)
	{
		int r= Utilities::RandomInt(0, opt.popSize - 1);
		POP_ITR it = pop.begin();
		std::advance(it, r);
		if (it->first < index->first)
			index = it;
	}
	return index;
}

void EA_Alg::GA::SinglePointCrossOver(const CHROMOSOME& p1, const CHROMOSOME& p2, CHROMOSOME& child)
{
	child.clear();
	int cutPoint = Utilities::RandomInt(0, std::min(p1.size(), p2.size())-2);
	int i = 0;
	for ( ;i < cutPoint; i++)
		child.push_back(p1[i]);
	int j = p2.size()-1;
	while ((j >= 0) && (child.size() < p2.size()))
	{
		int it_ch = 0;
		while ((it_ch < child.size()) && (child[it_ch].first != p2[j].first)) it_ch++;
		if (it_ch >= child.size()) child.push_back(p2[j]);
		j--;
	}
	if (child.size() < p1.size())
	{
		int j = p1.size()-1;
		while ((j >= 0) && (child.size() < p1.size()))
		{
			int it_ch = 0;
			while ((it_ch < child.size())&& (child[it_ch].first != p1[j].first)) it_ch++;
			if (it_ch >= child.size()) child.push_back(p1[j]);
			j--;
		}
	}
}

void EA_Alg::GA::Mutation(const CHROMOSOME& ch_in, CHROMOSOME& ch_out)
{
	ch_out = ch_in;
	ch_out.erase(ch_out.begin() + Utilities::RandomInt(0, ch_out.size() - 1));
	
	NODE_ID node = -1;
	while (node==-1)
	{
		node= Utilities::RandomInt(0, ga_gra.numNodes - 1);
		for (HALPA::SEED n : ch_out)
			if (n.first == node)
			{
				node = -1;
				break;
			}
	}
	HALPA::SEED s;
	s.first = node;
	s.second.first = GenerateRandomLabel();
	s.second.second = ga_gra.inflNodes[s.first];
	ch_out.push_back(s);
}

void EA_Alg::GA::Fitness(const CHROMOSOME& chr, FITNESS& outValue)
{
	MBR_CLUSTERING mbrCL;
	mbrCL.resize(ga_gra.numNodes);
	vector<MEMBERSHIPE> maxMbr;
	maxMbr.resize(ga_gra.numNodes, 0);
	HALPA::Run(ga_gra, chr, mbrCL,maxMbr);
	HALPA::Normalizing(mbrCL);
	outValue.value.clear();
	outValue.value.push_back(1-Utilities::OverlappingModularity(ga_gra, mbrCL));
	if (opt.Multi)
	{
		double sIN, sOUT;
		SumEdges(ga_gra, mbrCL, sIN, sOUT);
		outValue.value.push_back(sOUT/ (sIN+sOUT));
	}
}

LABEL_ID EA_Alg::GA::GenerateRandomLabel()
{
	return Utilities::RandomInt(0,ga_gra.numNodes-1);
}

bool EA_Alg::GA::Dominates(const FITNESS& f1, const FITNESS& f2)
{
	int nObj = f1.value.size();
	bool allCheck = true;
	bool anyCheck = false;

	for (int i = 0; i < nObj; i++)
	{
		if (f1.value[i] <= f2.value[i]) {
			if (f1.value[i] < f2.value[i])
				anyCheck = true;
		}
		else
		{
			allCheck = false;
			break;
		}
	}
	return (allCheck && anyCheck);
}

void EA_Alg::GA::NoneDominateedSorting(POPULATION& pop, vector<FRONT>& fs)
{
	int nPOP = pop.size();
	FRONT f;
	fs.clear();
	POP_ITR it_pop = pop.begin();
	for (; it_pop!=pop.end(); it_pop++)
	{
		it_pop->first.cd = 0.0;
		it_pop->first.np = 0;
		it_pop->first.rank = -1;
		it_pop->first.SP.clear();
	}
	POP_ITR it_pop1 = pop.begin();
	for (; it_pop1 != pop.end(); it_pop1++)
	{
		POP_ITR it_pop2 = it_pop1;
		it_pop2++;
		for (;it_pop2!=pop.end();it_pop2++)
		{
			if (Dominates(it_pop1->first, it_pop2->first))
			{
				it_pop1->first.SP.push_back(it_pop2);
				it_pop2->first.np++;
			}
			else
				if (Dominates(it_pop2->first, it_pop1->first))
				{
					it_pop2->first.SP.push_back(it_pop1);
					it_pop1->first.np++;
				}
		}
		if (it_pop1->first.np == 0)
		{
			f.push_back(it_pop1);
			it_pop1->first.rank = 1;
		}
	}
	fs.push_back(f);
	int k = 0;
	while (true)
	{
		FRONT tmpFront;
		FRONT::const_iterator p = fs[k].begin();
		for (; p != fs[k].end(); p++)
		{
			FRONT::const_iterator q = (* p)->first.SP.begin();
			for (; q != (*p)->first.SP.end(); q++)
			{
				(* q)->first.np--;
				if ((* q)->first.np == 0)
				{
					tmpFront.push_back(( * q));
					(* q)->first.rank = k + 1;
				}
			}
		}
		if (tmpFront.size() == 0)
			break;
		fs.push_back(tmpFront);
		k++;
		tmpFront.clear();
	}
}
struct CASTS { POP_ITR id; double value=0.0; };
bool CompCasts(CASTS c1, CASTS c2) { return (c1.value < c2.value); }
void EA_Alg::GA::CalcCrowdingDistance(POPULATION& pop, vector<FRONT>& fs)
{
	int nObj = pop[0].first.value.size();
	std::vector<FRONT>::const_iterator k = fs.begin();
	for (; k != fs.end(); k++)
	{
		for (int o = 0; o < nObj; o++)
		{
			std::vector<CASTS> Casts;
			FRONT::const_iterator F = k->begin();
			for (; F != k->end(); F++)
			{
				CASTS ca;
				ca.id = *F;
				ca.value = (*F)->first.value[o];
				Casts.push_back(ca);
			}
			std::sort(Casts.begin(), Casts.end(), CompCasts);
			Casts[0].id->first.cd += INFINITY;
			Casts[Casts.size() - 1].id->first.cd += INFINITY;
			double p = std::abs(Casts[0].value - Casts[Casts.size() - 1].value);
			for (int i = 1; i < Casts.size() - 1; i++)
				Casts[i].id->first.cd += std::abs(Casts[i - 1].value - Casts[i + 1].value);
		}
	}
}

void EA_Alg::GA::SumEdges(const OWG& gra,  MBR_CLUSTERING& mbrCL, double& sumIN, double& sumOUT)
{
	sumIN  = 0.0;
	sumOUT = 0.0;
	for (NODE_ID u = 0; u < gra.edges.size(); u++)
	{
		for (pair<NODE_ID, WEIGHT> v:gra.edges[u])
		{
			if (u > v.first) continue;
			bool f = false;
			for(pair<LABEL_ID, MEMBERSHIPE> label:mbrCL[u])
				if (mbrCL[v.first].Find(label.first) != mbrCL[v.first].end())
				{
					f = true;
					sumIN += v.second;
					break;
				}
			if (!f)
				sumOUT += v.second;
		}
	}
}

void EA_Alg::GA::CreatePOP()
{
	// Fill gSeeds
	for (int s = 0; s < ga_gra.numNodes; s++)
	{
		if (ga_gra.capecity[s] == 0) continue;
		map<NODE_ID, WEIGHT>::const_iterator it = ga_gra.edges[s].begin();
		while ((it != ga_gra.edges[s].end()))
		{
			if ((it->second > 0.0)&&(ga_gra.capecity[s] >= it->second) && (it->second >= (ga_gra.capecity[it->first] / 2.0)))
			{
				gSeeds.push_back(s);
				break;
			}		
			it++;
		}
	}
	//calculate cumulative probability
	cumP.clear();
	std::multimap<double, NODE_ID,std::greater<double>> influ;
	
	//for (int n = 0; n < ga_gra.numNodes; n++)
	//	if(ga_gra.inflNodes[n]>0)
	//		influ.insert(std::make_pair(ga_gra.inflNodes[n], n));
	vector<NODE_ID>::iterator it_gSeed = gSeeds.begin();
	for(;it_gSeed!=gSeeds.end();it_gSeed++)
		influ.insert(std::make_pair(ga_gra.inflNodes[*it_gSeed], *it_gSeed));

	double sumINFL = 0,maxInflu=-1;
	std::multimap<double, NODE_ID, std::greater<double>>::iterator it_influ = influ.begin();
	for (;it_influ != influ.end();it_influ++)
	{
		sumINFL += it_influ->first;
		if (it_influ->first > maxInflu) maxInflu = it_influ->first;
	}
	double totalP = 0;
	it_influ = influ.begin();
	for (;it_influ != influ.end();it_influ++)
	{
		totalP += static_cast<double>(it_influ->first) / sumINFL;//((static_cast<double>(it_influ->first) / maxInflu)>0.1) ? static_cast<double>(it_influ->first) / sumINFL : 0.0;
		cumP.push_back(std::make_pair(totalP,it_influ->second));
	}

	
	pop.resize(opt.popSize);
	int i = 0;
	while(i < opt.popSize)
	{
		//CreateChromosome(tmpPOP[i].second, tmpPOP[i].first);
		vector<thread> th_pool;
		for (int j = 0; (j < opt.maxThread)&&(i< opt.popSize); j++)
		{
			th_pool.push_back(thread(&GA::CreateChromosome2, this, std::ref(pop[i].second), std::ref(pop[i].first)));
			i++;
		}

		for (thread& th : th_pool)
			if (th.joinable())
				th.join();
		//std::cout << "Create pop " << i << "\n";
	}
	
	//std::sort(pop.begin(), pop.end(), CMP_Fit);
	// NonDominate sorting on POP
	if (opt.Multi)
	{
		std::vector<FRONT> fs;
		NoneDominateedSorting(pop, fs);
		CalcCrowdingDistance(pop, fs);
	}
}


void EA_Alg::GA::ManiLoop()
{
	int gen = 0;
	
	while (gen < opt.MaxGenCount)
	{
		//POPULATION tmpPOP; // This variable is for new chromosomes that obtain from cross over and mutation
		

		if (opt.ShowReprot)
			std::cout << "Start Gen: " << gen + 1 ;
		//Cross over
		if (opt.ShowReprot)
			std::cout << "\t\tCross over statrt \n";
		int numCross = floor(opt.pCrossOver * opt.popSize);
		
		
		//tmpPOP.resize(numCross);
		for (int i = 0; i < numCross; i++)
		{
			POP_ITR p1= TournamentSelection(opt.TournamentSize);
			POP_ITR p2 = TournamentSelection(opt.TournamentSize);
			CHROMOSOME chr;
			SinglePointCrossOver(p1->second, p2->second,chr);
			pop.push_back(std::make_pair(FITNESS(),chr));
			
		}
		
		if (opt.ShowReprot)
			std::cout << "\t\tCross over end \n";
		//Mutation
		if (opt.ShowReprot)
			std::cout << "\t\tMutation start \n";
		//tmpPOP.clear();
		int numMute = floor(opt.pMutation * opt.popSize);
		
		for (int j = 0; j < numMute; j++)
		{
			POP_ITR it_pop = pop.begin();
			std::advance(it_pop, Utilities::RandomInt(0, opt.popSize - (j+1)));
			CHROMOSOME tmpChr;
			Mutation(it_pop->second, tmpChr);
			pop.erase(it_pop);
			pop.push_back(std::make_pair(FITNESS(), tmpChr));
		}
		if (opt.ShowReprot)
			std::cout << "\t\tMutation End \n";
		
		// survival selection
		if (opt.ShowReprot)
			std::cout << "\t\tsurvival selection start \n";
		// Add new chromosomes to pop
				
		unsigned int ii = pop.size()-(numCross+numMute)-1;
		while (ii < pop.size())
		{
			vector<thread> th_pool;
			for (int j = 0; (j < opt.maxThread) && (ii < pop.size()); j++)
			{
				//Fitness(tmpPOP[ii].second, tmpPOP[ii].first);
				th_pool.push_back(thread(&GA::Fitness, this, std::ref(pop[ii].second), std::ref(pop[ii].first))); // Start thread
				ii++;
			}
			for (thread& th : th_pool)
				if (th.joinable())
					th.join(); // End thread
		}
		
				
		// NonDominate sorting on POP
		if (opt.Multi)
		{
			std::vector<FRONT> fs;
			NoneDominateedSorting(pop, fs);
			CalcCrowdingDistance(pop, fs);
			std::sort(pop.begin(), pop.end(), CMP_CrowdingDistance);
			std::sort(pop.begin(), pop.end(), CMP_Rank);
		}
		else
			std::sort(pop.begin(), pop.end(), CMP_Fit);

		// Eliteism
		
		while ((pop.size()>opt.popSize))
		{
			int it_p = pop.size() - 1;
			while ((pop.size() > opt.popSize)&&(it_p>=0))
			{
				if (Utilities::RandomDouble(0.0, 1.0) > 0.5)
					pop.erase(pop.begin() + it_p);
				it_p--;
			}
		}
		/*POPULATION::reverse_iterator it_pop = pop.rbegin();
		while ((pop.size() > opt.popSize)&&(it_pop!=pop.rend()))
		{
			if (Utilities::RandomDouble(0.0, 1.0) > 0.5)
				it_pop = std::make_reverse_iterator( pop.erase( (it_pop+1).base()));
			else
				it_pop++;
		}*/
		pop.erase(pop.begin() + opt.popSize, pop.end());

		if (opt.ShowReprot)
			std::cout << "\t\tsurvival selection End \n";

		
		genReport += "Best ( " + std::to_string(gen) + ") ";
		for (double v : pop.begin()->first.value)
			genReport += "\t" + std::to_string(v) + "\t";
		genReport += "\n\n";
		if (opt.ShowReprot)
		{
			POP_ITR itt = pop.begin();
			std::cout << "\t\t Best: " << 1 - itt->first.value[0] << "\n";
		}
		gen++;
	}
	if (opt.ShowReprot)
		std::cout<<POPulationToStr(pop);
}

void EA_Alg::GA::SaveResult(std::string fName)
{
	std::ofstream outFile{ fName,std::ios::out };
	if (!outFile) {
		std::cerr << "File could not open " << fName << std::endl;
		exit(EXIT_FAILURE);
	}
	outFile << opt.ToString() << "\n" << genReport << "</Gen>" << "\n" << POPulationToStr(pop) << std::endl;
	outFile.close();
}

bool EA_Alg::CMP_Rank(FIT_CHRO& s1, FIT_CHRO& s2)
{return (s1.first.rank < s2.first.rank); }

bool EA_Alg::CMP_CrowdingDistance(FIT_CHRO& s1, FIT_CHRO& s2)
{return (s1.first.cd > s2.first.cd);}

bool EA_Alg::CMP_Fit(FIT_CHRO& f1, FIT_CHRO& f2)
{return(f1.first.value[0] < f2.first.value[0]); }



double EA_Alg::RandomP() {return Utilities::RandomDouble(0.0, 1.0); }

string EA_Alg::POPulationToStr(const POPULATION& pop)
{
	string str = "\n<POP>";
	vector<pair< FITNESS, CHROMOSOME>>::const_iterator it = pop.begin();
	for (; it != pop.end(); it++)
	{
		str += "<f>"+std::to_string( it->first.value[0]) + "</f>\t";
		str += "<chr>";
		for (HALPA::SEED s : it->second)
			str += "<v>" + std::to_string(s.first) + "</v>  " + "<L>" + std::to_string(s.second.first) + "</L>\t\t";
		str += "</chr>\n";
	}
	str += "</POP>";
	return str;
}
