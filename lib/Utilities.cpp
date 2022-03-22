#include"Utilities.h"
void Utilities::PrintClustering(const CLUSTERING& cl)
{
	for (CLUSTER cluster : cl)
	{
		for (NODE_ID node : cluster)
			std::cout << node << "\t";
		std::cout << "\n";
	}
}
void Utilities::PrintClustering(const CLUSTERING& cl, const MAP_NAME_ID& nodeNames)
{
	std::cout << "*****Print clusterint***********\n";
	for (CLUSTER cluster : cl)
	{
		for (NODE_ID nID : cluster)
		{
			NODE_NAME name;
			for (std::pair<NODE_NAME, NODE_ID> item : nodeNames)
				if (item.second == nID)
				{
					name = item.first;
					break;
				}
			std::cout << name << " ";
		}

		std::cout << "\n\n\n";
	}
	std::cout << "****************\n";
}
void Utilities::MbrClusteringToClustering(const MBR_CLUSTERING& mbrCL, CLUSTERING& cl)
{
	cl.clear();
	cl.resize(1);
	for (int n = 0; n < mbrCL.size(); n++)
	{
		NodeLabeles::const_iterator it = mbrCL[n].begin();
		for (; it != mbrCL[n].end(); it++)
		{
			if (it->first >= cl.size())
				cl.resize(it->first + 1);
			cl[it->first].insert(n);
		}
	}

	// for pure 
	CLUSTERING::iterator it_cl = cl.begin();
	for (; it_cl != cl.end();)
		if (it_cl->size() == 0)
			it_cl = cl.erase(it_cl);
		else
			it_cl++;
}
void Utilities::ClusteringToDeque(const CLUSTERING& cl, std::deque<std::deque<int>>& dCLU)
{
	dCLU.clear();
	for (CLUSTER cluster : cl)
	{
		std::deque<int> d1;
		for (NODE_ID n : cluster)
			d1.push_back(n);
		dCLU.push_back(d1);
	}
}
double Utilities::OverlappingModularity(const OWG& gra, const MBR_CLUSTERING& mbrCL)
{
	double sum = 0.0;
	const double m = gra.numEdges;
	CLUSTERING cl;
	Utilities::MbrClusteringToClustering(mbrCL, cl);
	for (CLUSTER cluster : cl)
	{
		for (NODE_ID u : cluster)
		{
			int du = gra.edges[u].size();

			for (NODE_ID v : cluster)
			{

				int dv = gra.edges[v].size();
				int a = 0;
				if (gra.edges[u].find(v) != gra.edges[u].end())
					a = 1;

				sum += (a - (((double)du * dv) / (2.0 * m))) * (1.0 / static_cast<double>(mbrCL[u].size())) * (1.0 / static_cast<double> (mbrCL[v].size()));
			}
		}
	}
	return sum / static_cast<double>(2.0 * m); //This function considers U to V edge and do not consider V to U 
}

void Utilities::UpdateClusteringWithNodeName(CLUSTERING& cl, const MAP_NAME_ID& nodeNames)
{
	CLUSTERING tmpClustering;
	for (CLUSTER cluster : cl)
	{
		CLUSTER tmpCluster;
		for (NODE_ID n : cluster)
		{
			for(std::pair<NODE_NAME,NODE_ID> item:nodeNames)
				if (item.second==n)
				{
					tmpCluster.insert(item.first);
					break;
				}
		}
		tmpClustering.push_back(tmpCluster);
	}
	cl = tmpClustering;
}

void Utilities::IntersectionTwoCluster(const CLUSTER& cl1, const CLUSTER& cl2, CLUSTER& interCL)
{
	CLUSTER::const_iterator f1, f2, l1, l2;
	interCL.clear();
	f1 = cl1.begin(); l1 = cl1.end();
	f2 = cl2.begin(); l2 = cl2.end();

	while (f1 != l1 && f2 != l2)
	{
		if (*f1 < *f2) ++f1;
		else if (*f2 < *f1) ++f2;
		else {
			interCL.insert(*f1);
			 ++f1; ++f2;
		}
	}
}

double Utilities::F1Score(const CLUSTERING& trueComms, const CLUSTERING& preComms)
{

	CLUSTERING::const_iterator it1 = preComms.begin();

	std::vector<double> F1;
	for (; it1 != preComms.end(); it1++)
	{
		int sizeIntersection = 0, sizeTrueComm = 1;
		CLUSTERING::const_iterator it2 = trueComms.begin();
		for (; it2 != trueComms.end(); it2++)
		{
			CLUSTER tmp;// = MyTooles::intersection(*it1, *it2);
			IntersectionTwoCluster(*it1, *it2, tmp);
			if (sizeIntersection < tmp.size()) {
				sizeIntersection = tmp.size();
				sizeTrueComm = it2->size();
			}
		}
		double r = (double)sizeIntersection / it1->size();
		double p = (double)sizeIntersection / sizeTrueComm;

		if ((p == 0) && (r == 0))
			F1.push_back(0.0);
		else
			F1.push_back(2 * ((r * p) / (r + p)));
	}
	double sum = 0.0;
	for (auto f1 : F1)
		sum += f1;
	return sum / (double)F1.size();
}

void Utilities::PrintInfoClustering(const MBR_CLUSTERING& mbrCl)
{
	CLUSTERING cl;
	MbrClusteringToClustering(mbrCl, cl);
	int numOV = 0;
	for (int i = 0; i < mbrCl.size(); i++)
		if (mbrCl[i].size() > 1)
			numOV++;
	
	int numComm1=0;
	for (CLUSTER cluster : cl)
		if (cluster.size() == 1)
			numComm1++;

	std::cout << "Num overlapping nodes:\t" << numOV << "\n";
	std::cout << "Num comm :\t" << cl.size() << "\n";
	std::cout << "Num comm 1:\t" << numComm1 << "\n";
}

void Utilities::CalCriteria(const OWG& gra,const MBR_CLUSTERING& mbrCL, const CLUSTERING& gtCL, double& nmi, double& q, double& f1)
{
	std::deque<std::deque<int>> nmi_CL, nmi_gtCL;
	ClusteringToDeque(gtCL, nmi_gtCL);
	CLUSTERING cl;
	Utilities::MbrClusteringToClustering(mbrCL, cl);
	Utilities::ClusteringToDeque(cl, nmi_CL);
	NMI::Mutal mu;
	nmi = mu.mutual3(nmi_CL, nmi_gtCL);

	// Q_ov
	q = Utilities::OverlappingModularity(gra, mbrCL);
	//F1-score
	f1 = Utilities::F1Score(gtCL, cl);
}

double Utilities::RandomDouble(double s, double e)
{
	// construct a trivial random generator engine from a time-based seed:
	unsigned seedRand = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seedRand);
	std::uniform_real_distribution<double> distribution(s, e);
	return distribution(generator);
}

int Utilities::RandomInt(int s, int e)
{

	// construct a trivial random generator engine from a time-based seed:
	unsigned seedRand = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seedRand);
	std::uniform_int_distribution<int> distribution(s, e);
	return distribution(generator);  // generates number in the range s..e 
	
}

