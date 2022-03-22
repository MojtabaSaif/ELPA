#pragma once
#include<string>
#include <iostream>
#include <fstream>
#include <thread>
#include <mutex>
#include <map>
#include<algorithm>
#include<random>
#include<thread>
#include"Type.h"
#include"Alg.h"
#include"Utilities.h"
namespace EA_Alg {

	typedef  vector<HALPA::SEED> CHROMOSOME;
	//typedef  double FITNESS;
	struct FITNESS;
	typedef pair<FITNESS, CHROMOSOME> FIT_CHRO;
	typedef  vector<FIT_CHRO> POPULATION;
	typedef  POPULATION::iterator POP_ITR;
	bool CMP_Rank(FIT_CHRO& s1, FIT_CHRO& s2) ;
	bool CMP_CrowdingDistance(FIT_CHRO& s1, FIT_CHRO& s2) ;
	bool CMP_Fit(FIT_CHRO& f1, FIT_CHRO& f2);
	typedef vector<POP_ITR> FRONT;
	struct FITNESS {
		vector<double> value;
		//Add field for multi-obj
		std::vector<POP_ITR> SP;      // Dominate set
		unsigned int np = 0;				  // Dominate count
		double cd = 0.0;           //Crowding distance
		unsigned int rank = -1;

		bool operator<(const FITNESS& t) const
		{
			if (value.size() == 1)
				return (this->value[0] < t.value[0]);
			else
				return (this->rank<t.rank);
		}
	};
	struct OptionGA {
		int			lenSol = 2;                    // Define the length of each solution
		double	     pSeeds = 0.5;               // Determind probelity of selected node from the seeds in Init pop
		double	     pCrossOver = 0.8;         // Cross over rate
		double	    pMutation = 0.0;          //Mutation rate
		int			popSize = 100;              // Define the size of the population
		//int		    numComm = lenSol;    // Number of the community in graph
		int          MaxGenCount = 50;     // 
		int          TournamentSize = 4;   //
		bool         ShowReprot = true;
		bool		 Multi = true;
		int          maxThread = 1;
		std::string  ToString();
		void         setFromFile(std::string fName);
		OptionGA() { maxThread = std::thread::hardware_concurrency(); };
		void         SetMaxThread(int v);
	};
	using std::thread;
	double RandomP();
	string POPulationToStr(const POPULATION& pop);
	class GA
	{
	private:
		POPULATION pop;
		OptionGA opt;
		const OWG& ga_gra;
		string genReport = "<Gen>\nGen\tBestFit\tAvg Fit\n";
		vector<pair<double,NODE_ID>> cumP;
		vector<NODE_ID> gSeeds;
		string genReprot;
		/// <summary>
		/// Create a chromosome and return the fitness of it
		/// </summary>
		/// <param name="chr">Chromosome</param>
		/// <returns>Chromosome fitness</returns>
		void CreateChromosome(CHROMOSOME& chr,FITNESS& f);
		void CreateChromosome2(CHROMOSOME& chr, FITNESS& f);
		void CreateChromosome3(CHROMOSOME& chr, FITNESS& fit);
		void CreateChromosome4(CHROMOSOME& chr, FITNESS& fit);
		POP_ITR TournamentSelection(int size);
		void SinglePointCrossOver(const CHROMOSOME& p1, const CHROMOSOME& p2, CHROMOSOME& child);
		void Mutation(const CHROMOSOME& ch_in,CHROMOSOME& ch_out);
		void Fitness(const CHROMOSOME& chr, FITNESS& outValue);
		LABEL_ID GenerateRandomLabel();
		bool Dominates(const FITNESS& f1, const FITNESS& f2);
		void NoneDominateedSorting(POPULATION& pop, vector<FRONT>& fs);
		void CalcCrowdingDistance(POPULATION& pop, vector<FRONT>& fs);
		void SumEdges(const OWG&,  MBR_CLUSTERING& mbrCL, double& sumIN, double& sumOUT);
	public:
		GA(const OWG& gra, OptionGA& option) :ga_gra(gra), opt(option) { };
		void CreatePOP();
		void ManiLoop();
		POPULATION GetPOP() { return pop; };
		void SaveResult(std::string fName);
	};
}

