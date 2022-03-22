// ELPA.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include<string>
#include<set>
#include<algorithm>
#include <iomanip>
#include<thread>
#include"..\lib\Type.h"
#include"..\lib\File.h"
#include"..\lib\Utilities.h"
#include"..\lib\Alg.h"
#include "..\lib\GA.h"
#include"..\NMI_Package\mutual.h"
#include"..\CommLineParser\cxxopts.hpp"
using std::string;
void ChromosomeEvaluation(const OWG& gra,const EA_Alg::CHROMOSOME& chr,const CLUSTERING& gtCL,double& nmi, double& q, double& f1)
{
    MBR_CLUSTERING mbrCL;
    vector<MEMBERSHIPE> maxMbr;
    mbrCL.resize(gra.numNodes);
    maxMbr.resize(gra.numNodes);
    HALPA::Run(gra, chr, mbrCL, maxMbr);
    HALPA::Normalizing(mbrCL);
   
    Utilities::CalCriteria(gra, mbrCL, gtCL, nmi, q, f1);
}
CLUSTERING ConvertCHR_Clustering(const OWG& gra,EA_Alg::CHROMOSOME chr)
{
    MBR_CLUSTERING mbrCL;
    vector<MEMBERSHIPE> maxMbr;
    mbrCL.resize(gra.numNodes);
    maxMbr.resize(gra.numNodes);
    HALPA::Run(gra, chr, mbrCL, maxMbr);
    HALPA::Normalizing(mbrCL);
    CLUSTERING cl;
    Utilities::MbrClusteringToClustering(mbrCL, cl);
    return cl;
}
int main(int argc, char** argv)
{
    int numRun = 1;
    string name =  "karate";
    string folder = "";// "Revise_LFR\\N2\\";

    string gFile = "..\\Data\\"+folder + name + ".txt";//Data\\LFR2\\pureLFR2K
    string gtFile = "..\\Data\\"+folder+"gt\\gt_" + name + ".txt";
    string outFile = "None";
    string gaSetting_File = "None";
    std::string gaResult_File = "None";
    
    cxxopts::Options commLine("flag", "A novel method for overlapping communities detection in large scale complex net");
    commLine.allow_unrecognised_options();
    commLine.add_options()
        ("f,graph", "Edge list of graph file name", cxxopts::value<std::string>()->default_value(gFile))
        ("t,gt", "partitions file name:", cxxopts::value<std::string>()->default_value(gtFile))
        ("o,out", "partitions file name:", cxxopts::value<std::string>()->default_value(outFile))
        ("gao", "GA setting file name", cxxopts::value<std::string>()->default_value(gaSetting_File))
        ("gaf", "GA result file name", cxxopts::value<std::string>()->default_value(gaResult_File))
        ("run", "The number of runs", cxxopts::value<int>()->default_value(std::to_string(numRun)))
        ("h,help", "Print usage");

    auto resultCommLine = commLine.parse(argc, argv);

    if (resultCommLine.count("run"))
        numRun = resultCommLine["run"].as<int>();
    if (resultCommLine.count("graph"))
        gFile = resultCommLine["graph"].as<std::string>();
    if (resultCommLine.count("gt"))
        gtFile = resultCommLine["gt"].as<std::string>();
    if (resultCommLine.count("o"))
        outFile = resultCommLine["o"].as<std::string>();
    if (resultCommLine.count("gao"))
        gaSetting_File=resultCommLine["gao"].as<std::string>();
    if (resultCommLine.count("gaf"))
        gaResult_File = resultCommLine["gaf"].as<std::string>();

    if (resultCommLine.count("h"))
    {
        std::cout << commLine.help();
        return 0;
    }

   
    OWG gra;
    if (!File::LoadFileEdgeList(gFile, gra)) return -1;
    std::cout << "name: " << gFile << "\n";
    std::cout << "Num of nodes: " << gra.numNodes << "\n";
    std::cout << "Num of edges: " << gra.numEdges << "\n";
    
    time_t st = time(0);
    HALPA::Weighting(gra);
    time_t en = time(0);
    std::cout << "Weighting Time : " << difftime(en, st) << "\n";

    EA_Alg::OptionGA opt;
    if (gaSetting_File == "None")
    {
        std::cout << "Use defult setting for EA \n";
        opt.setFromFile("GA.txt");
    }
    else
    {
        opt.setFromFile(gaSetting_File);
    }
    std::cout << opt.ToString()<<"\n";
    for (int count = 0; count < numRun; count++)
    {
        std::cout << "=============================";
        std::cout << " number\t" << count+1;
        std::cout <<"===============================\n\n";
        // opt.maxThread = 1;//  opt.popSize;
        EA_Alg::GA ga(gra, opt);
        st = time(0);
        ga.CreatePOP();
        en = time(0);
        std::cout << "Initial population Time : " << difftime(en, st) << "\n";

        st = time(0);
        ga.ManiLoop();
        en = time(0);
        std::cout << "Main loop Time : " << difftime(en, st) << "\n";

        CLUSTERING gtCL;
        std::deque<std::deque<int>> nmi_CL, nmi_gtCL;
        File::LoadClustering(gtFile, gra.nodeNames, gtCL);
        
        Utilities::ClusteringToDeque(gtCL, nmi_gtCL);

        EA_Alg::POPULATION pop = ga.GetPOP();

        int    i_NMI = -1, i_Q = -1, i_F1 = -1;
        CLUSTERING clNMI, clQ, clF1;
        double maxNMI = -2, maxQ = -2, maxF1 = -2;
        double sNMI = 0.0, sQ = 0.0, sF1 = 0.0;
        vector<double> nmi(pop.size()), q(pop.size()), f1(pop.size());
        int i = 0;
        while (i < pop.size())
        {
            vector<std::thread> th_pool;
            for (int j = 0; (j < opt.maxThread) && (i < opt.popSize); j++)
            {
                th_pool.push_back(std::thread(ChromosomeEvaluation, std::ref(gra), std::ref(pop[i].second), std::ref(gtCL), std::ref(nmi[i]), std::ref(q[i]), std::ref(f1[i])));
                i++;
            }

            for (std::thread& th : th_pool)
                if (th.joinable())
                    th.join();
        }

        
        for (int j = 0; j < pop.size(); j++)
        {
            sNMI += nmi[j];
            if (nmi[j] > maxNMI)
            {
                maxNMI = nmi[j];
                i_NMI = j;
                //clNMI = cl;
            }

            
            sQ += q[j];
            if (q[j] > maxQ)
            {
                maxQ = q[j];
                i_Q = j;
                //clQ = cl;
            }

            
            sF1 += f1[j];
            if (f1[j] > maxF1)
            {
                maxF1 = f1[j];
                i_F1 = j;
                //clF1 = cl;
            }
        }

        std::cout << std::setw(10)<<" "<< std::setw(10) << "NMI" << std::setw(10) << "Q_ov" << std::setw(10) << "F1-score" << "\n";
        std::cout << std::setw(10) << "Best: " << std::setw(10) << maxNMI << std::setw(10) << maxQ << std::setw(10) << maxF1 << "\n";
        std::cout << std::setw(10) << "Ave: " << std::setw(10) << sNMI / static_cast<double>(pop.size()) << std::setw(10) << sQ / static_cast<double>(pop.size()) << std::setw(10) << sF1 / static_cast<double>(pop.size()) << "\n";
        
        

        
        if (outFile != "None")
        {
            clNMI = ConvertCHR_Clustering(gra,pop[i_NMI].second);
            clQ = ConvertCHR_Clustering(gra, pop[i_Q].second);
            clF1 = ConvertCHR_Clustering(gra, pop[i_F1].second);
            if(numRun==1)
            {
                File::SaveClustering(outFile + "_NMI", gra.nodeNames, clNMI);
                File::SaveClustering(outFile + "_Q", gra.nodeNames, clQ);
                File::SaveClustering(outFile + "_F1", gra.nodeNames, clF1);
            }
            else
            {
                File::SaveClustering(outFile +"_" + std::to_string(count) + "_NMI", gra.nodeNames, clNMI);
                File::SaveClustering(outFile +"_" + std::to_string(count) + "_Q", gra.nodeNames, clQ);
                File::SaveClustering(outFile +"_" + std::to_string(count) + "_F1", gra.nodeNames, clF1);
            }
            
        }
        if (gaResult_File != "None")
        {
            if(numRun==1)
                ga.SaveResult(gaResult_File);
            else
                ga.SaveResult(gaResult_File +"_" + std::to_string(count));
        }
        std::cout <<"***********************************\n";
    }
   
	return 0;
}
