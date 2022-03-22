// Original_LPA.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include<string>
#include<set>
#include<algorithm>
#include"..\lib\Type.h"
#include"..\lib\File.h"
#include"..\lib\Utilities.h"
#include"..\lib\Alg.h"
#include"..\NMI_Package\mutual.h"
#include"..\CommLineParser\cxxopts.hpp"
using std::string;
int main(int argc, char** argv)
{
    string name = "football";

    string gFile = "..\\Data\\" + name + ".txt";//Data\\LFR2\\pureLFR2K
    string gtFile = "..\\Data\\gt\\gt_" + name + ".txt";
    string outFile = "None";
    cxxopts::Options commLine("flag", "A novel method for overlapping communities detection in large scale complex net");
    commLine.allow_unrecognised_options();
    commLine.add_options()
        ("f,graph", "Edge list of graph file name", cxxopts::value<std::string>()->default_value(gFile))
        ("t,gt", "partitions file name:", cxxopts::value<std::string>()->default_value(gtFile))
        ("o,out", "partitions file name:", cxxopts::value<std::string>()->default_value(outFile))
        ("h,help", "Print usage");

    auto resultCommLine = commLine.parse(argc, argv);

    if (resultCommLine.count("graph"))
        gFile = resultCommLine["graph"].as<std::string>();
    if (resultCommLine.count("gt"))
        gtFile = resultCommLine["gt"].as<std::string>();
    if (resultCommLine.count("o"))
        outFile = resultCommLine["o"].as<std::string>();


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
    CLUSTERING cl;
    LPA::Run(gra, cl);
    //convert CLUSTERING result to MBR_CLUSTERING
    MBR_CLUSTERING mbrCL(gra.numNodes);
    int l = 0;
    for (CLUSTER cluster : cl)
    {
        for (NODE_ID n : cluster)
            mbrCL[n].push_back(std::make_pair(l, 1));
        l++;
    }

        
    //Utilities::PrintClustering(cl, gra.nodeNames);
    CLUSTERING gtCL;
    std::deque<std::deque<int>> nmi_CL, nmi_gtCL;

    File::LoadClustering(gtFile, gra.nodeNames, gtCL);
    Utilities::ClusteringToDeque(cl, nmi_CL);
    Utilities::ClusteringToDeque(gtCL, nmi_gtCL);
    NMI::Mutal mu;
    double nmi = mu.mutual3(nmi_CL, nmi_gtCL);
    std::cout << "NMI: " << nmi << std::endl;
    std::cout << "Q_ov: " << Utilities::OverlappingModularity(gra, mbrCL) << std::endl;
    std::cout << "F1: " << Utilities::F1Score(gtCL, cl) << "\n";
    if (outFile != "None")
    {
        if (File::SaveClustering(outFile, gra.nodeNames, cl))
            std::cout << "Save result in " << outFile << "\n";
    }

}

