

#ifndef SUBGRAPHMATCHING_MATCHINGCOMMAND_H
#define SUBGRAPHMATCHING_MATCHINGCOMMAND_H

#include "utility/commandparser.h"
#include <map>
#include <iostream>
enum OptionKeyword {
    Algorithm = 0,          // -a, The algorithm name, compulsive parameter
    QueryGraphFile = 1,     // -q, The query graph file path, compulsive parameter
    DataGraphFile = 2,      // -d, The data graph file path, compulsive parameter
    ThreadCount = 3,        // -n, The number of thread, optional parameter
    DepthThreshold = 4,     // -d0,The threshold to control the depth for splitting task, optional parameter
    WidthThreshold = 5,     // -w0,The threshold to control the width for splitting task, optional parameter
    IndexType = 6,          // -i, The type of index, vertex centric or edge centric
    Filter = 7,             // -filter, The strategy of filtering
    Order = 8,              // -order, The strategy of ordering
    Engine = 9,             // -engine, The computation engine
    MaxOutputEmbeddingNum = 10, // -num, The maximum output embedding num
    SpectrumAnalysisTimeLimit = 11, // -time_limit, The time limit for executing a query in seconds
    SpectrumAnalysisOrderNum = 12, // -order_num, The number of matching orders generated
    DistributionFilePath = 13,          // -dis_file, The output path of the distribution array
    CSRFilePath = 14,                    // -csr, The input csr file path
    IsEigenCheck =15,                //-eigen, Boolean, choose whether to enhance the filter with eigen check or not
    TopSEigen = 16,                  //-tops,  How many eigenvalues use for eigen check
    Dataset =17,                    // name of the dataset
    QuerySize = 18,
    QueryProperty = 19,
    QueryNumber = 20,
    alphaPar=21, //-alpha
    betaPar=22,  //-beta
};

class MatchingCommand : public CommandParser{
private:
    std::map<OptionKeyword, std::string> options_key;
    std::map<OptionKeyword, std::string> options_value;

private:
    void processOptions();

public:
    MatchingCommand(int argc, char **argv);

    std::string getQueryNumber(){
        return options_value[OptionKeyword::QueryNumber];
    }

    std::string getQuerySize() {
        return options_value[OptionKeyword::QuerySize];
    }

    std::string getQueryProperty() {
        return options_value[OptionKeyword::QueryProperty];
    }

    std::string getDatasetName() {
        return options_value[OptionKeyword::Dataset];
    }

    std::string getEigenOrNot() {
        return options_value[OptionKeyword::IsEigenCheck];
    }

    std::string getHowManyEigen(){
        return options_value[OptionKeyword::TopSEigen];
    }

    std::string getDataGraphFilePath() {
        return options_value[OptionKeyword::DataGraphFile];
    }

    std::string getQueryGraphFilePath() {
        return options_value[OptionKeyword::QueryGraphFile];
    }

    std::string getAlgorithm() {
        return options_value[OptionKeyword::Algorithm];
    }

    std::string getIndexType() {
        return options_value[OptionKeyword::IndexType] == "" ? "VertexCentric" : options_value[OptionKeyword::IndexType];
    }
    std::string getThreadCount() {
        return options_value[OptionKeyword::ThreadCount] == "" ? "1" : options_value[OptionKeyword::ThreadCount];
    }

    std::string getDepthThreshold() {
        return options_value[OptionKeyword::DepthThreshold] == "" ? "0" : options_value[OptionKeyword::DepthThreshold];
    }

    std::string getWidthThreshold() {
        return options_value[OptionKeyword::WidthThreshold] == "" ? "1" : options_value[OptionKeyword::WidthThreshold];
    }

    std::string getFilterType() {
        return options_value[OptionKeyword::Filter] == "" ? "CFL" : options_value[OptionKeyword::Filter];
    }

    std::string getOrderType() {
        return options_value[OptionKeyword::Order] == "" ? "GQL" : options_value[OptionKeyword::Order];
    }

    std::string getEngineType() {
        return options_value[OptionKeyword::Engine] == "" ? "LFTJ" : options_value[OptionKeyword::Engine];
    }

    std::string getMaximumEmbeddingNum() {
        return options_value[OptionKeyword::MaxOutputEmbeddingNum] == "" ? "MAX" : options_value[OptionKeyword::MaxOutputEmbeddingNum];
    }

    std::string getTimeLimit() {
        return options_value[OptionKeyword::SpectrumAnalysisTimeLimit] == "" ? "60" : options_value[OptionKeyword::SpectrumAnalysisTimeLimit];
    }

    std::string getOrderNum() {
        return options_value[OptionKeyword::SpectrumAnalysisOrderNum] == "" ? "100" : options_value[OptionKeyword::SpectrumAnalysisOrderNum];
    }

    std::string getDistributionFilePath() {
        return options_value[OptionKeyword::DistributionFilePath] == "" ? "temp.distribution" : options_value[OptionKeyword::DistributionFilePath];
    }

    std::string getCSRFilePath() {
        return options_value[OptionKeyword::CSRFilePath] == "" ? "" : options_value[OptionKeyword::CSRFilePath];
    }
    std::string getalpha() {
        return options_value[OptionKeyword::alphaPar] == "" ? "25" : options_value[OptionKeyword::alphaPar];
    }
        std::string getbeta() {
        return options_value[OptionKeyword::betaPar] == "" ? "500" : options_value[OptionKeyword::betaPar];
    }

};


#endif //SUBGRAPHMATCHING_MATCHINGCOMMAND_H
