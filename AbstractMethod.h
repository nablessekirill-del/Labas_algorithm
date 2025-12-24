#ifndef ABSTRACTMETHOD_H
#define ABSTRACTMETHOD_H

#include <vector>
#include <string>
#include <QString>

struct GraphSeriesData {
    std::string name;
    std::vector<double> x;
    std::vector<double> y;
    bool isScatter = false;
};

class AbstractMethod {
public:
    virtual ~AbstractMethod() {}
    virtual QString calculate(const std::vector<double>& data, const std::vector<int>& cens) = 0;
    virtual bool hasGraph() { return false; } // Без const
    virtual std::vector<GraphSeriesData> getGraphData() = 0;
};

#endif
