#ifndef NELDERMEAD_H
#define NELDERMEAD_H

#include <vector>
#include <functional>

int neldermead(
    std::vector<double>& x0,
    double eps,
    std::function<double(std::vector<double>)> func
    );

#endif // NELDERMEAD_H
