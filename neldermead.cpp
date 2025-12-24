#include "neldermead.h"
#include <cmath>
#include <algorithm>

using namespace std;

struct SimplexVertex {
    vector<double> x;
    double fx;
    SimplexVertex(int n) : x(n), fx(0.0) {}
    bool operator<(const SimplexVertex& other) const { return fx < other.fx; }
};

// Вспомогательная функция для генерации точки
void get_point(const vector<double>& x1, const vector<double>& x2, double coeff, int n, vector<double>& result) {
    for (int j = 0; j < n; ++j) {
        result[j] = (1.0 + coeff) * x1[j] - coeff * x2[j];
    }
}

int neldermead(vector<double>& x0, double eps, function<double(vector<double>)> func) {
    int n = x0.size();
    const double alpha = 1.0, beta = 0.5, gamma = 2.0;

    vector<SimplexVertex> simplex(n + 1, SimplexVertex(n));

    // 1. Построение начального симплекса
    simplex[0].x = x0;
    simplex[0].fx = func(x0);

    for (int i = 1; i <= n; ++i) {
        simplex[i].x = x0;
        simplex[i].x[i - 1] += (x0[i - 1] != 0) ? 0.05 * x0[i - 1] : 0.00025;
        simplex[i].fx = func(simplex[i].x);
    }

    int iterations = 0;
    vector<double> centroid(n), xr(n), xe(n), xc(n);

    while (iterations < 1000) {
        iterations++;
        sort(simplex.begin(), simplex.end());

        // Проверка сходимости
        double std_dev = 0;
        for (int i = 0; i <= n; ++i) std_dev += pow(simplex[i].fx - simplex[0].fx, 2);
        if (sqrt(std_dev / (n + 1)) < eps) break;

        // Центроид
        fill(centroid.begin(), centroid.end(), 0.0);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) centroid[j] += simplex[i].x[j] / n;

        // Отражение
        get_point(centroid, simplex[n].x, alpha, n, xr);
        double fxr = func(xr);

        if (simplex[0].fx <= fxr && fxr < simplex[n - 1].fx) {
            simplex[n].x = xr; simplex[n].fx = fxr;
            continue;
        }

        // Растяжение
        if (fxr < simplex[0].fx) {
            get_point(centroid, simplex[n].x, gamma, n, xe);
            double fxe = func(xe);
            if (fxe < fxr) { simplex[n].x = xe; simplex[n].fx = fxe; }
            else { simplex[n].x = xr; simplex[n].fx = fxr; }
            continue;
        }

        // Сжатие
        get_point(centroid, simplex[n].x, beta, n, xc);
        double fxc = func(xc);
        if (fxc < simplex[n].fx) {
            simplex[n].x = xc; simplex[n].fx = fxc;
            continue;
        }

        // Редукция (если ничего не помогло)
        for (int i = 1; i <= n; ++i) {
            for (int j = 0; j < n; ++j) {
                simplex[i].x[j] = simplex[0].x[j] + 0.5 * (simplex[i].x[j] - simplex[0].x[j]);
            }
            simplex[i].fx = func(simplex[i].x);
        }
    }

    x0 = simplex[0].x;
    return iterations;
}
