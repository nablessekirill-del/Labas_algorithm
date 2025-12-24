#ifndef METHOD_SHAPIROWILK_H
#define METHOD_SHAPIROWILK_H

#include "AbstractMethod.h"
#include <vector>
#include <cmath>
#include <algorithm>
#include <QString>
#include <boost/math/distributions/normal.hpp>

class Method_ShapiroWilk : public AbstractMethod {
public:
    bool hasGraph() override { return false; }

    QString calculate(const std::vector<double>& data, const std::vector<int>& cens) override {
        Q_UNUSED(cens);
        if (data.size() < 4) return "Ошибка: Недостаточно данных (n должно быть >= 3)";

        int n = static_cast<int>(data[0]);
        std::vector<double> sample;
        for (int i = 1; i <= n; ++i) {
            if (i < data.size()) sample.push_back(data[i]);
        }
        std::sort(sample.begin(), sample.end());

        double sum = std::accumulate(sample.begin(), sample.end(), 0.0);
        double mean = sum / n;
        double sSquared = 0;
        for (double x : sample) sSquared += (x - mean) * (x - mean);
        double stdDev = std::sqrt(sSquared / (n - 1));

        std::vector<double> a(n);
        std::vector<double> m(n);
        boost::math::normal_distribution<double> dist(0.0, 1.0);

        double m_norm = 0;
        for (int i = 1; i <= n; ++i) {
            m[i-1] = boost::math::quantile(dist, (i - 0.375) / (n + 0.25));
            m_norm += m[i-1] * m[i-1];
        }
        m_norm = std::sqrt(m_norm);
        for (int i = 0; i < n; ++i) a[i] = m[i] / m_norm;

        double b = 0;
        for (int i = 0; i < n; ++i) {
            b += a[i] * sample[i];
        }
        double W_obs = (b * b) / sSquared;

        double W_crit = getWcrit05(n);

        QString res = "Критерий Шапиро-Уилка\n";
        res += "Уровень значимости alpha = 0.05\n\n";
        res += QString("Размер выборки n = %1\n").arg(n);
        res += QString("Выборочное среднее x̄ = %1\n").arg(mean, 0, 'f', 6);
        res += QString("Выборочное СКО s = %1\n\n").arg(stdDev, 0, 'f', 7);

        res += QString("Сумма квадратов отклонений s^2 = %1\n").arg(sSquared, 0, 'f', 7);
        res += QString("Коэффициент b = %1\n").arg(std::abs(b), 0, 'f', 6);
        res += QString("Наблюдаемое значение статистики Wнабл = %1\n").arg(W_obs, 0, 'f', 4);
        res += QString("Критическое значение Wкр = %1\n\n").arg(W_crit, 0, 'f', 3);

        if (W_obs >= W_crit) {
            res += "Вывод: Wнабл > Wкр, нет оснований отвергать нулевую гипотезу.\n";
            res += "Распределение можно считать нормальным.";
        } else {
            res += "Вывод: Wнабл < Wкр, нулевая гипотеза отвергается.\n";
            res += "Распределение нельзя считать нормальным.";
        }

        return res;
    }

private:
    double getWcrit05(int n) {
        static const double Wcrit05[] = {
            0.0, 0.0, 0.0, 0.767, 0.748, 0.762, 0.788, 0.803, 0.818, 0.829,
            0.842, 0.850, 0.859, 0.866, 0.874, 0.881, 0.887, 0.892, 0.897, 0.901,
            0.905, 0.908, 0.911, 0.914, 0.916, 0.918, 0.920, 0.923, 0.924, 0.926,
            0.927, 0.929, 0.930, 0.931, 0.933, 0.934, 0.935, 0.936, 0.938, 0.939,
            0.940, 0.941, 0.942, 0.943, 0.944, 0.945, 0.945, 0.946, 0.947, 0.947, 0.947
        };
        if (n < 3) return 0;
        if (n > 50) return 0.947;
        return Wcrit05[n];
    }

    std::vector<GraphSeriesData> getGraphData() override { return {}; }
};

#endif
