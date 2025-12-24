#ifndef METHOD_GRUBBS_H
#define METHOD_GRUBBS_H

#include "AbstractMethod.h"
#include <vector>
#include <algorithm>
#include <cmath>
#include <QString>
#include <boost/math/distributions/students_t.hpp>

class Method_Grubbs : public AbstractMethod {
public:

    bool hasGraph() override { return false; }

    QString calculate(const std::vector<double>& data, const std::vector<int>& cens) override {
        Q_UNUSED(cens);


        if (data.size() < 6) return "Ошибка: Недостаточно данных для формата Граббса";

        int n_val = static_cast<int>(data[0]);
        int state = static_cast<int>(data[1]);
        double alpha = data[2];
        double a_param = data[3];
        double s_param = data[4];

        std::vector<double> sample;
        for (size_t i = 5; i < data.size(); ++i) {
            sample.push_back(data[i]);
        }

        double n = static_cast<double>(sample.size());
        if (n < 3) return "Ошибка: выборка слишком мала";

        double sum = 0;
        for (double x : sample) sum += x;
        double mean = sum / n;

        double sq_sum = 0;
        double minVal = sample[0];
        double maxVal = sample[0];
        for (double x : sample) {
            sq_sum += (x - mean) * (x - mean);
            if (x < minVal) minVal = x;
            if (x > maxVal) maxVal = x;
        }
        double stdDev = std::sqrt(sq_sum / (n - 1.0));

        double u_max = (maxVal - mean) / stdDev;
        double u_min = (mean - minVal) / stdDev;
        double u_obs = (state == 0) ? std::max(u_max, u_min) : (state == 1 ? u_max : u_min);

        double u_crit = calculateUCrit(sample.size(), alpha, (state != 0));

        QString res;
        res += "Критерий Граббса для нормального распределения\n";
        res += "---------------------------------------------\n";
        res += QString("Размер выборки n      = %1\n").arg(n_val);
        res += QString("state                 = %1\n").arg(state);
        res += QString("Уровень значимости α  = %1\n\n").arg(alpha, 0, 'f', 6);

        res += "Истинные параметры распределения (использовались при генерации):\n";
        res += QString("a (мат. ожидание)     = %1\n").arg(a_param, 0, 'f', 6);
        res += QString("s (СКО)                = %1\n\n").arg(s_param, 0, 'f', 6);

        res += "Выборка:\n";
        for (size_t i = 0; i < sample.size(); ++i) {
            res += QString("x[%1] = %2\n").arg(i).arg(sample[i], 0, 'f', 6);
        }

        res += "\nОценки по выборке:\n";
        res += QString("Среднее значение           = %1\n").arg(mean, 0, 'f', 6);
        res += QString("Стандартное отклонение     = %1\n").arg(stdDev, 0, 'f', 6);
        res += QString("Наблюдаемая статистика u   = %1\n").arg(u_obs, 0, 'f', 6);
        res += QString("Критическое значение u_alpha   = %1\n\n").arg(u_crit, 0, 'f', 6);

        res += "Вывод: ";
        if (u_obs > u_crit) {
            res += "u > u_α, нулевая гипотеза отвергается (в выборке есть подозрительный выброс).";
        } else {
            res += "u <= u_α, нулевая гипотеза НЕ отвергается (подозрительных выбросов нет).";
        }

        return res;
    }


    std::vector<GraphSeriesData> getGraphData() override {
        return std::vector<GraphSeriesData>();
    }



private:
    double calculateUCrit(int n, double alpha, bool oneSided) {
        using namespace boost::math;
        try {
            students_t dist(n - 2);
            double p = oneSided ? (alpha / n) : (alpha / (2.0 * n));
            double t = quantile(complement(dist, p));
            return ((n - 1.0) / std::sqrt(static_cast<double>(n))) * std::sqrt((t * t) / (n - 2.0 + t * t));
        } catch (...) { return 0.0; }
    }
};

#endif
