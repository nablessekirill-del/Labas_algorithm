#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "Method_MLE_Normal.h"
#include "Method_MLE_Weibull.h"
#include "Method_MLS_Normal.h"
#include "Method_MLS_Weibull.h"
#include "Method_Grubbs.h"
#include "Method_FisherStudent.h"
#include "Method_Anova.h"
#include "Method_ShapiroWilk.h"
#include "Method_Wilcoxon.h"
#include <QtCharts/QLineSeries>
#include <QtCharts/QScatterSeries>
#include <QtCharts/QValueAxis>
#include <QtCharts/QLogValueAxis>
#include <QRegularExpression>
#include <QFileDialog>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <QDateTime>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow), chart(new QChart())
{
    ui->setupUi(this);
    registerMethods();

    if (ui->chartContainer) {
        ui->chartContainer->setChart(chart);
        ui->chartContainer->setRenderHint(QPainter::Antialiasing);
        ui->chartContainer->setRubberBand(QChartView::RectangleRubberBand);
    }

    // Соединяем сигналы и слоты
    connect(ui->treeWidget, &QTreeWidget::itemClicked, this, &MainWindow::handleMenuClick);
    connect(ui->btn_1, &QPushButton::clicked, this, &MainWindow::on_btn_1_clicked);
    connect(ui->btn_load, &QPushButton::clicked, this, &MainWindow::onLoadInpFile);
    //connect(ui->btn_2, &QPushButton::clicked, this, &MainWindow::on_btn_2_clicked);
}

void MainWindow::registerMethods() {
    methodsMap["Нормальное распределение"] = new Method_MLE_Normal();
    methodsMap["Распределение Вейбулла-Гнеденко"] = new Method_MLE_Weibull();
    methodsMap["Нормальное распределение MLS"] = new Method_MLS_Normal();
    methodsMap["Распределение Вейбулла-Гнеденко MLS"] = new Method_MLS_Weibull();
    methodsMap["Критерий Граббса"] = new Method_Grubbs();
    methodsMap["Критерий Фишера-Стьюдента"] = new Method_FisherStudent();
    methodsMap["Однофакторный дисперсионный анализ (ANOVA)"] = new Method_Anova();
    methodsMap["Критерий Шапило-Уилка (W-критерий)"] = new Method_ShapiroWilk();
    methodsMap["Двухвыборочный критерий Уилкоксона"] = new Method_Wilcoxon();
}

MainWindow::~MainWindow() {
    qDeleteAll(methodsMap);
    delete ui;
}

void MainWindow::handleMenuClick(QTreeWidgetItem *item, int column) {
    Q_UNUSED(column);
    if (!item) return;
    QString methodName = item->text(0);
    if (!methodsMap.contains(methodName)) return;
}

void MainWindow::on_btn_1_clicked() {
    QString inputStr = ui->textEdit_input->toPlainText();
    if (inputStr.isEmpty()) {
        ui->textEdit_output->setText("Ошибка: Введите данные или загрузите файл!");
        return;
    }

    std::vector<double> data;
    std::vector<int> cens;

    QStringList lines = inputStr.split('\n', Qt::SkipEmptyParts);
    for (int i = 0; i < lines.size(); ++i) {
        QString line = lines[i].trimmed();
        if (line == "Data" && i + 1 < lines.size()) {
            QStringList vals = lines[i+1].trimmed().split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
            for(const QString &v : vals) data.push_back(v.toDouble());
        }
        else if (line == "Censorizes" && i + 1 < lines.size()) {
            QStringList vals = lines[i+1].trimmed().split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
            for(const QString &v : vals) cens.push_back(v.toInt());
        }
    }

    if (data.empty()) {
        QStringList parts = inputStr.split(QRegularExpression("\\s+"), Qt::SkipEmptyParts);
        for(const QString& p : parts) {
            bool ok;
            double val = p.toDouble(&ok);
            if(ok) data.push_back(val);
        }
    }

    if (data.empty()) {
        ui->textEdit_output->setText("Ошибка: Не удалось распознать числа в блоке Data!");
        return;
    }

    QTreeWidgetItem* item = ui->treeWidget->currentItem();
    if (!item) {
        ui->textEdit_output->setText("Ошибка: Сначала выберите распределение в списке слева!");
        return;
    }

    QString methodName = item->text(0);
    if (methodsMap.contains(methodName)) {
        AbstractMethod* method = methodsMap[methodName];

        if (cens.size() != data.size()) cens.assign(data.size(), 0);

        QString report = method->calculate(data, cens);
        ui->textEdit_output->setText(report);

        if (method->hasGraph()) {
            plotGraph(method->getGraphData());
        }
    }
}

void MainWindow::onLoadInpFile() {
    QString fileName = QFileDialog::getOpenFileName(this, "Открыть .inp файл", "", "Input Files (*.inp);;All Files (*)");
    if (fileName.isEmpty()) return;

    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return;

    QTextStream in(&file);
    QString fullContent = in.readAll();
    file.close();

    ui->textEdit_input->setText(fullContent);
    on_btn_1_clicked();
}

void MainWindow::plotGraph(const std::vector<GraphSeriesData>& seriesList) {
    chart->removeAllSeries();

    // Очищаем старые оси
    for(auto axis : chart->axes()) {
        chart->removeAxis(axis);
        delete axis;
    }

    if (seriesList.empty()) return;

    // 1. ОПРЕДЕЛЯЕМ ТИП ОСИ
    bool useLogX = false;
    if (ui->treeWidget->currentItem()) {
        QString currentMethod = ui->treeWidget->currentItem()->text(0);
        if (currentMethod.contains("Вейбулл")) {
            useLogX = true;
        }
    }

    // 2. СОЗДАЕМ ОСИ
    QAbstractAxis *axisX;
    if (useLogX) {
        QLogValueAxis *logAxis = new QLogValueAxis();
        logAxis->setBase(10);
        logAxis->setLabelFormat("%g"); // Красивые числа: 1000, 2000...
        logAxis->setTitleText("Время (лог. шкала)");
        axisX = logAxis;
    } else {
        QValueAxis *linAxis = new QValueAxis();
        linAxis->setTitleText("Время");
        axisX = linAxis;
    }

    QValueAxis *axisY = new QValueAxis();
    axisY->setTitleText("Квантиль");

    chart->addAxis(axisX, Qt::AlignBottom);
    chart->addAxis(axisY, Qt::AlignLeft);

    // 3. НАСТРОЙКА ГРАНИЦ
    double xMin = 1e18, xMax = -1e18;
    double yMin = 1e18, yMax = -1e18;

    // Сначала найдем экстремумы
    for (const auto& sData : seriesList) {
        for (size_t i = 0; i < sData.x.size(); ++i) {
            if(sData.x[i] < xMin) xMin = sData.x[i];
            if(sData.x[i] > xMax) xMax = sData.x[i];
            if(sData.y[i] < yMin) yMin = sData.y[i];
            if(sData.y[i] > yMax) yMax = sData.y[i];
        }
    }

    // Установка диапазонов с отступами
    if (useLogX) {
        // Для логарифмической оси мин. значение не может быть <= 0
        axisX->setRange(std::max(0.1, xMin * 0.7), xMax * 1.3);
    } else {
        double xDelta = (xMax - xMin) * 0.1;
        axisX->setRange(xMin - xDelta, xMax + xDelta);
    }
    axisY->setRange(yMin - 0.5, yMax + 0.5);

    // 4. ДОБАВЛЕНИЕ СЕРИЙ
    for (const auto& sData : seriesList) {
        QXYSeries *series;
        QString sName = QString::fromStdString(sData.name);

        if (sData.isScatter) {
            auto *scatter = new QScatterSeries();
            scatter->setMarkerSize(10);
            if (sName == "Censored") {
                scatter->setBrush(Qt::NoBrush);
                scatter->setPen(QPen(Qt::red, 2));
                scatter->setMarkerShape(QScatterSeries::MarkerShapeRectangle);
            } else {
                scatter->setColor(Qt::blue);
                scatter->setMarkerShape(QScatterSeries::MarkerShapeCircle);
            }
            series = scatter;
        } else {
            auto *line = new QLineSeries();
            QPen pen;
            if (sName == "MLE Линия") pen = QPen(Qt::darkBlue, 3);
            else pen = QPen(Qt::gray, 1, Qt::DashLine);
            line->setPen(pen);
            series = line;
        }

        series->setName(sName);
        for (size_t i = 0; i < sData.x.size(); ++i) {
            series->append(sData.x[i], sData.y[i]);
        }

        chart->addSeries(series);
        series->attachAxis(axisX);
        series->attachAxis(axisY);
    }
    chart->legend()->setVisible(true);
}

void MainWindow::on_btn_2_clicked() {
    // Проверяем текущий вывод из textEdit_output
    QString outputText = ui->textEdit_output->toPlainText();

    if (outputText.isEmpty()) {
        QMessageBox::warning(this, "Предупреждение",
                             "Нет данных для сохранения. Сначала выполните расчет.");
        return;
    }

    // Проверяем, не является ли вывод ошибкой
    if (outputText.contains("Ошибка:") || outputText.startsWith("Ошибка")) {
        QMessageBox::warning(this, "Предупреждение",
                             "Нельзя сохранить результат с ошибкой. Сначала исправьте ошибки.");
        return;
    }

    // Получаем путь для сохранения файла
    QString fileName = QFileDialog::getSaveFileName(this,
                                                    "Сохранить результаты",
                                                    QDir::homePath() + "/analysis_results.out",
                                                    "Output Files (*.out);;All Files (*)");

    if (fileName.isEmpty()) {
        return;  // Пользователь отменил сохранение
    }

    // Убедимся, что у файла правильное расширение
    if (!fileName.endsWith(".out", Qt::CaseInsensitive)) {
        fileName += ".out";
    }

    saveOutputToFile(fileName, outputText);
}

void MainWindow::saveOutputToFile(const QString& filePath, const QString& content) {
    QFile file(filePath);

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        QMessageBox::critical(this, "Ошибка",
                              "Не удалось открыть файл для записи: " + file.errorString());
        return;
    }

    QTextStream out(&file);

    // Добавляем заголовок с информацией о методе и дате
    out << "================================================================================\n";
    out << "Анализ данных - Результаты расчета\n";

    QTreeWidgetItem* item = ui->treeWidget->currentItem();
    if (item) {
        out << "Метод: " << item->text(0) << "\n";
    }

    out << "Дата создания: " << QDateTime::currentDateTime().toString("dd.MM.yyyy hh:mm:ss") << "\n";
    out << "================================================================================\n\n";

    // Сохраняем основные результаты
    out << content;

    file.close();

    QMessageBox::information(this, "Успешно",
                             "Результаты сохранены в файл:\n" + filePath);
}
