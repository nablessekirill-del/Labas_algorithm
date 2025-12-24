#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTreeWidget>
#include <QMap>
#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include "AbstractMethod.h"

QT_USE_NAMESPACE

    QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow {
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void handleMenuClick(QTreeWidgetItem *item, int column);
    void onLoadInpFile();
    void on_btn_1_clicked();
    void on_btn_2_clicked();  // Добавляем слот для сохранения

private:
    Ui::MainWindow *ui;
    QMap<QString, AbstractMethod*> methodsMap;
    void registerMethods();
    void plotGraph(const std::vector<GraphSeriesData>& seriesList);
    void saveOutputToFile(const QString& filePath, const QString& content);  // Обновленный метод

    QChart *chart;
    // currentOutputText больше не нужна, так как мы берем данные напрямую из textEdit_output
};

#endif
