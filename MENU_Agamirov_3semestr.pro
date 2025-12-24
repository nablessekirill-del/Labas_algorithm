QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++17


# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    analysis.cpp \
    main.cpp \
    mainwindow.cpp \
    neldermead.cpp

HEADERS += \
    AbstractMethod.h \
    Method_Anova.h \
    Method_FisherStudent.h \
    Method_Grubbs.h \
    Method_MLE_Normal.h \
    Method_MLE_Weibull.h \
    Method_MLS_Normal.h \
    Method_MLS_Weibull.h \
    Method_ShapiroWilk.h \
    Method_Wilcoxon.h \
    analysis.h \
    mainwindow.h \
    neldermead.h

FORMS += \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

# Путь к заголовочным файлам
INCLUDEPATH += /opt/homebrew/Cellar/boost/1.89.0_1/include

# Путь к скомпилированным библиотекам
LIBS += -L/opt/homebrew/Cellar/boost/1.89.0_1/lib

QT += charts
