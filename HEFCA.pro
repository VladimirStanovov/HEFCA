TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    sample.cpp

include(deployment.pri)
qtcAddDeployment()

HEADERS += \
    random_numbers.h \
    sample.h

