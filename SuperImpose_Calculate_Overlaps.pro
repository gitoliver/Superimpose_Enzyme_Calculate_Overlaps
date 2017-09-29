TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11
QMAKE_CXXFLAGS += -std=c++0x
SOURCES += main.cpp \
    io.cpp

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../gems/gmml/bin/release/ -lgmml  -pthread -lpthread -Wl,--no-as-needed
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../gems/gmml/bin/debug/ -lgmml -pthread -lpthread -Wl,--no-as-needed
else:unix: LIBS += -L$$PWD/../../gems/gmml/bin/ -lgmml -pthread -lpthread -Wl,--no-as-needed

INCLUDEPATH += $$PWD/../../gems/gmml/bin
DEPENDPATH += $$PWD/../../gems/gmml/bin

HEADERS += \
    io.h

TARGET = SuperImpose_Calculate_Overlaps
