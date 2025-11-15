#include "MainWindow.h"
#include <QtWidgets/QApplication>
#include <QGuiApplication>
#include <QTranslator>  
#include <QLibraryInfo> 
#include <QLocale>    
#include <qicon.h>

int main(int argc, char * argv[])
{

    QApplication app(argc, argv);

    app.setWindowIcon(QIcon(":/res/favicon2.ico"));

    QTranslator qtTranslator;

    if(qtTranslator.load("qt_" + QLocale::system().name(),
        QLibraryInfo::path(QLibraryInfo::TranslationsPath)))
    {
        app.installTranslator(&qtTranslator);
    }

    QTranslator qtBaseTranslator;
    if(qtBaseTranslator.load("qtbase_" + QLocale::system().name(),
        QLibraryInfo::path(QLibraryInfo::TranslationsPath)))
    {
        app.installTranslator(&qtBaseTranslator);
    }


    MainWindow window;
    window.show();
    return app.exec();
}