#pragma once

#include <QtWidgets/QMainWindow>
#include "matrix/matrix.h" 
#include "matrix/rational.h"

QT_BEGIN_NAMESPACE
class QTableWidget;
class QPushButton;
class QSpinBox;
class QPlainTextEdit;
class QGridLayout;
class QComboBox;
QT_END_NAMESPACE

using MatrixD = matrix<double>;

using MatrixR = matrix<rational>;

class MainWindow: public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget * parent = nullptr);
    ~MainWindow();

private slots:

    void updateMatrixADim();
    void updateMatrixBDim();

    void onAdd();
    void onSubtract();
    void onMultiply();

    void onDetA();
    void onInvA();
    void onRefA();
    void onRrefA();
    void onEigenvaluesA();

    void onAdjA();
    void onRankA();
    void onTraceA();

    void onSolveHomogeneous();
    void onSolveNonHomogeneous();

private:
    void setupUi(); 

    MatrixD tableToMatrix_Double(QTableWidget * table);
    MatrixR tableToMatrix_Rational(QTableWidget * table);

    void matrixToTable(const MatrixD & m, QTableWidget * table);

    void displayResult(const QString & title, const QString & content);
    void displayError(const std::exception & e);

    QGridLayout * mainLayout;
    QTableWidget * matrixATable;
    QTableWidget * matrixBTable;
    QPlainTextEdit * resultTextEdit;

    QSpinBox * rowsASpin;
    QSpinBox * colsASpin;
    QSpinBox * rowsBSpin;
    QSpinBox * colsBSpin;

    QComboBox * typeSelector;

};