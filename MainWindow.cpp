#include "MainWindow.h"
#include <QTableWidget>
#include <QPushButton>
#include <QSpinBox>
#include <QPlainTextEdit>
#include <QGridLayout>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QLabel>
#include <QHeaderView>
#include <QString>
#include <string>
#include <QComboBox> 
#include <QIcon>
#include <qicon.h>

QString formatVector(const std::vector<double> & vec)
{
    QString result = "[ ";
    for(size_t i = 0; i < vec.size(); ++i)
    {
        result += QString::number(vec[i], 'f', 4);
        if(i < vec.size() - 1)
        {
            result += ", ";
        }
    }
    result += " ]ᵀ";
    return result;
}

QString formatHomogeneousSolution(const std::vector<std::vector<double>> & basis)
{
    if(basis.empty())
    {
        return "计算出错或返回为空。";
    }

    bool is_trivial_solution = true;
    if(basis.size() == 1)
    {
        for(double val : basis[0])
        {
            if(!is_zero(val))
            {
                is_trivial_solution = false;
                break;
            }
        }
    }
    else
    {
        is_trivial_solution = false;
    }

    if(is_trivial_solution)
    {
        return QString("方程组只有唯一的零解:\nx = %1").arg(formatVector(basis[0]));
    }

    QString result = "方程组有无穷多解，通解为:\nx = ";
    for(size_t i = 0; i < basis.size(); ++i)
    {
        result += QString("c%1 * ").arg(i + 1) + formatVector(basis[i]);
        if(i < basis.size() - 1)
        {
            result += " + \n    "; 
        }
    }
    result += "\n(其中 c1, c2, ... 为任意常数)";
    return result;
}


QString formatNonHomogeneousSolution(std::vector<std::vector<double>> result_data)
{
    if(result_data.empty())
    {
        return "方程组无解。";
    }

    std::vector<double> particular_solution = result_data.back();
    result_data.pop_back();

    if(result_data.empty())
    {
        return QString("方程组有唯一解:\nx = %1").arg(formatVector(particular_solution));
    }

    QString result = "方程组有无穷多解，通解为:\nx = ";
    result += formatVector(particular_solution) + " + \n    ";

    for(size_t i = 0; i < result_data.size(); ++i)
    {
        result += QString("c%1 * ").arg(i + 1) + formatVector(result_data[i]);
        if(i < result_data.size() - 1)
        {
            result += " + \n    ";
        }
    }
    result += "\n(其中 c1, c2, ... 为任意常数)";
    return result;
}

QString formatVector_Rational(const std::vector<rational> & vec)
{
    QString result = "[ ";
    for(size_t i = 0; i < vec.size(); ++i)
    {
        result += QString::fromStdString(std::format("{}", vec[i]));
        if(i < vec.size() - 1)
        {
            result += ", ";
        }
    }
    result += " ]ᵀ";
    return result;
}

QString formatHomogeneousSolution_Rational(const std::vector<std::vector<rational>> & basis)
{
    if(basis.empty()) return "计算出错或返回为空。";

    bool is_trivial_solution = true;
    if(basis.size() == 1)
    {
        for(const auto & val : basis[0])
        {
            if(val.num != 0) { is_trivial_solution = false; break; }
        }
    }
    else
    {
        is_trivial_solution = false;
    }

    if(is_trivial_solution)
    {
        return QString("方程组只有唯一的零解:\nx = %1").arg(formatVector_Rational(basis[0]));
    }

    QString result = "方程组有无穷多解，通解为:\nx = ";
    for(size_t i = 0; i < basis.size(); ++i)
    {
        result += QString("c%1 * ").arg(i + 1) + formatVector_Rational(basis[i]);
        if(i < basis.size() - 1) result += " + \n    ";
    }
    result += "\n(其中 c1, c2, ... 为任意常数)";
    return result;
}

QString formatNonHomogeneousSolution_Rational(std::vector<std::vector<rational>> result_data)
{
    if(result_data.empty()) return "方程组无解。";

    std::vector<rational> particular_solution = result_data.back();
    result_data.pop_back();

    if(result_data.empty())
    {
        return QString("方程组有唯一解:\nx = %1").arg(formatVector_Rational(particular_solution));
    }

    QString result = "方程组有无穷多解，通解为:\nx = ";
    result += formatVector_Rational(particular_solution) + " + \n    ";
    for(size_t i = 0; i < result_data.size(); ++i)
    {
        result += QString("c%1 * ").arg(i + 1) + formatVector_Rational(result_data[i]);
        if(i < result_data.size() - 1) result += " + \n    ";
    }
    result += "\n(其中 c1, c2, ... 为任意常数)";
    return result;
}


MainWindow::MainWindow(QWidget * parent)
    : QMainWindow(parent)
{
    setWindowIcon(QIcon(":/res/favicon2.ico"));

    setupUi();
    updateMatrixADim();
    updateMatrixBDim();
}

MainWindow::~MainWindow()
{
}

void MainWindow::setupUi()
{
    QWidget * centralWidget = new QWidget(this);
    this->setCentralWidget(centralWidget);
    this->setWindowTitle("矩阵计算器 (By Iviesever)");

    // --- 左侧：矩阵A ---
    QGroupBox * groupA = new QGroupBox("矩阵 A / 方程组系数矩阵 Ax=b", this);
    QVBoxLayout * layoutA = new QVBoxLayout();
    QHBoxLayout * dimLayoutA = new QHBoxLayout();
    rowsASpin = new QSpinBox(this);
    rowsASpin->setRange(1, 20);
    colsASpin = new QSpinBox(this);
    colsASpin->setRange(1, 20);
    dimLayoutA->addWidget(new QLabel("行:", this));
    dimLayoutA->addWidget(rowsASpin);
    dimLayoutA->addWidget(new QLabel("列:", this));
    dimLayoutA->addWidget(colsASpin);
    matrixATable = new QTableWidget(this);
    layoutA->addLayout(dimLayoutA);
    layoutA->addWidget(matrixATable);
    groupA->setLayout(layoutA);

    // --- 右侧：矩阵B ---
    QGroupBox * groupB = new QGroupBox("矩阵 B / 非齐次项 b", this);
    QVBoxLayout * layoutB = new QVBoxLayout();
    QHBoxLayout * dimLayoutB = new QHBoxLayout();
    rowsBSpin = new QSpinBox(this);
    rowsBSpin->setRange(1, 20);
    colsBSpin = new QSpinBox(this);
    colsBSpin->setRange(1, 1); // 默认b为列向量

    // colsBSpin->setReadOnly(true);  // 删除或注释掉这一行
    colsBSpin->setRange(1, 20);      // 确保用户可以设置范围
    colsBSpin->setValue(1);          // 将默认值设为1，方便解方程
    connect(colsBSpin, &QSpinBox::valueChanged, this, &MainWindow::updateMatrixBDim);

    dimLayoutB->addWidget(new QLabel("行:", this));
    dimLayoutB->addWidget(rowsBSpin);
    dimLayoutB->addWidget(new QLabel("列:", this));
    dimLayoutB->addWidget(colsBSpin);
    matrixBTable = new QTableWidget(this);
    layoutB->addLayout(dimLayoutB);
    layoutB->addWidget(matrixBTable);
    groupB->setLayout(layoutB);

    // 在主布局创建之后，添加类型选择器
    QHBoxLayout * topControlsLayout = new QHBoxLayout();
    typeSelector = new QComboBox(this);
    typeSelector->addItems({ "浮点数 (double)", "分数 (rational)" });
    topControlsLayout->addWidget(new QLabel("计算类型:", this));
    topControlsLayout->addWidget(typeSelector);
    topControlsLayout->addStretch();

    // --- 中间：操作按钮 ---
    QGroupBox * groupOps = new QGroupBox("操作", this);
    QGridLayout * opsLayout = new QGridLayout();

    // 创建按钮
    QPushButton * btnAdd = new QPushButton("A + B", this);
    QPushButton * btnSub = new QPushButton("A - B", this);
    QPushButton * btnMul = new QPushButton("A * B", this);
    QPushButton * btnDetA = new QPushButton("det(A)", this);
    QPushButton * btnInvA = new QPushButton("inv(A)", this);
    QPushButton * btnRefA = new QPushButton("ref(A)", this);
    QPushButton * btnRrefA = new QPushButton("rref(A)", this);
    QPushButton * btnEigenA = new QPushButton("lambda(A)", this);

    QPushButton * btnAdjA = new QPushButton("adj(A)", this);
    QPushButton * btnRankA = new QPushButton("rank(A)", this);
    QPushButton * btnTraceA = new QPushButton("trace(A)", this);

    QPushButton * btnSolveHom = new QPushButton("解 Ax = 0", this);
    QPushButton * btnSolveNonHom = new QPushButton("解 Ax = b", this);

    opsLayout->addWidget(btnAdd, 0, 0);
    opsLayout->addWidget(btnSub, 0, 1);
    opsLayout->addWidget(btnMul, 0, 2);
    opsLayout->addWidget(btnDetA, 1, 0);
    opsLayout->addWidget(btnInvA, 1, 1);
    opsLayout->addWidget(btnEigenA, 1, 2);
    opsLayout->addWidget(btnRefA, 2, 0);
    opsLayout->addWidget(btnRrefA, 2, 1);
    opsLayout->addWidget(btnAdjA, 2, 2);

    opsLayout->addWidget(btnRankA, 3, 0);
    opsLayout->addWidget(btnTraceA, 3, 1);

    opsLayout->addWidget(btnSolveHom, 4, 0, 1, 3);
    opsLayout->addWidget(btnSolveNonHom, 5, 0, 1, 3);
    groupOps->setLayout(opsLayout);

    // --- 底部：结果显示 ---
    QGroupBox * groupResult = new QGroupBox("结果", this);
    QVBoxLayout * resultLayout = new QVBoxLayout();
    resultTextEdit = new QPlainTextEdit(this);
    resultTextEdit->setReadOnly(true);
    resultTextEdit->setFont(QFont("Courier New", 10));
    resultLayout->addWidget(resultTextEdit);
    groupResult->setLayout(resultLayout);

    // --- 主布局 ---
    QGridLayout * mainLayout = new QGridLayout(centralWidget);
    mainLayout->addLayout(topControlsLayout, 0, 0, 1, 2);
    mainLayout->addWidget(groupA, 1, 0);
    mainLayout->addWidget(groupB, 1, 1);
    mainLayout->addWidget(groupOps, 2, 0, 1, 2);
    mainLayout->addWidget(groupResult, 3, 0, 1, 2);

    // --- 连接信号和槽 ---
    connect(rowsASpin, &QSpinBox::valueChanged, this, &MainWindow::updateMatrixADim);
    connect(colsASpin, &QSpinBox::valueChanged, this, &MainWindow::updateMatrixADim);
    // 当A的行数改变时，B的行数也自动改变，以满足 Ax=b 的要求
    connect(rowsASpin, &QSpinBox::valueChanged, rowsBSpin, &QSpinBox::setValue);
    connect(rowsBSpin, &QSpinBox::valueChanged, this, &MainWindow::updateMatrixBDim);

    connect(btnAdd, &QPushButton::clicked, this, &MainWindow::onAdd);
    connect(btnSub, &QPushButton::clicked, this, &MainWindow::onSubtract);
    connect(btnMul, &QPushButton::clicked, this, &MainWindow::onMultiply);
    connect(btnDetA, &QPushButton::clicked, this, &MainWindow::onDetA);
    connect(btnInvA, &QPushButton::clicked, this, &MainWindow::onInvA);
    connect(btnRefA, &QPushButton::clicked, this, &MainWindow::onRefA);
    connect(btnRrefA, &QPushButton::clicked, this, &MainWindow::onRrefA);
    connect(btnEigenA, &QPushButton::clicked, this, &MainWindow::onEigenvaluesA);

    connect(btnAdjA, &QPushButton::clicked, this, &MainWindow::onAdjA);
    connect(btnRankA, &QPushButton::clicked, this, &MainWindow::onRankA);
    connect(btnTraceA, &QPushButton::clicked, this, &MainWindow::onTraceA);

    connect(btnSolveHom, &QPushButton::clicked, this, &MainWindow::onSolveHomogeneous);
    connect(btnSolveNonHom, &QPushButton::clicked, this, &MainWindow::onSolveNonHomogeneous);
}

// 将表格解析为分数矩阵的函数
MatrixR MainWindow::tableToMatrix_Rational(QTableWidget * table)
{
    MatrixR m(table->rowCount(), table->columnCount());
    for(int r = 0; r < table->rowCount(); ++r)
    {
        for(int c = 0; c < table->columnCount(); ++c)
        {
            QTableWidgetItem * item = table->item(r, c);
            if(item && !item->text().isEmpty())
            {
                QString text = item->text();
                long long num = 0, den = 1;

                if(text.contains('/'))
                {
                    QStringList parts = text.split('/');
                    if(parts.size() == 2)
                    {
                        num = parts[0].toLongLong();
                        den = parts[1].toLongLong();
                    }
                }
                else
                {
                    num = text.toLongLong();
                }
                m[r, c] = rational(num, den);
            }
            else
            {
                m[r, c] = rational(0, 1); // 空单元格视为0
            }
        }
    }
    return m;
}

void MainWindow::onSolveHomogeneous()
{
    try
    {
        if(typeSelector->currentText() == "浮点数 (double)")
        {
            MatrixD matA = tableToMatrix_Double(matrixATable);
            auto result = matA.solve_homogeneous();
            QString content = formatHomogeneousSolution(result);
            displayResult("解 Ax = 0:", content);
        }
        else
        {
            MatrixR matA = tableToMatrix_Rational(matrixATable);
            auto result = matA.solve_homogeneous();
            QString content = formatHomogeneousSolution_Rational(result);
            displayResult("解 Ax = 0:", content);
        }
    }
    catch(const std::exception & e) { displayError(e); }
}

void MainWindow::onSolveNonHomogeneous()
{
    try
    {
        if(typeSelector->currentText() == "浮点数 (double)")
        {
            MatrixD matA = tableToMatrix_Double(matrixATable);
            MatrixD matB = tableToMatrix_Double(matrixBTable);
            auto result = matA.solve_non_homogeneous(matB);
            QString content = formatNonHomogeneousSolution(result);
            displayResult("解 Ax = b:", content);
        }
        else
        {
            MatrixR matA = tableToMatrix_Rational(matrixATable);
            MatrixR matB = tableToMatrix_Rational(matrixBTable);
            auto result = matA.solve_non_homogeneous(matB);
            QString content = formatNonHomogeneousSolution_Rational(result);
            displayResult("解 Ax = b:", content);
        }
    }
    catch(const std::exception & e) { displayError(e); }
}



void MainWindow::updateMatrixADim()
{
    matrixATable->setRowCount(rowsASpin->value());
    matrixATable->setColumnCount(colsASpin->value());
    matrixATable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
    // 自动同步b的行数
    rowsBSpin->setValue(rowsASpin->value());
}

void MainWindow::updateMatrixBDim()
{
    matrixBTable->setRowCount(rowsBSpin->value());
    matrixBTable->setColumnCount(colsBSpin->value());
    matrixBTable->horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);
}

MatrixD MainWindow::tableToMatrix_Double(QTableWidget * table)
{
    MatrixD m(table->rowCount(), table->columnCount());
    for(int r = 0; r < table->rowCount(); ++r)
    {
        for(int c = 0; c < table->columnCount(); ++c)
        {
            QTableWidgetItem * item = table->item(r, c);
            if(item && !item->text().isEmpty())
            {
                bool ok;
                double val = item->text().toDouble(&ok);
                m[r, c] = ok ? val : 0.0;
            }
            else
            {
                m[r, c] = 0.0;
            }
        }
    }
    return m;
}

void MainWindow::matrixToTable(const MatrixD & m, QTableWidget * table)
{
    table->setRowCount(m.rows());
    table->setColumnCount(m.cols());
    for(size_t r = 0; r < m.rows(); ++r)
    {
        for(size_t c = 0; c < m.cols(); ++c)
        {
            table->setItem(r, c, new QTableWidgetItem(QString::number(m[r, c])));
        }
    }
}

void MainWindow::displayResult(const QString & title, const QString & content)
{
    resultTextEdit->setPlainText(title + "\n" + content);
}

void MainWindow::displayError(const std::exception & e)
{
    resultTextEdit->setPlainText("错误:\n" + QString::fromStdString(e.what()));
}

void MainWindow::onAdjA()
{
    try
    {
        if(typeSelector->currentText() == "浮点数 (double)")
        {
            MatrixD matA = tableToMatrix_Double(matrixATable);
            MatrixD result = matA.adj();
            displayResult("adj(A):", QString::fromStdString(std::format("{}", result)));
        }
        else
        {
            MatrixR matA = tableToMatrix_Rational(matrixATable);
            MatrixR result = matA.adj();
            displayResult("adj(A):", QString::fromStdString(std::format("{}", result)));
        }
    }
    catch(const std::exception & e) { displayError(e); }
}

void MainWindow::onRankA()
{
    try
    {
        if(typeSelector->currentText() == "浮点数 (double)")
        {
            MatrixD matA = tableToMatrix_Double(matrixATable);
            size_t result = matA.r();
            displayResult("rank(A):", QString::number(result));
        }
        else
        {
            MatrixR matA = tableToMatrix_Rational(matrixATable);
            size_t result = matA.r();
            displayResult("rank(A):", QString::number(result));
        }
    }
    catch(const std::exception & e) { displayError(e); }
}

void MainWindow::onTraceA()
{
    try
    {
        if(typeSelector->currentText() == "浮点数 (double)")
        {
            MatrixD matA = tableToMatrix_Double(matrixATable);
            double result = matA.tr();
            displayResult("trace(A):", QString::number(result));
        }
        else
        {
            MatrixR matA = tableToMatrix_Rational(matrixATable);
            rational result = matA.tr(); 
            displayResult("trace(A):", QString::fromStdString(std::format("{}", result)));
        }
    }
    catch(const std::exception & e) { displayError(e); }
}

void MainWindow::onAdd()
{
    try
    {
        if(typeSelector->currentText() == "浮点数 (double)")
        {
            MatrixD matA = tableToMatrix_Double(matrixATable);
            MatrixD matB = tableToMatrix_Double(matrixBTable);
            MatrixD result = matA + matB;
            displayResult("A + B:", QString::fromStdString(std::format("{}", result)));
        }
        else
        {
            MatrixR matA = tableToMatrix_Rational(matrixATable);
            MatrixR matB = tableToMatrix_Rational(matrixBTable);
            MatrixR result = matA + matB;
            displayResult("A + B:", QString::fromStdString(std::format("{}", result)));

        }
    }
    catch(const std::exception & e) { displayError(e); }
}

void MainWindow::onSubtract()
{
    try
    {
        if(typeSelector->currentText() == "浮点数 (double)")
        {
            MatrixD matA = tableToMatrix_Double(matrixATable);
            MatrixD matB = tableToMatrix_Double(matrixBTable);
            MatrixD result = matA - matB;
            displayResult("A - B:", QString::fromStdString(std::format("{}", result)));
        }
        else
        {
            MatrixR matA = tableToMatrix_Rational(matrixATable);
            MatrixR matB = tableToMatrix_Rational(matrixBTable);
            MatrixR result = matA - matB;
            displayResult("A - B:", QString::fromStdString(std::format("{}", result)));
        }
    }
    catch(const std::exception & e) { displayError(e); }
}

void MainWindow::onMultiply()
{
    try
    {
        if(typeSelector->currentText() == "浮点数 (double)")
        {
            MatrixD matA = tableToMatrix_Double(matrixATable);
            MatrixD matB = tableToMatrix_Double(matrixBTable);
            MatrixD result = matA * matB;
            displayResult("A * B:", QString::fromStdString(std::format("{}", result)));
        }
        else
        {
            MatrixR matA = tableToMatrix_Rational(matrixATable);
            MatrixR matB = tableToMatrix_Rational(matrixBTable);
            MatrixR result = matA * matB;
            displayResult("A * B:", QString::fromStdString(std::format("{}", result)));
        }
    }
    catch(const std::exception & e) { displayError(e); }
}

void MainWindow::onDetA()
{
    try
    {
        if(typeSelector->currentText() == "浮点数 (double)")
        {
            MatrixD matA = tableToMatrix_Double(matrixATable);
            double result = matA.det();
            displayResult("det(A):", QString::number(result));
        }
        else
        {
            MatrixR matA = tableToMatrix_Rational(matrixATable);
            rational result = matA.det();
            displayResult("det(A):", QString::fromStdString(std::format("{}", result)));
        }
    }
    catch(const std::exception & e) { displayError(e); }
}

void MainWindow::onInvA()
{
    try
    {
        if(typeSelector->currentText() == "浮点数 (double)")
        {
            MatrixD matA = tableToMatrix_Double(matrixATable);
            MatrixD result = matA.inv();
            displayResult("inv(A):", QString::fromStdString(std::format("{}", result)));
        }
        else
        {
            MatrixR matA = tableToMatrix_Rational(matrixATable);
            MatrixR result = matA.inv();
            displayResult("inv(A):", QString::fromStdString(std::format("{}", result)));
        }
    }
    catch(const std::exception & e) { displayError(e); }
}

void MainWindow::onRefA()
{
    try
    {
        if(typeSelector->currentText() == "浮点数 (double)")
        {
            MatrixD matA = tableToMatrix_Double(matrixATable);
            MatrixD result = matA.ref();
            displayResult("REF(A):", QString::fromStdString(std::format("{}", result)));
        }
        else
        {
            MatrixR matA = tableToMatrix_Rational(matrixATable);
            MatrixR result = matA.ref();
            displayResult("REF(A):", QString::fromStdString(std::format("{}", result)));
        }
    }
    catch(const std::exception & e) { displayError(e); }
}

void MainWindow::onRrefA()
{
    try
    {
        if(typeSelector->currentText() == "浮点数 (double)")
        {
            MatrixD matA = tableToMatrix_Double(matrixATable);
            MatrixD result = matA.rref();
            displayResult("RREF(A):", QString::fromStdString(std::format("{}", result)));
        }
        else
        {
            MatrixR matA = tableToMatrix_Rational(matrixATable);
            MatrixR result = matA.rref();
            displayResult("RREF(A):", QString::fromStdString(std::format("{}", result)));
        }
    }
    catch(const std::exception & e) { displayError(e); }
}

void MainWindow::onEigenvaluesA()
{
    try
    {
        if(typeSelector->currentText() == "浮点数 (double)")
        {
            MatrixD matA = tableToMatrix_Double(matrixATable);
            auto result = matA.lambda();
            QString content = "[\n";
            for(const auto & val : result)
            {
                content += "  " + QString::fromStdString(std::format("{}", val)) + "\n";
            }
            content += "]";
            displayResult("特征值(A):", content);
        }
        else
        {
            // 在分数模式下，求特征值通常没有解析解，需要数值方法，所以切换到double计算
            displayResult("提示:", "特征值计算通常在浮点数域进行。\n请在“浮点数 (double)”模式下使用此功能。");
        }
    }
    catch(const std::exception & e) { displayError(e); }
}