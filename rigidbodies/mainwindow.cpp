#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "simparameters.h"
#include "controller.h"

MainWindow::MainWindow(Controller &cont, int fps, QWidget *parent) :
    QMainWindow(parent),
    cont_(cont),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->GLWidget->setController(&cont);
    simRunning_ = false;
    connect(&renderTimer_, SIGNAL(timeout()), this, SLOT(updateGL()));
    renderTimer_.start(1000/fps);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_actionExit_triggered()
{
    close();
}

void MainWindow::setParametersFromUI()
{
    SimParameters params;

    params.simRunning = simRunning_;

    params.timeStep = ui->timeStepEdit->text().toDouble();
    params.NewtonTolerance = ui->newtonTolEdit->text().toDouble();
    params.NewtonMaxIters = ui->newtonMaxItersEdit->text().toInt();
    params.penaltyStiffness = ui->penaltyStiffnessEdit->text().toDouble();

    params.activeForces = 0;
    if(ui->gravityCheckBox->isChecked())
        params.activeForces |= SimParameters::F_GRAVITY;

    params.gravityG = ui->gravityGEdit->text().toDouble();

    params.bodyDensity = ui->densityEdit->text().toDouble();

    if(ui->sphereButton->isChecked())
        params.launchBody = SimParameters::R_SPHERE;
    else if(ui->twoByFourButton->isChecked())
        params.launchBody = SimParameters::R_2BY4;
    else if(ui->bunnyButton->isChecked())
        params.launchBody = SimParameters::R_BUNNY;
    else if(ui->customButton->isChecked())
        params.launchBody = SimParameters::R_CUSTOM;
    else if(ui->planeButton->isChecked())
        params.launchBody = SimParameters::R_PLANE;

    params.launchVel = ui->launchVelEdit->text().toDouble();
    params.randomLaunchAngVel = ui->randomAngularVelCheckBox->isChecked();
    params.randomLaunchOrientation = ui->randomOrienatationCheckBox->isChecked();
    params.randomLaunchVelMagnitude = ui->randomVelMagEdit->text().toDouble();

    if(ui->cageNever->isChecked())
        params.useCage = SimParameters::C_NEVER;
    else if(ui->cageAlways->isChecked())
        params.useCage = SimParameters::C_ALWAYS;
    else if(ui->cageBroad->isChecked())
        params.useCage = SimParameters::C_BROADPHASE;


    setUIFromParameters(params);
    QMetaObject::invokeMethod(&cont_, "updateParameters", Q_ARG(SimParameters, params));
}

void MainWindow::setUIFromParameters(const SimParameters &params)
{
    if(params.simRunning)
    {
        ui->startSimulationButton->setText(QString("Pause Simulation"));
        simRunning_ = true;
    }
    else
    {
        ui->startSimulationButton->setText(QString("Start Simulation"));
        simRunning_ = false;
    }

    ui->timeStepEdit->setText(QString::number(params.timeStep));
    ui->newtonTolEdit->setText(QString::number(params.NewtonTolerance));
    ui->newtonMaxItersEdit->setText(QString::number(params.NewtonMaxIters));
    ui->penaltyStiffnessEdit->setText(QString::number(params.penaltyStiffness));

    ui->gravityCheckBox->setChecked(params.activeForces & SimParameters::F_GRAVITY);
    ui->gravityGEdit->setText(QString::number(params.gravityG));       

    ui->densityEdit->setText(QString::number(params.bodyDensity));

    switch(params.launchBody)
    {
    case SimParameters::R_SPHERE:
        ui->sphereButton->setChecked(true);
        break;
    case SimParameters::R_2BY4:
        ui->twoByFourButton->setChecked(true);
        break;
    case SimParameters::R_BUNNY:
        ui->bunnyButton->setChecked(true);
        break;
    case SimParameters::R_CUSTOM:
        ui->customButton->setChecked(true);
        break;
    case SimParameters::R_PLANE:
        ui->planeButton->setChecked(true);
        break;
    }

    switch(params.useCage)
    {
    case SimParameters::C_NEVER:
        ui->cageNever->setChecked(true);
        break;
    case SimParameters::C_ALWAYS:
        ui->cageAlways->setChecked(true);
        break;
    case SimParameters::C_BROADPHASE:
        ui->cageBroad->setChecked(true);
    }

    ui->launchVelEdit->setText(QString::number(params.launchVel));
    ui->randomOrienatationCheckBox->setChecked(params.randomLaunchOrientation);
    ui->randomAngularVelCheckBox->setChecked(params.randomLaunchAngVel);
    ui->randomVelMagEdit->setText(QString::number(params.randomLaunchVelMagnitude));
}

void MainWindow::updateGL()
{
    ui->GLWidget->tick();
    ui->GLWidget->update();
}

void MainWindow::on_actionReset_Everything_triggered()
{
    QMetaObject::invokeMethod(&cont_, "reset");
}

void MainWindow::on_actionReset_triggered()
{
    QMetaObject::invokeMethod(&cont_, "clearScene");
}

void MainWindow::on_startSimulationButton_clicked()
{
    simRunning_ = !simRunning_;
    setParametersFromUI();
}

void MainWindow::on_timeStepEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_newtonTolEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_newtonMaxItersEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_gravityCheckBox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_gravityGEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_sphereButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_twoByFourButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_bunnyButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_customButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_launchVelEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_randomOrienatationCheckBox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_randomAngularVelCheckBox_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_randomVelMagEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_densityEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_penaltyStiffnessEdit_editingFinished()
{
    setParametersFromUI();
}

void MainWindow::on_planeButton_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_cageNever_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_cageAlways_clicked()
{
    setParametersFromUI();
}

void MainWindow::on_cageBroad_clicked()
{
    setParametersFromUI();
}
