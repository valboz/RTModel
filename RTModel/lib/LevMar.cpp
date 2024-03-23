// LevMar.cpp : main project file.
// This program performs Levenberg-Marquardt fitting on a given initial condition

#include "LevMarFit.h"

int main(int argc, char *argv[])
{
	LevMar *MyLevMar= new LevMar(argc,argv);
	MyLevMar->Run();
	delete MyLevMar;
    return 0;
}
