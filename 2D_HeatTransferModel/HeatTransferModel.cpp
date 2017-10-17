#include <time.h>
#include <iostream>
#include <math.h>
#include <stdlib.h> 
#include <fstream>
#include <random>
#include <classstatement.h>
using namespace std;
int main()
{
	const int section = 12, coolsection = 8, moldsection = 4, measurednumb = 8, predictstep = 5;
	float ccml[section + 1] = { 0.0f,0.2f,0.4f,0.6f,0.8f,1.0925f,2.27f,4.29f,5.831f,9.6065f,13.6090f,19.87014f,28.599f };
	float measuredpoistion[measurednumb] = { 0.9463f, 1.6812f, 3.28f, 5.0605f, 7.7188f, 11.6077f, 16.7395f, 24.235f };
	float hinit[section] = { 1380.0f,1170.0f,980.0f,800.0f,1223.16f,735.05f,424.32f,392.83f,328.94f,281.64f,246.16f,160.96f };
	float taimmeantemperature[coolsection] = { 966.149841f, 925.864746f, 952.322083f, 932.175537f, 914.607117f, 890.494263f, 870.804443f, 890.595825f };
	float measuredtemperature[measurednumb] = { 937.141f, 930.948f, 960.807f, 932.294f, 916.781f, 892.089f, 872.358f, 899.282f };
	float *htemp = new float[section];
	float *hresult = new float[section];
	float m_lx = 0.25f, m_ly = 1.79f, m_lz = 0.25f, m_tf = 1400.0f, m_vcast = -0.02f, rangeh = 50.0f, T_Cast = 1558.0f, dh = 1.0f, m_tao = 0.25f;
	int m_nx = 25, m_ny = 179, m_nz = 170, m_tnpts = 1501, sim_tnpts, maxiter_time = 100;
	int popsize = 30;
	float **allmeantemperature;
	float *staticmeantemperature;
	ContinuousCaster CasterOne = ContinuousCaster(section, coolsection, moldsection, ccml);
	Steel steel;
	Gradientbasedalgorithm Gradientmethod = Gradientbasedalgorithm(CasterOne, measuredtemperature);
	PSOalgorithm PSO = PSOalgorithm(measurednumb, popsize, rangeh, hinit, CasterOne);

	/*Temperature1d SteelTemperature1d = Temperature1d(m_nx, m_tao, m_lx, m_vcast, CasterOne, steel);
	Temperature1d SteelTemperature1dtemp = Temperature1d(m_nx, m_tao, m_lx, m_vcast, CasterOne, steel);

	allmeantemperature = new float*[CasterOne.coolsection];
	for (int i = 0; i < CasterOne.coolsection; i++)
		allmeantemperature[i] = new float[measurednumb];
	staticmeantemperature = new float[measurednumb];

	SteelTemperature1d.initcondition1d();
	SteelTemperature1dtemp.initcondition1d();

	for (int i = 0; i < CasterOne.section; i++)
		htemp[i] = hinit[i] + rangeh * rand() / RAND_MAX;

	for (int iter_time = 0; iter_time < maxiter_time; iter_time++)
	{
		PSO.fitnessevaluation(measuredtemperature, measuredpoistion, SteelTemperature1dtemp);
		PSO.findgbest();
		PSO.updatepoph();
		cout << endl;
	}*/

	
	Temperature1d SteelTemperature1d = Temperature1d(m_nx, m_tao, m_lx, m_vcast, CasterOne, steel);
	Temperature1d SteelTemperature1dtemp = Temperature1d(m_nx, m_tao, m_lx, m_vcast, CasterOne, steel);
	float noisemean = 0.0, noisestd = 1.0;
	Generatemeasuredtemperature gen_measuredtemperature = Generatemeasuredtemperature(CasterOne, measurednumb, noisemean, noisestd, hinit, measuredpoistion);
	gen_measuredtemperature.simulationtemperature(SteelTemperature1d, measuredtemperature);

	allmeantemperature = new float*[CasterOne.coolsection];
	for (int i = 0; i < CasterOne.coolsection; i++)
		allmeantemperature[i] = new float[measurednumb];
	staticmeantemperature = new float[measurednumb];

	SteelTemperature1d.initcondition1d();
	SteelTemperature1dtemp.initcondition1d();

	for (int i = 0; i < CasterOne.section; i++)
		htemp[i] = hinit[i] + rangeh * rand() / RAND_MAX;

	for (int iter_time = 0; iter_time < maxiter_time; iter_time++)
	{
		for (int i = 0; i < measurednumb + 1; i++)
		{
			SteelTemperature1dtemp.initcondition1d();
			if (i == measurednumb)
			{
				for (int j = 0; j < CasterOne.section; j++)
					hresult[j] = htemp[j];
				while (SteelTemperature1dtemp.tstep < SteelTemperature1dtemp.tnpts)
				{
					SteelTemperature1dtemp.differencecalculation1d(hresult);
				}
				SteelTemperature1dtemp.computetemperature1d(measuredpoistion, measurednumb);
				for (int j = 0; j < CasterOne.coolsection; j++)
					staticmeantemperature[j] = SteelTemperature1dtemp.computetemperature[j];
			}
			else
			{
				for (int j = 0; j < CasterOne.section; j++)
					hresult[j] = htemp[j];
				hresult[moldsection + i] = hresult[moldsection + i] + dh;
				while (SteelTemperature1dtemp.tstep < SteelTemperature1dtemp.tnpts)
				{
					SteelTemperature1dtemp.differencecalculation1d(hresult);
				} 
				SteelTemperature1dtemp.computetemperature1d(measuredpoistion, measurednumb);
				for (int j = 0; j < measurednumb; j++)
					allmeantemperature[i][j] = SteelTemperature1dtemp.computetemperature[j];
			}
		}

		Gradientmethod.init(allmeantemperature, staticmeantemperature, dh);
		Gradientmethod.gradientcalculation();
		Gradientmethod.linesearch();
		Gradientmethod.updateh(htemp);
		Gradientmethod.print();

		cout << "htemp = " << endl;
		for (int i = 0; i < CasterOne.section; i++)
			cout << htemp[i] << ", ";
	}


	/*clock_t t_start = clock();
	Temperature2d SteelTemperature2d = Temperature2d(m_nx, m_ny, m_tao, m_lx, m_ly, m_vcast, CasterOne, steel);
	Temperature2d SteelTemperature2dtemp = Temperature2d(m_nx, m_ny, m_tao, m_lx, m_ly, m_vcast, CasterOne, steel);

	allmeantemperature = new float*[CasterOne.coolsection];
	for (int i = 0; i < CasterOne.coolsection; i++)
	allmeantemperature[i] = new float[measurednumb];
	staticmeantemperature = new float[measurednumb];

	SteelTemperature2d.initcondition2d();
	SteelTemperature2dtemp.initcondition2d();

	for (int i = 0; i < CasterOne.section; i++)
	htemp[i] = hinit[i] + rangeh * rand() / RAND_MAX;

	for (int iter_time = 0; iter_time < maxiter_time; iter_time++)
	{
	for (int i = 0; i < measurednumb + 1; i++)
	{
	SteelTemperature2dtemp.initcondition2d();
	if (i == measurednumb)
	{
	for (int j = 0; j < CasterOne.section; j++)
	hresult[j] = htemp[j];
	while (SteelTemperature2dtemp.tstep < SteelTemperature2dtemp.tnpts)
	{
	SteelTemperature2dtemp.differencecalculation2d(hresult);
	}
	SteelTemperature2dtemp.computetemperature2d(measuredpoistion, measurednumb);
	for (int j = 0; j < CasterOne.coolsection; j++)
	staticmeantemperature[j] = SteelTemperature2dtemp.computetemperature[j];
	}
	else
	{
	for (int j = 0; j < CasterOne.section; j++)
	hresult[j] = htemp[j];
	hresult[moldsection + i] = hresult[moldsection + i] + dh;
	while (SteelTemperature2dtemp.tstep < SteelTemperature2dtemp.tnpts)
	{
	SteelTemperature2dtemp.differencecalculation2d(hresult);
	}
	SteelTemperature2dtemp.computetemperature2d(measuredpoistion, measurednumb);
	for (int j = 0; j < measurednumb; j++)
	allmeantemperature[i][j] = SteelTemperature2dtemp.computetemperature[j];
	}
	}

	Gradientmethod.init(allmeantemperature, staticmeantemperature, dh);
	Gradientmethod.gradientcalculation();
	Gradientmethod.linesearch();
	Gradientmethod.updateh(htemp);
	Gradientmethod.print();

	cout << "htemp = " << endl;
	for (int i = 0; i < CasterOne.section; i++)
	    cout << htemp[i] << ", ";
	cout << endl;
	}
	clock_t t_end = clock(); 
    cout << "The running time is " << (t_end - t_start) << " (ms)"<< endl;*/

	/*Temperature2d SteelTemperature2d = Temperature2d(m_nx, m_ny, m_tao, m_lx, m_ly, m_vcast, CasterOne, steel);
	SteelTemperature2d.initcondition2d();

	while (SteelTemperature2d.tstep < SteelTemperature2d.tnpts)
	{
		SteelTemperature2d.print2d();
		SteelTemperature2d.differencecalculation2d(hinit);
	}
	SteelTemperature2d.computetemperature2d(measuredpoistion, measurednumb);*/

	/*clock_t t_start = clock();
    m_nx = 25;
	m_ny = 3001;
	m_nz = 25;
	m_lx = 0.25;
	m_ly = 28.599f;
	m_lz = 0.25;
	m_tnpts = 10001;
	m_tf = 2000.0f;
	sim_tnpts = 0;
	Temperature SteelTemperature3dplant = Temperature(m_nx, m_ny, m_nz, m_tnpts, m_tf, m_lx, m_ly, m_lz, m_vcast, CasterOne, steel);
	Temperature SteelTemperature3dmodel = Temperature(m_nx, m_ny, m_nz, m_tnpts, m_tf, m_lx, m_ly, m_lz, m_vcast, CasterOne, steel);
	Temperature SteelTemperature3dtemp = Temperature(m_nx, m_ny, m_nz, m_tnpts, m_tf, m_lx, m_ly, m_lz, m_vcast, CasterOne, steel);
	
	allmeantemperature = new float* [CasterOne.coolsection];
	for (int i = 0; i < CasterOne.coolsection; i++)
		allmeantemperature[i] = new float[CasterOne.coolsection];
	staticmeantemperature = new float[CasterOne.coolsection];
		
	SteelTemperature3dplant.initcondition3d(T_Cast);	
	while (SteelTemperature3dplant.tstep <= sim_tnpts)
	{
		if (SteelTemperature3dplant.tstep % 100 == 0)
		{
			SteelTemperature3dplant.computemeantemperature3d();
			SteelTemperature3dplant.print3d();
		}
		SteelTemperature3dplant.differencecalculation3d(hinit);
	}

	SteelTemperature3dtemp = SteelTemperature3dplant;
	SteelTemperature3dmodel = SteelTemperature3dplant;

	while (SteelTemperature3dplant.tstep < m_tnpts && SteelTemperature3dplant.tstep > sim_tnpts)
	{	
		if (SteelTemperature3dplant.tstep % 1 == 0)
		{
			for (int i = 0; i < CasterOne.coolsection + 1; i++)
			{
				SteelTemperature3dtemp = SteelTemperature3dmodel;
				SteelTemperature3dtemp.computemeantemperature3d();

				if (i == CasterOne.coolsection)
				{
					for (int j = 0; j < CasterOne.section; j++)
						htemp[j] = hinit[j];
					for (int p = 0; p < predictstep; p++)
					    SteelTemperature3dtemp.differencecalculation3d(htemp);
					SteelTemperature3dtemp.computemeantemperature3d();
					for (int j = 0; j < CasterOne.coolsection; j++)
					    staticmeantemperature[j] = SteelTemperature3dtemp.meantemperature[j];
				}
				else 
				{
					for (int j = 0; j < CasterOne.section; j++)
						htemp[j] = hinit[j];
					htemp[moldsection + i] = htemp[moldsection + i] + dh;
					for (int p = 0; p < predictstep; p++)
					    SteelTemperature3dtemp.differencecalculation3d(htemp);
					SteelTemperature3dtemp.computemeantemperature3d();
					for (int j = 0; j < CasterOne.coolsection; j++)
						allmeantemperature[i][j] = SteelTemperature3dtemp.meantemperature[j];
				}
			}
			Gradientmethod.init(allmeantemperature, staticmeantemperature, dh);
			Gradientmethod.gradientcalculation();
			Gradientmethod.linesearch();
			Gradientmethod.updateh(hinit);
			//Gradientmethod.print();
		}
		SteelTemperature3dplant.differencecalculation3d(hinit);
		SteelTemperature3dmodel.differencecalculation3d(hinit);
		SteelTemperature3dplant.computemeantemperature3d();
		SteelTemperature3dmodel.computemeantemperature3d();
		if (SteelTemperature3dplant.tstep % 10 == 0)
		{
			cout << "model timestep = " << SteelTemperature3dmodel.tstep << " " << endl;
			SteelTemperature3dmodel.print3d();
		}
		Gradientmethod.outputdata(SteelTemperature3dplant.tstep, hinit);
	}
   clock_t t_end = clock();*/
   //cout << "The running time is " << (t_end - t_start) << " (ms)"<< endl;
}


