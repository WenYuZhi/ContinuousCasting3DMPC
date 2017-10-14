#include <time.h>
#include <iostream>
#include <math.h>
#include <stdlib.h> 
#include <fstream>
using namespace std;
class ContinuousCaster
{
    public:
		int section, coolsection, moldsection;
		float *ccml;
		ContinuousCaster(int, int, int ,float*);
		~ContinuousCaster();
		void print();
};

ContinuousCaster::ContinuousCaster(int parasection, int paracoolsection, int paramoldsection, float* paraccml)
{
	section = parasection;
	coolsection = paracoolsection;
	moldsection = paramoldsection;
	ccml = new float [section + 1];
	for (int i = 0; i < section + 1; i++)
		ccml[i] = paraccml[i];
}

ContinuousCaster::~ContinuousCaster()
{
	delete ccml;
}

void ContinuousCaster:: print()
{
	cout << "section = " << section << " ";
	cout << "coolsection = " << coolsection << " ";
	cout << "moldsection = " << moldsection << " " << endl;
	for (int i = 0; i < section; i++)
		cout << ccml[i] << ", ";
}

class Steel
{
	public:
	    float pho;
	    float ce;
	    float lamda;
		void physicalpara(float);
};

class Temperature
{
    private:
		float vcast;
		float h;
		float* T_New;
		float* T_Last;
		float* T_Surface;
		bool disout;
		int nx, ny, nz, tnpts;
		float dx, dy, dz, tao, tf, lx, ly, lz;
		ContinuousCaster* mCasterOne;
		Steel* steel;	
    public:
		int tstep;
		float* meantemperature;
		float* computetemperature;
		Temperature(int, int, int, int, float, float, float, float, float, ContinuousCaster &, Steel &);
		Temperature(const Temperature &);
		~Temperature();
		void differencecalculation3d(float*);
		void boundarycondition3d(ContinuousCaster & CasterOne, float*, int);
		void initcondition3d(float);
		void initcondition3d(float*);
		void print3d(int);
		void print3d();
		void computetemperature3d(float *, int);
		void computemeantemperature3d();
		void operator=(const Temperature &);
};

class Temperature2d
{
    private:
	    float vcast;
	    float h;
	    float* T_New;
	    float* T_Last;
	    float* T_Surface;
	    bool disout;
	    int nx, ny, tnpts;
	    float dx, dy, tao, tf, lx, ly;
	    ContinuousCaster* mCasterOne;
	    Steel* steel;
    public:
		int tstep;
		float* computetemperature;
	    Temperature2d(int, int, int, float, float, float, float, ContinuousCaster &, Steel &);
		~Temperature2d();
		void differencecalculation2d(float*);
		void boundarycondition2d(ContinuousCaster & CasterOne, float*);
		void initcondition2d();
		void computetemperature2d(float *, int);
};

Temperature2d::Temperature2d(int m_nx, int m_ny, int m_tnpts, float m_tf, float m_lx, float m_ly, float m_vcast, ContinuousCaster & m_CasterOne, Steel & m_steel)
{
	nx = m_nx;
	ny = m_ny;
	tnpts = m_tnpts;
	tf = m_tf;
	lx = m_lx;
	ly = m_ly;
	dx = m_lx / float(m_nx - 1);
	dy = m_ly / float(m_ny - 1);
	tao = m_tf / float(m_tnpts - 1);
	T_New = new float [nx * ny];
	T_Last = new float [nx * ny];
	T_Surface = new float[tnpts];
	vcast = m_vcast;
	tstep = 0;
	disout = true;
}

Temperature2d::~Temperature2d()
{
	delete [] T_New;
	delete [] T_Last;
	delete [] T_Surface;
	delete [] computetemperature;
}

void Temperature2d::differencecalculation2d(float* hinit)  
{
	float a, Tw = 30.0, T_Up, T_Down, T_Right, T_Left, T_Middle;
	this->boundarycondition2d(*mCasterOne, hinit);
	if (disout == 0)
	{
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				steel->physicalpara(T_Last[i * ny + j]);
				a = steel->lamda / (steel->pho * steel->ce);
				if (i == 0 && j != 0 && j != ny - 1)  //1
				{
					T_Up = T_Last[(i + 1)*ny + j];
					T_Down = T_Last[(i + 1)*ny + j] - 2 * dx * h * (T_Last[i*ny + j] - Tw) / steel->lamda;
					T_Right = T_Last[i*ny + j + 1];
					T_Left = T_Last[i*ny + j - 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == nx - 1 && j != 0 && j != ny - 1)//2
				{
					T_Up = T_Last[(i - 1)*ny + j] - 2 * dx*h*(T_Last[i*ny + j] - Tw) / steel->lamda;
					T_Down = T_Last[(i - 1)*ny + j];
					T_Right = T_Last[i*ny + j + 1];
					T_Left = T_Last[i*ny + j - 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (j == 0 && i != 0 && i != nx - 1)//3
				{
					T_Up = T_Last[(i + 1)*ny + j];
					T_Down = T_Last[(i - 1)*ny + j];
					T_Right = T_Last[i*ny + j + 1];
					T_Left = T_Last[i*ny + j + 1] - 2 * dy*h*(T_Last[i*ny + j] - Tw) / steel->lamda;
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (j == ny - 1 && i != 0 && i != nx - 1)//4
				{
					T_Up = T_Last[(i + 1)*ny + j];
					T_Down = T_Last[(i - 1)*ny + j];
					T_Right = T_Last[i*ny + j - 1] - 2 * dy*h*(T_Last[i*ny + j] - Tw) / steel->lamda;
					T_Left = T_Last[i*ny + j - 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == 0 && j == 0)//5
				{
					T_Up = T_Last[(i + 1)*ny + j];
					T_Down = T_Last[(i + 1)*ny + j];
					T_Right = T_Last[i*ny + j + 1];
					T_Left = T_Last[i*ny + j + 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == 0 && j == ny - 1)//6
				{
					T_Up = T_Last[(i + 1)*ny + j];
					T_Down = T_Last[(i + 1)*ny + j];
					T_Right = T_Last[i*ny + j - 1];
					T_Left = T_Last[i*ny + j - 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == nx - 1 && j == 0)//7
				{
					T_Up = T_Last[(i - 1)*ny + j];
					T_Down = T_Last[(i - 1)*ny + j];
					T_Right = T_Last[i*ny + j + 1];
					T_Left = T_Last[i*ny + j + 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == nx - 1 && j == ny - 1)//8
				{
					T_Up = T_Last[(i - 1)*ny + j];
					T_Down = T_Last[(i - 1)*ny + j];
					T_Right = T_Last[i*ny + j - 1];
					T_Left = T_Last[i*ny + j - 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else//9
				{
					T_Up = T_Last[(i + 1)*ny + j];
					T_Down = T_Last[(i - 1)*ny + j];
					T_Right = T_Last[i*ny + j + 1];
					T_Left = T_Last[i*ny + j - 1];
					T_Middle = T_Last[i*ny + j];
					T_New[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
			}
		}
	}
	else
	{
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < ny; j++)
			{
				steel->physicalpara(T_New[i * ny + j]);
				a = steel->lamda / (steel->pho * steel->ce);
				if (i == 0 && j != 0 && j != ny - 1)  //1
				{
					T_Up = T_New[(i + 1)*ny + j];
					T_Down = T_New[(i + 1)*ny + j] - 2 * dx*h*(T_New[i*ny + j] - Tw) / steel->lamda;
					T_Right = T_New[i*ny + j + 1];
					T_Left = T_New[i*ny + j - 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == nx - 1 && j != 0 && j != ny - 1)//2
				{
					T_Up = T_New[(i - 1)*ny + j] - 2 * dx*h*(T_New[i*ny + j] - Tw) / steel->lamda;
					T_Down = T_New[(i - 1)*ny + j];
					T_Right = T_New[i*ny + j + 1];
					T_Left = T_New[i*ny + j - 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (j == 0 && i != 0 && i != nx - 1)//3
				{
					T_Up = T_New[(i + 1)*ny + j];
					T_Down = T_New[(i - 1)*ny + j];
					T_Right = T_New[i*ny + j + 1];
					T_Left = T_New[i*ny + j + 1] - 2 * dy*h*(T_New[i*ny + j] - Tw) / steel->lamda;
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (j == ny - 1 && i != 0 && i != nx - 1)//4
				{
					T_Up = T_New[(i + 1)*ny + j];
					T_Down = T_New[(i - 1)*ny + j];
					T_Right = T_New[i*ny + j - 1] - 2 * dy*h*(T_New[i*ny + j] - Tw) / steel->lamda;
					T_Left = T_New[i*ny + j - 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == 0 && j == 0)//5
				{
					T_Up = T_New[(i + 1)*ny + j];
					T_Down = T_New[(i + 1)*ny + j];
					T_Right = T_New[i*ny + j + 1];
					T_Left = T_New[i*ny + j + 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == 0 && j == ny - 1)//6
				{
					T_Up = T_New[(i + 1)*ny + j];
					T_Down = T_New[(i + 1)*ny + j];
					T_Right = T_New[i*ny + j - 1];
					T_Left = T_New[i*ny + j - 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == nx - 1 && j == 0)//7
				{
					T_Up = T_New[(i - 1)*ny + j];
					T_Down = T_New[(i - 1)*ny + j];
					T_Right = T_New[i*ny + j + 1];
					T_Left = T_New[i*ny + j + 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else if (i == nx - 1 && j == ny - 1)//8
				{
					T_Up = T_New[(i - 1)*ny + j];
					T_Down = T_New[(i - 1)*ny + j];
					T_Right = T_New[i*ny + j - 1];
					T_Left = T_New[i*ny + j - 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
				else//9
				{
					T_Up = T_New[(i + 1)*ny + j];
					T_Down = T_New[(i - 1)*ny + j];
					T_Right = T_New[i*ny + j + 1];
					T_Left = T_New[i*ny + j - 1];
					T_Middle = T_New[i*ny + j];
					T_Last[i * ny + j] = (a / (dx*dx))*T_Up + (a / (dx*dx))*T_Down + (1 - 2 * a / (dy*dy) - 2 * a / (dx*dx))*T_Middle + (a / (dy*dy))*T_Right + (a / (dy*dy))*T_Left;
				}
			}
		}
	}
	disout = !disout;
	tstep++;
	T_Surface[tstep] = T_Last[int((ny - 1) / 2)];
}

void Temperature2d::boundarycondition2d(ContinuousCaster & CasterOne, float *hinit)
{
	float zposition = tstep * tao * vcast;
	//cout << "z = " << zposition;
	for (int i = 0; i < CasterOne.section; i++)
		if (zposition >= *(CasterOne.ccml + i) && zposition <= *(CasterOne.ccml + i + 1))
			h = *(hinit + i);
	//cout << "h = " << h << endl;
}

void Temperature2d::initcondition2d()
{
	float T_Cast = 1558.0f;
	tstep = 0;
	for (int i = 0; i < nx; i++)
		for (int j = 0; j < ny; j++)
		{
			T_Last[ny * i + j] = T_Cast;
			T_New[ny * i + j] = T_Cast;
		}
	disout = 0;
}

void Temperature2d::computetemperature2d(float *measuredpoistion, int measurednumb)
{
	computetemperature = new float[measurednumb];
	for (int i = 0; i < measurednumb; i++)
		for (int j = 0; j < tnpts; j++)
		{
			if ((fabs(j * vcast * tao - *(measuredpoistion + i))) <= tao)
				computetemperature[i] = T_Surface[j];
		}
	cout << endl << "computetemperature = " << endl;
	for (int i = 0; i < measurednumb; i++)
		cout << computetemperature[i] << ", ";
	cout << endl;
}

Temperature::Temperature(int m_nx, int m_ny, int m_nz, int m_tnpts, float m_tf, float m_lx, float m_ly, float m_lz, float m_vcast, ContinuousCaster & m_CasterOne, Steel & m_steel)
{
	mCasterOne = &m_CasterOne;
	steel = &m_steel;
	nx = m_nx;
	ny = m_ny;
	nz = m_nz;
	tnpts = m_tnpts;
	tf = m_tf;
	lx = m_lx;
	ly = m_ly;
	lz = m_lz;
	dx = m_lx / float(m_nx - 1);
	dy = m_ly / float(m_ny - 1);
	dz = m_lz / float(m_nz - 1);
	tao = m_tf / float(m_tnpts - 1);
	T_New = new float[nx * ny * nz];
	T_Last = new float[nx * ny * nz];
	T_Surface = new float[ny];
	meantemperature = new float[mCasterOne->coolsection];
	vcast = m_vcast;
	tstep = 0;
	disout = true;
}

Temperature::Temperature(const Temperature & m_SteelTemperature)
{
	mCasterOne = m_SteelTemperature.mCasterOne;
	steel = m_SteelTemperature.steel;
	nx = m_SteelTemperature.nx;
	ny = m_SteelTemperature.ny;
	nz = m_SteelTemperature.nz;
	tnpts = m_SteelTemperature.tnpts;
	tf = m_SteelTemperature.tf;
	lx = m_SteelTemperature.lx;
	ly = m_SteelTemperature.ly;
	lz = m_SteelTemperature.lz;
	dx = m_SteelTemperature.lx / float(m_SteelTemperature.nx - 1);
	dy = m_SteelTemperature.ly / float(m_SteelTemperature.ny - 1);
	dz = m_SteelTemperature.lz / float(m_SteelTemperature.nz - 1);
	tao = m_SteelTemperature.tf / float(m_SteelTemperature.tnpts - 1);

	T_New = new float[nx * ny * nz];
	T_Last = new float[nx * ny * nz];
	for (int j = 0; j < ny; j++)
		for (int i = 0; i < nx; i++)
			for (int k = 0; k < nz; k++)
			{
				T_Last[nx * nz * j + nz * i + k] = m_SteelTemperature.T_Last[nx * nz * j + nz * i + k];
				T_New[nx * nz * j + nz * i + k] = m_SteelTemperature.T_New[nx * nz * j + nz * i + k];
			}
	T_Surface = new float[ny];
	for (int j = 0; j < ny; j++)
		T_Surface[j] = m_SteelTemperature.T_Surface[j];

	meantemperature = new float[mCasterOne->coolsection];
	for (int i = 0; i < mCasterOne->coolsection; i++)
		meantemperature[i] = m_SteelTemperature.meantemperature[i];

	vcast = m_SteelTemperature.vcast;
	tstep = m_SteelTemperature.tstep;
	disout = m_SteelTemperature.disout;
	h = m_SteelTemperature.h;
}

Temperature::~Temperature()
{
	delete [] T_New;
	delete [] T_Last;
	delete [] T_Surface;
	delete [] meantemperature;
}


void Temperature::differencecalculation3d(float *hinit)
{
	float a, Tw = 30.0, T_Up, T_Down, T_Right, T_Left, T_Forw, T_Back, T_Middle, T_Cast = 1558.0f;
	if (disout)
	{
		for (int j = 0; j < ny; j++)
		{
			this->boundarycondition3d(*mCasterOne, hinit, j);
			for (int i = 0; i < nx; i++)
				for (int m = 0; m < nz; m++)
				{
					steel->physicalpara(T_Last[nx * nz * j + nz * i + m]);
					a = steel->lamda / (steel->pho * steel->ce);
					if (j == 0 && i != 0 && i != (nx - 1) && m != 0 && m != (nz - 1)) //1
					{
						T_New[nx * nz * j + nz * i + m] = T_Cast;
					}

					else if (j == 0 && i == 0 && m != 0 && m != (nz - 1)) //2
					{
						T_New[nx * nz * j + nz * i + m] = T_Cast;
					}

					else if (j == 0 && i == (nx - 1) && m != 0 && m != (nz - 1))//3
					{
						T_New[nx * nz * j + nz * i + m] = T_Cast;
					}

					else if (j == 0 && i != 0 && i != (nx - 1) && m == 0) //4
					{
						T_New[nx * nz * j + nz * i + m] = T_Cast;
					}

					else if (j == 0 && i != 0 && i != (nx - 1) && m == (nz - 1)) //5
					{
						T_New[nx * nz * j + nz * i + m] = T_Cast;
					}

					else if (j == 0 && i == 0 && m == 0)  //6
					{
						T_New[nx * nz * j + nz * i + m] = T_Cast;
					}

					else if (j == 0 && i == 0 && m == (nz - 1))  //7
					{
						T_New[nx * nz * j + nz * i + m] = T_Cast;
					}

					else if (j == 0 && i == (nx - 1) && m == 0)  //8
					{
						T_New[nx * nz * j + nz * i + m] = T_Cast;
					}

					else if (j == 0 && i == (nx - 1) && m == (nz - 1)) //9
					{
						T_New[nx * nz * j + nz * i + m] = T_Cast;
					}

					else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m != 0 && m != (nz - 1)) //10
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
						T_Back = T_Last[nx * nz * j + nz * i + m - 1];
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i == 0 && m != 0 && m != (nz - 1)) //11
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
						T_Back = T_Last[nx * nz * j + nz * i + m - 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i == (nx - 1) && m != 0 && m != (nz - 1)) //12
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
						T_Back = T_Last[nx * nz * j + nz * i + m - 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m == 0)  //13
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
						T_Back = T_Last[nx * nz * j + nz * i + m + 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m == (nz - 1))  //14
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m - 1];
						T_Back = T_Last[nx * nz * j + nz * i + m - 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i == 0 && m == 0)  //15
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
						T_Back = T_Last[nx * nz * j + nz * i + m + 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i == 0 && m == (nz - 1))  //16
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m - 1];
						T_Back = T_Last[nx * nz * j + nz * i + m - 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i == (nx - 1) && m == 0)  //17
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
						T_Back = T_Last[nx * nz * j + nz * i + m + 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j == (ny - 1) && i == (nx - 1) && m == (nz - 1))  //18
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Right = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m - 1];
						T_Back = T_Last[nx * nz * j + nz * i + m - 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i != 0 && i != (nx - 1) && m == 0)  //19
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
						T_Back = T_Last[nx * nz * j + nz * i + m + 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i != 0 && i != (nx - 1) && m == (nz - 1))  //20
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m - 1] - 2 * dz * h * (T_Middle - Tw) / steel->lamda;
						T_Back = T_Last[nx * nz * j + nz * i + m - 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i == 0 && m == 0) //21
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
						T_Back = T_Last[nx * nz * j + nz * i + m + 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m == 0)  //22
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
						T_Back = T_Last[nx * nz * j + nz * i + m + 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i == 0 && m == (nz - 1)) //23
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m - 1];
						T_Back = T_Last[nx * nz * j + nz * i + m - 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m == (nz - 1)) //24
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m - 1];
						T_Back = T_Last[nx * nz * j + nz * i + m - 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i == 0 && m != 0 && m != (nz - 1))  //25
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
						T_Back = T_Last[nx * nz * j + nz * i + m - 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m != 0 && m != (nz - 1)) //26
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i - 1) + m] - 2 * dx * h * (T_Middle - Tw) / steel->lamda;
						T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
						T_Back = T_Last[nx * nz * j + nz * i + m - 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}

					else  //27
					{
						T_Middle = T_Last[nx * nz * j + nz * i + m];
						T_Up = T_Last[nx * nz * j + nz * (i + 1) + m];
						T_Down = T_Last[nx * nz * j + nz * (i - 1) + m];
						T_Right = T_Last[nx * nz * (j + 1) + nz * i + m];
						T_Left = T_Last[nx * nz * (j - 1) + nz * i + m];
						T_Forw = T_Last[nx * nz * j + nz * i + m + 1];
						T_Back = T_Last[nx * nz * j + nz * i + m - 1];
						T_New[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
							+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
					}
				}
		}
		for (int k = 0; k < ny; k++)
			//T_Surface[k] = 1558.0f;
			T_Surface[k] = T_New[nx * nz * k + nz * int((nx - 1)/2) + nz - 1];
	}

	else
		{
			for (int j = 0; j < ny; j++)
			{
				this->boundarycondition3d(*mCasterOne, hinit, j);
				for (int i = 0; i < nx; i++)
					for (int m = 0; m < nz; m++)
					{
						steel->physicalpara(T_Last[nx * nz * j + nz * i + m]);
						a = steel->lamda / (steel->pho * steel->ce);
						if (j == 0 && i != 0 && i != (nx - 1) && m != 0 && m != (nz - 1)) //1
						{
							T_Last[nx * nz * j + nz * i + m] = T_Cast;
						}

						else if (j == 0 && i == 0 && m != 0 && m != (nz - 1)) //2
						{
							T_Last[nx * nz * j + nz * i + m] = T_Cast;
						}

						else if (j == 0 && i == (nx - 1) && m != 0 && m != (nz - 1))//3
						{
							T_Last[nx * nz * j + nz * i + m] = T_Cast;
						}

						else if (j == 0 && i != 0 && i != (nx - 1) && m == 0) //4
						{
							T_Last[nx * nz * j + nz * i + m] = T_Cast;
						}

						else if (j == 0 && i != 0 && i != (nx - 1) && m == (nz - 1)) //5
						{
							T_Last[nx * nz * j + nz * i + m] = T_Cast;
						}

						else if (j == 0 && i == 0 && m == 0)  //6
						{
							T_Last[nx * nz * j + nz * i + m] = T_Cast;
						}

						else if (j == 0 && i == 0 && m == (nz - 1))  //7
						{
							T_Last[nx * nz * j + nz * i + m] = T_Cast;
						}

						else if (j == 0 && i == (nx - 1) && m == 0)  //8
						{
							T_Last[nx * nz * j + nz * i + m] = T_Cast;
						}

						else if (j == 0 && i == (nx - 1) && m == (nz - 1)) //9
						{
							T_Last[nx * nz * j + nz * i + m] = T_Cast;
						}

						else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m != 0 && m != (nz - 1)) //10
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == 0 && m != 0 && m != (nz - 1)) //11
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == (nx - 1) && m != 0 && m != (nz - 1)) //12
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m == 0)  //13
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m + 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i != 0 && i != (nx - 1) && m == (nz - 1))  //14
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m - 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == 0 && m == 0)  //15
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m + 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == 0 && m == (nz - 1))  //16
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m - 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == (nx - 1) && m == 0)  //17
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m + 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j == (ny - 1) && i == (nx - 1) && m == (nz - 1))  //18
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m - 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i != 0 && i != (nx - 1) && m == 0)  //19
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m + 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i != 0 && i != (nx - 1) && m == (nz - 1))  //20
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m - 1] - 2 * dz * h * (T_Middle - Tw) / steel->lamda;
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == 0 && m == 0) //21
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m + 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m == 0)  //22
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m + 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == 0 && m == (nz - 1)) //23
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m - 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m == (nz - 1)) //24
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m - 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == 0 && m != 0 && m != (nz - 1))  //25
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else if (j != 0 && j != (ny - 1) && i == (nx - 1) && m != 0 && m != (nz - 1)) //26
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i - 1) + m] - 2 * dx * h * (T_Middle - Tw) / steel->lamda;
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}

						else  //27
						{
							T_Middle = T_New[nx * nz * j + nz * i + m];
							T_Up = T_New[nx * nz * j + nz * (i + 1) + m];
							T_Down = T_New[nx * nz * j + nz * (i - 1) + m];
							T_Right = T_New[nx * nz * (j + 1) + nz * i + m];
							T_Left = T_New[nx * nz * (j - 1) + nz * i + m];
							T_Forw = T_New[nx * nz * j + nz * i + m + 1];
							T_Back = T_New[nx * nz * j + nz * i + m - 1];
							T_Last[nx * nz * j + nz * i + m] = (a*tao / (dx*dx))*T_Up + a*(tao / (dx*dx))*T_Down + ((1 - 2 * a*tao / (dx*dx) - 2 * a*tao / (dy*dy) - 2 * a*tao / (dz*dz) + tao*vcast / dy))*T_Middle
								+ (a*tao / (dy*dy))*T_Right + (a*tao / (dy*dy) - tao*vcast / dy)*T_Left + (a*tao / (dz*dz))*T_Forw + (a*tao / (dz*dz))*T_Back;
						}
					}
			}
			for (int k = 0; k < ny; k++)
				//T_Surface[k] = 1558.0f;
				T_Surface[k] = T_Last[nx * nz * k + nz * int((nx - 1) / 2) + nz - 1];
		}
	disout = !disout;
	tstep++;
}

void Temperature::boundarycondition3d(ContinuousCaster & CasterOne, float *hinit, int j)
{
	float yposition = dy * j;
	for (int i = 0; i < CasterOne.section; i++)
		if (yposition >= *(CasterOne.ccml + i) && yposition <= *(CasterOne.ccml + i + 1))
			h = *(hinit + i);
}

void Temperature::initcondition3d(float T_Cast)
{
	tstep = 0;
	for (int j = 0; j < ny; j++)
	  for (int i = 0; i < nx; i++)
		for (int k = 0; k < nz; k++)
		{
			 T_Last[nx * nz * j + nz * i + k] = T_Cast;
			 T_New[nx * nz * j + nz * i + k] = T_Cast;
		}
	disout = true;
}

void Temperature::initcondition3d(float *T_Cast)
{
	tstep = 0;
	for (int j = 0; j < ny; j++)
		for (int i = 0; i < nx; i++)
			for (int k = 0; k < nz; k++)
			{
				T_Last[nx * nz * j + nz * i + k] = T_Cast[nx * nz * j + nz * i + k];
				T_New[nx * nz * j + nz * i + k] = T_Cast[nx * nz * j + nz * i + k];
			}
	disout = true;
}


void Temperature::computetemperature3d(float *measuredpoistion, int measurednumb)
{
	computetemperature = new float[measurednumb];
	for (int i = 0; i < measurednumb; i++)
		for (int j = 0; j < ny; j++)
			if (fabs(j * dy - *(measuredpoistion + i)) <= dy)
				computetemperature[i] = T_Surface[j];
}

void Temperature::computemeantemperature3d()
{
	float y;
	int count = 0;
	for (int i = 0; i < mCasterOne->coolsection; i++)
	{
		meantemperature[i] = 0.0;
		for (int j = 0; j < ny; j++)
		{
			y = j * dy;
			if (y > *((mCasterOne->ccml) + i + mCasterOne->moldsection) && y <= *((mCasterOne->ccml) + i + 1 + +mCasterOne->moldsection))
			{
				meantemperature[i] += T_Surface[j];
				count++;
			}
		}
		meantemperature[i] = meantemperature[i] / count;
		count = 0;
	}

}

void Temperature::print3d(int measurednumb)
{
	if (tstep == 0)
	{
		cout << "lx = " << lx << ", " << "nx = " << nx << ", ";
		cout << "ly = " << ly << ", " << "ny = " << ny << ", ";
		cout << "lz = " << lz << ", " << "nz = " << nz << ", ";
		cout << "casting speed = " << vcast << ", " << endl;
		cout << "dx = " << dx << ", ";
		cout << "dy = " << dy << ", ";
		cout << "dz = " << dz << ", ";
		cout << "time step = " << tao << ", " << endl;
	}
	else
	{
		cout << "computetemperature = " << endl;
		for (int i = 0; i < measurednumb; i++)
			cout << computetemperature[i] << ", ";
		cout << endl;
	}
}

void Temperature::print3d()
{
	if (tstep == 0)
	{
		cout << "lx = " << lx << ", " << "nx = " << nx << ", ";
		cout << "ly = " << ly << ", " << "ny = " << ny << ", ";
		cout << "lz = " << lz << ", " << "nz = " << nz << ", ";
		cout << "casting speed = " << vcast << ", " << endl;
		cout << "dx = " << dx << ", ";
		cout << "dy = " << dy << ", ";
		cout << "dz = " << dz << ", ";
		cout << "time step = " << tao << ", " << endl;
	}
	else
	{
		cout << "tstep = "<< tstep << endl;
		cout << "meantemperature = " << endl;
		for (int i = 0; i < mCasterOne->coolsection; i++)
			cout << meantemperature[i] << ", ";
		cout << endl;
	}
}

void Temperature::operator=(const Temperature & m_SteelTemperature)
{
	mCasterOne = m_SteelTemperature.mCasterOne;
	steel = m_SteelTemperature.steel;
	nx = m_SteelTemperature.nx;
	ny = m_SteelTemperature.ny;
	nz = m_SteelTemperature.nz;
	tnpts = m_SteelTemperature.tnpts;
	tf = m_SteelTemperature.tf;
	lx = m_SteelTemperature.lx;
	ly = m_SteelTemperature.ly;
	lz = m_SteelTemperature.lz;
	dx = m_SteelTemperature.lx / float(m_SteelTemperature.nx - 1);
	dy = m_SteelTemperature.ly / float(m_SteelTemperature.ny - 1);
	dz = m_SteelTemperature.lz / float(m_SteelTemperature.nz - 1);
	tao = m_SteelTemperature.tf / float(m_SteelTemperature.tnpts - 1);
	delete[] T_New;
	delete[] T_Last;
	delete[] T_Surface;
	delete[] meantemperature;

	T_New = new float[nx * ny * nz];
	T_Last = new float[nx * ny * nz];
	for (int j = 0; j < ny; j++)
		for (int i = 0; i < nx; i++)
			for (int k = 0; k < nz; k++)
			{
				T_Last[nx * nz * j + nz * i + k] = m_SteelTemperature.T_Last[nx * nz * j + nz * i + k];
				T_New[nx * nz * j + nz * i + k] = m_SteelTemperature.T_New[nx * nz * j + nz * i + k];
			}
	T_Surface = new float[ny];
	for (int j = 0; j < ny; j++)
		T_Surface[j] = m_SteelTemperature.T_Surface[j];

	meantemperature = new float[mCasterOne->coolsection];
	for (int i = 0; i < mCasterOne->coolsection; i++)
		meantemperature[i] = m_SteelTemperature.meantemperature[i];

	vcast = m_SteelTemperature.vcast;
	tstep = m_SteelTemperature.tstep;
	disout = m_SteelTemperature.disout;
}

void Steel::physicalpara(float T)
{
	float Ts = 1462.0, Tl = 1518.0, lamds = 30, lamdl = 50, phos = 7000, phol = 7500, cel = 540.0, L = 265600.0, fs = 0.0;
	if (T < Ts)
	{
		fs = 0;
		pho = phos;
		lamda = lamds;
		ce = cel;
	}

	if (T >= Ts && T <= Tl)
	{
		fs = (T - Ts) / (Tl - Ts);
		pho = fs * phos + (1 - fs) * phol;
		lamda = fs*lamds + (1 - fs) * lamdl;
		ce = cel + L / (Tl - Ts);
	}

	if (T > Tl)
	{
		fs = 1;
		pho = phol;
		lamda = lamdl;
		ce = cel;
	}
}

class Optimizationalgorithm
{
    private:
	    int coolsection;
    public:
		float** allmeantemperature;
		float* staticmeantemperature;
		float* taimmeantemperature;
		float** Jacobian;
		float* gradient;
		float dh;
		float step;
		float costvalue;
		ContinuousCaster* mCasterOne;
		Optimizationalgorithm(ContinuousCaster &, float*);
		void gradientcalculation();
		void init(float**, float*, float);
		void linesearch();
		void updateh(float*);
		void print();
		void outputdata(int, float*);
};

Optimizationalgorithm::Optimizationalgorithm(ContinuousCaster & CasterOne, float* m_taimmeantemperature)
{
	mCasterOne = &CasterOne;
	coolsection = mCasterOne->coolsection;
	allmeantemperature = new float*[mCasterOne->coolsection];
	for (int i = 0; i < mCasterOne->coolsection; i++)
		allmeantemperature[i] = new float[coolsection];
	Jacobian = new float*[mCasterOne->coolsection];
	for (int i = 0; i < mCasterOne->coolsection; i++)
		Jacobian[i] = new float[coolsection];
	staticmeantemperature = new float[coolsection];
	gradient = new float[coolsection];
	taimmeantemperature = new float[coolsection];
	for (int i = 0; i < mCasterOne->coolsection; i++)
		taimmeantemperature[i] = m_taimmeantemperature[i];
}

void Optimizationalgorithm:: gradientcalculation()
{
	for (int i = 0; i < coolsection; i++)
	   for (int j = 0; j < coolsection; j++)
		   Jacobian[i][j] = (allmeantemperature[i][j] - staticmeantemperature[j]) / dh;
	for (int i = 0; i < coolsection; i++)
	{
		gradient[i] = 0.0;
		for (int j = 0; j < coolsection; j++)
			gradient[i] = gradient[i] + (taimmeantemperature[j] - staticmeantemperature[j]) * Jacobian[i][j];
	}
}

void::Optimizationalgorithm::linesearch()
{
	float step1 = 0.0, step2 = 0.0, eps = 1.0;
	for (int i = 0 ; i < coolsection; i++)
		for (int j = 0; j < coolsection; j++)
		{
			step1 += (staticmeantemperature[i] - taimmeantemperature[i])*Jacobian[i][j];
			step2 += Jacobian[i][j] * gradient[j] * Jacobian[i][j] * gradient[j];
		}
	step = fabs(step1 / (step2 + eps));
}

void::Optimizationalgorithm::updateh(float*hresult)
{
	//cout << " h = ";
	for (int i = 0; i < coolsection; i++)
	{
		hresult[i + mCasterOne->moldsection] = hresult[i + mCasterOne->moldsection] + step * gradient[i];
		//cout  << hresult[i + mCasterOne->moldsection] << " ";
	}
		
	costvalue = 0.0;
	for (int i = 0; i < coolsection; i++)
		costvalue += (staticmeantemperature[i] - taimmeantemperature[i]) * (staticmeantemperature[i] - taimmeantemperature[i]);
}

void Optimizationalgorithm::init(float**m_allmeantemperature, float *m_staticmeantemperature, float m_dh)
{
	dh = m_dh;
	for (int i = 0; i < coolsection; i++)
		for (int j = 0; j < coolsection; j++)
			allmeantemperature[i][j] = m_allmeantemperature[i][j];
	for (int i = 0; i < coolsection; i++)
		staticmeantemperature[i] = m_staticmeantemperature[i];
}

void Optimizationalgorithm::print()
{
	cout << endl;
	cout << "Jacobian = " << endl;
	for (int i = 0; i < coolsection; i++)
	{
		for (int j = 0; j < coolsection; j++)
			cout << Jacobian[i][j] << ",";
		cout << endl;
	}

	cout << "staticmeantemperature = " << endl;
	for (int i = 0; i < coolsection; i++)
		cout << staticmeantemperature[i] << ", ";
	cout << endl;

	/*cout << "allmeantemperature = " << endl;
	for (int i = 0; i < coolsection; i++)
	{
		for (int j = 0; j < coolsection; j++)
			cout << allmeantemperature[i][j] << ", ";
		cout << endl;
	}*/

	cout << "Gradient = ";
	for (int i = 0; i < coolsection; i++)
			cout << gradient[i] << ", ";
	cout << endl;
	cout << "step = " << step << endl;
	cout << "costvalue = " << costvalue << endl;
}

void Optimizationalgorithm::outputdata(int m_tstep, float*hinit)
{
	ofstream outputfile;
	if(m_tstep % 10 == 0)
	{
		outputfile.open("C:\\SteelTemperatureData.txt", ios::app);
		outputfile << m_tstep << endl;
		outputfile << step << ", " << costvalue << endl;
		for (int i = 0; i < coolsection; i++)
			outputfile << staticmeantemperature[i] << ", ";
		outputfile << endl;
		for (int i = 0; i < coolsection; i++)
			outputfile << (staticmeantemperature[i] - taimmeantemperature[i]) << ", ";
		outputfile << endl;
		for (int i = 0; i < coolsection; i++)
			outputfile << gradient[i] << ", ";
		outputfile << endl;
		for (int i = 0; i < coolsection; i++)
			outputfile << hinit[i + mCasterOne->moldsection] << ", ";
		outputfile << endl;
		outputfile << endl;
		outputfile.close();
	}
	
}

class PSOalgorithm
{
    public:
	int popsize;
	int measurednumb;
	int coolsection, moldsection, section;
	int labelgbest;
	float *fitnessvalue;
	float **poph;
	float **popv;
	float *pbest;
	float *gbest;
	float c1, c2, omga, vmax;
	float gbestvalue;
	PSOalgorithm::PSOalgorithm(int, int, float, float*, ContinuousCaster &);
	void fitnessevaluation(float*, float*, Temperature &, Temperature &);
	void findgbest();
	void updatepoph(float*);
	void initpoph(float*, float);
};

PSOalgorithm::PSOalgorithm(int m_measurednumb, int m_popsize, float rangeh, float *hinit, ContinuousCaster & CasterOne)
{
	coolsection = CasterOne.coolsection;
	moldsection = CasterOne.moldsection;
	section = CasterOne.section;
	measurednumb = m_measurednumb;
	popsize = m_popsize;
	c1 = 0.1f;
	c2 = 0.1f;
	vmax = 5.0;
	fitnessvalue = new float[popsize];
	popv = new float*[popsize];
	gbest = new float[section];
	for (int i = 0; i < popsize; i++)
		popv[i] = new float[section];
	poph = new float*[popsize];
	for (int i = 0; i < popsize; i++)
		poph[i] = new float[section];
	for (int i = 0; i < popsize; i++)
		for (int j = 0; j < section; j++)
		{
			if (j < moldsection)
				poph[i][j] = hinit[j];
			else
				poph[i][j] = hinit[j] + rangeh * rand() / float(RAND_MAX);
		}
	for (int i = 0; i < popsize; i++)
		for (int j = 0; j < section; j++)
			popv[i][j] = 0.0;
}

void PSOalgorithm::fitnessevaluation(float* measuredtemperature, float* measuredpoistion, Temperature & Temperature3dmodellast, Temperature & Temperature3dtemp)
{
	for (int i = 0; i < popsize; i++)
	{
		Temperature3dtemp = Temperature3dmodellast;
		Temperature3dtemp.differencecalculation3d(poph[i]);
		Temperature3dtemp.computetemperature3d(measuredpoistion, measurednumb);
		fitnessvalue[i] = 0.0;
		for (int j = 0; j < measurednumb; j++)
			fitnessvalue[i] +=(Temperature3dtemp.computetemperature[j] - measuredtemperature[j]) * (Temperature3dtemp.computetemperature[j] - measuredtemperature[j]);
	}
}

void PSOalgorithm::findgbest()
{
	gbestvalue = 10e10;
	for (int i = 0; i < popsize; i++)
	{
		if (fitnessvalue[i] < gbestvalue)
		{
			gbestvalue = fitnessvalue[i];
			labelgbest = i;
		}
	}
	for (int i = 0; i < section; i++)
		gbest[i] = poph[labelgbest][i];
}

void PSOalgorithm::initpoph(float*hinit, float rangeh)
{
	for (int i = 0; i < popsize; i++)
		for (int j = 0; j < section; j++)
		{
			if (j < moldsection)
				poph[i][j] = hinit[i];
			else
				poph[i][j] = hinit[i] + rangeh * rand() / float(RAND_MAX);
		}
}

void PSOalgorithm::updatepoph(float*hresult)
{
	for (int i = 0; i < popsize; i++)
		for (int j = moldsection; j < section; j++)
		{
			popv[i][j] = omga * popv[i][j] + c1* rand() * (gbest[j] - poph[i][j]);
			if(fabs(popv[i][j]) < vmax)
			    poph[i][j] = poph[i][j] + popv[i][j];
			else
				poph[i][j] = poph[i][j] + popv[i][j] * vmax / fabs(popv[i][j]);
		}
	for (int j = moldsection; j < section; j++)
		hresult[j] = pbest[j];
}

int main()
{
	const int section = 12, coolsection = 8, moldsection = 4, measurednumb = 8, predictstep = 5;
	float ccml[section + 1] = { 0.0f,0.2f,0.4f,0.6f,0.8f,1.0925f,2.27f,4.29f,5.831f,9.6065f,13.6090f,19.87014f,28.599f };
	float measuredpoistion[measurednumb] = { 0.9463f, 1.6812f, 3.28f, 5.0605f, 7.7188f, 11.6077f, 16.7395f, 24.235f };
	float *measuredtemperature = new float[measurednumb];
	float hinit[section] = { 1380.0f,1170.0f,980.0f,800.0f,1223.16f,735.05f,424.32f,392.83f,328.94f,281.64f,246.16f,160.96f };
	float taimmeantemperature[coolsection] = { 966.149841f, 925.864746f, 952.322083f, 932.175537f, 914.607117f, 890.494263f, 870.804443f, 890.595825f };
	float *htemp = new float[section];
	float *hresult = new float[section];
	float m_lx = 0.25f, m_ly = 0.25f, m_lz = 0.25f, m_tf = 1500.0f, m_vcast = -0.02f, rangeh = 50.0f, T_Cast = 1558.0f, dh = 1.0f;
	int m_nx = 21, m_ny = 21, m_nz = 21, m_tnpts = 3001, sim_tnpts;
	int popsize = 10;
	float **allmeantemperature;
	float *staticmeantemperature;
	ContinuousCaster CasterOne = ContinuousCaster(section, coolsection, moldsection, ccml);
	Steel steel;
	Optimizationalgorithm Gradientmethod = Optimizationalgorithm(CasterOne, taimmeantemperature);

	Temperature2d SteelTemperature2d = Temperature2d(m_nx, m_ny, m_tnpts, m_tf, m_lx, m_ly, m_vcast, CasterOne, steel);
	SteelTemperature2d.initcondition2d();

	while (SteelTemperature2d.tstep < m_tnpts)
	{
		SteelTemperature2d.differencecalculation2d(hinit);
		SteelTemperature2d.computetemperature2d(measuredpoistion, measurednumb);
	}

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


