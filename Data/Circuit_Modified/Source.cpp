#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <iomanip>
#include <Eigen\Eigen>

using namespace std;
using namespace Eigen;

double Pi = 4 * atan(1);
double Pha(complex<double> &C);

struct Component
{
	string Type;
	complex <double> Impedance;
	complex <double> Admittance;
	complex <double> Current;
	complex <double> Voltage;
	int Nodes[2];
};

int main()
{

	vector <Component> VSources;
	vector <Component> CSources;
	vector <Component> Elements;
	vector <Component> ALL;

	
	ifstream File;
	string Name;
	
	Ahmed:
	cout << "Enter File Name " << endl;
	cin >> Name;
	File.open(Name+".txt");
	if (File.is_open())
	{
		string Type;
		int Nodes[2];
		double Value;   double Omega;
		double Mod;     double Phase;

		while (File >> Type)
		{
			if (Type == "w")
				File >> Omega;

			if (Type[0] == 'V') {
				File >> Nodes[0];    File >> Nodes[1];
				File >> Mod;          File >> Phase;
				complex <double> Temp(polar(Mod, (Phase*Pi / 180)));

				Component V;
				V.Voltage = Temp;
				V.Nodes[0] = Nodes[0];
				V.Nodes[1] = Nodes[1];
				V.Type = Type;

				ALL.push_back(V);
				VSources.push_back(V);
			}

			else if (Type[0] == 'I') {
				File >> Nodes[0];    File >> Nodes[1];
				File >> Mod;          File >> Phase;
				complex <double> Temp(polar(Mod, (Phase*Pi / 180)));

				Component I;
				I.Current = Temp;
				I.Nodes[0] = Nodes[0];
				I.Nodes[1] = Nodes[1];
				I.Type = Type;

				ALL.push_back(I);
				CSources.push_back(I);
			}

			else if (Type[0] == 'L') {
				File >> Nodes[0];    File >> Nodes[1];
				File >> Value;
				complex <double> Temp1(0, (Omega*Value));
				complex <double> Temp2(0, (-1/(Omega*Value)));

				Component L;
				L.Impedance = Temp1;
				L.Admittance = Temp2;
				L.Nodes[0] = Nodes[0];
				L.Nodes[1] = Nodes[1];
				L.Type = Type;

				ALL.push_back(L);
				Elements.push_back(L);
			}

			else if (Type[0] == 'C') {
				File >> Nodes[0];    File >> Nodes[1];
				File >> Value;
				complex <double> Temp1(0, (-1 / (Omega*Value)));
				complex <double> Temp2(0, (Omega*Value));

				Component C;
				C.Impedance = Temp1;
				C.Admittance = Temp2;
				C.Nodes[0] = Nodes[0];
				C.Nodes[1] = Nodes[1];
				C.Type = Type;

				ALL.push_back(C);
				Elements.push_back(C);
			}

			else if (Type[0] == 'R') {
				File >> Nodes[0];    File >> Nodes[1];
				File >> Value;
				complex <double> Temp1(Value, 0);
				complex <double> Temp2((1/Value), 0);

				Component R;
				R.Impedance = Temp1;
				R.Admittance = Temp2;
				R.Nodes[0] = Nodes[0];
				R.Nodes[1] = Nodes[1];
				R.Type = Type;

				ALL.push_back(R);
				Elements.push_back(R);
			}
		}

		//Number of Nudes
		int NodesNum = 0;
		for (int i = 0;i < ALL.size(); i++)
		{
			if (ALL[i].Nodes[0] > NodesNum)
				NodesNum = ALL[i].Nodes[0];
			if (ALL[i].Nodes[1] > NodesNum)
				NodesNum = ALL[i].Nodes[1];
		}
		

		//Nodes Voltages
		vector <complex <double>> VNodes;
		for (int i = 0;i <= NodesNum ; i++)
		{
			complex<double> Temp(0, 0);
			VNodes.push_back(Temp);
		}

		int VNUM = VSources.size();

		Eigen::MatrixXcd A(NodesNum + VNUM , NodesNum + VNUM);
		A = MatrixXd::Constant(NodesNum + VNUM, NodesNum + VNUM, 0);
		Eigen::MatrixXcd G(NodesNum , NodesNum);
		G = MatrixXd::Constant(NodesNum , NodesNum , 0);
		Eigen::MatrixXcd B(NodesNum, VNUM);
		B = MatrixXd::Constant(NodesNum ,VNUM, 0);
		Eigen::MatrixXcd C(VNUM , NodesNum);
		C = MatrixXd::Constant(VNUM, NodesNum, 0);
		Eigen::MatrixXcd D(VNUM , VNUM);
		D = MatrixXd::Constant(VNUM, VNUM, 0);
		Eigen::MatrixXcd X(NodesNum + VNUM, 1);
		X = MatrixXd::Constant(NodesNum + VNUM, 1, 0);
		Eigen::MatrixXcd Z(NodesNum + VNUM, 1);
		Z = MatrixXd::Constant(NodesNum + VNUM, 1, 0);
		Eigen::MatrixXcd I(NodesNum, 1);
		I = MatrixXd::Constant(NodesNum , 1, 0);
		Eigen::MatrixXcd E(VNUM, 1);
		E = MatrixXd::Constant(VNUM, 1, 0);

		//G
		for (int i = 0;i < Elements.size();i++)
		{
			if (Elements[i].Nodes[0] != 0 && Elements[i].Nodes[1] != 0)
			{
				G(Elements[i].Nodes[0] - 1, Elements[i].Nodes[0] - 1) += Elements[i].Admittance;
				G(Elements[i].Nodes[0] - 1, Elements[i].Nodes[1] - 1) -= Elements[i].Admittance;
				G(Elements[i].Nodes[1] - 1, Elements[i].Nodes[1] - 1) += Elements[i].Admittance;
				G(Elements[i].Nodes[1] - 1, Elements[i].Nodes[0] - 1) -= Elements[i].Admittance;

			}
			else if (Elements[i].Nodes[0] == 0)
				G(Elements[i].Nodes[1] - 1, Elements[i].Nodes[1] - 1) += Elements[i].Admittance;
			else if (Elements[i].Nodes[1] == 0)
				G(Elements[i].Nodes[0] - 1, Elements[i].Nodes[0] - 1) += Elements[i].Admittance;
		}

		//B,C,E
		for (int i = 0;i < VSources.size();i++)
		{
			E(i, 0) = VSources[i].Voltage;
			if (VSources[i].Nodes[0] != 0 && VSources[i].Nodes[1] != 0)
			{
				B(VSources[i].Nodes[0] - 1, i) += 1;
				C(i, VSources[i].Nodes[0] - 1) += 1;
				B(VSources[i].Nodes[1] - 1, i) -= 1;
				C(i, VSources[i].Nodes[1] - 1) -= 1;
			}
			else if (VSources[i].Nodes[0] == 0)
			{
				B(VSources[i].Nodes[1] - 1, i) -= 1;
				C(i, VSources[i].Nodes[1] - 1) -= 1;
			}
			else if (VSources[i].Nodes[1] == 0)
			{
				B(VSources[i].Nodes[0] - 1, i) += 1;
				C(i, VSources[i].Nodes[0] - 1) += 1;
			}	
		}

		//A
		for (int i = 0;i < NodesNum;i++)
		{
			for (int j = 0;j < NodesNum;j++)
			{
				A(i, j) += G(i, j);
			}
		}
		for (int i = NodesNum;i < NodesNum + VNUM ;i++)
		{
			for (int j = 0;j < NodesNum;j++)
			{
				A(i, j) += C(i- NodesNum, j);
			}
		}
		for (int i = 0;i < NodesNum ;i++)
		{
			for (int j = NodesNum;j < NodesNum + VNUM;j++)
			{
				A(i, j) += B(i, j-NodesNum);
			}
		}
		for (int i = NodesNum ;i < NodesNum + VNUM;i++)
		{
			for (int j = NodesNum + 1;j < NodesNum + VNUM;j++)
			{
				A(i, j) += D(i-NodesNum, j-NodesNum);
			}
		}
		
		//I
		for (int i = 0;i < CSources.size();i++)
		{
			if (CSources[i].Nodes[0] != 0 && CSources[i].Nodes[1] != 0)
			{
				I(CSources[i].Nodes[0] - 1, 0) += CSources[i].Current;
				I(CSources[i].Nodes[1] - 1, 0) -= CSources[i].Current;
			}
			else if (CSources[i].Nodes[0] == 0)
				I(CSources[i].Nodes[1] - 1, 0) -= CSources[i].Current;
			else if (CSources[i].Nodes[1] == 0)
				I(CSources[i].Nodes[0] - 1, 0) += CSources[i].Current;
				
		}

		//Z
		for (int i = 0;i < NodesNum;i++)
		{
			Z(i, 0) = I(i, 0);
		}
		for (int i = NodesNum;i < NodesNum+ VNUM;i++)
		{
			Z(i, 0) = E(i-NodesNum, 0);
		}

		X += A.colPivHouseholderQr().solve(Z);

		for (int i = NodesNum;i < VNUM + NodesNum;i++)
		{
			complex<double>Temp(-X(i, 0).real(), -X(i, 0).imag());
			X(i, 0) = Temp;
		}

		for (int i = 0;i < NodesNum;i++)
		{
			VNodes[i + 1] = X(i, 0);
		}

		//Approximate
		for (int i = 1;i <= NodesNum;i++)
		{
			if (abs(VNodes[i].real()) < 0.000001)
			{
				complex <double> Temp(0,VNodes[i].imag());
				VNodes[i] = Temp;
			}
			if(abs(VNodes[i].imag()) < 0.000001)
			{
				complex <double> Temp(VNodes[i].real() , 0);
				VNodes[i] = Temp;
			}
		}

		for (int i = NodesNum;i < VNUM+ NodesNum;i++)
		{
			if (abs(X(i, 0).real()) < 0.000001)
			{
				complex <double> Temp(0, X(i , 0).imag());
				X(i, 0) = Temp;
			}
			if (abs(X(i, 0).imag()) < 0.000001)
			{
				complex <double> Temp(X(i, 0).real(), 0);
				X(i, 0) = Temp;
			}
		}

		for (int i = 0;i < ALL.size();i++)
		{
			if(ALL[i].Type[0] != 'V' && ALL[i].Type[0] != 'I')
				ALL[i].Current = (VNodes[ALL[i].Nodes[0]] - VNodes[ALL[i].Nodes[1]]) / ALL[i].Impedance;
		}

		for (int i = 0;i < ALL.size();i++)
		{
			if (ALL[i].Current.imag() == -0)
			{
				complex <double> Temp(ALL[i].Current.real(), 0);
				ALL[i].Current = Temp;
			}
			if (ALL[i].Current.real() == -0)
			{
				complex <double> Temp(0 , ALL[i].Current.imag());
				ALL[i].Current = Temp;
			}
		}

		//V on The Screen
		cout << endl;
		for (int i = 1;i <= NodesNum;i++)
		{
			string INFO = "V" + to_string(i);
			cout << INFO << "\t" << abs(VNodes[i]) << "   " <<  Pha(VNodes[i]) <<endl;
		}

		cout << endl;
		
		//I on The Screen
		int k = 0;
		for(int i=0;i<ALL.size();i++)
		{
			if(ALL[i].Type[0] == 'V')
			{
				string INFO = "I"  + to_string(ALL[i].Nodes[1]) + to_string(ALL[i].Nodes[0])+ " "  + ALL[i].Type;;
				cout << INFO << "\t" << abs(X(NodesNum+k , 0)) <<"   " << Pha(X(NodesNum + k, 0)) <<endl;
				k++;
			}
			else if (ALL[i].Type[0] == 'I')
			{
				string INFO = "I" + to_string(ALL[i].Nodes[1]) + to_string(ALL[i].Nodes[0]) + " "  + ALL[i].Type;
				cout  << INFO  << "\t" << abs(ALL[i].Current) << "   " << Pha(ALL[i].Current) << endl;
			}
			else
			{
				string INFO = "I"  + to_string(ALL[i].Nodes[0]) + to_string(ALL[i].Nodes[1]) + " "  + ALL[i].Type;
				cout << INFO  << "\t" << abs(ALL[i].Current) << "   " << Pha(ALL[i].Current) << endl;
			}
		}
	}
	else 
	{
	 cout << "File Not Found" <<endl;
	system("pause");
	}
	system("pause");
}

double Pha(complex <double> &C)
{
	double phase = 0;
	if (C.imag() == 0 && C.real() > 0)
	{
		phase = 0;
		return phase;
	}
	if (C.imag() == 0 && C.real() < 0)
	{
		phase = 180;
		return phase;
	}
	if (C.real() == 0 && C.imag() > 0)
	{
		phase = 90;
		return phase;
	}
	if (C.real() == 0 && C.imag() < 0)
	{
		phase = -90;
		return phase;
	}

	phase = atan(C.imag() / C.real()) * 180 / Pi;
	if (C.imag() < 0 && C.real() < 0)
		phase -= 180;
	else if (C.imag() > 0 && C.real() < 0)
		phase += 180;
	return phase;
}

