#include "Course2Header.h"

//typedef double Type;

Type Alpha(Type*& P) {
	Type r = sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2]);
	//return r*r;
	return   1;
}
Type Lambda(Type*& P) {
	return 1;
	
	Type r = sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2]);
	return r;
}
Type Q(Type*& P) {
	
	if (Cond2) {
		Type x = P[0];
		return (-3 * pow(x, 3. / 2) + 4 * sqrt(x) + sqrt(3)) / (4 * sqrt(x)) * exp(-pow(x, 3. / 2) / sqrt(3));

		return  2 * cos(P[0]);
	}
	else {
		Type r = sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2]);

		//return 1 * sqrt(3) / 2 * ((3 * 3. / 2) * sqrt(r) + (2. / sqrt(3) - sqrt(3) / 2) * r * r) * exp(-1. / sqrt(3) * pow(r, 3. / 2));
		
		Type omega = sqrt(3) / 2 * sqrt(Alpha(P) / Lambda(P));
	    return  1*(-Lambda(P)*(omega * omega - 2.*omega / r) + Alpha(P)) * exp(-omega * r);
	}

}
Type Sol(Type*& P) {

	if (Cond2)
		return exp(-pow(P[0], 3. / 2) / sqrt(3));// cos(P[0]);
	else {
		Type r = sqrt(P[0] * P[0] + P[1] * P[1] + P[2] * P[2]);
		//return 1 * exp(-pow(r, 3. / 2) / sqrt(3));

		return 1 * exp(-r * sqrt(3) / 2 * sqrt(Alpha(P)));
	}
}


Type Grad(Type*& P, Type*& n) {
	Type x = P[0];
	return (-sqrt(3) / 2 * sqrt(x) * exp(-pow(x, 3. / 2) / sqrt(3)))* n[0];
	//if(abs(n[0])>0.1)cout<<  P[0]<<"  " << n[0] << "\n";
	 return -sin(P[0]) * n[0];
}

int main(int argc, char* argv[]) {


	{
		/*std::ifstream ifile;
		std::ofstream ofile2;

		ifile.open("C:\\Users\\Artem\\Desktop\\TetraGrid\\Out55.txt");
		ofile2.open("C:\\Users\\Artem\\Desktop\\TetraGrid\\Out555.txt");

		Type r, U;

		while (ifile >> r >> U)
			ofile2 << r <<' '<< U + 0.005 << '\n';


		ifile.close();
		ofile2.close();
		return 0;*/
	}
	//SetFromMeshToVTK("C:\\Users\\Artem\\Desktop\\TetraGrid\\meshCube\\", "C:\\Users\\Artem\\Desktop\\TetraGrid\\meshCube\\Sphere0", 1,0, 0.1);
	//return 0;

	
	const char* SettingFile = "D:\\Desktop\\TetraGrid\\CourseSetFile.txt";
	//файл с сеткой
	std::string NameFileVTK;
	std::string NameFileDataGrid;
	std::string NameFileOut;
	//чтение начальных данных
	ReadSettings(SettingFile, NameFileVTK, NameFileOut, NameFileDataGrid);
	
	//сетка
	vtkSmartPointer<vtkUnstructuredGrid> unstructuredgrid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	//FullConvergence(NameFileVTK, NameFileDataGrid, NameFileOut, 3, unstructuredgrid, Sol);
	//return 0;
	
	clock_t start_clock = clock();
	if (ReadFileVTK(NameFileVTK, unstructuredgrid, true)) {
		std::cout << " Reading error!\n";
		return 1;
	}
	clock_t end_clock = clock();
	std::cout << "\n Read file VTK: " << (Type)(end_clock - start_clock) / CLOCKS_PER_SEC << "\n\n";

	Type foo;
	WriteFileNeighbors(NameFileDataGrid, NameFileDataGrid+"dV.txt",foo, unstructuredgrid);
	//return 0;


	//число ячеек в сетке
	size_t SizeGrid = unstructuredgrid->GetNumberOfCells();
	// матрица значений
	MyPair<Type>** MatrixEnergy;	
	//вектор правой части
	Type* VectorQ;
	//вектор решений
	std::vector<MySolution<Type>> VectorSol(SizeGrid);
	
	SetMemory(SizeGrid, VectorQ);
	SetMemory(SizeGrid, 5, MatrixEnergy); //5 --- 5 не нулевых элементов в строке
	const size_t N = SizeGrid * 5;
	
	
	Type t = Test_h_min * Test_h_min;// *My_h_min;  // 0.01 * Test_h_min; Sphere_h_min

	start_clock = clock();
	MakeMatrix(NameFileDataGrid, unstructuredgrid, VectorSol, VectorQ, MatrixEnergy, t);
	end_clock = clock();
	std::cout << "\n Make full Matrix: " << (Type)(end_clock - start_clock) / CLOCKS_PER_SEC << "\n\n";

	Type eps = 0.00001;
	size_t MaxIter = 100000;

	{
	/*	for (size_t i = 0; i < SizeGrid; i++)
		{
			for (size_t j = 0; j < 5; j++)
			{
				cout << MatrixEnergy[i][j] << ' ';
				
			}
			cout <<"   Q_i:"<< VectorQ[i];
			cout << '\n';
		}*/
	}

	start_clock = clock();

	

	SolveSparseTimeInDependentSystemEigenBiCGSTAB(SizeGrid, MaxIter, VectorSol, VectorQ, MatrixEnergy);
	{
		/*std::ofstream ofile;
		ofile.open(NameFileOut + "Q.txt");
		for (size_t i = 0; i < SizeGrid; i++)
			ofile << VectorQ[i] << '\n';
		ofile.close();*/
	}
    //WriteFileSolution(NameFileOut + "5.txt", SizeGrid, VectorSol);
	WriteFileSolution(NameFileOut+"test5.vtk", SizeGrid, VectorSol, unstructuredgrid);

	end_clock = clock();
	std::cout << "\n Solving system: " << (Type)(end_clock - start_clock) / CLOCKS_PER_SEC << "\n\n";
	


	cout << "Normmax: " << NormMax(SizeGrid, VectorSol, Sol) << '\n';
	//cout << "Norm2: " << NormGridF(NameFileDataGrid + "dV.txt", SizeGrid, VectorSol, Sol) << '\n';
	
	
	ClearMemory(SizeGrid, MatrixEnergy);
	ClearMemory(VectorQ);
	VectorSol.clear();

	system("pause");
	return 0;
}