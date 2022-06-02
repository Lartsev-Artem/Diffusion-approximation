#pragma once
typedef double Type;

/*Константы*/
#define Test_h_min 0.00135333 // минимальный  радиус вписанной сферы на тестовой сетке 
#define My_h_min 0.0413735 // минимальный  радиус вписанной сферы на VTK2 сетке 
#define Sphere_h_min 0.00654584 // минимальный  радиус вписанной сферы на Sphere.vtk сетке
#define PI 3.14159265358979323846
#define Velocity 29979245800.0
#define X_0 388165.0*pow(10,5)
#define U_0 pow(10,-10)
#define T_orb 4578.0 // 76.3*60 


/* константы условий*/
#define Cond2 0 // Гран. условия 2 рода true/false
#define Ellips 1  //Основная задача или нет
#define My_Zero pow(10,-50)

/*Стандартные заголовки*/
#include <iostream>
#include<string>
#include<algorithm>
#include<vector>
#include <cmath>


/*Eigen*/
#include<eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/IterativeLinearSolvers>

/*VTK*/
#include<vtk-9.0/vtkGenericDataObjectReader.h>
#include<vtk-9.0/vtkGenericDataObjectWriter.h>
#include<vtk-9.0/vtkUnstructuredGrid.h>
#include<vtk-9.0/vtkSmartPointer.h>
#include <vtk-9.0\vtkCellData.h>
#include <vtk-9.0\vtkDoubleArray.h>
#include <vtk-9.0\vtkXMLUnstructuredGridReader.h>


/*Свои заголовки*/
#include"MemoryPack.h"
#include "GMRES.h"
size_t MyMakeGrid3D();


Type Q(Type*& P);
Type Alpha(Type*& P);
Type Lambda(Type*& P);
Type Grad(Type*& P, Type*& n);

/*Типы данных*/

template <class Type>
struct MyPair{
	// пара номер столбца== id ячейки + значение
	int Id = -1;
	Type val = 0;

	friend std::ostream& operator<< (std::ostream& out, const MyPair& pair);
};
template <class Type>
struct MySolution{
	// точка в пространстве + значение в ней

	Type value;
	Type* Point;
	
	size_t SetData (Type*& P, Type val = 0) {
		value = val;
		Point = new Type[3];
		for (size_t i = 0; i < 3; i++)
			Point[i] = P[i];
		return 0;
	}

	MySolution() {
		value = 0;
		Point = NULL;
	}
	~MySolution() {
		delete[] Point;
	}

private:
	MySolution(const MySolution<Type>& copy);
};

template <class Type>
struct MyPoint {
	// точка в пространстве

	Type* Point;

	size_t SetData(Type*& P) {
		Point = new Type[3];
		for (size_t i = 0; i < 3; i++)
			Point[i] = P[i];
		return 0;
	}
	size_t SetData(Type a, Type b, Type c) {
		Point = new Type[3];
		Point[0] = a;
		Point[1] = b;
		Point[2] = c;
		return 0;
	}


	MyPoint() {
		Point = new Type[3]; //NULL;
	}
	~MyPoint() {
		delete[] Point;
	}


	MyPoint(const MyPoint<Type>& copy) {
		Point = new Type[3];
		for (size_t i = 0; i < 3; i++){
			Point[i] = copy.Point[i];
		}
	}
};

template <class Type>
struct MyCell {
	// точка в пространстве

	size_t* idCell;

	size_t SetData(size_t*& id) {
		idCell = new size_t[4];
		for (size_t i = 0; i < 4; i++)
			idCell[i] = id[i];
		return 0;
	}
	size_t SetData(MyCell& id) {
		idCell = new size_t[4];
		for (size_t i = 0; i < 4; i++)
			idCell[i] = id.idCell[i];
		return 0;
	}

	size_t SetData(size_t a, size_t b, size_t c,size_t d) {
		idCell = new size_t[4];
		idCell[0] = a;
		idCell[1] = b;
		idCell[2] = c;
		idCell[3] = d;
		return 0;
	}
	size_t ReSetData(size_t a, size_t b, size_t c, size_t d) {

		idCell[0] = a;
		idCell[1] = b;
		idCell[2] = c;
		idCell[3] = d;
		return 0;
	}


	MyCell() {
		idCell = new size_t[4];
	}
	~MyCell() {
		delete[] idCell;
	}


MyCell(const MyCell<Type>& copy) {
	idCell = new size_t[4];
	for (size_t i = 0; i < 4; i++){
		idCell[i] = copy.idCell[i];
	}
	}
};


template <class Type>
struct MyFace {
	// точка в пространстве

	size_t* idCell;
	size_t idTetra;

	size_t SetData(size_t*& id) {
		idCell = new size_t[3];
		for (size_t i = 0; i < 3; i++)
			idCell[i] = id[i];
		return 0;
	}
	size_t SetData(MyFace& id) {
		idCell = new size_t[3];
		for (size_t i = 0; i < 3; i++)
			idCell[i] = id.idCell[i];
		return 0;
	}

	size_t SetData(size_t a, size_t b, size_t c) {
		idCell = new size_t[3];
		idCell[0] = a;
		idCell[1] = b;
		idCell[2] = c;

		return 0;
	}
	size_t SetData(size_t a, size_t b, size_t c, size_t Id) {
		idCell = new size_t[3];
		idCell[0] = a;
		idCell[1] = b;
		idCell[2] = c;
		idTetra = Id;
		return 0;
	}

	size_t ReSetData(size_t a, size_t b, size_t c) {

		idCell[0] = a;
		idCell[1] = b;
		idCell[2] = c;
		return 0;
	}
	size_t ReSetData(size_t a, size_t b, size_t c, size_t Id) {

		idCell[0] = a;
		idCell[1] = b;
		idCell[2] = c;
		idTetra = Id;
		return 0;
	}



	MyFace() {
		idCell = new size_t[3];
		idTetra = -1;
	}
	~MyFace() {
		delete[] idCell;
	}


	MyFace(const MyFace<Type>& copy) {
		idCell = new size_t[3];
		for (size_t i = 0; i < 3; i++) {
			idCell[i] = copy.idCell[i];
		}
		idTetra = copy.idTetra;
	}
};

/*Ввод/Вывод данных*/

size_t WriteFileNeighbors(const std::string NameFileOut, const std::string NameFile_dv, Type& averageH, vtkSmartPointer<vtkUnstructuredGrid>& UGrid);
size_t ReadFileNeighbors(std::ifstream& ifile, size_t& curCell, int*& IdNeighbors, Type& VolumeCell,
	Type*& CenterCell, Type*& SquareFaces, std::vector<Type*>& NormalsCurCell, std::vector<Type*>& CenterCellNeighbors);

size_t ReadSettings(const char* NameFileSet, std::string& NameFileGrid, std::string& NameFileOut, std::string& NameFileDataGrid);
size_t ReadFileVTK(const std::string NameFileVTK, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, const bool print);
size_t WriteFileVTK(const std::string NameFileOut, vtkSmartPointer<vtkUnstructuredGrid>& UGrid);
size_t WriteFileSolution(const std::string NameFileOut, const size_t n, std::vector<MySolution<Type>>& VectorSol, size_t NumCoord=1);
size_t WriteFileSolution(const std::string NameFileOut, const size_t n, std::vector<MySolution<Type>>& VectorSol,
	vtkSmartPointer<vtkUnstructuredGrid>& UGrid);

/*Операции*/
Type NormMax(const size_t n, std::vector<MySolution<Type>>& x, Type(*f)(Type*&));
Type NormMax(const size_t n, Type*& x1, Type*& x2);
Type Norm2(const size_t n, Type*& x1, Type*& x2);
Type Norm2(const size_t n, std::vector<MySolution<Type>>& x1, Type*& x2);
Type Norm(Type*& v1);
size_t Normalize(Type*& v1);
Type Distance(const Type* h1, const Type* h2);

/*Преобразование данных*/

size_t FromMyToEigen(const size_t n, MyPair<Type>**& MatrixEnergy, Eigen::MatrixXd& matrix, const size_t numNonZeroСolumns = 5);
size_t FromMyToEigen(const size_t n, MyPair<Type>**& MatrixEnergy, Eigen::MatrixXd& matrix, Type*& MyVector, Eigen::VectorXd& EigenVector,
	const size_t numNonZeroСolumns = 5);
size_t FromEigenToMySol(const size_t n, Eigen::VectorXd& EigenVector, std::vector<MySolution<Type>>& VectorSol);

/*Вычисления компонент разностной схемы*/

Type SquareFace(size_t NumberCell, size_t NumberFace, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid);
Type GetVolumeCell(size_t NumberCell, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid);
size_t CenterOfTetra(size_t NumberCell, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Type*& PointInTetra);
size_t NormalToFace(size_t NumberCell, size_t NumberFace, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Type*& n);

size_t NormalAndSquareFace(size_t NumberCell, size_t NumberFace, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Type*& n, Type& S);
Type ProjectionOnNormal(Type*& n, Type*& Point1, Type*& Point2);


/*Создание матрицы*/
size_t MakeOneRow(const size_t NumCell, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid,
	MyPair<Type>**& MatrixEnergy, std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ);

size_t MakeOneRow(const size_t NumCell, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid,
	std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ, std::vector <Eigen::Triplet <Type>>& tripletList, MyPair<Type>**& MatrixEnergy);

size_t MakeMatrix(const std::string NameFileDataGrid, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid,
	std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ, MyPair<Type>**& MatrixEnergy, const Type t);


/*Решение системы*/

size_t SolutionSystemSeidel(const Type eps, const size_t MaxIter, const size_t SizeGrid,
	MyPair<Type>**& MatrixEnergy, std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ);

size_t SolveFullSystemEigenQR(const size_t SizeGrid, std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ, MyPair<Type>**& MatrixEnergy);

size_t SolveSparseTimeDependentInDirectSystemGMRES(const Type eps, const size_t SizeGrid, const Type t, const size_t MaxIter,
	std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ, MyPair<Type>**& MatrixEnergy, const size_t numNonZeroСolumns = 5);

size_t SolveSparseTimeDependentInDirectSystemEigenBiCGSTAB(const Type eps, const size_t SizeGrid, const Type t, const size_t MaxIter,
	std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ, MyPair<Type>**& MatrixEnergy,
	const size_t numNonZeroСolumns = 5);

size_t SolveSparseTimeInDependentSystemEigenBiCGSTAB(const size_t SizeGrid, const size_t MaxIter,
	std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ, MyPair<Type>**& MatrixEnergy,
	const size_t numNonZeroСolumns = 5);

size_t SolveFullTimeDependentSystemEigenQR(const Type eps, const size_t SizeGrid, const size_t MaxIter, const Type t,
	std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ, MyPair<Type>**& MatrixEnergy,
	const size_t numNonZeroСolumns = 5);


/*Функции*/

template<typename FileName, typename Type>
size_t FromMeshToVTK(FileName NameMesh, FileName NameVTK, Type foo = 0.1) {

	std::ifstream ifile;
	std::ofstream ofile;
	std::string str;
	ifile.open(NameMesh);
	ofile.open(NameVTK);

	ofile << "# vtk DataFile Version 3.0\n";
	ofile << "Elemensts Volumes\n";
	ofile << "ASCII\n";
	ofile << "DATASET UNSTRUCTURED_GRID\n";

	size_t NumPoint = 0;
	ifile >> NumPoint;
	std::cout << "POINTS " << NumPoint << '\n';
	ofile << "POINTS " << NumPoint << " double\n";

	Type x, y, z;

	for (size_t i = 0; i < NumPoint; i++) {
		ifile >> x >> y >> z;
		ofile << x << ' ' << y << ' ' << z << '\n';
	}

	size_t NumCells = 0;
	ifile >> NumCells;
	std::cout << "CELLS " << NumCells << '\n';
	ofile << "CELLS " << NumCells << ' ' << NumCells * 5 << '\n';

	Type a, b, c, d, drop;
	for (size_t i = 0; i < NumCells; i++) {
		ifile >> drop >> a >> b >> c >> d;
		ofile << 4 << ' ' << a - 1 << ' ' << b - 1 << ' ' << c - 1 << ' ' << d - 1 << '\n';
	}

	ofile << "CELL_TYPES " << NumCells << '\n';
	for (size_t i = 0; i < NumCells; i++) {
		ofile << 10 << '\n';
	}


	ifile.close();
	ofile.close();

	return 0;
}

template<typename FileName, typename Type>
size_t SetFromMeshToVTK(FileName MeshFile, FileName VTKFile, const size_t NumFile, const size_t Numstart = 0, Type foo = 0.1) {

	std::string NameMesh = MeshFile;
	std::string NameVTK = VTKFile;

	std::ifstream ifile;
	std::ofstream ofile;
	std::string str;

	for (size_t i = Numstart; i < NumFile; i++) {
		
		ifile.open(NameMesh + std::to_string(i) + ".txt");
		ofile.open(NameVTK + std::to_string(i) + ".vtk");



		ofile << "# vtk DataFile Version 3.0\n";
		ofile << "Elemensts Volumes\n";
		ofile << "ASCII\n";
		ofile << "DATASET UNSTRUCTURED_GRID\n";

		size_t NumPoint = 0;
		ifile >> NumPoint;
		std::cout << "POINTS " << NumPoint << '\n';
		ofile << "POINTS " << NumPoint << " double\n";

		Type x, y, z;

		for (size_t i = 0; i < NumPoint; i++) {
			ifile >> x >> y >> z;
			ofile << x << ' ' << y << ' ' << z << '\n';
		}

		size_t NumCells = 0;
		ifile >> NumCells;
		std::cout << "CELLS " << NumCells << '\n';
		ofile << "CELLS " << NumCells << ' ' << NumCells * 5 << '\n';

		size_t a, b, c, d, drop;
		for (size_t i = 0; i < NumCells; i++) {
			ifile >> drop >> a >> b >> c >> d;
			ofile << 4 << ' ' << a - 1 << ' ' << b - 1 << ' ' << c - 1 << ' ' << d - 1 << '\n';
		}

		ofile << "CELL_TYPES " << NumCells << '\n';
		for (size_t i = 0; i < NumCells; i++) {
			ofile << 10 << '\n';
		}


		ifile.close();
		ofile.close();
	}

	return 0;
}

template<typename FileName, typename Type>
Type NormMax(FileName NameSol, Type(*f)(Type)) {
	
	std::ifstream ifile;
	std::string str;
	ifile.open(NameSol);
	Type max = -1.;
	Type buf;
	Type r, val;
	while (ifile >> r >> val) {
		buf = abs(f(r) - val);
		if (buf > max);
		max = buf;
	}
	return max;

	ifile.close();

	return 0;
}

template<typename FileName, typename Type>
Type NormGridF(FileName Data, const size_t n, std::vector<MySolution<Type>>& x, Type(*f)(Type*&)) {
	Type sum = 0;
	std::ifstream ifile;
	ifile.open(Data);
	Type dV;
	//"C:\\Users\\Artem\\Desktop\\TetraGrid\\mesh\\Radius.txt"
	for (size_t i = 0; i < n; i++) {
		ifile >> dV;
		sum += pow(f(x[i].Point) - x[i].value, 2) * dV;
	}

	ifile.close();
	return sqrt(sum);
}

template<typename FileName, typename Type>
Type OutGrid(FileName NameSol, FileName NameSol2, std::vector<MyPoint<Type>>& Points, std::vector<MyCell<Type>>& Cells,
	std::vector<MyFace<Type>>& Faces) {

	std::ofstream ofile;
	std::string str;
	ofile.open(NameSol);

	ofile << "# vtk DataFile Version 3.0\n";
	ofile << "Elemensts Volumes\n";
	ofile << "ASCII\n";
	ofile << "DATASET UNSTRUCTURED_GRID\n";

	size_t NumPoint = Points.size();
	std::cout << "POINTS " << NumPoint << '\n';
	ofile << "POINTS " << NumPoint << " double\n";

	Type x, y, z;

	for (size_t i = 0; i < NumPoint; i++) {
		ofile << Points[i].Point[0] << ' ' << Points[i].Point[1] << ' ' << Points[i].Point[2] << '\n';
	}

	size_t NumCells = Cells.size();
	std::cout << "CELLS " << NumCells << '\n';
	ofile << "CELLS " << NumCells << ' ' << NumCells * 5 << '\n';

	Type a, b, c, d, drop;
	for (size_t i = 0; i < NumCells; i++) {
		ofile << 4 << ' ' << Cells[i].idCell[0] << ' ' << Cells[i].idCell[1] << ' ' << Cells[i].idCell[2] << ' ' << Cells[i].idCell[3] << '\n';
	}

	ofile << "CELL_TYPES " << NumCells << '\n';
	for (size_t i = 0; i < NumCells; i++) {
		ofile << 10 << '\n';
	}


	ofile.close();


	ofile.open(NameSol2);
	for (size_t i = 0; i < Faces.size(); i++) {
		ofile << Faces[i].idCell[0] << ' ' << Faces[i].idCell[1] << ' ' << Faces[i].idCell[2] << '\n';
	}
	ofile.close();

	return 0;
}

template<typename Type>
size_t FromVTKToMySet(std::string NameVTK, std::vector<MyPoint<Type>>& Points, std::vector<MyCell<Type>>& Cells) {

	Points.clear();
	Cells.clear();
	vtkSmartPointer<vtkUnstructuredGrid> unstructuredgrid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
	ReadFileVTK(NameVTK, unstructuredgrid, true);


	size_t NumPoint = unstructuredgrid->GetNumberOfPoints();

	Type P[3];
	MyPoint<Type> point;

	for (size_t i = 0; i < NumPoint; i++) {
		unstructuredgrid->GetPoint(i, P);
		point.Point[0] = P[0];
		point.Point[1] = P[1];
		point.Point[2] = P[2];
		Points.push_back(point);
	}

	size_t NumCells = unstructuredgrid->GetNumberOfCells();


	MyCell<Type> cell;
	for (size_t i = 0; i < NumCells; i++) {
		vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(i)->GetPointIds();
		cell.idCell[0] = idp->GetId(0);
		cell.idCell[1] = idp->GetId(1);
		cell.idCell[2] = idp->GetId(2);
		cell.idCell[3] = idp->GetId(3);
		Cells.push_back(cell);
	}

	std::cout << "NumPoint: " << NumPoint << '\n';
	std::cout << "NumCells: " << NumCells << '\n';

	return 0;
}

template<typename FileName, typename Type>
Type OutGrid(FileName NameSol, std::vector<MyPoint<Type>>& Points, std::vector<MyCell<Type>>& Cells) {

	std::ofstream ofile;
	std::string str;
	ofile.open(NameSol);

	ofile << "# vtk DataFile Version 3.0\n";
	ofile << "Elemensts Volumes\n";
	ofile << "ASCII\n";
	ofile << "DATASET UNSTRUCTURED_GRID\n";

	size_t NumPoint = Points.size();
	std::cout << "POINTS " << NumPoint << '\n';
	ofile << "POINTS " << NumPoint << " double\n";

	Type x, y, z;

	for (size_t i = 0; i < NumPoint; i++) {
		ofile << Points[i].Point[0] << ' ' << Points[i].Point[1] << ' ' << Points[i].Point[2] << '\n';
	}

	size_t NumCells = Cells.size();
	std::cout << "CELLS " << NumCells << '\n';
	ofile << "CELLS " << NumCells << ' ' << NumCells * 5 << '\n';

	Type a, b, c, d, drop;
	for (size_t i = 0; i < NumCells; i++) {
		ofile << 4 << ' ' << Cells[i].idCell[0] << ' ' << Cells[i].idCell[1] << ' ' << Cells[i].idCell[2] << ' ' << Cells[i].idCell[3] << '\n';
	}

	ofile << "CELL_TYPES " << NumCells << '\n';
	for (size_t i = 0; i < NumCells; i++) {
		ofile << 10 << '\n';
	}


	ofile.close();


	return 0;
}

template<typename Type>
size_t FullConvergence(std::string NameFileVTK, const std::string NameFileDataGrid, std::string NameFileOut, const size_t CountGrid,
	vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Type(*Sol)(Type*&)) {

	std::ofstream ofile;
	ofile.open(NameFileOut);
	ofile << "Normmax" << "     " << "Norm2\n";
	ofile.close();
	for (size_t H = 0; H < CountGrid; H++)
	{
		clock_t start_clock = clock();
		if (ReadFileVTK(NameFileVTK + std::to_string(H) + ".vtk", unstructuredgrid, true)) {
			std::cout << " Reading error!\n";
			return 1;
		}
		clock_t end_clock = clock();
		std::cout << "\n Read file VTK" << H << " :" << (Type)(end_clock - start_clock) / CLOCKS_PER_SEC << "\n\n";

		Type averageH;
		WriteFileNeighbors(NameFileDataGrid, NameFileDataGrid + "dV.txt", averageH, unstructuredgrid);

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

		start_clock = clock();
		MakeMatrix(NameFileDataGrid, unstructuredgrid, VectorSol, VectorQ, MatrixEnergy, 0);
		end_clock = clock();
		std::cout << "\n Make full Matrix: " << (Type)(end_clock - start_clock) / CLOCKS_PER_SEC << "\n\n";

		Type eps = 0.00001;
		size_t MaxIter = 100000;

		start_clock = clock();
		SolveSparseTimeInDependentSystemEigenBiCGSTAB(SizeGrid, MaxIter, VectorSol, VectorQ, MatrixEnergy);
		end_clock = clock();
		std::cout << "\n Solving system: " << (Type)(end_clock - start_clock) / CLOCKS_PER_SEC << "\n\n";

		ofile.open(NameFileOut, std::ios::app);

		ofile << NormMax(SizeGrid, VectorSol, Sol) << "    "
			<< NormGridF(NameFileDataGrid + "dV.txt", SizeGrid, VectorSol, Sol) <<
			"   " << averageH << '\n';
		ofile.close();


		ClearMemory(SizeGrid, MatrixEnergy);
		ClearMemory(VectorQ);
		VectorSol.clear();
	}

	return 0;
}




