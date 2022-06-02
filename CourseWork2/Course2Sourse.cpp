#include "Course2Header.h"
size_t Intersection(const Type* A, const Type* B, const Type* C, Type*& X0,  Type*& n, Type*& res) {

	Type a, b, c, d;
	Type t;

	a = A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2]);
	b = A[0] * (C[2] - B[2]) + B[0] * (A[2] - C[2]) + C[0] * (B[2] - A[2]);
	c = A[0] * (B[1] - C[1]) + B[0] * (C[1] - A[1]) + C[0] * (A[1] - B[1]);
	d = A[0] * (C[1] * B[2] - B[1] * C[2]) + B[0] * (A[1] * C[2] - C[1] * A[2]) + C[0] * (B[1] * A[2] - A[1] * B[2]);

	t = -(a * X0[0] + b * X0[1] + c * X0[2] + d) / (a * n[0] + b * n[1] + c * n[2]);

	for (size_t i = 0; i < 3; i++)
		res[i] = (n[i] * t + X0[i]);

	return 0;
}



/*¬вод/¬ывод данных*/

std::ostream& operator<< (std::ostream& out, const MyPair<Type>& pair)
{
	out << "(" << pair.Id << ", " << pair.val << ")";
	return out;
}
std::ofstream& operator<< (std::ofstream& out, const MySolution<Type>& Sol) {
	out << Sol.Point[0] << " " << Sol.value << '\n';
	return out;
}

size_t ReadSettings(const char* NameFileSet, std::string& NameFileGrid, std::string& NameFileOut, std::string& NameFileDataGrid) {

	std::ifstream ifile;
	std::string str;
	ifile.open(NameFileSet);

	getline(ifile, NameFileGrid);
	getline(ifile, NameFileOut);
	getline(ifile, NameFileDataGrid);


	ifile.close();
	return 0;
}
size_t ReadFileVTK(const std::string NameFileVTK, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, const bool print) {

	/*„тение исходного файла и запись в vtkUnstructuredGrid*/
	vtkSmartPointer<vtkGenericDataObjectReader> vtkreader =
		vtkSmartPointer<vtkGenericDataObjectReader>::New();
	vtkreader->ReadAllScalarsOn();
	vtkreader->SetFileName(NameFileVTK.c_str());
	vtkreader->Update();

	if (vtkreader->IsFileUnstructuredGrid()) {
		unstructuredgrid = vtkreader->GetUnstructuredGridOutput();
		unstructuredgrid->Modified();
	}
	else {
		std::cout << "Error read file\n";
		return 1;
	}

	if (print) {
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfPoints() << " points." << std::endl;
		std::cout << "Grid has " << unstructuredgrid->GetNumberOfCells() << " cells." << std::endl;
	}

	vtkreader->GetOutput()->GlobalReleaseDataFlagOn();
	return 0;
}
size_t WriteFileVTK(const std::string NameFileOut, vtkSmartPointer<vtkUnstructuredGrid>& UGrid) {

	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileName(NameFileOut.c_str());
	writer->SetInputData(UGrid);
	writer->Write();

	return 0;
}
size_t WriteFileSolution(const std::string NameFileOut, const size_t n, std::vector<MySolution<Type>>& VectorSol,size_t NumCoord) {
	std::ofstream ofile;
	ofile.open(NameFileOut.c_str());
	if (!ofile.is_open()) {
		std::cout << "File wasn't opened\n";
		return 1;
	}
	if (NumCoord > 3) {
		cout << "Error Size\n";
		return 1;
	}

	if (Cond2) {
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < NumCoord; j++)
				ofile << VectorSol[i].Point[j] << ' ';
			ofile << VectorSol[i].value << '\n';
		}
	}
	else {
		for (size_t i = 0; i < n; i++) {
			Type r = sqrt(pow(VectorSol[i].Point[0], 2) + pow(VectorSol[i].Point[1], 2) + pow(VectorSol[i].Point[2], 2));
			ofile << r << ' ';
			ofile << VectorSol[i].value << '\n';
		}
	}

	ofile.close();
	return 0;
}
size_t WriteFileNeighbors(const std::string NameFileOut, const std::string NameFile_dv, Type& averageH, vtkSmartPointer<vtkUnstructuredGrid>& UGrid) {
	std::ofstream ofile;
	ofile.open(NameFileOut.c_str(),std::ios::binary);
	if (!ofile.is_open()) {
		std::cout << "File wasn't opened\n";
		return 1;
	}


	std::ofstream ofile2;
	ofile2.open(NameFile_dv);

	ofile2 << std::setprecision(16);
	ofile << std::setprecision(16);

	const size_t NumberOfCells = UGrid->GetNumberOfCells();
	vtkSmartPointer<vtkIdList> idc =
		vtkSmartPointer<vtkIdList>::New();
	Type SquareFaces[4];
	//нормали к гран€м
	std::vector<Type*> NormalsCurCell;
	NormalsCurCell.resize(4);
	for (size_t i = 0; i < 4; i++)
		NormalsCurCell[i] = new Type[3];
	// центры соседних €чеек
	std::vector<Type*> CenterCellNeighbors;
	CenterCellNeighbors.resize(4);
	for (size_t i = 0; i < 4; i++)
		CenterCellNeighbors[i] = new Type[3];

	std::vector<int> IdNeighbors(4);
	Type VolumeCell = 0;
	Type* CenterCell = new Type[3];
	
	static Type Minh = 100;
	 averageH = 0;

	for (size_t curCell = 0; curCell < NumberOfCells; ++curCell){
		
		
		// 4 нормали и площади к текущей €чейке
		for (size_t i = 0; i < 4; i++)
			NormalAndSquareFace(curCell, i, UGrid, NormalsCurCell[i], SquareFaces[i]);
		// объем текущей €чейки
		Type CurVolume = GetVolumeCell(curCell, UGrid);

		CenterOfTetra(curCell, UGrid, CenterCell);
		// центры, Id соседей
		for (size_t i = 0; i < 4; i++) {
			// Id €чейки(idc) к i-ой грани
			UGrid->GetCellNeighbors(curCell, UGrid->GetCell(curCell)->GetFace(i)->GetPointIds(), idc);

			// вычисление Id,центра,значение соседа i-ой грани
			if (idc->GetNumberOfIds() == 1) {
				IdNeighbors[i] = idc->GetId(0);
				CenterOfTetra(idc->GetId(0), UGrid, CenterCellNeighbors[i]);
			}
			else if (idc->GetNumberOfIds() == 0) {
				IdNeighbors[i] = -1;
				//cout << "There aren't Neighbors => border\n";
			}
			else {
				IdNeighbors[i] = -2;
				cout << "More then 1 neighbors???\n";
			}

		
		}

		VolumeCell=GetVolumeCell(curCell, UGrid);
		CenterOfTetra(curCell, UGrid, CenterCell);

		Type curH = 3 * VolumeCell / (SquareFaces[0] + SquareFaces[1] + SquareFaces[2] + SquareFaces[3]);
		ofile2 << VolumeCell << '\n';
		averageH += curH;
		if (curH < Minh)
			Minh = curH;
		
		ofile << curCell << ' ';
		for (size_t j = 0; j < 4; j++)
			ofile << IdNeighbors[j] << ' ';
		ofile << '\n';

		for (size_t j = 0; j < 4; j++)
			ofile << SquareFaces[j] << ' ';
		ofile << VolumeCell << '\n';

		for (size_t j = 0; j < 3; j++)
			ofile << CenterCell[j] << ' ';
		ofile << '\n';

		for (size_t i = 0; i < 4; i++) {
			for (size_t j = 0; j < 3; j++)
				ofile << CenterCellNeighbors[i][j]<<' ';
			ofile << '\n';
		}
		for (size_t i = 0; i < 4; i++) {
			for (size_t j = 0; j < 3; j++)
				ofile << NormalsCurCell[i][j] << ' ';
			ofile << '\n';
		}
		//ofile << "\n\n";
	}

	ofile.close();
	ofile2.close();
	cout << "h: " << Minh << '\n';

	for (size_t i = 0; i < 4; i++) {
		delete[] NormalsCurCell[i];
		delete[] CenterCellNeighbors[i];
	}
	ClearMemory(CenterCell);
	averageH /= NumberOfCells;
	return 0;
}
size_t WriteFileSolution(const std::string NameFileOut, const size_t n, std::vector<MySolution<Type>>& VectorSol, 
	vtkSmartPointer<vtkUnstructuredGrid>& UGrid) {

	vtkSmartPointer<vtkDoubleArray> IllumArray =
		vtkSmartPointer<vtkDoubleArray>::New();
	
	if (Ellips)
		for (size_t i = 0; i < n; i++)
			IllumArray->InsertNextTuple1(VectorSol[i].value * U_0);// *AA / Velocity);
	else
		for (size_t i = 0; i < n; i++)
			IllumArray->InsertNextTuple1(VectorSol[i].value);

	

	vtkSmartPointer<vtkUnstructuredGrid> ungrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	ungrid = UGrid;
	ungrid->GetCellData()->SetActiveScalars("energy");
	ungrid->GetCellData()->SetScalars(IllumArray);
	

	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileName(NameFileOut.c_str());
	writer->SetInputData(ungrid);
	writer->Write();
	return 0;
}

size_t ReadFileNeighbors(std::ifstream& ifile, size_t& curCell,  int*& IdNeighbors, Type& VolumeCell,
	Type*& CenterCell,  Type*& SquareFaces, std::vector<Type*>& NormalsCurCell, std::vector<Type*>& CenterCellNeighbors) {
	std::string str;

	ifile >> curCell;
//	cout << "CurCell: " << curCell;

//	cout << "\nIds: ";
	for (size_t j = 0; j < 4; j++) {
		ifile >> IdNeighbors[j];
//		cout<< IdNeighbors[j] << ' ';
	}
	getline(ifile, str);

//	cout << "\nSquare: ";
	for (size_t j = 0; j < 4; j++) {
		ifile >> SquareFaces[j];
	//	cout << SquareFaces[j] << ' ';
	}
	ifile >> VolumeCell;
//	cout << "\nVolume " << VolumeCell;

	getline(ifile, str);

//	cout << "\nCenterCell: ";
	for (size_t j = 0; j < 3; j++) {
		ifile >> CenterCell[j];
//		cout <<  CenterCell[j] << ' ';
	}

	getline(ifile, str);
//	cout << "\nCenterCellNeighbors: ";
	for (size_t i = 0; i < 4; i++) {
		for (size_t j = 0; j < 3; j++) {
			ifile >> CenterCellNeighbors[i][j];
	//		cout << CenterCellNeighbors[i][j] << ' ';
		}
		getline(ifile, str);
//		cout << '\n';
	}

	//cout << "\nNormalsCurCell: ";
	for (size_t i = 0; i < 4; i++) {
		for (size_t j = 0; j < 3; j++) {
			ifile >> NormalsCurCell[i][j];
	//		cout << NormalsCurCell[i][j]<<' ';
		}
			getline(ifile, str);
	//		cout << '\n';
	}


	return 0;
}

/*ќперации*/


Type NormMax(const size_t n, std::vector<MySolution<Type>>& x, Type(*f)(Type*&))
{// норма разности

	Type max = -1;
	size_t index = -1;
	Type buf;
	for (size_t i = 0; i < n; i++) {
		buf = abs(f(x[i].Point) - x[i].value);
		if (buf > max){
			max = buf;
			index = i;
		}
	}
	return max;
}
Type Norm(Type*& v1) {
	Type sum = 0;
	for (size_t i = 0; i < 3; i++)
		sum += pow(v1[i], 2);
	return sqrt(sum);
}
Type Norm2(const size_t n, Type*& x1, Type*& x2)
{// норма разности
	Type sum = 0;
	for (size_t i = 0; i < n; i++)
		sum += abs(x1[i] - x2[i]);
	return sum;
}
Type NormMax(const size_t n, Type*& x1, Type*& x2)
{// норма разности
	Type max = -1;
	Type buf;
	for (size_t i = 0; i < n; i++) {
		buf = abs(x1[i] - x2[i]);
		if (buf>max);
		max = buf;
	}
	return max;
}

Type Norm2(const size_t n, std::vector<MySolution<Type>>& x1, Type*& x2)
{// норма разности
	Type sum = 0;
	for (size_t i = 0; i < n; i++)
		sum += abs(x1[i].Point[0] - x2[i]);
	return sum;
}
size_t Normalize(Type*& v1) {

	Type norm = Norm(v1);
	for (int i = 0; i < 3; i++)
		v1[i] /= norm;

	return 0;
}
Type Distance(const Type* h1, const Type* h2) {
	Type s = 0;
	for (size_t i = 0; i < 3; i++)
		s += pow((h2[i] - h1[i]), 2);

	return sqrt(s);
}

/*ѕреобразование данных*/

size_t FromMyToEigen(const size_t n, MyPair<Type>**& MatrixEnergy, Eigen::MatrixXd& matrix, const size_t numNonZero—olumns ) {
	// перевод начальной матрицы в полную eigen матрицу

	matrix = Eigen::MatrixXd::Zero(n, n);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < numNonZero—olumns; j++)
			if (MatrixEnergy[i][j].Id >= 0)
				// если индекс отрицательный, значит это граница области.
				// в MatrixEnergy  мусор, в Eigen нули ( из MatrixXd::Zero)
				matrix(i, MatrixEnergy[i][j].Id) = MatrixEnergy[i][j].val;
	}

	return 0;
}
size_t FromMyToEigen(const size_t n, MyPair<Type>**& MatrixEnergy, Eigen::MatrixXd& matrix, Type*& MyVector, Eigen::VectorXd& EigenVector,
	const size_t numNonZero—olumns ) {
	// перевод начальной матрицы и вектора в полную eigen матрицу и вектор
	matrix = Eigen::MatrixXd::Zero(n, n);
	EigenVector = Eigen::VectorXd::Zero(n);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < numNonZero—olumns; j++)
			if (MatrixEnergy[i][j].Id >= 0)
				// если индекс отрицательный, значит это граница области.
				// в MatrixEnergy  мусор, в Eigen нули ( из MatrixXd::Zero)
				matrix(i, MatrixEnergy[i][j].Id) = MatrixEnergy[i][j].val;
		EigenVector(i) = MyVector[i];
	}

	return 0;
}
size_t FromEigenToMySol(const size_t n, Eigen::VectorXd& EigenVector, std::vector<MySolution<Type>>& VectorSol) {
	// присвоение результата из eigen в свой формат

	for (size_t i = 0; i < n; i++)
		VectorSol[i].value = EigenVector(i);
	return 0;
}

/*¬ычислени€ компонент разностной схемы*/

Type SquareFace(size_t NumberCell, size_t NumberFace, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid) {
	// 1/2 нормали
	Type S = 0;
	vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(NumberCell)->GetFace(NumberFace)->GetPointIds();

	Type P0[3], P1[3], P2[3];
	unstructuredgrid->GetPoint(idp->GetId(0), P0);
	unstructuredgrid->GetPoint(idp->GetId(1), P1);
	unstructuredgrid->GetPoint(idp->GetId(2), P2);

	Type a[3], b[3];
	for (size_t i = 0; i < 3; i++) {
		a[i] = P1[i] - P0[i];
		b[i] = P2[i] - P0[i];
	}

	S = pow(a[1] * b[2] - a[2] * b[1], 2) + pow(a[0] * b[2] - a[2] * b[0], 2) + pow(a[0] * b[1] - a[1] * b[0], 2);
	return 0.5 * sqrt(S);
}
Type GetVolumeCell(size_t NumberCell, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid) {
	Type V = 0;
	Type P0[3], P1[3], P2[3], P3[3];

	vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(NumberCell)->GetPointIds();
	unstructuredgrid->GetPoint(idp->GetId(0), P0);
	unstructuredgrid->GetPoint(idp->GetId(1), P1);
	unstructuredgrid->GetPoint(idp->GetId(2), P2);
	unstructuredgrid->GetPoint(idp->GetId(3), P3);

	Type a[3], b[3], c[3];
	for (size_t i = 0; i < 3; i++) {
		a[i] = P1[i] - P0[i];
		b[i] = P2[i] - P0[i];
		c[i] = P3[i] - P0[i];
	}

	V = a[0] * (b[1] * c[2] - c[1] * b[2]) - a[1] * (b[0] * c[2] - c[0] * b[2]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
	return abs(V) / 6;
}
size_t CenterOfTetra(size_t NumberCell, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Type*& PointInTetra, bool old) {

	Type P0[3], P1[3], P2[3], P3[3];
	vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(NumberCell)->GetPointIds();

	unstructuredgrid->GetPoint(idp->GetId(0), P0);
	unstructuredgrid->GetPoint(idp->GetId(1), P1);
	unstructuredgrid->GetPoint(idp->GetId(2), P2);
	unstructuredgrid->GetPoint(idp->GetId(3), P3);

	Type c = Distance(P1, P2);
	Type b = Distance(P1, P3);
	Type a = Distance(P2, P3);
	Type Sum = a + b + c;
	Type Center[3] = { (a * P1[0] + b * P2[0] + c * P3[0]) / Sum,
					   (a * P1[1] + b * P2[1] + c * P3[1]) / Sum,
					   (a * P1[2] + b * P2[2] + c * P3[2]) / Sum, }; //центр треугольник 123

	for (size_t i = 0; i < 3; i++) {
		PointInTetra[i] = (Center[i] + P0[i]) / 2;
	}

	return 0;
}
size_t CenterOfTetra(size_t NumberCell, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid,Type*& PointInTetra) {

	auto MakeS{ [](Type* P0,Type* P1,Type* P2) {
		Type Sum = 0;
		Type a[3], b[3];
		for (size_t i = 0; i < 3; i++) {
			a[i] = P1[i] - P0[i];
			b[i] = P2[i] - P0[i];
		}

		Sum = pow(a[1] * b[2] - a[2] * b[1], 2) + pow(a[0] * b[2] - a[2] * b[0], 2) + pow(a[0] * b[1] - a[1] * b[0], 2);
		return 0.5 * sqrt(Sum);
} };

	Type P0[3], P1[3], P2[3], P3[3];
	vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(NumberCell)->GetPointIds();

	unstructuredgrid->GetPoint(idp->GetId(0), P0);
	unstructuredgrid->GetPoint(idp->GetId(1), P1);
	unstructuredgrid->GetPoint(idp->GetId(2), P2);
	unstructuredgrid->GetPoint(idp->GetId(3), P3);

	Type Squr[4] = { MakeS(P1,P2,P3),MakeS(P0,P2,P3), MakeS(P0,P1,P3),MakeS(P0,P1,P2) };
	

	Type Sum = Squr[0] + Squr[1] + Squr[2] + Squr[3];
	for (size_t i = 0; i < 3; i++) {
		PointInTetra[i] = (Squr[0] * P0[i] + Squr[1] * P1[i] + Squr[2] * P2[i] + Squr[3] * P3[i]) / Sum;
	}
	return 0;
}
size_t NormalToFace(size_t NumberCell, size_t NumberFace, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Type*& n) {
	vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(NumberCell)->GetFace(NumberFace)->GetPointIds();

	Type P0[3], P1[3], P2[3];
	unstructuredgrid->GetPoint(idp->GetId(0), P0);
	unstructuredgrid->GetPoint(idp->GetId(1), P1);
	unstructuredgrid->GetPoint(idp->GetId(2), P2);

	Type a[3], b[3];
	for (size_t i = 0; i < 3; i++) {
		a[i] = P1[i] - P0[i];
		b[i] = P2[i] - P0[i];
	}
	n[0] = a[1] * b[2] - a[2] * b[1];
	n[1] = -a[0] * b[2] + a[2] * b[0];
	n[2] = a[0] * b[1] - a[1] * b[0];

	Normalize(n);
	return 0;
}

size_t NormalAndSquareFace(size_t NumberCell, size_t NumberFace, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, Type*& n, Type& S) {
	vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(NumberCell)->GetFace(NumberFace)->GetPointIds();

	Type P0[3], P1[3], P2[3];
	unstructuredgrid->GetPoint(idp->GetId(0), P0);
	unstructuredgrid->GetPoint(idp->GetId(1), P1);
	unstructuredgrid->GetPoint(idp->GetId(2), P2);

	Type a[3], b[3];
	for (size_t i = 0; i < 3; i++) {
		a[i] = P1[i] - P0[i];
		b[i] = P2[i] - P0[i];
	}
	n[0] = a[1] * b[2] - a[2] * b[1];
	n[1] = -a[0] * b[2] + a[2] * b[0];
	n[2] = a[0] * b[1] - a[1] * b[0];

	S = sqrt(pow(n[0], 2) + pow(n[1], 2) + pow(n[2], 2)) / 2;
	Normalize(n);

	vtkSmartPointer<vtkIdList> idp2 = unstructuredgrid->GetCell(NumberCell)->GetPointIds();

	size_t id;
	for (size_t i = 0; i < 4; i++){
		int count = 0;
		for (size_t j = 0; j < 3; j++)
			if (idp2->GetId(i) != idp->GetId(j))
				count++;
		if (count == 3) {
			id = i;
			break;
		}

	}

	Type sum = 0;
	Type P3[3];
	unstructuredgrid->GetPoint(idp2->GetId(id), P3);
	/*for (size_t i = 0; i < 3; i++){
		sum += n[i] * (P3[i] - P0[i]);
	}*/

	sum = P1[0] * (P2[1] - P3[1]) * P0[2] + P0[0] * (P3[1] - P2[1]) * P1[2] +
		P0[0] * (P1[1] - P3[1]) * P2[2] + P2[2] * (P1[0] * P3[1] - P1[0] * P0[1]) +
		P3[0] * (P0[2] * (P1[1] - P2[1]) + P1[2] * (P2[1] - P0[1]) + P2[2] * (P0[1] - P1[1]))
		+ P3[2] * (P1[0] * (P0[1] - P2[1]) + P0[0] * (P2[1] - P1[1])) +
		P2[0] * (P0[2] * (P3[1] - P1[1]) + P1[2] * (P0[1] - P3[1]) + P3[2] * (P1[1] - P0[1]));

	if (sum < 0)
		for (size_t i = 0; i < 3; i++)
			n[i] *= -1;
	return 0;
}
Type ProjectionOnNormal(Type*& n, Type*& Point1, Type*& Point2) {
	/*Type r1 = 0;
	Type r2 = 0;
	for (size_t i = 0; i < 3; i++){
		r1 += pow(Point1[i], 2);
		r2 += pow(Point2[i], 2);
	}
	return ((sqrt(r2) - sqrt(r1))*n[0] );*/

	Type Projection = 0;
	for (size_t i = 0; i < 3; i++)
		Projection += (Point2[i] - Point1[i]) * n[i]; // Norm(n)==1
	return abs(Projection);
}

/*—оздание матрицы*/

size_t MakeOneRow(const size_t NumCell, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid, 
	MyPair<Type>**& MatrixEnergy, std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ) {

	//площади граней
	Type SquareFaces[4];
	//нормали к гран€м
	std::vector<Type*> NormalsCurCell;
	NormalsCurCell.resize(4);
	for (size_t i = 0; i < 4; i++)
		NormalsCurCell[i] = new Type[3];
	// центры соседних €чеек
	std::vector<Type*> CenterCellNeighbors;
	CenterCellNeighbors.resize(4);
	for (size_t i = 0; i < 4; i++)
		CenterCellNeighbors[i] = new Type[3];
	//центр текущей €чейки
	Type* CenterCurCell = new Type[3];
	// Id соседних €чеек
	vtkSmartPointer<vtkIdList> idc =
		vtkSmartPointer<vtkIdList>::New();

	// id всех соседей
	int* IdNeighbors = new int[4];

	//alpha в текущей €чейки
	Type CurAlpha;
	//lambda в текущей €чейки
	Type CurLambda;
	// отношени€ Lambda/v*S/l
	Type Frac[4];


	//¬џ„»—Ћ≈Ќ»я:

	// 4 нормали и площади к текущей €чейке
	for (size_t i = 0; i < 4; i++)
		NormalAndSquareFace(NumCell, i, unstructuredgrid, NormalsCurCell[i], SquareFaces[i]);
	// объем текущей €чейки
	Type CurVolume = GetVolumeCell(NumCell, unstructuredgrid);

	CenterOfTetra(NumCell, unstructuredgrid, CenterCurCell);
	// центры, Id соседей
	for (size_t i = 0; i < 4; i++) {
		// Id €чейки(idc) к i-ой грани
		unstructuredgrid->GetCellNeighbors(NumCell, unstructuredgrid->GetCell(NumCell)->GetFace(i)->GetPointIds(), idc);

		// вычисление Id,центра,значение соседа i-ой грани
		if (idc->GetNumberOfIds() == 1) {
			IdNeighbors[i] = idc->GetId(0);
			CenterOfTetra(idc->GetId(0), unstructuredgrid, CenterCellNeighbors[i]);
		}
		else if (idc->GetNumberOfIds() == 0) {
			IdNeighbors[i] = -1;
			//cout << "There aren't Neighbors => border\n";
		}
		else {
			IdNeighbors[i] = -2;
			cout << "More then 1 neighbors???\n";
		}
	}

	VectorSol[NumCell].SetData(CenterCurCell);
	//значени€ в €чейке
	VectorQ[NumCell] = Q(CenterCurCell);
	CurAlpha = Alpha(CenterCurCell);
	CurLambda = Lambda(CenterCurCell);

	//значение отношений Li*Sj/(Vi*(Xi-Xj)nj)
	for (size_t i = 0; i < 4; i++) {
		if (IdNeighbors[i] != -1) {
			Type Proj = ProjectionOnNormal(NormalsCurCell[i], CenterCellNeighbors[i], CenterCurCell);
			Frac[i] = -(CurLambda * SquareFaces[i]) / (Proj * CurVolume);
		}
		else Frac[i] = 0;
	}

	//заполнеие матрицы
	MatrixEnergy[NumCell][0].Id = NumCell;
	MatrixEnergy[NumCell][0].val = (Frac[0] + Frac[1] + Frac[2] + Frac[3]) + CurAlpha;
	for (size_t j = 0; j < 4; ++j) {
		MatrixEnergy[NumCell][j + 1].Id = IdNeighbors[j];
		MatrixEnergy[NumCell][j + 1].val = -Frac[j];
	}

	// удаление пам€ти
	for (size_t i = 0; i < 4; i++) {
		delete[] NormalsCurCell[i];
		delete[] CenterCellNeighbors[i];
	}
	ClearMemory(CenterCurCell, IdNeighbors);
	return 0;
}
size_t MakeOneRow(const size_t NumCell, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid,
	std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ, std::vector <Eigen::Triplet <Type>>& tripletList, MyPair<Type>**& MatrixEnergy) {

	//площади граней
	Type* SquareFaces = new Type[4];
	//нормали к гран€м
	std::vector<Type*> NormalsCurCell;
	NormalsCurCell.resize(4);
	for (size_t i = 0; i < 4; i++)
		NormalsCurCell[i] = new Type[3];
	// центры соседних €чеек
	std::vector<Type*> CenterCellNeighbors;
	CenterCellNeighbors.resize(4);
	for (size_t i = 0; i < 4; i++)
		CenterCellNeighbors[i] = new Type[3];
	//центр текущей €чейки
	Type* CenterCurCell = new Type[3];
	// Id соседних €чеек
	vtkSmartPointer<vtkIdList> idc =
		vtkSmartPointer<vtkIdList>::New();

	// id всех соседей
	int* IdNeighbors = new int[4];// = new int[4];


	//alpha в текущей €чейки
	Type CurAlpha;
	//lambda в текущей €чейки
	Type CurLambda;
	// отношени€ Lambda/v*S/l
	Type Frac[4];


	//¬џ„»—Ћ≈Ќ»я:

	// 4 нормали и площади к текущей €чейке
	for (size_t i = 0; i < 4; i++)
		NormalAndSquareFace(NumCell, i, unstructuredgrid, NormalsCurCell[i], SquareFaces[i]);
	// объем текущей €чейки
	Type CurVolume = GetVolumeCell(NumCell, unstructuredgrid);

	CenterOfTetra(NumCell, unstructuredgrid, CenterCurCell);
	// центры, Id соседей
	for (size_t i = 0; i < 4; i++) {
		// Id €чейки(idc) к i-ой грани
		unstructuredgrid->GetCellNeighbors(NumCell, unstructuredgrid->GetCell(NumCell)->GetFace(i)->GetPointIds(), idc);

		// вычисление Id,центра,значение соседа i-ой грани
		if (idc->GetNumberOfIds() == 1) {
			IdNeighbors[i] = idc->GetId(0);
			CenterOfTetra(idc->GetId(0), unstructuredgrid, CenterCellNeighbors[i]);
		}
		else if (idc->GetNumberOfIds() == 0) {
			IdNeighbors[i] = -1;
			//cout << "There aren't Neighbors => border\n";
		}
		else {
			IdNeighbors[i] = -2;
			cout << "More then 1 neighbors???\n";
		}
	}

	VectorSol[NumCell].SetData(CenterCurCell);
	//значени€ в €чейке
	VectorQ[NumCell] = Q(CenterCurCell);
	CurAlpha = Alpha(CenterCurCell);
	CurLambda = Lambda(CenterCurCell);

	//значение отношений Li*Sj/(Vi*(Xi-Xj)nj)
	for (size_t i = 0; i < 4; i++) {
		if (IdNeighbors[i] != -1) {
			Type Proj = ProjectionOnNormal(NormalsCurCell[i], CenterCellNeighbors[i], CenterCurCell);
			Frac[i] = (CurLambda * SquareFaces[i]) / (Proj * CurVolume);
		}
		else Frac[i] = 0;
	}

	//заполнеие матрицы
	MatrixEnergy[NumCell][0].Id = NumCell;
	MatrixEnergy[NumCell][0].val = (Frac[0] + Frac[1] + Frac[2] + Frac[3]) + CurAlpha;
	tripletList.push_back(Eigen::Triplet<Type>(NumCell, NumCell, (Frac[0] + Frac[1] + Frac[2] + Frac[3]) + CurAlpha));
	for (size_t j = 0; j < 4; ++j) {
		MatrixEnergy[NumCell][j + 1].Id = IdNeighbors[j];
		MatrixEnergy[NumCell][j + 1].val = -Frac[j];
		if (IdNeighbors[j] >= 0)
			tripletList.push_back(Eigen::Triplet<Type>(NumCell, IdNeighbors[j], -Frac[j]));
	}



	// удаление пам€ти
	for (size_t i = 0; i < 4; i++) {
		delete[] NormalsCurCell[i];
		delete[] CenterCellNeighbors[i];
	}
	ClearMemory(CenterCurCell);// , IdNeighbors);


	return 0;
}
size_t MakeMatrix(const std::string NameFileDataGrid, vtkSmartPointer<vtkUnstructuredGrid>& unstructuredgrid,
	std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ, MyPair<Type>**& MatrixEnergy, const Type t) {

	size_t SizeGrid = unstructuredgrid->GetNumberOfCells();
	vtkDataArray* density;
	vtkDataArray* AbsorpCoef;
	vtkDataArray* QArray;

	if (Ellips) {
		density = unstructuredgrid->GetCellData()->GetScalars("density");
		AbsorpCoef = unstructuredgrid->GetCellData()->GetScalars("AbsorpCoef");
		QArray = unstructuredgrid->GetCellData()->GetScalars("radEnLooseRate");
	
		/*density = unstructuredgrid->GetCellData()->GetScalars("alpha");
		AbsorpCoef = unstructuredgrid->GetCellData()->GetScalars("alpha");
		QArray = unstructuredgrid->GetCellData()->GetScalars("Q");*/

		std::cout << "AbsorpCoef_Size: " << AbsorpCoef->GetSize() << std::endl;
		std::cout << "Q_Size: " << QArray->GetSize() << std::endl;
	}


	std::ifstream ifile;
	ifile.open(NameFileDataGrid.c_str());
	if (!ifile.is_open()) {
		std::cout << "WTF???\n";
		return -1;
	}

	//площади граней
	Type* SquareFaces = new Type[4];
	//нормали к гран€м
	std::vector<Type*> NormalsCurCell;
	NormalsCurCell.resize(4);
	for (size_t i = 0; i < 4; i++)
		NormalsCurCell[i] = new Type[3];
	// центры соседних €чеек
	std::vector<Type*> CenterCellNeighbors;
	CenterCellNeighbors.resize(4);
	for (size_t i = 0; i < 4; i++)
		CenterCellNeighbors[i] = new Type[3];
	//центр текущей €чейки
	Type* CenterCurCell = new Type[3];

	// id всех соседей
	int* IdNeighbors = new int[4];// = new int[4];
	Type Volume;


	//alpha в текущей €чейки
	Type CurAlpha;
	//lambda в текущей €чейки
	Type CurLambda=1;
	// отношени€ Lambda/v*S/l
	Type Frac[4];
	size_t NumCell;

	Type kappa;
	for (size_t i = 0; i < SizeGrid; i++) {

		ReadFileNeighbors(ifile, NumCell, IdNeighbors, Volume, CenterCurCell, SquareFaces, NormalsCurCell, CenterCellNeighbors);
		VectorSol[NumCell].SetData(CenterCurCell);

		if (Ellips) {
			//значени€ в €чейке
			VectorQ[NumCell] = QArray->GetTuple1(NumCell) * T_orb / U_0;
			kappa = density->GetTuple1(NumCell);// *AbsorpCoef->GetTuple1(NumCell);
			if (kappa < My_Zero)
				kappa = My_Zero;
			CurAlpha = kappa * Velocity * T_orb;
		}
		else {
			//значени€ в €чейке
			VectorQ[NumCell] = Q(CenterCurCell); 

			CurAlpha = Alpha(CenterCurCell);
		}
		

		//значение отношений Li*Sj/(Vi*(Xi-Xj)nj)
		for (size_t i = 0; i < 4; i++) {
			if (IdNeighbors[i] != -1) 
			{
				Type Proj = ProjectionOnNormal(NormalsCurCell[i], CenterCurCell, CenterCellNeighbors[i]);
				if (Ellips) 
				{
					kappa = density->GetTuple1(IdNeighbors[i]);// *AbsorpCoef->GetTuple1(IdNeighbors[i]);
					if (kappa < My_Zero) 
					{
						kappa = My_Zero;
					}
					Frac[i] = ((Velocity * T_orb / (X_0 * X_0) / kappa / 3) * SquareFaces[i]) / (Proj * Volume);
				}
				else
					Frac[i] = (Lambda(CenterCellNeighbors[i]) * SquareFaces[i]) / (Proj * Volume); 
				
			}
			else {
			
				if (Cond2) {
					Frac[i] = 0;


					Type* PointFace = new Type[3];
					Type A[3], B[3], D[3];
					vtkSmartPointer<vtkIdList> idp = unstructuredgrid->GetCell(NumCell)->GetFace(i)->GetPointIds();
					unstructuredgrid->GetPoint(idp->GetId(0), A);
					unstructuredgrid->GetPoint(idp->GetId(1), B);
					unstructuredgrid->GetPoint(idp->GetId(2), D);
					//Intersection(A, B, D, CenterCurCell, NormalsCurCell[i], PointFace);

					{
						Type a =0, b =0, c = 0;
						for (size_t i = 0; i < 3; i++) {
							a += pow(A[i] - B[i], 2);
							b+= pow(A[i] - D[i], 2);
							c+= pow(D[i] - B[i], 2);
						}
						Type sum = sqrt(a) + sqrt(b) + sqrt(c);
						for (size_t i = 0; i < 3; i++)
							PointFace[i] = (sqrt(c) * A[i] + sqrt(b) * B[i] + sqrt(a) * D[i]) / sum;
						
					}
					//cout << "PointFace: " << PointFace[0] << ' ' << PointFace[1] << ' ' << PointFace[2] << '\n';

					VectorQ[NumCell] += SquareFaces[i] * Grad(PointFace, NormalsCurCell[i]) / Volume * Lambda(PointFace);
					//Frac[i]= -SquareFaces[i] * Grad(PointFace, NormalsCurCell[i]) / Volume * Lambda(PointFace);
					delete[] PointFace;
				}
				else {
					if (Ellips)
						//Frac[i] = (1. / 2) * (SquareFaces[i] / Volume);
						Frac[i] = (Velocity*T_orb/(2*X_0)) * (SquareFaces[i] / Volume);
					else
						Frac[i] = (sqrt(3) / 2) * sqrt(CurAlpha * CurLambda) * (SquareFaces[i] / Volume);
				}
			}
		}

		MatrixEnergy[NumCell][0].Id = NumCell;
		//MatrixEnergy[NumCell][0].val = -(Frac[0] + Frac[1] + Frac[2] + Frac[3]) - CurAlpha + 1. / t;
		MatrixEnergy[NumCell][0].val = (Frac[0] + Frac[1] + Frac[2] + Frac[3]) + CurAlpha;
		for (size_t j = 0; j < 4; ++j) {
			if (IdNeighbors[j] != -1) {
				MatrixEnergy[NumCell][j + 1].Id = IdNeighbors[j];
				MatrixEnergy[NumCell][j + 1].val = -Frac[j];
			}
			else {
				MatrixEnergy[NumCell][j + 1].Id = IdNeighbors[j];
				MatrixEnergy[NumCell][j + 1].val = 0;
			}
		}

	}
	ifile.close();
	return 0;
}

/*–ешение системы*/

size_t SolutionSystemSeidel(const Type eps, const size_t MaxIter, const size_t SizeGrid,
	MyPair<Type>**& MatrixEnergy, std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ) {
	Type* previousVariableValues;
	Type* currentVariableValues;
	SetMemory(SizeGrid, previousVariableValues, currentVariableValues);

	for (size_t i = 0; i < SizeGrid; i++)
		previousVariableValues[i] = VectorQ[i];
	size_t count = 0;
	do {
		// cчитаем значени€ неизвестных на текущей итерации
		for (int i = 0; i < SizeGrid; i++) {
			// »нициализируем i-ую неизвестную значением  свободного члена i-ой строки матрицы
			currentVariableValues[i] = VectorQ[i]; //f

			// ¬ычитаем сумму по всем отличным от i-ой неизвестным
			for (int j = 1; j < 5; j++) {
				if (MatrixEnergy[i][j].Id < 0) continue;
				if (MatrixEnergy[i][j].Id < i)
					currentVariableValues[i] -= MatrixEnergy[i][j].val * currentVariableValues[MatrixEnergy[i][j].Id]; //f-Ax^k
				if (MatrixEnergy[i][j].Id > i)
					currentVariableValues[i] -= MatrixEnergy[i][j].val * previousVariableValues[MatrixEnergy[i][j].Id]; 
			}
			// ƒелим на коэффициент при i-ой неизвестной
			currentVariableValues[i] /= MatrixEnergy[i][0].val; //x^k+1/a[i,i]
		}
		std::swap(previousVariableValues, currentVariableValues);
		if (count % 10) {
			if (Norm2(SizeGrid, previousVariableValues, currentVariableValues) <eps) {
				std::cout << " Seidel converged\n";
				break;
			}
			else if (count >MaxIter) {
				std::cout << " Non Solution\n";
				break;
			}
		}
		count++;
	} while (true);

	for (size_t i = 0; i < SizeGrid; i++)
		VectorSol[i].value = currentVariableValues[i];
	ClearMemory(previousVariableValues, currentVariableValues);
	return 0;
}

size_t SolveFullSystemEigenQR(const size_t SizeGrid, std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ, MyPair<Type>**& MatrixEnergy) {
	//ѕр€мое решение стационарной задачи
	// сходитс€ при больших альфа

	Eigen::MatrixXd matrix;
	Eigen::VectorXd Vec;

	//инициализаци€ Eigen(полных) данных
	FromMyToEigen(SizeGrid, MatrixEnergy, matrix, VectorQ, Vec);

	// решение системы полным QR разложением
	Eigen::VectorXd sol = matrix.fullPivHouseholderQr().solve(Vec);
	
	// присвоение результата из eigen в свой формат
	FromEigenToMySol(SizeGrid, sol, VectorSol);

	std::cout << "EigenQR end\n";
	return 0;
}

size_t SolveSparseTimeDependentInDirectSystemGMRES(const Type eps, const size_t SizeGrid, const Type t, const size_t MaxIter,
	std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ, MyPair<Type>**& MatrixEnergy, const size_t numNonZero—olumns) {
	//от времени не€вна€ разреженна€ в формате triplet (3 вектора)

	const size_t N = SizeGrid * 5; // пор€док «јƒј„»
	Type* A;  //вектор данных
	int* row; // вектор индексов строк
	int* col; // вектор индексов столбцов
	Type* curx; // текущее приблежение
	Type* prevx; // приближение на предыдущем шаге
	Type* b;  // вектор правой части
	size_t count = 0; // число ненулевых элементов 

	SetMemory(N, A,curx, prevx, b);
	SetMemory(N,  row, col);

	// заполнение векторв значени€ми
	for (size_t i = 0; i < SizeGrid; i++) {
		for (size_t j = 0; j < numNonZero—olumns; j++)
			if (MatrixEnergy[i][j].Id >= 0) {
				if (i != MatrixEnergy[i][j].Id) {
					A[count] = MatrixEnergy[i][j].val;
					col[count] = MatrixEnergy[i][j].Id;
					row[count++] = i;
				}
				else {
					A[count] = MatrixEnergy[i][j].val;// +1 / t; // смотреть разностную схему
					col[count] = MatrixEnergy[i][j].Id;
					row[count++] = i;
				}
			}
	}

	// задание начальных приблежений
	for (size_t i = 0; i < SizeGrid; i++) {
		curx[i] = prevx[i] = VectorQ[i];
		b[i] = VectorQ[i];
	}
	{
		mgmres_st(SizeGrid, count, row, col, A, curx, b, MaxIter, 50, eps, eps);

		for (size_t i = 0; i < SizeGrid; i++) {
			VectorSol[i].value = curx[i];
			
		}

		ClearMemory(A, row, col, curx, prevx, b);
		return 0;
	}

	size_t iterCount = 0;
	Type norm = 1000;
	do {
		std::swap(curx, prevx);

		//вектор правой части
		for (size_t i = 0; i < SizeGrid; i++)
			b[i] = VectorQ[i] + prevx[i] / t; // см. разностную схему

		mgmres_st(SizeGrid, count, row, col, A, curx, b, MaxIter/100, SizeGrid / 1000, eps, eps);

		if (NormMax(SizeGrid, curx, prevx) / t < eps) {
			std::cout << " GMRES converged\n";
			break;
		}
		else  if (iterCount++ > MaxIter) {
			std::cout << "GMRES iterations more max\n";
			break;
		}

	} while (true);

	for (size_t i = 0; i < SizeGrid; i++)
		VectorSol[i].value = curx[i];

	ClearMemory(A, row, col, curx, prevx, b);
	return 0;

}

size_t SolveSparseTimeDependentInDirectSystemEigenBiCGSTAB(const Type eps, const size_t SizeGrid, const Type t, const size_t MaxIter,
	std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ, MyPair<Type>**& MatrixEnergy,
	const size_t numNonZero—olumns) {
	//от времени не€вна€ разреженна€ в формате triplet, EigenBiCGSTAB

	std::vector <Eigen::Triplet <Type>> tripletList;
	tripletList.reserve(SizeGrid * 5);

	// заполнение триплетов
	for (size_t i = 0; i < SizeGrid; i++) {
		for (size_t j = 0; j < numNonZero—olumns; j++)
			if (MatrixEnergy[i][j].Id >= 0) {
				if (i != MatrixEnergy[i][j].Id)
					tripletList.push_back(Eigen::Triplet<Type>(i, MatrixEnergy[i][j].Id, MatrixEnergy[i][j].val));
				else
					tripletList.push_back(Eigen::Triplet<Type>(i, MatrixEnergy[i][j].Id, MatrixEnergy[i][j].val + 1 / t));
			}
	}

	Eigen::SparseMatrix<Type> mat(SizeGrid, SizeGrid); // разреженна€ матрица системы
	Eigen::VectorXd b(SizeGrid);  // вектор исходной правой части

	// заполнение матрица и вектора
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	for (size_t i = 0; i < SizeGrid; i++) {
		b[i] = VectorQ[i];
	}

	Eigen::BiCGSTAB<Eigen::SparseMatrix<Type>, Eigen::IncompleteLUT<Type>> sol(mat);

	//cout << "\n\n" << mat << "\n\n";
	sol.compute(mat);
	Eigen::VectorXd cur =  Eigen::VectorXd::Zero(SizeGrid);
	Eigen::VectorXd prev = Eigen::VectorXd::Zero(SizeGrid);

	Type norm = 1000;
	size_t count = 0;
	do {
	//	cout << prev[0] << '\n';
		prev = cur;
		cur = sol.solve(b + prev / t);
		norm = (cur - prev).norm();

		if (norm / t < eps) {
			std::cout << " BiCGSTAB time converged\n";
			break;
		}
		else  if (count++ > MaxIter) {
			std::cout << "BiCGSTAB time iterations more max\n";
			break;
		}

	} while (true);

	FromEigenToMySol(SizeGrid, cur, VectorSol);
	return 0;
}

size_t SolveSparseTimeInDependentSystemEigenBiCGSTAB(const size_t SizeGrid, const size_t MaxIter,
	std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ, MyPair<Type>**& MatrixEnergy,
	const size_t numNonZero—olumns) {
	//стационарна€ разреженна€ в формате triplet, EigenBiCGSTAB

	std::vector <Eigen::Triplet <Type>> tripletList;
	tripletList.reserve(SizeGrid * 5);

	Eigen::SparseMatrix<Type> mat(SizeGrid, SizeGrid); // разреженна€ матрица системы
	Eigen::VectorXd b(SizeGrid);  // вектор исходной правой части

	// заполнение триплетов
	for (size_t i = 0; i < SizeGrid; i++) {
		for (size_t j = 0; j < numNonZero—olumns; j++)
			if (MatrixEnergy[i][j].Id >= 0) {
					tripletList.push_back(Eigen::Triplet<Type>(i, MatrixEnergy[i][j].Id, MatrixEnergy[i][j].val));
			}
	}

	

	// заполнение матрица и вектора
	mat.setFromTriplets(tripletList.begin(), tripletList.end());
	for (size_t i = 0; i < SizeGrid; i++) {
		b[i] = VectorQ[i];
	}


	Eigen::BiCGSTAB<::Eigen::SparseMatrix<Type>, Eigen::IncompleteLUT<Type>> solver;
	solver.setMaxIterations(MaxIter);
	solver.compute(mat);
	Eigen::VectorXd x = solver.solve(b);



	FromEigenToMySol(SizeGrid, x, VectorSol);

	if (solver.iterations()< MaxIter) {
		std::cout << " BiCGSTAB  converged\n";
	}
	else {
		std::cout << "BiCGSTAB  iterations more max\n";
	}

	return 0;
}

size_t SolveFullTimeDependentSystemEigenQR(const Type eps, const size_t SizeGrid, const size_t MaxIter, const Type t,
	std::vector<MySolution<Type>>& VectorSol, Type*& VectorQ, MyPair<Type>**& MatrixEnergy,
	const size_t numNonZero—olumns) {

	Eigen::MatrixXd matrix;
	Eigen::VectorXd b;

	//инициализаци€ Eigen(полных) данных
	FromMyToEigen(SizeGrid, MatrixEnergy, matrix, VectorQ, b, numNonZero—olumns);
	for (size_t i = 0; i < SizeGrid; i++)
		matrix(i, i) += 1 / t;

	Eigen::VectorXd cur = Eigen::VectorXd::Zero(16);
	Eigen::VectorXd prev = Eigen::VectorXd::Zero(16);

	Type norm = 1000;
	size_t count = 0;

	do {
		prev = cur;
		cur = matrix.fullPivHouseholderQr().solve(b + prev / t);
		norm = (cur - prev).norm();

		if (norm / t < eps) {
			std::cout << " QR time converged\n";
			break;
		}
		else  if (count++ > MaxIter) {
			std::cout << "QR time iterations more max\n";
			break;
		}

	} while (true);

	FromEigenToMySol(SizeGrid, cur, VectorSol);

	return 0;
}


/*разное*/
size_t From_XmlVtu_To_VTK() {
	/*vtkSmartPointer<vtkXMLUnstructuredGridReader> vtkreader =
		vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();*/

	/*vtkreader->SetFileName("C:\\Users\\Artem\\Desktop\\TetraGrid\\Sphere.vtu");
	vtkreader->Update();

	vtkSmartPointer<vtkUnstructuredGrid> unstructuredgrid = vtkSmartPointer<vtkUnstructuredGrid>::New();

	unstructuredgrid = vtkreader->GetOutput();

	vtkSmartPointer<vtkGenericDataObjectWriter> writer =
		vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileName("C:\\Users\\Artem\\Desktop\\TetraGrid\\Sphere.vtk");
	writer->SetInputData(unstructuredgrid);
	writer->Write();*/

	return 0;
}



Type MakeFace(MyPoint<Type>& A, MyPoint<Type>& B) {

	return sqrt(pow(A.Point[0] - B.Point[0], 2) + pow(A.Point[1] - B.Point[1], 2) + pow(A.Point[2] - B.Point[2], 2));;
}

size_t FracTetra2D(std::vector<MyPoint<Type>>& Points, MyCell<Type> cell) {
	MyPoint<Type> newPoint;
	Type a, b, c, d;
	Type x, y, z;
	a = MakeFace(Points[cell.idCell[0]], Points[cell.idCell[1]]);
	b = MakeFace(Points[cell.idCell[0]], Points[cell.idCell[2]]);
	c = MakeFace(Points[cell.idCell[1]], Points[cell.idCell[2]]);
	Type sum = a + b + c;
	
	x = (Points[cell.idCell[0]].Point[0] * c + Points[cell.idCell[1]].Point[0] * b + Points[cell.idCell[2]].Point[0] * a);
	z = (Points[cell.idCell[0]].Point[2] * c + Points[cell.idCell[1]].Point[2] * b + Points[cell.idCell[2]].Point[2] * a);
	y = (Points[cell.idCell[0]].Point[1] * c + Points[cell.idCell[1]].Point[1] * b + Points[cell.idCell[2]].Point[1] * a);

	newPoint.SetData(x, y, z);
	Points.push_back(newPoint);

	return 0;
}

size_t FracTetra(const size_t NumberCell, std::vector<MyPoint<Type>>& Points, std::vector<MyCell<Type>>& Cells) {


	MyPoint<Type> newPoint;
	auto MakeS{ [](Type*& P0,Type*& P1,Type*& P2) {
		Type Sum = 0;
		Type a[3], b[3];
		for (size_t i = 0; i < 3; i++) {
			a[i] = P1[i] - P0[i];
			b[i] = P2[i] - P0[i];
		}

		Sum = pow(a[1] * b[2] - a[2] * b[1], 2) + pow(a[0] * b[2] - a[2] * b[0], 2) + pow(a[0] * b[1] - a[1] * b[0], 2);
		return 0.5 * sqrt(Sum);
} };
	auto FromPointToType{ [&Points,&Cells,&NumberCell](Type*& P0,Type*& P1,Type*& P2,Type*& P3) {
		for (size_t i = 0; i < 3; i++)
		{
			P0[i] = Points[Cells[NumberCell].idCell[0]].Point[i];
			P1[i] = Points[Cells[NumberCell].idCell[1]].Point[i];
			P2[i] = Points[Cells[NumberCell].idCell[2]].Point[i];
			P3[i] = Points[Cells[NumberCell].idCell[3]].Point[i];
				
		}
		return 0;
	} };



	Type* P0;
	Type* P1;
	Type* P2;
	Type* P3;
	SetMemory(3, P0, P1, P2, P3);
	FromPointToType(P0, P1, P2, P3);


	Type Squr[4] = { MakeS(P1,P2,P3),MakeS(P0,P2,P3), MakeS(P0,P1,P3),MakeS(P0,P1,P2) };


	Type Sum = Squr[0] + Squr[1] + Squr[2] + Squr[3];
	for (size_t i = 0; i < 3; i++) {
		newPoint.Point[i] = (Squr[0] * P0[i] + Squr[1] * P1[i] + Squr[2] * P2[i] + Squr[3] * P3[i]) / Sum;
	}

	Points.push_back(newPoint);
	ClearMemory(P0, P1, P2, P3);


	size_t newid = Points.size() - 1;

	
	MyCell<Type> curCell;
	curCell.SetData(Cells[NumberCell]);

	MyCell<Type> cell;
	cell.SetData(curCell.idCell[0], curCell.idCell[1], curCell.idCell[2], newid);
	Cells.push_back(cell);

	cell.ReSetData(curCell.idCell[0], Cells[NumberCell].idCell[1], Cells[NumberCell].idCell[3], newid);
	Cells.push_back(cell);

	cell.ReSetData(Cells[NumberCell].idCell[1], Cells[NumberCell].idCell[2], Cells[NumberCell].idCell[3], newid);
	Cells.push_back(cell);

	Cells[NumberCell].ReSetData(Cells[NumberCell].idCell[0], Cells[NumberCell].idCell[2], Cells[NumberCell].idCell[3], newid);

	return 0;
}

size_t AddFace(const size_t NumFace, std::vector<MyPoint<Type>>& Points, std::vector<MyCell<Type>>& Cells,
	std::vector<MyFace<Type>>& Faces) {

	MyPoint<Type> newPoint;
	MyCell<Type> newCell;
	Type a, b, c, d;
	Type x, y, z;
	a = MakeFace(Points[Faces[NumFace].idCell[0]], Points[Faces[NumFace].idCell[1]]);
	b = MakeFace(Points[Faces[NumFace].idCell[0]], Points[Faces[NumFace].idCell[2]]);
	c = MakeFace(Points[Faces[NumFace].idCell[1]], Points[Faces[NumFace].idCell[2]]);
	Type sum = a + b + c;

	//центр грани
	x = (Points[Faces[NumFace].idCell[0]].Point[0] * c + Points[Faces[NumFace].idCell[1]].Point[0] * b + Points[Faces[NumFace].idCell[2]].Point[0] * a)/sum;
	z = (Points[Faces[NumFace].idCell[0]].Point[2] * c + Points[Faces[NumFace].idCell[1]].Point[2] * b + Points[Faces[NumFace].idCell[2]].Point[2] * a)/sum;
	y = (Points[Faces[NumFace].idCell[0]].Point[1] * c + Points[Faces[NumFace].idCell[1]].Point[1] * b + Points[Faces[NumFace].idCell[2]].Point[1] * a)/sum;

	Type sum2 = sqrt(x * x + y * y + z * z);

	newPoint.SetData(x / sum2, y / sum2, z / sum2);
	Points.push_back(newPoint);

	size_t newid = Points.size() -1;
	MyFace<Type> curFace;
	curFace.SetData(Faces[NumFace]);

	MyFace<Type> face;
	face.SetData(curFace.idCell[0], curFace.idCell[1],  newid);
	Faces.push_back(face);

	face.ReSetData(curFace.idCell[0], curFace.idCell[2], newid);
	Faces.push_back(face);

	Faces[NumFace].ReSetData(curFace.idCell[1], curFace.idCell[2], newid);


	newCell.SetData(curFace.idCell[0], curFace.idCell[1], curFace.idCell[2], newid);
	Cells.push_back(newCell);

	return 0;
}

size_t HardAddFace(const size_t NumFace, std::vector<MyPoint<Type>>& Points, std::vector<MyCell<Type>>& Cells,
	std::vector<MyFace<Type>>& Faces) {
	Type* n = new Type[3];
	{
		Type P0[3] = { Points[Faces[NumFace].idCell[0]].Point[0],Points[Faces[NumFace].idCell[0]].Point[1] ,Points[Faces[NumFace].idCell[0]].Point[2] };
		Type P1[3] = { Points[Faces[NumFace].idCell[1]].Point[0],Points[Faces[NumFace].idCell[1]].Point[1] ,Points[Faces[NumFace].idCell[1]].Point[2] };
		Type P2[3] = { Points[Faces[NumFace].idCell[2]].Point[0],Points[Faces[NumFace].idCell[2]].Point[1] ,Points[Faces[NumFace].idCell[2]].Point[2] };

		Type a[3], b[3];
		for (size_t i = 0; i < 3; i++) {
			a[i] = P1[i] - P0[i];
			b[i] = P2[i] - P0[i];
		}
		n[0] = a[1] * b[2] - a[2] * b[1];
		n[1] = -a[0] * b[2] + a[2] * b[0];
		n[2] = a[0] * b[1] - a[1] * b[0];

		Normalize(n);

		size_t id;
		for (size_t i = 0; i < 4; i++) {
			int count = 0;
			for (size_t j = 0; j < 3; j++)
				if (Cells[Faces[NumFace].idTetra].idCell[i] != Faces[NumFace].idCell[j])
					count++;
			if (count == 3) {
				id = i;
				break;
			}

		}

		Type sum = 0;
		Type P3[3] = { Points[id].Point[0],Points[id].Point[1] ,Points[id].Point[2] };

		sum = P1[0] * (P2[1] - P3[1]) * P0[2] + P0[0] * (P3[1] - P2[1]) * P1[2] +
			P0[0] * (P1[1] - P3[1]) * P2[2] + P2[2] * (P1[0] * P3[1] - P1[0] * P0[1]) +
			P3[0] * (P0[2] * (P1[1] - P2[1]) + P1[2] * (P2[1] - P0[1]) + P2[2] * (P0[1] - P1[1]))
			+ P3[2] * (P1[0] * (P0[1] - P2[1]) + P0[0] * (P2[1] - P1[1])) +
			P2[0] * (P0[2] * (P3[1] - P1[1]) + P1[2] * (P0[1] - P3[1]) + P3[2] * (P1[1] - P0[1]));

		if (sum < 0)
			for (size_t i = 0; i < 3; i++)
				n[i] *= -1;
	}

	MyPoint<Type> newPoint;
	MyCell<Type> newCell;
	Type a, b, c, d;
	Type x, y, z;
	a = MakeFace(Points[Faces[NumFace].idCell[0]], Points[Faces[NumFace].idCell[1]]);
	b = MakeFace(Points[Faces[NumFace].idCell[0]], Points[Faces[NumFace].idCell[2]]);
	c = MakeFace(Points[Faces[NumFace].idCell[1]], Points[Faces[NumFace].idCell[2]]);
	Type sum = a + b + c;

	//центр грани
	x = (Points[Faces[NumFace].idCell[0]].Point[0] * c + Points[Faces[NumFace].idCell[1]].Point[0] * b + Points[Faces[NumFace].idCell[2]].Point[0] * a) / sum;
	z = (Points[Faces[NumFace].idCell[0]].Point[2] * c + Points[Faces[NumFace].idCell[1]].Point[2] * b + Points[Faces[NumFace].idCell[2]].Point[2] * a) / sum;
	y = (Points[Faces[NumFace].idCell[0]].Point[1] * c + Points[Faces[NumFace].idCell[1]].Point[1] * b + Points[Faces[NumFace].idCell[2]].Point[1] * a) / sum;


	Type sum2 = n[0] * x + n[1] * y + n[2] * z;
	Type t = -sum2 + sqrt(sum2 * sum2 - (x * x + y * y + z * z - 1));

	newPoint.SetData(n[0] * t + x, n[1] * t + y, n[2] * t + z);
	Points.push_back(newPoint);


	size_t newid = Points.size() - 1;
	MyFace<Type> curFace;
	curFace.SetData(Faces[NumFace]);

	newCell.SetData(curFace.idCell[0], curFace.idCell[1], curFace.idCell[2], newid);
	Cells.push_back(newCell);
	size_t newCellid = Cells.size() - 1;


	MyFace<Type> face;
	face.SetData(curFace.idCell[0], curFace.idCell[1], newid, newCellid);
	Faces.push_back(face);

	face.ReSetData(curFace.idCell[0], curFace.idCell[2], newid, newCellid);
	Faces.push_back(face);

	Faces[NumFace].ReSetData(curFace.idCell[1], curFace.idCell[2], newid, newCellid);


	
	delete[] n;
	return 0;
}



size_t MyMakeGrid3D() {
	std::vector<MyPoint<Type>> Points;
	std::vector<MyCell<Type>> Cells;
	std::vector<MyFace<Type>> Faces;
	std::string FileVTK = "C:\\Users\\Artem\\Desktop\\TetraGrid\\mesh2\\Sphere1.vtk";
	FromVTKToMySet(FileVTK, Points, Cells);

	int count1 = 0;
	while (count1++ < 4) {
		size_t size = Cells.size();
		for (size_t i = 0; i < size; i++) {
			if(sqrt(pow(Points[Cells[i].idCell[0]].Point[0],2)+ 
				pow(Points[Cells[i].idCell[0]].Point[1], 2)+
				pow(Points[Cells[i].idCell[0]].Point[2], 2))<0.5)
			FracTetra(i, Points, Cells);

		}
	}



	OutGrid("C:\\Users\\Artem\\Desktop\\TetraGrid\\MySet\\Sphere0.vtk", Points, Cells);

	return 1;
	MyPoint<Type> A, B, C, D;
	//A.SetData(0, 0, 0);
	//B.SetData(0, 0, 1);
	//C.SetData(0, 1, 0);
	//D.SetData(1, 0, 0);
	A.SetData(sqrt(8./9), 0, -1./3);
	B.SetData(-sqrt(2./9), sqrt(2./3), -1./3);
	C.SetData(-sqrt(2. / 9), -sqrt(2. / 3), -1. / 3);
	D.SetData(0, 0, 1);
	
	Points.push_back(A);
	Points.push_back(B);
	Points.push_back(C);
	Points.push_back(D);

	MyCell<Type> cell;
	cell.SetData(0, 1, 2, 3);
	Cells.push_back(cell);

	MyFace<Type> face;
	face.SetData(0, 1, 2,0);
	Faces.push_back(face);
	face.SetData(0, 1, 3,0);
	Faces.push_back(face);
	face.SetData(1, 2, 3,0);
	Faces.push_back(face);
	face.SetData(0, 2, 3,0);
	Faces.push_back(face);

	std::vector<MyPoint<Type>> bufPoints;
	std::vector<MyCell<Type>> bufCells;
	std::vector<MyFace<Type>> bufFaces;
	
	int count = 0;
	while (count++ < 2) {
		size_t size = Faces.size();
		for (size_t i = 0; i < size; i++)
			HardAddFace(i, Points, Cells, Faces);
	}

	/*int count = 0;
	while (count++ < 2) {
		size_t size = Cells.size();
		for (size_t i = 0; i < size; i++)
	     	FracTetra(i, Points, Cells);
	}*/

	

	OutGrid("C:\\Users\\Artem\\Desktop\\MyGrid.vtk", "C:\\Users\\Artem\\Desktop\\\Faces.txt",Points, Cells,Faces);

		
	return 0;
}





