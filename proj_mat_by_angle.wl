(* ::Package:: *)

indexToRep[n_]:= indexToRep[n] = Module[{j, m},
	j = Floor[Sqrt[n - 1]];
	m = n - 1 - j^2 - j;
	{j, m}
]

repToIndex[j_, m_]:= repToIndex[j, m] = j^2 + m + j + 1;

indexToRepSortedBym[n_]:= 
	indexToRepSortedBym[n] = Sort[indexToRep[n], #1[[2]] < #2[[2]] &]
	
approximateByProjection[mat_, n_, a_, csvFilePath_, matFilePath_]:= Block[
	{u, w, v, res, index, fileKey},
	{u, w, v} = SingularValueDecomposition[mat, Min[Dimensions[mat]]];
	res = u . DiagonalMatrix[If[Re[#] > 0.5, 1., 0.] &/@ Diagonal[w]] . ConjugateTranspose[v];
	Export[csvFilePath, res];
	Export[matFilePath, res];
	res
]

symmetrizeLowerDiag[mat_] := Transpose[LowerTriangularize[mat,-1]] + mat

sphNormCoeffs[l_,m_]:=Sqrt[(2*l+1)*(l-m)!/(4*Pi*(l+m)!)]

projInfinitySparseList[r1_?ListQ, r2_?ListQ, angles_]:=Block[
	{\[Theta], \[Phi], 
	j1 = First[r1],
	m1 = Last[r1],
	j2 = First[r2],
	m2 = Last[r2], 
	integral},
	integral = 
		2*Pi*sphNormCoeffs[j2, m2] * sphNormCoeffs[j1, m1]* Integrate[LegendreP[j1, m1, x]LegendreP[j2, m2, x], x];
	Subtract @@ (integral /. x -> {1, Cos[angles]})
]


projSparseAllAngles[n_, a_, dir_]:=
	Block[
	{
		ind = indexToRep/@Range[n^2],
		angles = Table[ArcCos[1 - k / (a - 0.5)], {k, 0, a - 1}],
		indSorted, indexofjm0, indSortedUntiljmAre0, arraysAllAngles, symmetrized,
		indices, fileKeys, csvFilePaths, matFilePaths
	},
	indSorted = Sort[ind, #1[[2]] < #2[[2]]&];
	indexofjm0 = Total[Range[1, n]];
	indSortedUntiljmAre0 = indSorted[[1;;indexofjm0]];
	arraysAllAngles = Table[
		SparseArray[
		{{i_,j_}:>
		Evaluate[projInfinitySparseList[indSortedUntiljmAre0[[i]],indSortedUntiljmAre0[[j]],angles]][[k]]/;
		i >= j &&Last[indSortedUntiljmAre0[[i]]] == Last[indSortedUntiljmAre0[[j]]]},
		{Length[indSortedUntiljmAre0], Length[indSortedUntiljmAre0]}],
		{k, 1, Length[angles]}
	];
	symmetrized = symmetrizeLowerDiag /@ arraysAllAngles;
	indices = Range[a];
	fileKeys = "n" <> ToString[n] <> "_a_" <>ToString[a] <> "_i_" <> ToString[#] &/@indices;
	csvFilePaths = dir <> # <> ".csv" &/@ fileKeys;
	matFilePaths = dir <> # <> ".mat" &/@ fileKeys;
	Print[csvFilePaths];
	(*mat_, n_, a_, filePath_*)
	Return @ MapThread[approximateByProjection[#1, n, a, #2, #3]&, {symmetrized, csvFilePaths, matFilePaths}]
]


projSparseAllAnglesParallel[n_, a_, dir_]:=
	Block[{
		ind = indexToRep/@Range[n^2],
		angles = Table[ArcCos[1 - k / (a - 0.5)], {k, 0, a - 1}],
		indSorted, indexofjm0, indSortedUntiljmAre0, arraysAllAngles, symmetrized,
		indices, fileKeys, csvFilePaths, matFilePaths
	},
	indSorted = Sort[ind,#1[[2]]< #2[[2]]&];
	indexofjm0 = Total[Range[1, n]];
	indSortedUntiljmAre0 = indSorted[[1;;indexofjm0]];
	CloseKernels[]; LaunchKernels[];
	arraysAllAngles = Parallelize @ Table[
		SparseArray[{{i_,j_}:>
		Evaluate[projInfinitySparseList[indSortedUntiljmAre0[[i]],indSortedUntiljmAre0[[j]],angles]][[k]]/;i>=j &&Last[indSortedUntiljmAre0[[i]]]==Last[indSortedUntiljmAre0[[j]]]},{Length[indSortedUntiljmAre0],Length[indSortedUntiljmAre0]}],{k,1,Length[angles]}];
	symmetrized = symmetrizeLowerDiag /@ arraysAllAngles;
	indices = Range[a];
	fileKeys = "n_" <> ToString[n] <> "_a_" <>ToString[a] <> "_i_" <> ToString[#] &/@indices;
	csvFilePaths = dir <> # <> ".csv" &/@ fileKeys;
	matFilePaths = dir <> # <> ".mat" &/@ fileKeys;
	(*mat_, n_, a_, filePath_*)
	Parallelize @ MapThread[approximateByProjection[#1, n, a, #2, #3]&, {symmetrized, csvFilePaths, matFilePaths}]
]


(* ::Section:: *)
(*Test*)


n = 10;
a = 5;
dir = NotebookDirectory[] <> "data_test/"

test = AbsoluteTiming[projSparseAllAnglesParallelizeOnTableParallelizeMapThread[n, a, dir]][[1]]



