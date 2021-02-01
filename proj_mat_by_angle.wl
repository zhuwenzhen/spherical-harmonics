(* ::Package:: *)

IndexToRep[n_]:= IndexToRep[n] = Module[{j, m},
	j = Floor[Sqrt[n - 1]];
	m = n - 1 - j^2 - j;
	{j, m}
]

RepToIndex[j_, m_]:= RepToIndex[j, m] = j^2 + m + j + 1;

IndexToRepSortedBym[n_]:= 
	IndexToRepSortedBym[n] = Sort[IndexToRep[n], #1[[2]] < #2[[2]] &]
	
ApproximateByProjection[mat_]:= Module[{u,w,v},
	{u, w, v} = SingularValueDecomposition[mat];
	Return[u . DiagonalMatrix[If[Re[#] > 0.5, 1., 0.] &/@ Diagonal[w]] . ConjugateTranspose[v]];
]

ProjInfinity[r1_, r2_, angle_?NumericQ]:= Block[
	{\[Theta], \[Phi]},
	If[
		Last[r1]!=Last[r2],
		0.,
		NIntegrate[
			Sin[\[Theta]]SphericalHarmonicY[First[r2],Last[r2], \[Theta], \[Phi]]SphericalHarmonicY[First[r1],Last[r1], \[Theta], -\[Phi]],
			{\[Theta],0,angle},{\[Phi],0,2Pi}, 
			Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
			WorkingPrecision -> 30
		]
	]
]

ProjInfinityExact[r1_, r2_, angle_?NumericQ]:= Block[
	{\[Theta], \[Phi]},
	If[
		Last[r1]!=Last[r2],
		0.,
		FullSimplify[
			Integrate[
				Sin[\[Theta]]SphericalHarmonicY[First[r2],Last[r2], \[Theta], \[Phi]]SphericalHarmonicY[First[r1],Last[r1], \[Theta], -\[Phi]],
				{\[Theta], 0, angle}, {\[Phi], 0, 2Pi}
			]
		]
	]
]

ProjInfinityGKrule[r1_, r2_, angle_?NumericQ]:= Block[
	{\[Theta], \[Phi]},
	If[
		Last[r1]!=Last[r2],
		0.,
		NIntegrate[
			Sin[\[Theta]]SphericalHarmonicY[First[r2],Last[r2],\[Theta],\[Phi]]SphericalHarmonicY[First[r1],Last[r1],\[Theta],-\[Phi]],
			{\[Theta], 0, angle},{\[Phi], 0, 2Pi}, 
			Method -> {"GaussKronrodRule", "SymbolicProcessing" -> 0},
			WorkingPrecision -> 30]
	]
]

ProjInfinityC[r1_, r2_, angle_?NumericQ]:= Block[
	{\[Theta], \[Phi]},
	If[
		Last[r1]!=Last[r2],
		0.0,
		NIntegrate[
			Sin[\[Theta]]SphericalHarmonicY[First[r2],Last[r2],\[Theta],\[Phi]]SphericalHarmonicY[First[r1],Last[r1],\[Theta],-\[Phi]],
			{\[Theta], 0, angle}, {\[Phi], 0, 2\[Pi]}, 
			Method -> {"GlobalAdaptive", Method -> "GaussKronrodRule", "SingularityDepth" -> Infinity},
			WorkingPrecision -> 30
		]
	]
]


ProjSorted[jmax_,angle_]:= Block[
	{ind = IndexToRep/@ Range[(jmax+1)^2], indSorted, overlap, M},
	overlap = Function[{r1,r2}, ProjInfinity[r1, r2, angle]];
	indSorted = Sort[ind, #1[[2]]<#2[[2]]&];
	M = Outer[overlap, indSorted, indSorted, 1];
	Chop[Re[ApproximateByProjection[M]]]
]

ProjSortedShort[jmax_,angle_]:= Block[
	{ind=IndexToRep/@Range[(jmax+1.)^2], indSorted, indexofjm0, indSortedUntiljmAre0, overlap, M},
	overlap = Function[{r1,r2}, ProjInfinityGKrule[r1, r2, angle]];
	indSorted = Sort[ind, #1[[2]]<#2[[2]]&];
	indexofjm0 = Total[Range[1,jmax+1]];
	indSortedUntiljmAre0 = indSorted[[1;;indexofjm0]];
	M = Outer[overlap, indSortedUntiljmAre0, indSortedUntiljmAre0, 1];
	Chop[Re[ApproximateByProjection[M]]]
]

ProjSortedShortParallel[jmax_,angle_]:= Block[
	{ind=IndexToRep/@Range[(jmax+1.)^2], indSorted, indexofjm0, indSortedUntiljmAre0, overlap, M},
	overlap = Function[{r1,r2}, ProjInfinityGKrule[r1, r2, angle]];
	indSorted = Sort[ind, #1[[2]]<#2[[2]]&];
	indexofjm0 = Total[Range[1,jmax+1]];
	indSortedUntiljmAre0 = indSorted[[1;;indexofjm0]];
	M = Parallelize @ Outer[overlap, indSortedUntiljmAre0, indSortedUntiljmAre0, 1];
	Chop[Re[ApproximateByProjection[M]]]
]
