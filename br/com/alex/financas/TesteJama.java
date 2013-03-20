package br.com.alex.financas;

import Jama.LUDecomposition;
import Jama.Matrix;

public class TesteJama {
	private static Matrix buildShiftMatrix(int dim) {
		double s[][] = new double[dim][dim];
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				s[i][j] = 0;
			}
		}
		for (int i = 0; i < dim - 1; i++) {
			s[i + 1][i] = 1;
		}
		Matrix S = new Matrix(s);
		return S;
	}

	static Matrix buildIdentityMatrix(int dim) {
		double arrI[][] = new double[dim][dim];
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				arrI[i][j] = 0;
				if (i == j) {
					arrI[i][j] = 1;
				}
			}
		}
		Matrix mI = new Matrix(arrI);
		return mI;
	}

	static Matrix buildIdentityColumnVector(int dim) {
		double arrI[][] = new double[dim][1];
		for (int i = 0; i < dim; i++) {
			arrI[i][0] = 1;
		}
		Matrix mI = new Matrix(arrI);
		return mI;
	}

	static Matrix buildIdentityLineVector(int dim) {
		return buildIdentityColumnVector(dim).transpose();

	}

	static Matrix buildEpsilonVector(int dim) {
		double arr[][] = new double[dim][1];
		for (int i = 0; i < dim; i++) {
			arr[i][0] = Math.random() - 0.5;
		}
		Matrix eps = new Matrix(arr);

		return eps;
	}

	public static void main(String[] args) {
		//TODO: acertar ou estimar alfa e beta.
		double alfa = 0.5;
		double delta = 1;
		long Q = 1000;
		int dim = 100;
		double P1 = 500000;
		
		for (int i = 0; i < dim && Q > 0; i++) {
			Matrix E = buildEpsilonVector(dim-i);
			Matrix S = buildShiftMatrix(dim-i);
//			Matrix I = buildIdentityMatrix(dim-i);
			Matrix E1 = new Matrix(dim-i, 1);
			E1.set(0, 0, E.get(0, 0));
			Matrix A = (buildIdentityMatrix(dim-i).minus(S.times(alfa))).times(
					delta).inverse()
					.plus(buildIdentityMatrix(dim-i).times(delta));
			Matrix B = (buildIdentityMatrix(dim-i).minus(S.times(alfa)))
					.times(E1).times(P1);

			Matrix gradF = A.plus(A.transpose());
//			System.out.println("A+At:");
//			gradF.print(10, 10);
//			System.out.println();
			Matrix L = new Matrix(dim + 1-i, dim + 1-i);
			L.setMatrix(0, dim - 1-i, 0, dim - 1-i, gradF);
			L.setMatrix(0, dim - 1 -i, dim-i, dim-i, buildIdentityColumnVector(dim-i));
			L.setMatrix(dim-i, dim-i, 0, dim - 1-i, buildIdentityLineVector(dim-i));
			L.set(dim-i, dim-i, 0);
//			System.out.println("L:");
//			L.print(10, 10);

			Matrix R = new Matrix(dim + 1-i, 1);
			R.setMatrix(0, dim - 1-i, 0, 0, B.uminus());
			R.set(dim-i, 0, Q);

//			System.out.println("A:");
//			A.print(10, 10);
//			System.out.println("B:");
//			B.print(10, 10);
//			System.out.println("S:");
//			S.print(10, 10);
//			System.out.println("E:");
//			E.print(10, 10);
//			System.out.println("E1:");
//			E1.print(10, 10);
//			System.out.println("R:");
//			R.print(10, 10);

			LUDecomposition lu = new LUDecomposition(L);
//			System.out.println("Resultado:");
			Matrix result = lu.solve(R);
//			result.print(10, 10);
			Matrix qV = result.getMatrix(0, dim - 1-i, 0, 0);
//			System.out.println("Contratos");
//			qV.print(10, 10);			
			Matrix P = buildIdentityMatrix(dim-i).minus(S.times(alfa)).inverse()
					.times(S.times(qV).times(alfa).times(delta).plus(E)).plus(
							E1.times(P1));
//			System.out.println("Precos");
//			P.print(10, 10);
			P1 = P.get(0,0);
			long contratos = Math.round(qV.get(0, 0));
			Q = Q - contratos;
			System.out.println("Contratos :"+contratos);
			System.out.println("Preço :"+P1);
		}
		System.out.println("Compra o restante: " +Q);
	}

}
