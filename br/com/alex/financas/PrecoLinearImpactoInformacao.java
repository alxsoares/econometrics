package br.com.alex.financas;

import java.util.Random;

import Jama.LUDecomposition;
import Jama.Matrix;

/**
 * 
 * Implementação baseada no paper Optimal control of execution costs de Dimitris
 * Bertsimas, Andrew W. Lo
 * 
 * Modelo da seção 2.4: Linear price impact with information
 * 
 * @author alex
 * 
 */
public class PrecoLinearImpactoInformacao {
	/**
	 * Coeficiente de impacto no preço devido à uma compra de S contratos.
	 * Impacto = \theta*S
	 */
	private final static double theta = 5 * 10E-5;

	/**
	 * Medida da sensibilidade do modelo a mudanças de mercado.
	 */
	private final static double gamma = 5D;

	/**
	 * Coeficiente de ajuste do vetor de informações.
	 */
	private final static double rho = 0.5D;

	/**
	 * desvio padrão da sequência aleatória \epsilon
	 */
	private final static double stdDevEps = 0.125D;

	/**
	 * desvio padrão da sequência aleatória \eta
	 */
	private final static double stdDevEta = 0.001D;

	private static final double meanEps = 0;

	private static final double meanEta = 0;

	private static Random epsDistrib;

	private static Random etaDistrib;

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

	/**
	 * Gero os números pseudo-aleatórios que representam um shock aleatório no
	 * preço A sequência tem média 0 e desvio padrão stdDevEps definidos acima
	 * 
	 * @return
	 */
	private static double getEpsilon() {
		if (epsDistrib == null) {
			epsDistrib = new Random(System.currentTimeMillis());
		}
		return epsDistrib.nextGaussian() * stdDevEps + meanEps;
	}

	/**
	 * Gero os números pseudo-aleatórios que representam um shock aleatório no
	 * vetor X de informações de mercado. A sequência tem média 0 e desvio
	 * padrão stdDevNi definidos acima.
	 * 
	 * @return
	 */
	private static double getEta() {
		if (etaDistrib == null) {
			etaDistrib = new Random(System.currentTimeMillis());
		}
		return etaDistrib.nextGaussian() * stdDevEta + meanEta;
	}

	static Matrix buildEpsilonVector(int dim) {
		double arr[][] = new double[dim][1];
		for (int i = 0; i < dim; i++) {
			arr[i][0] = getEpsilon();
		}
		Matrix eps = new Matrix(arr);

		return eps;
	}

	static Matrix buildEtaVector(int dim) {
		double arr[][] = new double[dim][1];
		for (int i = 0; i < dim; i++) {
			arr[i][0] = getEta();
		}
		Matrix ni = new Matrix(arr);

		return ni;
	}

	static Matrix buildIndentyMinuxShiftInverse(int dim) {
		Matrix identity = buildIdentityMatrix(dim);
		Matrix shift = buildShiftMatrix(dim);
		Matrix m = identity.minus(shift).inverse();
		return m;
	}

	/**
	 * A = ((Identidade - Shift)^-1)*(\theta*Shift)
	 * 
	 * @param dim
	 * @return
	 */
	static Matrix buildMatrixA(int dim) {
		Matrix shift = buildShiftMatrix(dim);
		Matrix inverse = buildIndentyMinuxShiftInverse(dim);
		Matrix r = inverse.times(theta).times(shift);
		return r;
	}

	/**
	 * P1 = p1*(1,0,0....)
	 * 
	 * @param dim
	 * @param p1
	 * @return
	 */
	static Matrix buildOneElementVector(int dim, double p1) {
		Matrix column = new Matrix(dim, 1);
		column.set(0, 0, p1);
		return column;
	}

	/**
	 * Matrix B = (\gamma*X +\epsilon +(1 -Shift)^-1*p1)^T
	 * 
	 * @param dim
	 * @param x
	 * @param eps
	 * @param p1
	 * @return
	 */
	static Matrix buildMatrixB(int dim, Matrix x, Matrix eps, Matrix p1) {
		Matrix b1 = x.times(gamma).plus(eps);
		Matrix b2 = buildIndentyMinuxShiftInverse(dim).times(p1);
		return b1.plus(b2);
	}

	/**
	 * vetor de informações \gamma*X + \eta 
	 * @param dim
	 * @param x
	 * @param eta
	 * @return
	 */
	static Matrix autoregressiveX(int dim, Matrix x, Matrix eta) {
		Matrix newX = new Matrix(dim, 1);
		newX.set(0, 0, x.get(0, 0) * gamma);
		return newX.plus(eta);
	}

	/**
	 * Calcula o vetor de preços.
	 * @param p1
	 * @param q2
	 * @param x2
	 * @param eps2
	 * @return
	 */
	static Matrix autoregressivePrice(Matrix p1, Matrix q2, Matrix x2, Matrix eps2)
	{
		Matrix thetaQ2 = q2.times(theta);
		Matrix gammaX2 = x2.times(gamma);
		Matrix p2 = p1.plus(thetaQ2).plus(gammaX2).plus(eps2);		
		return p2;
	}
	
	
	 
	/**
	 * @param args
	 */
	public static void main(String[] args) {

		long Q = 100000;
		int dim = 20;
		double P1 = 5000000;
		Matrix epsilon = buildEpsilonVector(dim);
		Matrix eta = buildEtaVector(dim);
		double x0 = 0;
		double cost = 0;
		for (int i = 0; i < dim; i++) {
			Matrix epsi = buildOneElementVector(dim-i, epsilon.get(i, 0));
			Matrix etai = buildOneElementVector(dim -i, eta.get(i, 0));
			Matrix a = buildMatrixA(dim-i);
			Matrix x = autoregressiveX(dim -i, buildOneElementVector(dim-i, x0), etai);
			Matrix b = buildMatrixB(dim-i, x, epsi, buildOneElementVector(dim-i, P1));			
			Matrix aat = a.plus(a.transpose());
			Matrix oneColumm = buildIdentityColumnVector(dim-i);
			Matrix oneRow = oneColumm.transpose();
			/**
			 * Dx = E
			 */
			//D
			Matrix d = new Matrix(dim -i +1, dim-i+1);
			d.setMatrix(0, dim - 1-i, 0, dim - 1-i, aat);
			d.setMatrix(0, dim - 1 -i, dim-i, dim-i, oneColumm);
			d.setMatrix(dim-i, dim-i, 0, dim - 1-i, oneRow);
			d.set(dim-i, dim-i, 0);
            
			//E
			Matrix e = new Matrix(dim + 1-i, 1);
			e.setMatrix(0, dim - 1-i, 0, 0, b.uminus());
			e.set(dim-i, 0, Q);
			//Resolvo sistema
			
			LUDecomposition lu = new LUDecomposition(d);
			Matrix result = lu.solve(e);
			//Contratos
			Matrix q = result.getMatrix(0, dim - 1-i, 0, 0);
			Matrix p = autoregressivePrice(buildOneElementVector(dim-i, P1), q, x, epsi);
			double contracts = q.get(0, 0);			
			P1 = p.get(0, 0);			
			Q = Q - Math.round(contracts);	
			System.out.println("Preco : "+ P1);
			System.out.println("Comprados: "+ contracts);
			cost+=P1*contracts;
			
		}
		System.out.println("Resto : "+ Q);
		System.out.println("Custo :" + cost);

	}

}
