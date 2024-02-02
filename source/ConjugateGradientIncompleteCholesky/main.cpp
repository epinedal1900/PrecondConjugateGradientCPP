#include <Basic/Log.h>
#include <Basic/Float.h>
#include <Container/CSRMatrix.h>
#include <Container/Vector.h>
#include <Solver/TriangularSystem.h>


template <typename T>
int ConjugateGradientIncompleteCholesky(const CSRMatrix<T>& A, Vector<T>& x, const Vector<T>& b, T tolerance, int max_steps, const CSRMatrix<T>& L, const CSRMatrix<T>& Lt, int threads) noexcept 
{
	try
	{
		Log(1, "ConjugateGradientIncompleteCholesky:");
		Log(1, "-Tolerance:   %.5e", tolerance);
		Log(1, "-MaxSteps:    %i", max_steps);

		int n = x.size;

		Vector<T> g(n); // Gradient
		Vector<T> p(n); // Descent direcction
		Vector<T> w(n); // w = A*p
		Vector<T> q(n); // Preconditioner system solution
		Vector<T> r(n); // Preconditioner system intermediate solution

		T gg = 0.0;
		#pragma omp parallel for default(shared) reduction(+:gg) schedule(guided, 64) num_threads(threads)
		for (int i = 1; i <= n; ++i)
		{
			int* __restrict A_index_i = A.index[i];
			T* __restrict A_entry_i = A.entry[i];

			T sum = 0.0;
			int k_max = A.Count(i);
			for (register int k = 1; k <= k_max; ++k)
			{
				sum += A_entry_i[k]*x.entry[A_index_i[k]];
			}
			g.entry[i] = sum - b.entry[i]; // g = AX - b
			gg += g.entry[i]*g.entry[i]; // gg = g'*g
		}

		// Solve for q: MQ = g
		LowerTriangularSystem(L, r, g);
		UpperTriangularSystem(Lt, q, r);

		T gq = 0.0;
		#pragma omp parallel for default(shared) reduction(+:gq) schedule(guided, 64) num_threads(threads)
		for (int i = 1; i <= n; ++i)
		{
			p.entry[i] = -q.entry[i]; // p = -q
			gq += g.entry[i]*q.entry[i]; // gq = g'*q
		}

		Log(2, "-Step  r'*r");
		T epsilon = tolerance*tolerance;
		int step = 0;
		while (step < max_steps)
		{
			// Test for invalid value
			if (Float<T>::IsNaN(gg))
			{
				Log(1, "-[Error] ConjugateGradientIncompleteCholesky got NaN.");
				return -1;
			};

			// Test termination condition
			if (gg <= epsilon) // Norm(Gn) <= tolerance
			{
				break;
			}

			T pw = 0.0;
			#pragma omp parallel for default(shared) reduction(+:pw) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				int* __restrict A_index_i = A.index[i];
				T* __restrict A_entry_i = A.entry[i];

				T sum = 0.0;
				int k_max = A.Count(i);
				for (register int k = 1; k <= k_max; ++k)
				{
					sum += A_entry_i[k]*p.entry[A_index_i[k]];
				}
				w.entry[i] = sum; // w = AP
				pw += p.entry[i]*w.entry[i]; // pw = p'*w
			}

			T alpha = gq/pw; // alpha = (g'*q)/(p'*w)

			T gngn = 0.0;
			#pragma omp parallel for default(shared) reduction(+:gngn) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				x.entry[i] += alpha*p.entry[i]; // Xn = x + alpha*p
				g.entry[i] += alpha*w.entry[i]; // Gn = g + alpha*w
				gngn += g.entry[i]*g.entry[i]; // gngn = Gn'*Gn
			}

			// Solve for q: MQ = g
			LowerTriangularSystem(L, r, g);
			UpperTriangularSystem(Lt, q, r);

			T gnqn = 0.0;
			#pragma omp parallel for default(shared) reduction(+:gnqn) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				gnqn += g.entry[i]*q.entry[i]; // gnqn = Gn'*Qn
			}

			T beta = gnqn/gq; // beta = (Gn'*Gn)/(g'*g)

			#pragma omp parallel for default(shared) schedule(guided, 64) num_threads(threads)
			for (int i = 1; i <= n; ++i)
			{
				p.entry[i] = beta*p.entry[i] - q.entry[i]; // Pn = -q + beta*p
			}

			gg = gngn;
			gq = gnqn;

			Log(2, "%5i  %.5e", step, gg);
			++step;
		}
		Log(1, "-Total steps: %i", step);

		if (step >= max_steps)
		{
			Log(1, "-[Error] ConjugateGradientIncompleteCholesky did not converge in %i steps", max_steps);
			return -1;
		}

		return step;
	}
	catch (Exception&)
	{
		ReThrow();
	}
}
