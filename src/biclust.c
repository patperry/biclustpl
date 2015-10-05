//  Copyright 2015 Cheryl J. Flynn and Patrick O. Perry.
//
//  Licensed under the Apache License, Version 2.0 (the "License");
//  you may not use this file except in compliance with the License.
//  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.

#include <Rdefines.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>


static double entropy_gaussian(double size, double sum)
{
	double mean;

	if (size <= 0)
		return 0;

	mean = sum / size;
	return 0.5 * sum * mean;
}


static double entropy_binomial(double size, double sum)
{
	double csum, p, q;

	if (size <= 0)
		return 0;

	if (sum <= 0 || sum >= size)
		return 0;

	csum = size - sum;
	p = sum / size;
	q = csum / size;
	return sum * log(p) + csum * log(q);
}


static double entropy_poisson(double size, double sum)
{
	double mean;

	if (size <= 0)
		return 0;

	if (sum <= 0)
		return 0;

	mean = sum / size;
	return sum * (log(mean) - 1);
}


struct ploglik {
	int m, n;
	int K, L;

	const double *x;

	int *row_cl;
	int *col_cl;

	double obj_tot;
	double *obj;
	double *sum;
	double *row_sum;
	double *col_sum;

	double *size;
	double *row_size;
	double *col_size;

	double (*entropy)(double, double);
};


void ploglik_init(struct ploglik *pl, int m, int n, int K, int L,
		  const double *x, double (*entropy)(double, double))
{
	pl->x = x;
	pl->m = m;
	pl->n = n;
	pl->K = K;
	pl->L = L;
	pl->row_cl = (int *)R_alloc(m, sizeof(int));
	pl->col_cl = (int *)R_alloc(n, sizeof(int));
	pl->sum = (double *)R_alloc(K * L, sizeof(double));
	pl->obj = (double *)R_alloc(K * L, sizeof(double));
	pl->row_sum = (double *)R_alloc(m * L, sizeof(double));
	pl->col_sum = (double *)R_alloc(n * K, sizeof(double));
	pl->size = (double *)R_alloc(K * L, sizeof(double));
	pl->row_size = (double *)R_alloc(m * L, sizeof(double));
	pl->col_size = (double *)R_alloc(n * K, sizeof(double));
	pl->entropy = entropy;
}


void ploglik_assign(struct ploglik *pl, const struct ploglik *src)
{
	const int m = pl->m, n = pl->n, K = pl->K, L = pl->L;

	memcpy(pl->row_cl, src->row_cl, m * sizeof(int));
	memcpy(pl->col_cl, src->col_cl, n * sizeof(int));
	memcpy(pl->sum, src->sum, K * L * sizeof(double));
	memcpy(pl->obj, src->obj, K * L * sizeof(double));
	memcpy(pl->row_sum, src->row_sum, m * L * sizeof(double));
	memcpy(pl->col_sum, src->col_sum, n * K * sizeof(double));
	memcpy(pl->size, src->size, K * L * sizeof(double));
	memcpy(pl->row_size, src->row_size, m * L * sizeof(double));
	memcpy(pl->col_size, src->col_size, n * K * sizeof(double));
	pl->obj_tot = src->obj_tot;
}


double ploglik_eval(struct ploglik *pl)
{
	const int m = pl->m, n = pl->n, K = pl->K, L = pl->L;
	const double *x = pl->x;
	double obj_kl, n_kl, s_kl, x_ij;
	int i, j, k, l;

	// initialize sums and sizes to 0
	memset(pl->sum, 0, K * L * sizeof(double));
	memset(pl->size, 0, K * L * sizeof(double));
	memset(pl->row_sum, 0, m * L * sizeof(double));
	memset(pl->row_size, 0, m * L * sizeof(double));
	memset(pl->col_sum, 0, n * K * sizeof(double));
	memset(pl->col_size, 0, n * K * sizeof(double));

	// compute block sums and sizes
	for (j = 0; j < n; j++) {
		l = pl->col_cl[j];
		for (i = 0; i < m; i++) {
			x_ij = x[i + j * m];
			if (R_finite(x_ij)) {
				k = pl->row_cl[i];
				pl->sum[k + l * K] += x_ij;
				pl->size[k + l * K] += 1;
				pl->row_sum[i + l * m] += x_ij;
				pl->row_size[i + l * m] += 1;
				pl->col_sum[j + k * n] += x_ij;
				pl->col_size[j + k * n] += 1;
			}
		}
	}

	// compute objective function
	pl->obj_tot = 0;
	for (l = 0; l < L; l++) {
		for (k = 0; k < K; k++) {
			n_kl = pl->size[k + l * K];
			s_kl = pl->sum[k + l * K];
			obj_kl = pl->entropy(n_kl, s_kl);

			pl->obj[k + l * K] = obj_kl;
			pl->obj_tot += obj_kl;
		}
	}

	return pl->obj_tot;
}


enum update_type {
	NO_UPDATE,
	UPDATE
};

// compute how much the objective changes from adding weight w to
// row i in class k; (use weight = -1 to move i out of class k;
// use weight = +1 to move i into class k)
double ploglik_add_row(struct ploglik *pl, enum update_type type,
		       int i, int k, double w)
{
	const int m = pl->m, n = pl->n, K = pl->K, L = pl->L;
	const double *x = pl->x;
	double n_l, n_kl_new, obj_delta, obj_kl_old, obj_kl_new,
	       s_l, s_kl_new, x_ij;
	int j, l;

	// compute how much objective changes from moving i out of class k
	obj_delta = 0.0;
	for (l = 0; l < L; l++) {
		obj_kl_old = pl->obj[k + l * K];
		n_l = pl->row_size[i + l * m];
		s_l = pl->row_sum[i + l * m];
		n_kl_new = pl->size[k + l * K] + w * n_l;
		s_kl_new = pl->sum[k + l * K] + w * s_l;
		obj_kl_new = pl->entropy(n_kl_new, s_kl_new);

		// update block sums
		if (type == UPDATE) {
			pl->size[k + l * K] = n_kl_new;
			pl->sum[k + l * K] = s_kl_new;
			pl->obj[k + l * K] = obj_kl_new;
		}

		obj_delta += obj_kl_new - obj_kl_old;

	}

	//Rprintf("[obj_old = %lg]", pl->obj_tot);
	//Rprintf("[obj_delta = %lg]", obj_delta);

	// update column sums and objective function
	if (type == UPDATE) {
		for (j = 0; j < n; j++) {
			x_ij = x[i + j * m];

			if (R_finite(x_ij)) {
				pl->col_sum[j + k * n] += w * x_ij;
				pl->col_size[j + k * n] += w;
			}
		}

		pl->obj_tot += obj_delta;
	}

	//Rprintf("[obj_tot_new = %lg]", pl->obj_tot);

	return obj_delta;
}


double ploglik_add_col(struct ploglik *pl, enum update_type type,
		       int j, int l, double w)
{
	const int m = pl->m, n = pl->n, K = pl->K;
	const double *x = pl->x;
	double n_k, n_kl_new, obj_delta, obj_kl_old, obj_kl_new,
	       s_k, s_kl_new, x_ij;
	int i, k;

	// compute how much objective changes from moving i out of class l
	obj_delta = 0.0;
	for (k = 0; k < K; k++) {
		obj_kl_old = pl->obj[k + l * K];
		n_k = pl->col_size[j + k * n];
		s_k = pl->col_sum[j + k * n];
		n_kl_new = pl->size[k + l * K] + w * n_k;
		s_kl_new = pl->sum[k + l * K] + w * s_k;
		obj_kl_new = pl->entropy(n_kl_new, s_kl_new);

		// update block sums
		if (type == UPDATE) {
			pl->size[k + l * K] = n_kl_new;
			pl->sum[k + l * K] = s_kl_new;
			pl->obj[k + l * K] = obj_kl_new;
		}

		obj_delta += obj_kl_new - obj_kl_old;

	}

	// update row sums and objective function
	if (type == UPDATE) {
		for (i = 0; i < m; i++) {
			x_ij = x[i + j * m];

			if (R_finite(x_ij)) {
				pl->row_sum[i + l * m] += w * x_ij;
				pl->row_size[i + l * m] += w;
			}
		}

		pl->obj_tot += obj_delta;
	}

	return obj_delta;
}


// compute best improvement gotten from moving row i; return best
// new cluster in kbest
double ploglik_improve_row(const struct ploglik *pl, int i, int *kbest)
{
	struct ploglik *pl1 = (struct ploglik *)pl; // discard const
	const int K = pl->K, kold = pl->row_cl[i];
	double delta, delta_best, delta_out, delta_in;
	int k;

	*kbest = kold;
	delta_best = 0.0;

	if (K == 1)
		return delta_best;

	// compute how much objective changes from moving i out of
	// class kold
	delta_out = ploglik_add_row(pl1, NO_UPDATE, i, kold, -1);

	// compute how much objective changes from moving i
	// into class k
	for (k = 0; k < K; k++) {
		if (k == kold)
			continue;

		delta_in = ploglik_add_row(pl1, NO_UPDATE, i, k, +1);
		delta = delta_out + delta_in;

		if (delta > delta_best) {
			delta_best = delta;
			*kbest = k;
		}
	}

	//Rprintf("best for row %d is cl %d -> %d (%lg)\n", i,
	//	kold, *kbest, delta_best);

	return delta_best;
}


// compute best improvement gotten from moving column j; return best
// new cluster in lbest
double ploglik_improve_col(const struct ploglik *pl, int j, int *lbest)
{
	struct ploglik *pl1 = (struct ploglik *)pl; // discard const
	const int L = pl->L, lold = pl->col_cl[j];
	double delta, delta_best, delta_out, delta_in;
	int l;

	*lbest = lold;
	delta_best = 0.0;

	if (L == 1)
		return delta_best;

	// compute how much objective changes from moving j out of class lold
	delta_out = ploglik_add_col(pl1, NO_UPDATE, j, lold, -1);

	// compute how much objective changes from moving j into class l
	for (l = 0; l < L; l++) {
		if (l == lold)
			continue;

		delta_in = ploglik_add_col(pl1, NO_UPDATE, j, l, +1);
		delta = delta_out + delta_in;

		if (delta > delta_best) {
			delta_best = delta;
			*lbest = l;
		}
	}

	//Rprintf("best for col %d is cl %d -> %d (%lg)\n", j,
	//	lold, *lbest, delta_best);

	return delta_best;
}


enum move_type {
	ROW_MOVE,
	COL_MOVE
};


struct move {
	enum move_type type;
	int index;
	int cluster;
	double obj_delta;
};


struct plan {
	struct move *move;
	int size;
};


void plan_init(struct plan *plan, int m, int n)
{
	int size = m + n;

	plan->size = size;
	plan->move = (struct move *)R_alloc(size, sizeof(struct move));
}


static int move_compare(const void *pa, const void *pb)
{
	const struct move *a = pa, *b = pb;

	if (a->obj_delta > b->obj_delta) {
		return -1;
	} else if (a->obj_delta == b->obj_delta) {
		return 0;
	} else {
		return +1;
	}
}


double plan_compute(struct plan *plan, const struct ploglik *pl)
{
	double delta, delta_best;
	int cl_best, i, j, m = pl->m, n = pl->n;

	delta_best = 0;

	for (i = 0; i < m; i++) {
		delta = ploglik_improve_row(pl, i, &cl_best);
		plan->move[i].type = ROW_MOVE;
		plan->move[i].index = i;
		plan->move[i].cluster = cl_best;
		plan->move[i].obj_delta = delta;

		if (delta > delta_best)
			delta_best = delta;
	}

	for (j = 0; j < n; j++) {
		delta = ploglik_improve_col(pl, j, &cl_best);
		plan->move[m + j].type = COL_MOVE;
		plan->move[m + j].index = j;
		plan->move[m + j].cluster = cl_best;
		plan->move[m + j].obj_delta = delta;

		if (delta > delta_best)
			delta_best = delta;
	}

	if (delta_best > 0) {
		qsort(plan->move, m + n, sizeof(*plan->move), move_compare);
	}

	return delta_best;
}


int plan_eval_path(struct plan *plan, struct ploglik *pl)
{
	int step, size = plan->size;
	const struct move *move;
	int i, cl_new, cl_old;
	double obj_best = pl->obj_tot;
	int step_best = -1;

	for (step = 0; step < size; step++) {
		move = &plan->move[step];
		cl_new = move->cluster;
		i = move->index;
		switch (move->type) {
		case ROW_MOVE:
			cl_old = pl->row_cl[i];
			if (cl_new != cl_old) {
				ploglik_add_row(pl, UPDATE, i, cl_old, -1);
				ploglik_add_row(pl, UPDATE, i, cl_new, +1);
			}
			break;
		case COL_MOVE:
			cl_old = pl->col_cl[i];
			if (cl_new != cl_old) {
				ploglik_add_col(pl, UPDATE, i, cl_old, -1);
				ploglik_add_col(pl, UPDATE, i, cl_new, +1);
			}
			break;
		}
		if (pl->obj_tot > obj_best) {
			obj_best = pl->obj_tot;
			step_best = step;
		}
	}

	return step_best;
}


double plan_eval(struct plan *plan, struct ploglik *pl, int max_step)
{
	int step;
	const struct move *move;
	int i, cl;
	double obj;

	for (step = 0; step <= max_step; step++) {
		move = &plan->move[step];
		i = move->index;
		cl = move->cluster;
		switch (move->type) {
		case ROW_MOVE:
			pl->row_cl[i] = cl;
			break;
		case COL_MOVE:
			pl->col_cl[i] = cl;
			break;
		}
	}

	obj = ploglik_eval(pl);
	return obj;
}


SEXP biclust_dense(SEXP sx, SEXP srow_nclusters, SEXP srow_clusters0,
		   SEXP scol_nclusters, SEXP scol_clusters0,
		   SEXP sfamily, SEXP sepsilon, SEXP smaxit, SEXP strace)
{
	const double *x = REAL(sx);
	const int m = nrows(sx);
	const int n = ncols(sx);
	const int K = asInteger(srow_nclusters);
	const int L = asInteger(scol_nclusters);
	const int *row_cl0 = INTEGER(srow_clusters0);
	const int *col_cl0 = INTEGER(scol_clusters0);
	const int maxit = asInteger(smaxit);
	const double epsilon = asReal(sepsilon);
	const int trace = asLogical(strace);
	double (*entropy)(double, double);
	struct ploglik pl0, pl1;
	double val0, val1;
	struct plan plan;
	int best_step, converged, it, nmove;
	SEXP ans, names, sizes, sums, row_sizes, row_sums, row_clusters,
	     col_sizes, col_sums, col_clusters;

	if (strcmp(CHAR(STRING_ELT(sfamily, 0)), "binomial") == 0) {
		entropy = entropy_binomial;
	} else if (strcmp(CHAR(STRING_ELT(sfamily, 0)), "gaussian") == 0) {
		entropy = entropy_gaussian;
	} else if (strcmp(CHAR(STRING_ELT(sfamily, 0)), "poisson") == 0) {
		entropy = entropy_poisson;
	} else {
		error("invalid 'family' argument");
	}

	// allocate memory
	ploglik_init(&pl0, m, n, K, L, x, entropy);
	ploglik_init(&pl1, m, n, K, L, x, entropy);
	plan_init(&plan, m, n);

	// evaluate initial profile likelihood
	memcpy(pl0.row_cl, row_cl0, m * sizeof(int));
	memcpy(pl0.col_cl, col_cl0, n * sizeof(int));
	val0 = ploglik_eval(&pl0);

	converged = 0;
	nmove = 0;

	for (it = 0; it < maxit; it++) {
		if (trace) {
			Rprintf("iter: %d; loglik: %lg\n", it, val0);
		}

		if (plan_compute(&plan, &pl0) == 0) {
			// no possible local improvement
			converged = 1;
			break;
		}

		ploglik_assign(&pl1, &pl0);
		best_step = plan_eval_path(&plan, &pl1);

		val1 = plan_eval(&plan, &pl0, best_step);
		nmove += best_step + 1;

		if (val1 - val0 <= epsilon) {
			converged = 1;
			break;
		}

		val0 = val1;
	}

	if (!converged)
		warning("algorithm did not converge");

	sizes = PROTECT(allocMatrix(REALSXP, K, L));
	memcpy(REAL(sizes), pl0.size, K * L * sizeof(double));

	sums = PROTECT(allocMatrix(REALSXP, K, L));
	memcpy(REAL(sums), pl0.sum, K * L * sizeof(double));

	row_sizes = PROTECT(allocMatrix(REALSXP, m, L));
	memcpy(REAL(row_sizes), pl0.row_size, m * L * sizeof(double));

	row_sums = PROTECT(allocMatrix(REALSXP, m, L));
	memcpy(REAL(row_sums), pl0.row_sum, m * L * sizeof(double));

	row_clusters = PROTECT(allocVector(INTSXP, m));
	memcpy(INTEGER(row_clusters), pl0.row_cl, m * sizeof(int));

	col_sizes = PROTECT(allocMatrix(REALSXP, n, K));
	memcpy(REAL(col_sizes), pl0.col_size, n * K * sizeof(double));

	col_sums = PROTECT(allocMatrix(REALSXP, n, K));
	memcpy(REAL(col_sums), pl0.col_sum, n * K * sizeof(double));

	col_clusters = PROTECT(allocVector(INTSXP, n));
	memcpy(INTEGER(col_clusters), pl0.col_cl, n * sizeof(int));

	ans = PROTECT(allocVector(VECSXP, 12));
	names = PROTECT(allocVector(STRSXP, 12));
	SET_VECTOR_ELT(ans, 0, ScalarReal(pl0.obj_tot));
	SET_STRING_ELT(names, 0, mkChar("loglik"));
	SET_VECTOR_ELT(ans, 1, row_clusters);
	SET_STRING_ELT(names, 1, mkChar("row_clusters"));
	SET_VECTOR_ELT(ans, 2, col_clusters);
	SET_STRING_ELT(names, 2, mkChar("col_clusters"));
	SET_VECTOR_ELT(ans, 3, sizes);
	SET_STRING_ELT(names, 3, mkChar("sizes"));
	SET_VECTOR_ELT(ans, 4, sums);
	SET_STRING_ELT(names, 4, mkChar("sums"));
	SET_VECTOR_ELT(ans, 5, row_sizes);
	SET_STRING_ELT(names, 5, mkChar("row_sizes"));
	SET_VECTOR_ELT(ans, 6, row_sums);
	SET_STRING_ELT(names, 6, mkChar("row_sums"));
	SET_VECTOR_ELT(ans, 7, col_sizes);
	SET_STRING_ELT(names, 7, mkChar("col_sizes"));
	SET_VECTOR_ELT(ans, 8, col_sums);
	SET_STRING_ELT(names, 8, mkChar("col_sums"));
	SET_VECTOR_ELT(ans, 9, ScalarInteger(it));
	SET_STRING_ELT(names, 9, mkChar("niter"));
	SET_VECTOR_ELT(ans, 10, ScalarInteger(nmove));
	SET_STRING_ELT(names, 10, mkChar("nmove"));
	SET_VECTOR_ELT(ans, 11, ScalarLogical(converged));
	SET_STRING_ELT(names, 11, mkChar("conv"));
	SET_NAMES(ans, names);

	UNPROTECT(10);
	return ans;
}
