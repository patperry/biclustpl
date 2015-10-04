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

#include <stdlib.h>
#include <string.h>


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
};


void ploglik_init(struct ploglik *pl, int m, int n, int K, int L,
		  const double *x)
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
	pl->row_sum = (double *)R_alloc(L * m, sizeof(double));
	pl->col_sum = (double *)R_alloc(K * n, sizeof(double));
	pl->size = (double *)R_alloc(K * L, sizeof(double));
	pl->row_size = (double *)R_alloc(L * m, sizeof(double));
	pl->col_size = (double *)R_alloc(K * n, sizeof(double));
}


void ploglik_assign(struct ploglik *pl, const struct ploglik *src)
{
	const int m = pl->m, n = pl->n, K = pl->K, L = pl->L;

	memcpy(pl->row_cl, src->row_cl, m * sizeof(int));
	memcpy(pl->col_cl, src->col_cl, n * sizeof(int));
	memcpy(pl->sum, src->sum, K * L * sizeof(double));
	memcpy(pl->obj, src->obj, K * L * sizeof(double));
	memcpy(pl->row_sum, src->row_sum, L * m * sizeof(double));
	memcpy(pl->col_sum, src->col_sum, K * n * sizeof(double));
	memcpy(pl->size, src->size, K * L * sizeof(double));
	memcpy(pl->row_size, src->row_size, L * m * sizeof(double));
	memcpy(pl->col_size, src->col_size, K * n * sizeof(double));
	pl->obj_tot = src->obj_tot;
}


double ploglik_eval(struct ploglik *pl)
{
	const int m = pl->m, n = pl->n, K = pl->K, L = pl->L;
	const double *x = pl->x;
	double obj_kl, n_kl, m_kl, s_kl, x_ij;
	int i, j, k, l;

	// initialize sums and sizes to 0
	memset(pl->sum, 0, K * L * sizeof(double));
	memset(pl->size, 0, K * L * sizeof(double));
	memset(pl->row_sum, 0, L * m * sizeof(double));
	memset(pl->row_size, 0, L * m * sizeof(double));
	memset(pl->col_sum, 0, K * n * sizeof(double));
	memset(pl->col_size, 0, K * n * sizeof(double));

	// compute block sums and sizes
	for (j = 0; j < n; j++) {
		l = pl->col_cl[j];
		for (i = 0; i < m; i++) {
			x_ij = x[i + j * m];
			if (R_finite(x_ij)) {
				k = pl->row_cl[i];
				pl->sum[k + l * K] += x_ij;
				pl->size[k + l * K] += 1;
				pl->row_sum[l + i * L] += x_ij;
				pl->row_size[l + i * L] += 1;
				pl->col_sum[k + j * K] += x_ij;
				pl->col_size[k + j * K] += 1;
			}
		}
	}

	// compute objective function
	memset(pl->obj, 0, K * L * sizeof(double));
	pl->obj_tot = 0;
	for (l = 0; l < L; l++) {
		for (k = 0; k < K; k++) {
			n_kl = pl->size[k + l * K];
			s_kl = pl->sum[k + l * K];

			if (n_kl > 0) {
				// 0.5 * n_kl (m_kl)^2
				m_kl = s_kl / n_kl;
				obj_kl = 0.5 * s_kl * m_kl;
				pl->obj[k + l * K] = obj_kl;
				pl->obj_tot += obj_kl;
			}
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
	double m_kl_new, n_l, n_kl_new, obj_delta, obj_kl_old,
	       obj_kl_new, s_l, s_kl_new, x_ij;
	int j, l;

	// compute how much objective changes from moving i out of
	// class kold
	obj_delta = 0.0;
	for (l = 0; l < L; l++) {
		obj_kl_old = pl->obj[k + l * K];
		n_l = pl->row_size[l + i * L];
		s_l = pl->row_sum[l + i * L];
		n_kl_new = pl->size[k + l * K] + w * n_l;
		s_kl_new = pl->sum[k + l * K] + w * s_l;

		if (n_kl_new > 0) {
			m_kl_new = s_kl_new / n_kl_new;
			obj_kl_new = 0.5 * s_kl_new * m_kl_new;
		} else {
			obj_kl_new = 0.0;
		}

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
				pl->col_sum[k + j * K] += w * x_ij;
				pl->col_size[k + j * K] += w;
			}
		}

		pl->obj_tot += obj_delta;
	}

	//Rprintf("[obj_tot_new = %lg]", pl->obj_tot);

	return obj_delta;
}


// compute best improvement gotten from moving row i; return best
// new cluster in kbest
double ploglik_improve_row(const struct ploglik *pl, int i, int *kbest)
{
	struct ploglik *pl1 = (struct ploglik *)pl; // discard const
	const int K = pl->K, kold = pl->row_cl[i];
	double obj, obj_best, obj_out, obj_in;
	int k;

	*kbest = kold;
	obj_best = 0.0;

	if (K == 1)
		return obj_best;

	// compute how much objective changes from moving i out of
	// class kold
	obj_out = ploglik_add_row(pl1, NO_UPDATE, i, kold, -1);

	// compute how much objective changes from moving i
	// into class k
	for (k = 0; k < K; k++) {
		if (k == kold)
			continue;

		obj_in = ploglik_add_row(pl1, NO_UPDATE, i, k, +1);
		obj = obj_out + obj_in;

		if (obj > obj_best) {
			obj_best = obj;
			*kbest = k;
		}
	}

	//Rprintf("best for row %d is cl %d -> %d (%lg)\n", i,
	//	kold, *kbest, obj_best);

	return obj_best;
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
	double *path;
	int size;
};


void plan_init(struct plan *plan, int m, int n)
{
	int size = m + n;

	plan->size = size;
	plan->move = (struct move *)R_alloc(size, sizeof(struct move));
	plan->path = (double *)R_alloc(size, sizeof(double));
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
	double obj, obj_best;
	int cl_best, i, j, m = pl->m, n = pl->n;

	obj_best = 0;

	for (i = 0; i < m; i++) {
		obj = ploglik_improve_row(pl, i, &cl_best);
		plan->move[i].type = ROW_MOVE;
		plan->move[i].index = i;
		plan->move[i].cluster = cl_best;
		plan->move[i].obj_delta = obj;

		if (obj > obj_best)
			obj_best = obj;
	}

	for (j = 0; j < n; j++) {
		obj = 0.0;
		cl_best = pl->col_cl[j];
		plan->move[m + j].type = COL_MOVE;
		plan->move[m + j].index = j;
		plan->move[m + j].cluster = cl_best;
		plan->move[m + j].obj_delta = obj;

		if (obj > obj_best)
			obj_best = obj;
	}

	if (obj_best > 0) {
		qsort(plan->move, m + n, sizeof(*plan->move), move_compare);
	}

	return obj_best;
}


int plan_eval_path(struct plan *plan, struct ploglik *pl)
{
	int step, size = plan->size;
	const struct move *move;
	int i, cl_new, cl_old;
	double obj_best = 0.0;
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
			break;
		}
		if (pl->obj_tot > obj_best) {
			obj_best = pl->obj_tot;
			step_best = step;
		}

		//Rprintf("step %d (%lg)\n", step, pl->obj_tot);
	}

	//Rprintf("step_best: %d (%lg)\n", step_best, obj_best);

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
		cl = move->cluster;
		i = move->index;
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


SEXP biclust_dense(SEXP sx, SEXP srow_nclusters, SEXP srow_clusters,
		   SEXP scol_nclusters, SEXP scol_clusters)
{
	int iter, max_iter = 10000;
	const double *x = REAL(sx);
	const int m = nrows(sx);
	const int n = ncols(sx);
	const int K = asInteger(srow_nclusters);
	const int L = asInteger(scol_nclusters);
	const int *row_cl0 = INTEGER(srow_clusters);
	const int *col_cl0 = INTEGER(scol_clusters);
	int converged;
	double tol = 1e-6;

	struct ploglik pl0, pl1;
	double val0, val1;
	struct plan plan;
	int best_step;

	// allocate memory
	ploglik_init(&pl0, m, n, K, L, x);
	ploglik_init(&pl1, m, n, K, L, x);
	plan_init(&plan, m, n);

	// evaluate initial profile likelihood
	memcpy(pl0.row_cl, row_cl0, m * sizeof(int));
	memcpy(pl0.col_cl, col_cl0, n * sizeof(int));
	val0 = ploglik_eval(&pl0);

	converged = 0;

	for (iter = 0; iter < max_iter; iter++) {
		Rprintf("iter %d; value %lg\n", iter, val0);

		if (plan_compute(&plan, &pl0) == 0) {
			Rprintf("no local moves\n");
			// no possible local improvement
			converged = 1;
			break;
		}

		ploglik_assign(&pl1, &pl0);
		best_step = plan_eval_path(&plan, &pl1);

		val1 = plan_eval(&plan, &pl0, best_step);

		if (val1 - val0 <= tol) {
			converged = 1;
			break;
		}

		val0 = val1;
	}

	return NULL_USER_OBJECT;
}
