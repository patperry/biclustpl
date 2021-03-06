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

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "biclust.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

SEXP biclust_dense(SEXP sx, SEXP srow_nclusters, SEXP srow_clusters,
		   SEXP scol_nclusters, SEXP scol_clusters,
		   SEXP sfamily, SEXP sepsilon, SEXP smaxit, SEXP strace);


static const R_CallMethodDef CallEntries[] = {
	CALLDEF(biclust_dense, 9),
	{NULL, NULL, 0}
};


void R_init_biclustpl(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
}


void R_unload_biclustpl(DllInfo *dll)
{
}
