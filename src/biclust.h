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

#ifndef BICLUST_H
#define BICLUST_H

SEXP biclust_dense(SEXP x, SEXP row_nclusters, SEXP row_clusters,
		   SEXP col_nclusters, SEXP col_clusters);

#endif // BICLUST_H
