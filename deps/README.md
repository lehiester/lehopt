lehopt - Copyright (c) Luke Hiester

This README describes the dependencies distributed along with lehopt and
their licensing.



# CSparse_modified

CSparse is Copyright (c) 2006-2022, Timothy A. Davis. All Rights Reserved.

Modifications copyright (c) Luke Hiester, as described below.

CSparse is licensed under the GNU Lesser General Public License 2.1 or
later. A copy is provided here in CSparse_modified/Doc/lesser.txt.
Accordingly, CSparse_modified is hereby released under the same terms.

CSparse is, at the time of writing, contained in the SuiteSparse monorepo.
Current version here is from SuiteSparse v7.8.3 (Oct. 22, 2024):
https://github.com/DrTimothyAldenDavis/SuiteSparse/releases/tag/v7.8.3

CSparse does *not* make use of BLAS routines, unlike UMFPACK.

The copy redistributed here is slightly modified, in the following ways:
- Changed `cs_lu` to return the singleton `CS_LU_RANK_DEFICIENT` to
  indicate rank deficiency as detected at cs_lu.c, line 61; declared and
  defined this singleton at bottom of Include/cs.h and Source/cs_lu.c
  respectively. All prior line numbers preserved (except for the closing
  `#endif` in cs.h).
- Removed the following directories and their contents:
  build, Demo, MATLAB, Matrix, Tcov
