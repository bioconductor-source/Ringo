/* compare two lists of genomic regions, the lists being x and y;
   return a matrix of pairwise-overlap of dimensions nx * ny

   J. Toedling, EMBL EBI, 2008
   adopted some code from Biobase's "matchpt.c" by O.Sklyar */

/* chromx, startx and endx are expected to be of same length,
   as are chromy, starty and endy */

#include <Rinternals.h>
#include <Rdefines.h>
#include <string.h>

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
SEXP overlap_xy(SEXP chromx, SEXP startx, SEXP endx, SEXP chromy, SEXP starty, SEXP endy) {

    int i, j, overlap, nprotect, nptx, npty;
    int xs, xe, ys, ye, ts, te;
    const char *chrx, *chry;
    SEXP res, dim;

    nprotect = 0;

    nptx = length(startx); // number of regions in list 1
    npty = length(starty); // number of regions in list 2

    PROTECT(res = allocVector(INTSXP, nptx*npty));
    nprotect++;

    // initialize result with all overlaps = 0
    for (i = 0; i<nptx; i++){
      for (j = 0; j<npty; j++){
	INTEGER(res)[j*nptx+i] = 0;
      }
    }

    // fill in overlaps
    for (i = 0; i<nptx; i++) {
      chrx = CHAR(STRING_ELT(chromx,i));
      for (j = 0; j<npty; j++) {
	chry = CHAR(STRING_ELT(chromy,j));
        // if on different chromosomes, overlap remains 0
	// strcmp returns 0, i.e. false, if strings equal
 	if (strcmp(chrx, chry)!=0) continue;
        // get start and end of x[i] and y[j]
	xs = INTEGER(startx)[i];
	xe = INTEGER(endx)[i];
	ys = INTEGER(starty)[j];
	ye = INTEGER(endy)[j];
	if (ys < xs){ 
          // if region y[j] is actually left of x[i], temporarily swap
	  ts = xs; te = xe;
          xs = ys; xe = ye;
          ys = ts; ye = te;
	}
	overlap = 0; 
        /* is the start of y[j] actually upstream of the end of x[i]?
           otherwise there's no overlap */
	if (ys <=  xe){ 
	  if (ye <= xe){ // is y[j] completely inside x[i]?
	    overlap = ye - ys + 1;
	  } else {
	    overlap = xe - ys + 1;
	  } 
	}//if (ys <=  xe)
	INTEGER(res)[j*nptx+i] = overlap;
      }
    }

    // return a matrix of dimension nptx * npty
    PROTECT(dim=allocVector(INTSXP, 2));
    nprotect++;
    INTEGER(dim)[0] = nptx;
    INTEGER(dim)[1] = npty;
    SET_DIM(res, dim);
  
    UNPROTECT(nprotect);
    return res;
}//overlap_xy
