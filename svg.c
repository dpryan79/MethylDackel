#include "MethylDackel.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

//Agresti-Coull confidence interval around a methylation value
//If which is 1, return the upper 99.9% CI, otherwise the return the lower 99.9% CI
double CI(uint32_t um, uint32_t m, int which) {
    double ZZ, Z, N_dot, P_dot, X, N;

    X = (double) m;
    N = (double) (m+um);
    ZZ = 10.8275661707; //qnorm(0.9995)^2
    Z = 3.2905267315; //qnorm(0.9995)
    N_dot = N + ZZ;
    P_dot = (1.0/N_dot)*(X+0.5*ZZ);
    if(which)
        return P_dot + Z*sqrt((P_dot/N_dot)*(1-P_dot));
    else
        return P_dot - Z*sqrt((P_dot/N_dot)*(1-P_dot));
}

double getMaxY(strandMeth *m) {
    double maximum = 0.0;
    double val;
    int i;

    for(i=0; i<m->l; i++) {
        if(m->meth1[i]+m->unmeth1[i]) {
            val = CI(m->unmeth1[i], m->meth1[i], 1);
            maximum = (val > maximum) ? val : maximum;
        }
        if(m->meth2[i]+m->unmeth2[i]) {
            val = CI(m->unmeth2[i], m->meth2[i], 1);
            maximum = (val > maximum) ? val : maximum;
        }
    }
    maximum += 0.03; //Allow some buffer

    //ceil() to nearest multiple of 0.05
    if(5*(((int)ceil(100*maximum))/5)-(int)ceil(100*maximum)) {
        maximum = (1+((int)ceil(100*maximum))/5)*0.05;
    } else {
        maximum = ((int)ceil(100*maximum)/5)*0.05;
    }

    if(maximum > 0.8) maximum = 1.0;
    assert(maximum > 0.0);
    return maximum;
}

double getMinY(strandMeth *m) {
    double minimum = 1.0;
    double val;
    int i;

    for(i=0; i<m->l; i++) {
        if(m->meth1[i]+m->unmeth1[i]) {
            val = CI(m->unmeth1[i], m->meth1[i], 0);
            minimum = (val < minimum) ? val : minimum;
        }
        if(m->meth2[i]+m->unmeth2[i]) {
            val = CI(m->unmeth2[i], m->meth2[i], 0);
            minimum = (val < minimum) ? val : minimum;
        }
    }
    minimum -= 0.03; //Allow some buffer

    //floor() to nearest multiple of 0.05
    minimum = 0.01*(5*(((int) (100*minimum))/5));

    if(minimum < 0.2) minimum = 0.0;
    assert(minimum<1.0);
    return minimum;
}

int getMinX(strandMeth *m, int which) {
    int i;
    for(i=0; i<m->l; i++) {
        if(which==1) {
            if(m->unmeth1[i]+m->meth1[i]) return i;
        } else {
            if(m->unmeth2[i]+m->meth2[i]) return i;
        }
    }
    return m->l;
}

int getMaxX(strandMeth *m) {
    int i;
    for(i=m->l; i>0; i--) {
        if(m->unmeth1[i-1]+m->meth1[i-1]) break;
        if(m->unmeth2[i-1]+m->meth2[i-1]) break;
    }
    if(i>=0) {
        //Round to the next multiple of 5
        if(i%5) i += 5-(i%5);
        return i;
    }
    fprintf(stderr, "[getMaxX] There were apparently no methylation calls for this strand!.\n");
    assert(1==0);
}

int *getXTicks(int maxX, int *n) {
    int *o, maxN = 7, i, span = 5;

    /***************************************************************************
     *
     * Make ticks in intervals of 5, 10, 25, 50, etc. such that there are no
     * more than maxN ticks.
     *
    ***************************************************************************/
    *n = maxX/5;
    if(*n > maxN) {
        span = 10;
        *n = maxX/span;
    } else if(*n > maxN) {
        span = 25;
        *n = maxX/span;
    } else if(*n > maxN) {
        span = 50;
        *n = maxX/span;
    } else if(*n > maxN) {
        span = 100;
        *n = maxX/span;
    } else if(*n > maxN) {
        span = 250;
        *n = maxX/span;
    } else if(*n > maxN) {
        span = 500;
        *n = maxX/span;
    } else if(*n > maxN) {
        span = 1000;
        *n = maxX/span;
    }

    o = malloc((*n) * sizeof(int));
    assert(o);

    for(i=0; i<*n; i++) o[i] = (i+1)*span;

    return o;
}

double *getYTicks(double minY, double maxY, int *n) {
    double *o, span = maxY-minY;
    int i;

    *n = (int) (1+ceil(span/0.05));
    if(span < 0.05) *n = 2; //Just show the bounds

    o = malloc(sizeof(double) * (*n));
    assert(o);

    for(i=0; i<*n; i++) o[i] = 0.05*i+minY;

    return o;
}

double remapY(double orig, double minY, double maxY, int buffer, int dim) {
    return buffer + dim - ((double) dim)*(orig-minY)/(maxY-minY);
}

double remapX(int orig, int maxX, int buffer, int dim) {
    return buffer + ((double) dim)*orig/((double) maxX);
}

void plotCI(FILE *of, int minX, int maxX, strandMeth *m, int which, char *col, int buffer, int dim, double minY, double maxY) {
    uint32_t *meth, *umeth;
    int32_t i;
    double val;

    if(which==1) {
        meth = m->meth1;
        umeth = m->unmeth1;
    } else {
        meth = m->meth2;
        umeth = m->unmeth2;
    }

    //Start the path
    val = CI(umeth[minX], meth[minX], 0);
    fprintf(of, "<path d=\"M %f %f\n", remapX(minX+1, maxX, buffer,dim), remapY(val, minY, maxY, buffer, dim));
    for(i=minX+1; i<=m->l; i++) {
        if(meth[i]||umeth[i]) {
            val = CI(umeth[i], meth[i], 0);
            fprintf(of, "  L %f %f\n", remapX(i+1, maxX, buffer,dim), remapY(val, minY, maxY, buffer, dim));
        }
    }
    for(i=m->l-1; i>=0; i--) {
        if(meth[i]||umeth[i]) {
            val = CI(umeth[i], meth[i], 1);
            fprintf(of, "  L %f %f\n", remapX(i+1, maxX, buffer,dim), remapY(val, minY, maxY, buffer, dim));
        }
    }
    fprintf(of, "Z\" fill=\"%s\" fill-opacity=\"0.2\"/>\n", col);
}

void plotVals(FILE *of, int minX, int maxX, strandMeth *m, int which, char *col, int buffer, int dim, double minY, double maxY) {
    uint32_t *meth, *umeth;
    int32_t i;
    double val;

    if(which==1) {
        meth = m->meth1;
        umeth = m->unmeth1;
    } else {
        meth = m->meth2;
        umeth = m->unmeth2;
    }
    assert(minX>=0);

    //Start the path
    val = meth[minX]/((double) (meth[minX]+umeth[minX]));
    fprintf(of, "<path d=\"M %f %f\n", remapX(minX+1, maxX, buffer,dim), remapY(val, minY, maxY, buffer, dim));
    for(i=minX+1; i<=m->l; i++) {
        if(meth[i]||umeth[i]) {
            val = meth[i]/((double) (meth[i]+umeth[i]));
            fprintf(of, "  L %f %f\n", remapX(i+1, maxX, buffer,dim), remapY(val, minY, maxY, buffer, dim));
        }
    }
    fprintf(of, "\" stroke=\"%s\" stroke-width=\"2\" fill-opacity=\"0\"/>\n", col);
}

//Determine inclusion threshold suggestions. The general algorithm is as follows:
// (1) Get the average methylation of the middle 60% of the array
// (2) Determine the min/max of the confidence intervals of the middle 60%
// (3) Flag any point for exclusion if 
//     (A) Its lower/upper CI is above/below #1
//     (B) Its value is above/below #2
//     (C) Its absolute difference from #1 is at least 0.05
// (4) Moving from the middle of the array out, stop at the first found position matching the above
// (5) If no position is found, the bounds are 0 (bounds are 1-based)
void getThresholds(strandMeth *m, int which, int *lthresh, int *rthresh) {
    uint32_t *meth, *umeth;
    int i, total = 0, middle = m->l/2;
    double average = 0.0, minCI = 1.0, maxCI = 0.0, tmp, tmp2;

    if(which==1) {
        meth = m->meth1;
        umeth = m->unmeth1;
    } else {
        meth = m->meth2;
        umeth = m->unmeth2;
    }

    //Step 1 & 2
    for(i=(int) (0.2*m->l); i<=(int) (0.8*m->l); i++) {
        if(meth[i] || umeth[i]) {
            total++;
            average += ((double) meth[i])/((double)(meth[i]+umeth[i]));
            tmp = CI(umeth[i], meth[i], 1);
            if(minCI > tmp) minCI = tmp;
            tmp = CI(umeth[i], meth[i], 0);
            if(maxCI < tmp) maxCI = tmp;
        }
    }
    if(total) average /= total;
    else {
        *lthresh = 0;
        *rthresh = 0;
         return;
    }

    //lthresh
    for(i=middle; i>=0; i--) {
        if(meth[i] || umeth[i]) {
            tmp = ((double) meth[i])/((double)(meth[i]+umeth[i]));
            tmp2 = CI(umeth[i], meth[i], 1);
            if(tmp2 < average && tmp < minCI && fabs(tmp-average) > 0.05) break;
            tmp2 = CI(umeth[i], meth[i], 0);
            if(tmp2 > average && tmp > maxCI && fabs(tmp-average) > 0.05) break;
        }
    }
    if(i>=0) *lthresh = i+2;
    else *lthresh = 0;

    //rthresh
    for(i=middle+1; i<m->l; i++) {
        if(meth[i] || umeth[i]) {
            tmp = ((double) meth[i])/((double)(meth[i]+umeth[i]));
            tmp2 = CI(umeth[i], meth[i], 1);
            if(tmp2 < average && tmp < minCI && fabs(tmp-average) > 0.05) break;
            tmp2 = CI(umeth[i], meth[i], 0);
            if(tmp2 > average && tmp > maxCI && fabs(tmp-average) > 0.05) break;
        }
    }
    if(i<m->l) *rthresh = i;
    else *rthresh = 0;
}

//"which" denotes where the methylation metrics came from:
// bit 0: CpG
// bit 1: CHG
// bit 2: CHH
void makeSVGs(char *opref, strandMeth **meths, int which) {
    double minY = 1.0, maxY = 0.0;
    int minX1 = -1, minX2 = -1, maxX = 0, hasRead1 = 0, hasRead2 = 0;
    int i, j, buffer = 80, dim = 500, nXTicks, nYTicks;
    char *oname = malloc(sizeof(char) *(strlen(opref)+strlen("_CTOT.svg ")));
    char *titles[4] = {"Original Top", "Original Bottom",
                       "Complementary to the Original Top", "Complementary to the Original Bottom"};
    char *abbrevs[4] = {"OT", "OB", "CTOT", "CTOB"};
    char *col1 = "rgb(248,118,109)";
    char *col2 = "rgb(0,191,196)";
    FILE *of;
    double *yTicks;
    int *xTicks, lthresh1, lthresh2, rthresh1, rthresh2;
    int alreadyPrinting = 0, doingLabel = 0;

    for(i=0; i<4; i++) {
        if(meths[i]->l) {
            //Get the plot bounds and tick positions
            minY = getMinY(meths[i]);
            maxY = getMaxY(meths[i]);
            minX1 = getMinX(meths[i], 1);
            minX2 = getMinX(meths[i], 2);
            maxX = getMaxX(meths[i]);
            xTicks = getXTicks(maxX, &nXTicks);
            yTicks = getYTicks(minY, maxY, &nYTicks);

            //Basic plot setup
            sprintf(oname, "%s_%s.svg", opref, abbrevs[i]);
            of = fopen(oname, "w");
            fprintf(of, "<svg height=\"%i\" width=\"%i\"\n", dim+2*buffer, dim+2*buffer);
            fprintf(of, "    xmlns=\"http://www.w3.org/2000/svg\"\n");
            fprintf(of, "    xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n");
            fprintf(of, "    xmlns:ev=\"http://www.w3.org/2001/xml-events\">\n");
            fprintf(of,"<title>%s Strand</title>\n", titles[i]);
            fprintf(of,"<rect x=\"0\" y=\"0\" width=\"%i\" height=\"%i\" fill=\"white\" />\n", dim+2*buffer, dim+2*buffer);
            fprintf(of,"<text x=\"%i\" y=\"%i\" text-anchor=\"middle\">%s Strand</text>\n", buffer+(dim>>1), 20, titles[i]);
            fprintf(of,"<line x1=\"%i\" y1=\"%i\" x2=\"%i\" y2=\"%i\" stroke=\"black\" />\n", buffer, buffer, buffer, buffer+dim);
            fprintf(of,"<line x1=\"%i\" y1=\"%i\" x2=\"%i\" y2=\"%i\" stroke=\"black\" />\n", buffer, buffer+dim, buffer+dim, buffer+dim);

            //Ticks and labels
            fprintf(of,"<text x=\"15\" y=\"%i\" transform=\"rotate(270 15, %i)\" text-anchor=\"middle\" dominant-baseline=\"text-before-edge\">", buffer+(dim>>1), buffer+(dim>>1));
            doingLabel = 0;
            if(which&1) {
                doingLabel = 1;
                fprintf(of, "CpG");
            }
            if(which&2) {
                if(doingLabel) fprintf(of, "/CHG");
                else fprintf(of, "CHG");
                doingLabel = 1;
            }
            if(which&4) {
                if(doingLabel) fprintf(of, "/CHH");
                else fprintf(of, "CHH");
                doingLabel = 1;
            }
            if(doingLabel) fprintf(of, " ");
            fprintf(of, "Methylation %%</text>\n");
            fprintf(of,"<text x=\"%i\" y=\"%i\" text-anchor=\"middle\">Position along mapped read (5'->3' of + strand)</text>\n", buffer+(dim>>1), buffer+dim+40);
            fprintf(of,"<line x1=\"%i\" y1=\"%i\" x2=\"%i\" y2=\"%i\" stroke=\"black\" />\n", \
                buffer, buffer+dim, buffer, buffer+dim+5);
            fprintf(of,"<text x=\"%i\" y=\"%i\" text-anchor=\"middle\">%i</text>\n", \
                buffer, buffer+dim+20, 0);
            for(j=0; j<nXTicks; j++) {
                fprintf(of,"<line x1=\"%f\" y1=\"%i\" x2=\"%f\" y2=\"%i\" stroke-dasharray=\"5 5\" stroke=\"grey\" />\n", \
                    remapX(xTicks[j], maxX, buffer, dim), buffer, remapX(xTicks[j], maxX, buffer, dim), buffer+dim);
                fprintf(of,"<line x1=\"%f\" y1=\"%i\" x2=\"%f\" y2=\"%i\" stroke=\"black\" />\n", \
                    remapX(xTicks[j], maxX, buffer, dim), buffer+dim, remapX(xTicks[j], maxX, buffer, dim), buffer+dim+5);
                fprintf(of,"<text x=\"%f\" y=\"%i\" text-anchor=\"middle\">%i</text>\n", \
                    remapX(xTicks[j], maxX, buffer, dim), buffer+dim+20, xTicks[j]);
            }
            for(j=0; j<nYTicks; j++) {
                fprintf(of,"<line x1=\"%i\" y1=\"%f\" x2=\"%i\" y2=\"%f\" stroke=\"black\" />\n", \
                    buffer, remapY(yTicks[j], minY, maxY, buffer, dim), buffer-5, remapY(yTicks[j], minY, maxY, buffer, dim));
                fprintf(of,"<text x=\"%i\" y=\"%f\" text-anchor=\"middle\" dominant-baseline=\"middle\">%4.2f</text>\n", \
                    buffer-25, remapY(yTicks[j], minY, maxY, buffer, dim), yTicks[j]);
            }

            //Values
            for(j=0; j<meths[i]->l; j++) {
                if(meths[i]->unmeth1[j]+meths[i]->meth1[j]) hasRead1 = 1;
                if(meths[i]->unmeth2[j]+meths[i]->meth2[j]) hasRead2 = 1;
                if(hasRead1 && hasRead2) break;
            }

            //Draw the lines
            if(hasRead1) plotCI(of, minX1, maxX, meths[i], 1, col1, buffer, dim, minY, maxY);
            if(hasRead2) plotCI(of, minX2, maxX, meths[i], 2, col2, buffer, dim, minY, maxY);
            if(hasRead1) plotVals(of, minX1, maxX, meths[i], 1, col1, buffer, dim, minY, maxY);
            if(hasRead2) plotVals(of, minX2, maxX, meths[i], 2, col2, buffer, dim, minY, maxY);

            //Get cutting threshold suggestions
            getThresholds(meths[i], 1, &lthresh1, &rthresh1);
            getThresholds(meths[i], 2, &lthresh2, &rthresh2);
            if(lthresh1+lthresh2+rthresh1+rthresh2) {
                fprintf(of, "<text x=\"%i\" y=\"%i\" text-anchor=\"end\">--%s %i,%i,%i,%i</text>\n", \
                    2*buffer+dim-10, 2*buffer+dim-10, abbrevs[i], lthresh1, rthresh1, lthresh2, rthresh2);
                if(lthresh1) fprintf(of, "<line x1=\"%f\" y1=\"%i\" x2=\"%f\" y2=\"%i\" stroke-dasharray=\"5 1\" stroke=\"%s\" stroke-width=\"1\" />\n", \
                    remapX(lthresh1, maxX, buffer, dim), dim+buffer, remapX(lthresh1, maxX, buffer, dim), buffer, col1);
                if(rthresh1) fprintf(of, "<line x1=\"%f\" y1=\"%i\" x2=\"%f\" y2=\"%i\" stroke-dasharray=\"5 1\" stroke=\"%s\" stroke-width=\"1\" />\n", \
                    remapX(rthresh1, maxX, buffer, dim), dim+buffer, remapX(rthresh1, maxX, buffer, dim), buffer, col1);
                if(lthresh2) fprintf(of, "<line x1=\"%f\" y1=\"%i\" x2=\"%f\" y2=\"%i\" stroke-dasharray=\"5 1\" stroke=\"%s\" stroke-width=\"1\" />\n", \
                    remapX(lthresh2, maxX, buffer, dim), dim+buffer, remapX(lthresh2, maxX, buffer, dim), buffer, col2);
                if(rthresh2) fprintf(of, "<line x1=\"%f\" y1=\"%i\" x2=\"%f\" y2=\"%i\" stroke-dasharray=\"5 1\" stroke=\"%s\" stroke-width=\"1\" />\n", \
                    remapX(rthresh2, maxX, buffer, dim), dim+buffer, remapX(rthresh2, maxX, buffer, dim), buffer, col2);
            }

            //Add some legend boxes on the right
            if(hasRead1) {
                fprintf(of, "<rect x=\"%i\" y=\"%i\" width=\"20\" height=\"20\" fill=\"%s\" />\n", dim+buffer+10, (dim>>1)+buffer-20, col1);
                fprintf(of, "<text x=\"%i\" y=\"%i\" text-anchor=\"start\" dominant-baseline=\"middle\">#1</text>\n", dim+buffer+35, (dim>>1)+buffer-10);
            }
            if(hasRead2) {
                fprintf(of, "<rect x=\"%i\" y=\"%i\" width=\"20\" height=\"20\" fill=\"%s\" />\n", dim+buffer+10, (dim>>1)+buffer, col2);
                fprintf(of, "<text x=\"%i\" y=\"%i\" text-anchor=\"start\" dominant-baseline=\"middle\">#2</text>\n", dim+buffer+35, (dim>>1)+buffer+10);
            }

            //Finish the image
            fprintf(of, "</svg>\n");

            //Print the trimming options to stderr if applicable
            if(lthresh1 + rthresh1 + lthresh2 + rthresh2) {
                if(!alreadyPrinting) fprintf(stderr, "Suggested inclusion options:");
                fprintf(stderr, " --%s %i,%i,%i,%i", abbrevs[i], lthresh1,rthresh1,lthresh2,rthresh2);
                alreadyPrinting=1;
            }

            //Clean up
            fclose(of);
            free(xTicks);
            free(yTicks);
            hasRead1 = 0;
            hasRead2 = 0;
        }
    }
    if(alreadyPrinting) fprintf(stderr, "\n");
    free(oname);
}

void makeTXT(strandMeth **m) {
    char *abbrevs[4] = {"OT", "OB", "CTOT", "CTOB"};
    int i, j;

    printf("Strand\tRead\tPosition\tnMethylated\tnUnmethylated\n");
    for(i=0; i<4; i++) {
        if(m[i]->l) {
            for(j=0; j<m[i]->l; j++) {
                if(m[i]->meth1[j] || m[i]->unmeth1[j])
                    printf("%s\t1\t%i\t%"PRIu32"\t%"PRIu32"\n", abbrevs[i], j+1, m[i]->meth1[j], m[i]->unmeth1[j]);
                if(m[i]->meth2[j] || m[i]->unmeth2[j])
                    printf("%s\t2\t%i\t%"PRIu32"\t%"PRIu32"\n", abbrevs[i], j+1, m[i]->meth2[j], m[i]->unmeth2[j]);
            }
        }
    }
}
