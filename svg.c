#include "PileOMeth.h"
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
            val = CI(m->unmeth1[i], m->meth1[i], 0);
            maximum = (val > maximum) ? val : maximum;
        }
        if(m->meth2[i]+m->unmeth2[i]) {
            val = CI(m->unmeth2[i], m->meth2[i], 0);
            maximum = (val > maximum) ? val : maximum;
        }
    }
    maximum *= 1.05; //Allow some buffer

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
            val = CI(m->unmeth1[i], m->meth1[i], 1);
            minimum = (val < minimum) ? val : minimum;
        }
        if(m->meth1[i]+m->unmeth1[i]) {
            val = CI(m->unmeth1[i], m->meth1[i], 1);
            minimum = (val < minimum) ? val : minimum;
        }
    }
    minimum -= 0.05*minimum; //Allow some buffer

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
    for(i=m->l-1; i>=0; i--) {
        if(m->unmeth1[i]+m->meth1[i]) return i+1;
        if(m->unmeth2[i]+m->meth2[i]) return i+1;
    }
    fprintf(stderr, "[getMaxX] There were apparently no methylation calls for this strand!.\n");
    assert(1==0);
}

int *getXTicks(int maxX, int *n) {
    int *o, maxN = 7, i;

    /***************************************************************************
     *
     * Make ticks in intervals of 5, 10, 25, 50, etc. such that there are no
     * more than maxN ticks.
     *
    ***************************************************************************/
    *n = maxX/5;
    if(*n > maxN) *n = maxX/10;
    if(*n > maxN) *n = maxX/25;
    if(*n > maxN) *n = maxX/50;
    if(*n > maxN) *n = maxX/100;
    if(*n > maxN) *n = maxX/250;
    if(*n > maxN) *n = maxX/500;
    if(*n > maxN) *n = maxX/1000;

    o = malloc((*n) * sizeof(int));
    assert(o);

    for(i=0; i<*n; i++) o[i] = (i+1)*maxX/(*n);

    return o;
}

double *getYTicks(double minY, double maxY, int *n) {
    double *o, span = maxY-minY;
    int i;

    *n = (int) (1+span/0.05);
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

void plotCI(FILE *of, int minX, strandMeth *m, int which, char *col, int buffer, int dim, double minY, double maxY) {
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
    fprintf(of, "<path d=\"M %f %f\n", remapX(minX+1, m->l, buffer,dim), remapY(val, minY, maxY, buffer, dim));
    for(i=minX+1; i<m->l; i++) {
        if(meth[i]||umeth[i]) {
            val = CI(umeth[i], meth[i], 0);
            fprintf(of, "  L %f %f\n", remapX(i+1, m->l, buffer,dim), remapY(val, minY, maxY, buffer, dim));
        }
    }
    for(i=m->l-1; i>=0; i--) {
        if(meth[i]||umeth[i]) {
            val = CI(umeth[i], meth[i], 1);
            fprintf(of, "  L %f %f\n", remapX(i+1, m->l, buffer,dim), remapY(val, minY, maxY, buffer, dim));
        }
    }
    fprintf(of, "Z\" fill=\"%s\" fill-opacity=\"0.5\"/>\n", col);
}

void plotVals(FILE *of, int minX, strandMeth *m, int which, char *col, int buffer, int dim, double minY, double maxY) {
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
    fprintf(of, "<path d=\"M %f %f\n", remapX(minX+1, m->l, buffer,dim), remapY(val, minY, maxY, buffer, dim));
    for(i=minX+1; i<m->l; i++) {
        if(meth[i]||umeth[i]) {
            val = meth[i]/((double) (meth[i]+umeth[i]));
            fprintf(of, "  L %f %f\n", remapX(i+1, m->l, buffer,dim), remapY(val, minY, maxY, buffer, dim));
        }
    }
    fprintf(of, "\" stroke=\"%s\" stroke-width=\"2\" fill-opacity=\"0\"/>\n", col);
}

void makeSVGs(char *opref, strandMeth **meths) {
    double minY = 1.0, maxY = 0.0;
    int minX1 = -1, minX2 = -1, maxX = 0, hasRead1 = 0, hasRead2 = 0;
    int i, j, buffer = 40, dim = 500, nXTicks, nYTicks;
    char *oname = malloc(sizeof(char) *(strlen(opref)+strlen("_CTOT.svg ")));
    char *titles[4] = {"Original Top", "Original Bottom",
                       "Complementary to the Original Top", "Complementary to the Original Bottom"};
    char *abbrevs[4] = {"OT", "OB", "CTOT", "CTOB"};
    FILE *of;
    double *yTicks;
    int *xTicks;

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
            fprintf(of,"<text x=\"10\" y=\"%i\" transform=\"rotate(270 10, %i)\" text-anchor=\"middle\" dominant-baseline=\"text-before-edge\">%% Methylation</text>\n", buffer+(dim>>1), buffer+(dim>>1));
            fprintf(of,"<text x=\"%i\" y=\"%i\" text-anchor=\"middle\">Position</text>\n", buffer+(dim>>1), 2*buffer+dim);
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
            fprintf(of,"<text x=\"%i\" y=\"%i\" text-anchor=\"middle\" dominant-baseline=\"middle\">%g</text>\n", \
                buffer-20, buffer+dim, minY);
            fprintf(of,"<text x=\"%i\" y=\"%i\" text-anchor=\"middle\" dominant-baseline=\"middle\">%g</text>\n", \
                buffer-20, buffer, maxY);
            for(j=0; j<nYTicks; j++) {
                fprintf(of,"<line x1=\"%i\" y1=\"%f\" x2=\"%i\" y2=\"%f\" stroke=\"black\" />\n", \
                    buffer, remapY(yTicks[j], minY, maxY, buffer, dim), buffer-5, remapY(yTicks[j], minY, maxY, buffer, dim));
            }

            //Values
            for(j=0; j<meths[i]->l; j++) {
                if(meths[i]->unmeth1[j]+meths[i]->meth1[j]) hasRead1 = 1;
                if(meths[i]->unmeth2[j]+meths[i]->meth2[j]) hasRead2 = 1;
                if(hasRead1 && hasRead2) break;
            }

            //Draw the lines
            if(hasRead1) plotCI(of, minX1, meths[i], 1, "rgb(248,118,109)", buffer, dim, minY, maxY);
            if(hasRead2) plotCI(of, minX2, meths[i], 2, "rgb(0,191,196)", buffer, dim, minY, maxY);
            if(hasRead1) plotVals(of, minX1, meths[i], 1, "rgb(248,118,109)", buffer, dim, minY, maxY);
            if(hasRead2) plotVals(of, minX2, meths[i], 2, "rgb(0,191,196)", buffer, dim, minY, maxY);

            //Finish the image
            fprintf(of, "</svg>\n");

            //Clean up
            fclose(of);
            free(xTicks);
            free(yTicks);
            hasRead1 = 0;
            hasRead2 = 0;
        }
    }
    free(oname);
}
