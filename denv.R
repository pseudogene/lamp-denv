#
# Copyright 2014-2016, Michaël Bekaert <michael.bekaert@stir.ac.uk>
#
# This file is part of lamp-denv.
#
# lamp-denv is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# lamp-denv is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with lamp-denvs. If not, see <http://www.gnu.org/licenses/>.
#
library(adegenet);
library(ape);
library(parallel);

process_dengue <- function( prefix, fasta, output ) {
    if (!dir.exists(output)) {
        dir.create(output);
    }
    denv <- fasta2genlight(fasta, chunk=10, saveNbAlleles=TRUE, quiet=TRUE, parallel=TRUE);
    save(denv, file=paste(output, "/", prefix, ".Rsave", sep=""));

    png(paste(output,"/location_of_the_SNPs.png", sep=""));
    temp <- density(position(denv), bw=10);
    plot(temp, type="n", xlab="Position in the alignment", main="Location of the SNPs", xlim=c(0,10761));
    polygon(c(temp$x,rev(temp$x)), c(temp$y, rep(0,length(temp$x))), col=transp("blue",.3));
    points(position(denv), rep(0, nLoc(denv)), pch="|", col="blue");
    dev.off();

    #png(paste(output,"/Missing_values.png", sep=""));
    #temp <- density(glNA(denv), bw=10);
    #plot(temp, type="n", xlab="Position in the alignment", main="Location of the missing values (NAs)", xlim=c(0,10761));
    #polygon(c(temp$x,rev(temp$x)), c(temp$y, rep(0,length(temp$x))), col=transp("blue",.3));
    #points(glNA(denv), rep(0, nLoc(denv)), pch="|", col="blue");
    #dev.off();

    #png(paste(output,"/Alleles_per_loci.png", sep=""));
    #temp <- table(unlist(other(denv)));
    #barplot(temp, main="Distribution of the number \nof alleles per loci", xlab="Number of alleles", ylab="Number of sites", col=heat.colors(4));
    #dev.off();

    #temp <- temp[-1];
    #temp <- 100*temp/sum(temp);
    #round(temp,1);

    #png(paste(output,"/2nd_Allele_frequencies.png", sep=""));
    #myFreq <- glMean(denv);
    #hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies", main="Distribution of (second) allele frequencies");
    #temp <- density(myFreq);
    #lines(temp$x, temp$y*1.8,lwd=3);
    #dev.off();

    #png(paste(output,"/Allele_frequencies.png", sep=""));
    #myFreq <- glMean(denv);
    #myFreq <- c(myFreq, 1-myFreq);
    #hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies", main="Distribution of allele frequencies", nclass=20);
    #temp <- density(myFreq, bw=.05);
    #lines(temp$x, temp$y*2,lwd=3);
    #dev.off();

    pca1 <- glPca(denv, nf=3);
    myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4);
    
    #png(paste(output,"/PCA_scatterplot.png", sep=""));
    #scatter(pca1, posi="topright");
    #title("PCA of the DENV data\n axes 1-2");
    #dev.off();

    pdf(paste(output,"/PCA_colourcatterplot.pdf", sep=""));
    abline(h=0,v=0, col="grey");
    add.scatter.eig(pca1$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3);
    title("PCA of the DENV data\n axes 1-2");
    dev.off();

    png(paste(output,"/PCA_colourcatterplot.png", sep=""));
    abline(h=0,v=0, col="grey");
    add.scatter.eig(pca1$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3);
    title("PCA of the DENV data\n axes 1-2");
    dev.off();

    write.table(pca1$scores,file=paste(output, "/", prefix, ".pca.tsv", sep=""),sep="\t");

    tre <- nj(dist(as.matrix(denv)));
    tre;
    pdf(paste(output,"/NJ_tree.pdf", sep=""));
    plot(tre, typ="fan", show.tip=FALSE)
    tiplabels(pch=20, col=myCol, cex=4);
    title("NJ tree of the DENV data");
    dev.off();
    
    pdf(paste(output,"/NJ_tree.named.pdf", sep=""));
    plot(tre, typ="fan", cex=0.7);
    title("NJ tree of the DENV data");
    dev.off();

    png(paste(output,"/NJ_tree.png", sep=""));
    plot(tre, typ="fan", show.tip=FALSE)
    tiplabels(pch=20, col=myCol, cex=4);
    title("NJ tree of the DENV data");
    dev.off();

    write.tree(tre, file = paste(output, "/", prefix, ".tree.tre", sep=""));

    save(denv, tre, pca1, file=paste(output, "/", prefix, ".Rsave", sep=""));
}

process_dengue("denv1","denv1.fasta.align.fa","denv1.results");
process_dengue("denv2","denv2.fasta.align.fa","denv2.results");
process_dengue("denv3","denv3.fasta.align.fa","denv3.results");
process_dengue("denv4","denv4.fasta.align.fa","denv4.results");
