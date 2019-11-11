GMS Introduction to GWAS
======================================


Suggested Software
----------------------

PLINK

R



Timetable
---------

Session 1 (9:30-10:30) Historical Perspective and Concepts

Session 2 (10:30-11:30) Steps involved in Genome-wide association studies

Break: 11:30-12:00

Session 3 (12:00-13:00): Post-GWAS Analyses: Characterising GWAS Loci

Lunch-break: 13:00-14:00

Session 4 (14:00-16:30): Practical

**Practical Requirements**
-------------------------

**A.**	Download the following files into your working directory:
1.	https://portals.broadinstitute.org/collaboration/giant/images/6/6e/Whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt.gz
2.	https://portals.broadinstitute.org/collaboration/giant/images/9/97/Whradjbmi.giant-ukbb.meta-analysis.females.23May2018.HapMap2_only.txt.gz
3.	https://portals.broadinstitute.org/collaboration/giant/images/7/71/Whradjbmi.giant-ukbb.meta-analysis.males.23May2018.HapMap2_only.txt.gz

4.	From GMS github (/GWAS_practical/):

    a.	regions.for.credible.sets.txt

    b.	crediblesets.R

    c.	sexdiff.ecf

    d.	snps.for.sexdiff.lst



**B.**	Check if the R library EasyStrata is available 

    In R:

      library(EasyStrata)

   If EasyStrata library is not installed, proceed as follows:
   
   1.	Download EasyStrata (https://homepages.uni-regensburg.de/~wit59712/easystrata/EasyStrata_8.6.tar.gz)
    
   2.	Install EasyStrata 

    In R: 

       install.packages("/path2tarball/EasyStrata_8.6.tar.gz")

