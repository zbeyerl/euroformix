#' @title getGlist 
#' @author Oyvind Bleka
#' @description Returns a list of genotypes with corresponding genotype probabilities for an unknown individual.
#' @details The function returns the list of all possible observed genotypes for each marker. A wrapper of the function calcGjoint. 
#' @param popFreq A list of allele frequencies for a given population. Each named element in the list must be a allele-named vector with allele frequencies. 
#' @param fst Assumed theta/fst-correction
#' @param refKlist contains a list of vector of alleles for the known typed reference profiles (a,b,c,d...)
#' @param refRlist contains a list of vector of alleles for a related reference profile (a,b)
#' @param ibd The identical by decent coefficients of the relationship denotation
#' @return Glist A list with genotypes and genotype probabilities for each locus.
#' @export 
 getGlist <- function(popFreq,fst=0,refKlist=NULL,refRlist=NULL,ibd=c(1,0,0)) {
  locs <- names(popFreq)
  Glist <- list()
  for (loc in locs) {
    fstMarker = fst #set to default (can be a vector)
    if(length(fst)>1) fstMarker = fst[names(fst)==loc] #extract fst to use for marker
    if(length(fstMarker)==0) stop("The locus name in fst vector was not recognized!")
    Glist[[loc]] = calcGjoint(freq=popFreq[[loc]],nU=1,fst=fstMarker,refK=refKlist[[loc]],refR=refRlist[[loc]],ibd=ibd) 
  } 
  return(Glist)
 }
