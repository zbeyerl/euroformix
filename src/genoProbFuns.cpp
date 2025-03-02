//helpfunction to get allele probability:
//#include <vector> //vector storage

getKit <- function(kit=NULL, what=NA, fileName = "kit.txt", folderName=NULL) {  
  .separator <- .Platform$file.sep # Platform dependent path separator. 
  
  if(is.null(folderName)) {
    packagePath <- path.package("euroformix", quiet = FALSE) # Get package path.
    folderName <- paste(packagePath,"extdata",sep=.separator) #get folder containing the filename
  }
  filePath <- paste(folderName, fileName, sep=.separator) #get full pathname of kit file
  .kitInfo <- read.delim(file=filePath, header = TRUE, sep = "\t", quote = "\"",dec = ".", fill = TRUE, stringsAsFactors=FALSE)
 
  # Available kits. Must match else if construct.
  kits<-unique(.kitInfo$Short.Name)
	if (is.null(kit)) {	# Check if NULL
		res<-kits
	} else {	# String provided.
		# Check if number or string.
		if (is.numeric(kit)) {
			index<-kit # Set index to number.
		} else {
			index<-match(toupper(kit),toupper(kits)) # Find matching kit index (case insensitive)
		}
		if (any(is.na(index))) { 		# No matching kit.
			return(NA)
		# Assign matching kit information.
		} else {
		  currentKit <- .kitInfo[.kitInfo$Short.Name==kits[index], ]
              res <- data.frame(Panel = currentKit$Panel,
                        Short.Name = currentKit$Short.Name,
                        Full.Name = currentKit$Full.Name,
                        Marker = currentKit$Marker,
                        Allele = currentKit$Allele,
                        Size = currentKit$Size,
                        Size.Min = currentKit$Size.Min,
                        Size.Max = currentKit$Size.Max,
                        Virtual = currentKit$Virtual,
                        Color = currentKit$Color,
                        Repeat = currentKit$Repeat,
                        Marker.Min = currentKit$Marker.Min,
                        Marker.Max = currentKit$Marker.Max,
                        Offset = currentKit$Offset,
                        Gender.Marker = currentKit$Gender.Marker,
                        stringsAsFactors = FALSE)
		  res$Marker <- factor(res$Marker, levels=unique(res$Marker))# Create useful factors. 
		} 
	}
 if (!is.null(kit)) {

    if(is.na(what)){  # Return all kit information.
      return(res)
 } else if (toupper(what) == "GENDER"){  # Return gender marker as string. 
      genderMarker <- as.character(unique(res$Marker[res$Gender.Marker == TRUE]))
      if(length(genderMarker) > 1){
        warning(paste("More than one gender marker returned for kit", kit))
      }
      return(genderMarker);


double prob_a(double Pa, double mm, double nn, double fst) {
	return( (mm*(fst) + (1-(fst))*Pa)/(1+(nn-1)*(fst)) ); 
}

double prob_relUnknown(int aindU, int bindU, int Ugind, double *Fvec, double fst, double *maTypedvec, double nTyped, int aindR, int bindR, int Rgind, double *ibd) {
	//Ugind = genotype index of unknown indiviudal 
	//Rgind = genotype index of related indiviudal 		
	//ibd = ibd vector (NOC long)
	
	double genoSum; //used to sum the genotype probability tfor unknowns		
	//int aindU = outG1vec[Ugind  ]; //get allele index of genotype g_1 (unknown)
	//int bindU = outG1vec[Ugind+1]; //get allele index of genotype g_2 (unknown)
	//int aindR = outG1vec[Rgind  ]; //get allele index of genotype g_1 (related)
	//int bindR = outG1vec[Rgind+1]; //get allele index of genotype g_2 (related)
	
	bool Uhom = aindU==bindU; //boolean of whether unknown genotype is homozygote
	
if(genderMarker) { 
		Uhom; }
	//First step: Calculate random match probability of unrelated situation:
	genoSum = prob_a(Fvec[aindU],maTypedvec[aindU],nTyped,fst); //init with prob 1st allele  
	if(Uhom) { //if unknown is homozygote					
		genoSum *= prob_a(Fvec[aindU],maTypedvec[aindU]+1,nTyped+1,fst); //calculate random match prob (always used) 					
	} else { //if unknown is heterozygote variant
          genoSum *= 2*prob_a(Fvec[bindU],maTypedvec[bindU],nTyped+1,fst); //calculate prob 2st allele  (and scale with 2)	}
}
	
	//Extension with kappa-coefficient: SEE FORMULAS IN TABLE A.3 in Book "A forensic practicioners guide...."
	if( Rgind != -1 ) { //if related is specified (not -1 index)
		genoSum *= ibd[0]; //multiply with kappa0
		if( Ugind==Rgind ) { //if Unknown and Related are same genotype
			genoSum+=ibd[2]; //sum with kappa2
			
			if(Uhom) { //if unknown is homozygote
				genoSum += prob_a(Fvec[aindU],maTypedvec[aindU],nTyped,fst)*ibd[1] ; //multiply with kappa1 
			} else { //if unknown is heterozygote variant
				genoSum += (prob_a(Fvec[aindU],maTypedvec[aindU],nTyped,fst)+prob_a(Fvec[bindU],maTypedvec[bindU],nTyped,fst))*ibd[1]/2 ; //multiply with kappa1						
			}	
			
		} else { //if not the same genotype we need to check overlap (a,b)~(c,d)
			bool A1eq = aindU==aindR; //check if a=c 
			bool A2eq = aindU==bindR; //check if a=d 
			bool B1eq = bindU==aindR; //check if b=c 
			bool B2eq = bindU==bindR; //check if b=d 

			if( A1eq || A2eq || B1eq || B2eq) { //if any overlap (one shared allele): THERE ARE 3 OUTCOME!!
				if(Uhom) { //if Unknown is homozygote we know what allele to use
					genoSum += prob_a(Fvec[aindU],maTypedvec[aindU],nTyped,fst)*ibd[1]/2; //multiply with kappa1 	
				} else { //if Unknown is heterozygote we need to found the non-overlapping allele 
					bool Rhom = aindR==bindR; //check if related genotype is homozygote
					int cind; //temporary variable to select non-overlapping allele
					if(A1eq || A2eq) { //Allele bindU is non-overlapping
						cind = bindU;
					} else {
						cind = aindU; //Allele aindU is non-overlapping
					}
					if(Rhom) { //should not divide by 2 if related geno is homozygote
						genoSum += prob_a(Fvec[cind],maTypedvec[cind],nTyped,fst)*ibd[1]; //multiply with kappa1 																			
					} else { //should divide by 2 if related geno is heterozygote
						genoSum += prob_a(Fvec[cind],maTypedvec[cind],nTyped,fst)*ibd[1]/2; //multiply with kappa1 																												
					}
				}
			}	//otherwise none shared allele (and we are done (kappa0 already included)												
		}
	}
	return( genoSum ); //return genotype prob
}
