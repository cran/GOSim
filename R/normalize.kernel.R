normalize.kernel = function(Ker, kerself1=NULL, kerself2=NULL, method="none"){
	if(method != "none"){
		if(is.null(kerself1) || is.null(kerself2)){
			if(method == "sqrt"){ # result between -1 and 1
				Kd<-sqrt(diag(Ker) + 1e-10)
				Ker<-Ker/(Kd%*%t(Kd))			
			}
			else if(method == "Lin"){ # result: diagonal = 1		
				Kd = diag(Ker)
				Ker = 2*Ker / outer(Kd, Kd, "+")
			}
			else if(method == "Tanimoto"){ 
				Kd = diag(Ker)
				Ker = Ker / (outer(Kd, Kd, "+") - Ker)
			}								
			else
				stop(paste("Unknown normalization method", method))
			diag(Ker) = 1
		}
		else{
			if(method == "sqrt")
				return(Ker / sqrt(kerself1 * kerself2))
			else if(method == "Lin")
				return(2*Ker / (kerself1 + kerself2))
			else if (method == "Tanimoto")
				return(Ker / (kerself1 + kerself2 - Ker))
			else
				stop(paste("Unknown normalization method", method))
		}
	}
	Ker
}

