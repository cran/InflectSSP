#' This function determines the 4 parameter or 3 parameter log fit for the proteome level curve.
#' @param Data_Quantified the median abundance values calculated in the Quantify function
#' @importFrom stats nls
#' @export
#' @return the curve fit parameters for the control and condition curves at the proteome level
#' @examples
#' \dontrun{
#' Data_CurveFit1Parameters<-CurveFit1(Data_Quantified)
#' }

CurveFit1<-function(Data_Quantified){
#Median abundance is fit using nls function. A four parameter log fit is attempted at first. If there is a convergence failure, the program changes to three parameter log fit.
  Abundance<-nls(Median_Abundance ~ B_Proteome+((T_Proteome-B_Proteome)/(1 + 10^(b_Proteome*(xmid_Proteome-log10(Temperature)))))^1, data = Data_Quantified, start = c(T_Proteome=1,B_Proteome=.1,xmid_Proteome=1.7,b_Proteome=-10),nls.control(maxiter = 500,warnOnly = TRUE))
  if(Abundance$convInfo$isConv==TRUE){
    Proteome_Abundance<-summary(Abundance)
    T_Proteome_Fit<-Proteome_Abundance$coefficients[1]
    B_Proteome_Fit<-Proteome_Abundance$coefficients[2]
    xmid_Proteome_Fit<-Proteome_Abundance$coefficients[3]
    b_Proteome_Fit<-Proteome_Abundance$coefficients[4]
   } else {
    Abundance<-nls(Median_Abundance ~ 0+((T_Proteome-0)/(1 + 10^(b_Proteome*(xmid_Proteome-log10(Temperature_C1)))))^1, data = Data_Quantified, start = c(T_Proteome=1,xmid_Proteome=1.7,b_Proteome=-10),nls.control(maxiter = 500,warnOnly = TRUE))
    Proteome_Abundance<-summary(Abundance)
    T_Proteome_Fit<-Proteome_Abundance$coefficients[1]
    B_Proteome_Fit<-0
    xmid_Proteome_Fit<-Proteome_Abundance$coefficients[3]
    b_Proteome_Fit<-Proteome_Abundance$coefficients[4]
  }

  Data_CurveFit1Parameters<-cbind(Data_Quantified,cbind(B_Proteome_Fit,cbind(T_Proteome_Fit,cbind(xmid_Proteome_Fit,cbind(b_Proteome_Fit)))))

  return(Data_CurveFit1Parameters)
}




