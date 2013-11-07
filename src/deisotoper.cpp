#include <Rcpp.h>
#include <base/ms/deisotoper.h>
#include <rcppdeisotoperenvelope.h>

RcppExport SEXP deisotoper_main(SEXP Dsexp_mZ_, SEXP Dsexp_intensity_, SEXP Dsexp_Z_, SEXP Dsexp_DF_, SEXP Dsexp_massError_)
{
    try {			// or use BEGIN_RCPP macro
	
    Rcpp::DataFrame DF = Rcpp::DataFrame(Dsexp_DF_);

	Rcpp::NumericVector mZ_ = Rcpp::NumericVector(Dsexp_mZ_);
	Rcpp::NumericVector intensity_ =
	    Rcpp::NumericVector(Dsexp_intensity_);

	Rcpp::IntegerVector Z = Rcpp::IntegerVector(Dsexp_Z_);

    double massError = Rcpp::as<double>(Dsexp_massError_);

    ralab::base::ms::RcppIsotopeenvelope ril(DF);
    ralab::base::ms::Deisotoper x(Z.begin(), Z.end(), massError);
    x.setIsotopPatternMap(& ril);

    x.computeIsotopChains(mZ_.begin() , mZ_.end(), intensity_.begin());
    x.assignIsotopIntensities(mZ_.begin() , mZ_.end(), intensity_.begin());

    /*
	Rcpp::List NDF =
	    Rcpp::List::create(Rcpp::Named("mZ") = mZ_,
				    Rcpp::Named("intensity") = intensity_,
                    Rcpp::Named("result") = x.deisotop(mZ_.begin() , mZ_.end(), intensity_.begin()));
    */

	Rcpp::List NDF =
	    Rcpp::List::create(Rcpp::Named("result") = x.getIsotopChainResults(), 
        Rcpp::Named("score") = x.getIsotopInnerProductsResults(),
        Rcpp::Named("score1") = x.getIsotopInnerProducts1Results(),
        Rcpp::Named("group") = x.getIsotopGroups());

	return (NDF);
    }
    catch(std::exception & ex) {	// or use END_RCPP macro
	forward_exception_to_r(ex);
    }
    catch( ...) {
	::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;		// -Wall
}
