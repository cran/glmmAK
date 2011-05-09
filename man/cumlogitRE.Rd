\name{cumlogitRE}
\alias{cumlogitRE}
\title{
  Logit and cumulative logit model with random effects 
}
\description{
  This function implements MCMC sampling for the logit model with binary
  response and the cumulative logit model for multinomial ordinal
  response. Details are given in Komárek and Lesaffre (2008). 
  On as many places as possible, the same notation as
  in this paper is used also in this manual page.

  In general, the following (cumulative) logit model for response \eqn{Y}
  is assumed:
        \tabular{rcl}{
          \eqn{\log\Bigl\{\frac{P(Y\geq 1)}{P(Y=0)}\Bigr\}}{log[P(Y>=1)/P(Y=0)]} \tab
	  \eqn{=}{=} \tab
	  \eqn{\eta^{(1)}}{eta[1]} \cr
          \eqn{\log\Bigl\{\frac{P(Y\geq 2)}{P(Y\leq 1)}\Bigr\}}{log[P(Y>=2)/P(Y<=1)]} \tab
	  \eqn{=}{=} \tab
	  \eqn{\eta^{(2)}}{eta[2]} \cr
	  \tab \eqn{\vdots}{...} \tab \cr
          \eqn{\log\Bigl\{\frac{P(Y=C)}{P(Y\leq C-1)}\Bigr\}}{log[P(Y=C)/P(Y<=C-1)]} \tab
	  \eqn{=}{=} \tab
	  \eqn{\eta^{(C)},}{eta[C],}
	}
	where the form of the linear predictors
  \eqn{\eta_1,\ldots,\eta_C}{eta[1],\ldots,eta[C]} depends on whether a
  hierarchical centering is used or not. In the following,
  \eqn{\beta}{beta} denotes fixed effects and \eqn{b}{b} random effects.

  \bold{No hierarchical centering (DEFAULT)}\cr
  The linear predictor for the \eqn{c}th logit has the following form
  \deqn{\eta^{(c)} = \beta^{(c)'} (v',\,v_b') + \beta^{*'}(x',\,x_b') + b'(v_b',\,x_b')\qquad(c=1,\dots,C),}{%
    eta[c] = beta[c]'(v', v(b)') + beta[*]'(x', x(b)') + b'(v(b)', x(b)')
    (c=1,\ldots,C),}
  where
  \eqn{\beta=(\beta^{(1)'},\ldots,\beta^{(C)'},\beta^{*'})'}{beta=(beta[1]',\ldots,beta[C]',beta[*])'}
  is the vector or the fixed-effects and \eqn{b} is a vector of random
  effects with zero location. 
  
  \bold{Hierarchical centering}\cr
  The linear predictor for the \eqn{c}th logit has the following form
  \deqn{\eta^{(c)} = \beta^{(c)'} v + \beta^{*'}x + b^{(c)'}v_b + b^{*'}x_b\qquad(c=1,\ldots,C),}{%
    eta[c] = beta[c]'v + beta[*]'x + b[c]'v(b) + b[*]x(b) (c=1,\ldots,C),}
  where
  \eqn{\beta=(\beta^{(1)'},\ldots,\beta^{(C)'},\beta^{*'})'}{beta=(beta[1]',\ldots,beta[C]',beta[*])'}
  is the vector or the fixed-effects and
  \eqn{b=(b^{(1)'},\ldots,b^{(C)'},b^{*'})'}{b=(b[1]',\ldots,b[C]',b[*]')'} is a vector of random
  effects with location
  \eqn{\alpha=(\alpha^{(1)'},\ldots,\alpha^{(C)'},\alpha^{*'})'.}{alpha=(alpha[1]',\ldots,alpha[C]',alpha[*])'.}  

  
  \bold{Normal random effects (\code{drandom="normal"})}\cr
  A vector of random effects is assumed to follow a (multivariate)
  normal distribution.

  That is, if there is \bold{no hierarchical
  centering} we assume for \eqn{b}:
  \deqn{b \sim \mbox{N}(\mathbf{0},\, \mathbf{D}_b),}{%
        b ~ N(0, D[b]),
      }
  where \eqn{\mathbf{D}_b}{D[b]} is their variance-covariance matrix.
      
  If the random effects are \bold{hierarchically centered} then we
  assume for \eqn{b=(b^{(1)'},\ldots,b^{(C)'},b^{*'})'}{b=(b[1]',\ldots,b[C]',b[*]')'}:
  \deqn{b \sim \mbox{N}(\alpha,\, \mathbf{D}_b),}{%
        b ~ N(alpha, D[b]),
      }
  where 
  \eqn{\alpha=(\alpha^{(1)'},\ldots,\alpha^{(C)'},\alpha^{*'})'}{alpha=(alpha[1]',\ldots,alpha[C]',alpha[*])'}  
  is a vector of random effect locations (means) and \eqn{\mathbf{D}_b}{D[b]} is their variance-covariance matrix.

    Further, in a Bayesian model,
    it is assumed that \eqn{\mathbf{D}_b^{-1}}{D[b]^(-1)} has
    a Wishart \eqn{W(\nu_b,\,S_b)}{W(nu[b], S[b])} prior with
    \eqn{\nu_b}{nu[b]}  degrees of freedom and scale
    matrix \eqn{S_b}{S[b]}. That is, a priori
    \deqn{\mathbf{D}_b^{-1} \sim \mbox{W}(\nu_b,\,S_b),}{%
          D[b]^(-1) ~ W(nu[b], S[b]),}
    \deqn{\mbox{E}(\mathbf{D}_b^{-1}) = \nu_b\,S_b.}{%
	  E(D[b]^(-1)) = nu[b]*S[b].}
    Note that \eqn{\nu_b}{nu[b]} must be higher than the number of random
    effects minus 1.

    Alternatively, when there is only a univariate random effect with
    the variance \eqn{d_b^2}{d[b]^2}, it is possible to specify
    a uniform prior for the standard deviation of the random
    effect. That is, a priori
    \deqn{d_b \sim \mbox{Unif}(0,\,S).}{d[b] ~ Unif(0, S).}

    \bold{G-spline distributed random effects
      (\code{drandom="gspline"})}\cr
    See \code{\link{gspline1}} and \code{\link{gspline2}} for
    description of the G-spline (penalized Gaussian mixture)
    distribution.
    Further, see Komárek and Lesaffre (2008)
    for description of the prior distribution on the G-spline. Brief
    description follows here as well.

    \bold{Univariate} G-spline \cr
    If there is only a univariate random effect in the model and its
    distribution is specified as a G-spline that it is assumed that
    \deqn{
      b \sim \alpha + \sum_{j=-K}^K w_j \mbox{N}\Bigl(\tau\mu_j,\,(\tau\sigma)^2\Bigr),
    }{%
      b ~ alpha + sum[j=-K]^K w[j] N(tau*mu[j], (tau*sigma)^2),
    }
    where \eqn{\alpha}{alpha} is a location parameter, \eqn{\tau}{tau}
    is a scale parameter and \eqn{w=(w_{-K}, \ldots,
      w_K)'}{w=(w[-K],\ldots,w[K])'} is a vector of G-spline
    weights. For \bold{hierarchically centered} random effects, the
    location parameter \eqn{\alpha}{alpha} is fixed to zero.

    Further,
    \eqn{M=(\mu_{-K},\ldots,\mu_{K})'}{M=(mu[-K],\ldots,mu[K])'} is a
    vector of fixed equidistant knots (component means), where
    \deqn{\mu_j = j\delta\qquad (j=-K,\ldots,K)}{mu[j] = j*delta
      (j=-K,\ldots,K)}
    and \eqn{\sigma}{sigma} is a fixed basis standard deviation.

    The constraints \eqn{0<w_j<1\;(j=-K,\ldots,K)}{0<w[j]<1
      (j=-K,\ldots,K)} and \eqn{\sum_{j=-K}^K w_j = 1}{sum[j=-K]^K w[j]
      = 1} are sufficient for G-spline to be a density. To avoid
    constraint estimation we will estimate transformed weights
    \eqn{a=(a_{-K},\ldots,a_K)'}{a=(a[-K],\ldots,a[K])'} instead which
    relates to the original weights by
    \deqn{w_j = \frac{\exp(a_j)}{\sum_{k=-K}^K\exp(a_k)}\qquad (j=-K,\ldots,K),}{%
      w[j] = exp(a[j])/sum[k=-K]^K exp(a[k]) (j=-K,\ldots,K).
    }
    \deqn{a_j = \log\frac{w_j}{w_0}\qquad (j=-K,\ldots,K),}{%
      a[j] = log(w[j]/w[0]) (j=-K,\ldots,K).
    }     
    In the estimation procedure a penalty on the \eqn{a} coefficients in
    the form of a Gaussian Markov random field prior is imposed.

    \bold{Bivariate} G-spline \cr
    If there is a bivariate random effect \eqn{b=(b_1,\,b_2)'}{b=(b[1], b[2])'} in the model and its
    distribution is specified as a G-spline that it is assumed that
    \deqn{
      b \sim (\alpha_1,\,\alpha_2)' + \sum_{j_1=-K_1}^{K_1}\sum_{j_2=-K_2}^{K_2}
        w_{j_1,j_2}
        \mbox{N}\Bigl((\tau_1\mu_{1,j_1},\,\tau_2\mu_{2,j_2})',\,\mbox{diag}\bigl((\tau_1\sigma_1)^2,\,(\tau_2\sigma_2)^2\bigr)\Bigr),
    }{%
      b ~ (alpha[1], alpha[2])' + sum[j1=-K1]^K1 sum[j2=-K2]^K2
        w[j1,j2] N((tau[1]*mu[1,j1], tau[2]*mu[2,j2])',
  	  diag((tau[1]*sigma[1])^2, (tau[2]*sigma[2])^2)),
    }
    where \eqn{\alpha_1,\,\alpha_2}{alpha[1], alpha[2]} are location
    parameters, \eqn{\tau_1,\,\tau_2}{tau[1], tau[2]}
    are scale parameters and \eqn{W=(w_{-K_1,-K_2}, \ldots,
      w_{K_1,K_2})'}{W=(w[-K1,-K2],\ldots,w[K1,K2])'} is a matrix of G-spline
    weights. For \bold{hierarchically centered} random effects, the
    location parameters \eqn{\alpha_1,\,\alpha_2}{alpha[1], alpha[2]} are
    fixed to zero.

    Further,
    \eqn{M_1=(\mu_{1,-K_1},\ldots,\mu_{1,K_1})'}{M[1]=(mu[1,-K1],\ldots,mu[1,K1])'} is a
    vector of fixed equidistant knots (component means) in the first margin, where
    \deqn{\mu_{1,j_1} = j_1\delta_1\qquad (j_1=-K_1,\ldots,K_1)}{mu[1,j1] = j1*delta[1]
      (j1=-K1,\ldots,K1)}
    and \eqn{\sigma_1}{sigma_1} is a fixed basis standard deviation in
    the first margin. Similarly for the second margin.

    Similar reparametrization of the G-spline weights as in the
    univariate case is used to avoid constrained estimation.
}
\usage{
cumlogitRE(y, v, x, vb, xb, cluster,                       
      intcpt.random=FALSE,
      hierar.center=FALSE,                       
      drandom=c("normal", "gspline"),
      C=1,
      logit.order=c("decreasing", "increasing"),
      prior.fixed,
      prior.random,
      prior.gspline,
      init.fixed,
      init.random,
      init.gspline,                 
      nsimul = list(niter=10, nthin=1, nburn=0, nwrite=10),
      store = list(prob=FALSE, b=FALSE, alloc=FALSE, acoef=FALSE),
      dir=getwd(),
      precision=8)
}
\arguments{
  \item{y}{response vector taking values \eqn{0, 1, \dots,
      C}{0,1,...,C}.}
  \item{v}{vector, matrix or data.frame with covarites for \bold{fixed}
    effects whose effect does not necessarily satisfy proportional odds
    assumption.   

    If the argument \code{intcpt.random} is set to \code{FALSE} then the
    fixed intercept is included by default in the model. The intercept
    column should be included neither in \code{x}, nor in \code{v}.
  }  
  \item{x}{vector, matrix or data.frame with covarites for \bold{fixed}
    effects whose effect is assumed to satisfy proportional odds
    assumption.
  }
  \item{vb}{vector, matrix or data.frame with covarites for
    \bold{random} effects whose effect does not necessarily
    satisfy proportional odds assumption.

    If you want to include \bold{random intercept}, do it by setting the argument
    \code{intcpt.random} to \code{TRUE}. The intercept column should be included neither in
    \code{xb}, nor in \code{vb}. 
  }  
  \item{xb}{vector, matrix or data.frame with covarites for
    \bold{random} effects whose effect is assumed to satisfy
    proportional odds assumption.
  }
  \item{cluster}{vector which determines clusters. Needed only when
    there are any random effects in the model.
  }  
  \item{intcpt.random}{logical indicating whether a \bold{random}
    intercept should be included in the model.
  }
  \item{hierar.center}{logical indicating whether a hierarchical
    centering of random effects should be used or not.
  }  
  \item{drandom}{character indicating assumed distribution of random
    effects (if there are any).
  }
  \item{C}{number of response categories minus 1.}
  \item{logit.order}{either "decreasing" or "increasing" indicating in
    which direction the logits are formed. See the same argument in
    \code{\link{cumlogit}} for more details.

    Currently, only "decreasing" is implemented for \code{cumlogitRE}.
  }  
  \item{prior.fixed}{list specifying the prior distribution for the
    fixed effects (regression coefficients). 
    \describe{
      \item{mean}{vector giving the prior mean for each regression
	coefficient. It can be a single number only in which case it is
	recycled.
      }
      \item{var}{vector giving the prior variance for each regression
	coefficient. It can be a single number only in which case it is
	recycled.
      }
    }
    In the case that prior mean and/or variances are different for different
    regression coefficients then they should be given in the following
    order:
    intercept for the first logit,
    \eqn{\beta}{beta}(\code{v}) for the first logit, \dots, 
    intercept for the last logit,
    \eqn{\beta}{beta}(\code{v}) for the last logit,    
    \eqn{\beta}{beta}(\code{x}),
  }
  \item{prior.random}{list specifying the prior distribution for the
    parameters of the distribution of random effects. Composition of
    this list depends on the chosen distribution of random effects
    (\bold{normal} or \bold{Gspline}).
        
    \describe{
      \item{Mdistrib}{character specifying the prior distribution of the
	location of random effects. It is ignored if
	\code{hierar.center}=\code{FALSE}, in which case the location of
	the random effects is fixed to zero.
	It can be one of the following.

        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{fixed}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
            Locations of random effects are assumed
	    to be fixed and are not updated.

        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{normal}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
 	    Locations of random effects are assumed to be
            apriori normally distributed. Parameters of the normal
 	    distribution are specified further by items
	    \code{Mmean} and \code{Mvar}.

            This is also a default choice when \code{Mdistrib} is
  	    not specified.
      }	 
      \item{Mmean}{vector giving the prior means for the locations of
	random effects. It can be a single number only in which case it is
	recycled.

	It is ignored when \code{hierar.center} is \code{FALSE} in which
	case all random effects have zero location.
      }
      \item{Mvar}{vector giving the prior variances for the locations of
	random effects. It can be a single number only in which case it is
	recycled.

        In the case that prior means and/or variances are different for different
        locations of random effects then they should be given in the similar
        order as specified above for argument \code{prior.fixed}.

	It is ignored when \code{hierar.center} is \code{FALSE} in which
	case all random effects have zero location.	
      }
      \item{Ddistrib}{character specifying the prior distribution of the
        covariance matrix of random effects. It can be one of the
	following.

        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{fixed}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
            covariance matrix of random effects is assumed to
  	    be fixed and is not updated.
	
        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{wishart}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
            inverse of the covariance matrix of
	    \bold{normally} distributed random
	    effects is assumed to have a priori Wishart
	    distribution. Parameters of the Wishart distribution are
	    specified further by items \code{Ddf} and \code{DinvScale}.

	    This is also a default choice when \code{Ddistrib} is not
	    specified and random effects are normally distributed.

            This option is not allowed when the random effects
	    distribution is modelled using \bold{G-splines}.
	
        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{sduniform}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
            when random effects are \bold{normally}
	    distributed then this option may be used when there is only a
	    univariate random effect in the model. Its standard
	    deviation is then assumed to follow a uniform distribution
	    on the interval \eqn{(0,S)}{(0,S)}. Value of \eqn{S}{S} is
	    specified further by an item \code{Dupper}.

            When random effects distribution is modelled using
	    \bold{G-splines} then \code{sduniform} can be used also for
	    multivariate random effects. It is then assumed that the
	    overall standard deviation \eqn{\tau_m}{tau[m]} in the
	    \eqn{m}{m}th margin follows a uniform distribution on the
	    interval \eqn{(0,S_m)}{(0,S[m])}.	    

        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{gamma}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
            when random effects are \bold{normally}
	    distributed then this option may be used when there is only a
	    univariate random effect in the model. Its inverse variance
	    is then assumed to have a Gamma prior with the shape specified
	    further by the item \code{Dshape} and the rate (inverse
	    scale) specified by the item \code{DinvScale}.

            When random effects distribution is modelled using
	    \bold{G-splines} then \code{gamma} can be used also for
	    multivariate random effects. It is then assumed that the
	    overall inverse variance \eqn{\tau_m^{-2}}{tau[m]^[-2]} in the
	    \eqn{m}{m}th margin follows	a Gamma prior with the shapes
	    specified further by the item \code{Dshape} and the rates
	    (inverse scales) specified by the item \code{DinvScale}.
      }
      \item{Ddf}{degrees of freedom \eqn{\nu_b}{nu[b]} for the Wishart prior distribution
	of the inverse covariance matrix of \bold{normally} distributed random effects. 
      }
      \item{Dshape}{number or vector with shape parameters for the gamma
	priors of the variance components of the random effects.

	If it is a single number and random effects are multivariate it
	may be recycled.
      }	
      \item{DinvScale}{number or matrix determining the inverse scale
	matrix \eqn{S_b^{-1}}{S[b]^(-1)} for the Wishart prior distribution of the inverse
	covariance matrix of \bold{normally} distributed random effects.

        Or number or vector giving the rate (inverse scale) parameter(s)
	for the gamma priors of the variance components.

	If it is a single number and \code{Ddistrib} is \bold{wishart}
	then it is assumed that \eqn{S_b^{-1}}{S[b]^(-1)}
	is diagonal with that single number on a diagonal.
      }
      \item{Dupper}{upper limit for the uniform distribution of the
	standard deviation of the random effect(s) when
	\code{Ddistrib} is equal to \bold{sduniform}.

	If it is a single number and random effects are multivariate it
	may be recycled.
      }	
    }
  }
  \item{prior.gspline}{%%% prior.gspline
    list specifying the G-spline distribution of
    random effects and prior distribution of the G-spline
    parameters. This argument is required only when \code{drandom} is equal
    to \bold{gspline}. In the following let \eqn{q} denote the
    number (dimension) of random effects.

    The list \code{prior.gspline} can have the following components.
    \describe{%%% components of prior.gspline
      \item{K}{vector of length \eqn{q} or a number (it is
	recycled) which specifies, for each marginal G-spline, the
	number of knots on each side of the zero knot. That is, the
	\eqn{m}-th marginal G-spline has \eqn{2K_m+1}{2K[m]+1} knots.

        It is set to 15 if not explicitely specified.
      }
      \item{delta}{vector of length \eqn{q} or a number (it is
	recycled) which specifies the distance between two consecutive
	knots for each marginal G-spline. That is, the \eqn{m}{m}-th
	marginal G-spline has the following knots
        \deqn{\mu_{m,j} = j\,\delta_m,\; j=-K_m,\dots,K_m.}{%
          mu[m,j] = j*delta[m], j=-K[m],...,K[m].}

	It is set to \eqn{0.3}{0.3} if not explicitely specified.
      }
      \item{sigma}{vector of length \eqn{q} or a number (it is
	recycled) which specifies the basis standard deviation of each
	marginal G-spline.

	\eqn{\sigma_m}{sigma[m]} is set to
	\eqn{(2/3)\delta_m}{(2/3)*delta[m]} if not explicitely specified.
      }
      \item{CARorder}{vector of length \eqn{q} or a number (it is
	recycled) giving the order of the intrinsic
	conditional autoregression in the Gaussian Markov random field
	prior for the transformed G-spline weights in each margin.

        It does not need to be specified when \code{neighbor.system} is
	different from \code{uniCAR}.
	
	It is set to 3 if not explicitely specified.
      }
      \item{neighbor.system}{
	character specifying the type of the
	Gaussian Markov random field in the prior for the transformed
	G-spline weights \eqn{a}{a} of a \bold{bivariate} G-spline. It does not
	have to be specified for univariate G-splines.	
	It can be one of the following.
	
        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{uniCAR}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
            univariate (in each margin) conditional
	    autoregression as described in
	    Komárek and Lesaffre (2008). That is a priori
	    \deqn{
            p(a\,|\,\lambda) \propto
            \exp\Biggl[
            -\biggl\{
            \frac{\lambda_1}{2}\sum_{j_1}\cdots\sum_{j_q}\Bigl(\Delta_1^d a_{j_1,\ldots,j_q}\Bigr)^2
            + \cdots + 
            \frac{\lambda_q}{2}\sum_{j_1}\cdots\sum_{j_q}\Bigl(\Delta_q^d a_{j_1,\ldots,j_q}\Bigr)^2
            \biggr\}
            \Biggr],
	    }{%
            p(a|lambda) propto
            \exp[
            -{
            lambda[1]/2 \sum[j1]\dots\sum[jq](\Delta[1]^d a[j1,\ldots,jq])^2
            + \dots + 
            lambda[q]/2 \sum[j1]\dots\sum[jq](\Delta[q]^d a[j1,\ldots,jq])^2
            }
            ],
	    }
            where \eqn{\Delta_m^d}{Delta[m]^d} is a difference operator of order \eqn{d} in the \eqn{m}th margin,
            e.g.,
	    \deqn{\Delta_1^3 a_{j_1,j_2,\ldots,j_q} =
	      a_{j_1,j_2,\ldots,j_q} - 3a_{j_1-1,j_2,\ldots,j_q} +
	      3a_{j_1-2,j_2,\ldots,j_q} - a_{j_1-3,j_2,\ldots,j_q},}{%
              Delta[1]^3 a[j1,j2,\ldots,jq] =
	      a[j1,j2,\ldots,jq] - 3a[j1-1,j2,\ldots,jq] +
	      3a[j1-2,j2,\ldots,jq] - a[j1-3,j2,\ldots,jq],
	    }  
            and \eqn{\lambda = (\lambda_1,\ldots,\lambda_q)'}{
	       lambda = (lambda[1], \ldots, lambda[q])'} are smoothing hyperparameters. 

            This is also a default choice when \code{neighbor.system} is not specified.	    
	
        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{eight.neighbors}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
            this prior is applicable for \bold{bivariate} G-splines only
	    and is based on eight nearest neighbors in a spatial
	    meaning.  That is, except on edges, each full conditional of
	    \eqn{a}{a} depends only on eight nearest neighbors and local quadratic smoothing.
	    The prior is then defined as
	    \deqn{
              p(a\,|\,\lambda) \propto \exp
	      \Bigl\{-\frac{\lambda}{2}\sum_{j_1=-K_1}^{K_1-1}\sum_{j_2=-K_2}^{K_2-1}
	      \bigl(\Delta a_{j_1,j_2} \bigr)^2\Bigr\},
	    }{%
             p(a|lambda) propto exp(-lambda/2 *
	     sum[j1=-K1][K1-1]sum[j2=-K2][K2-1] (Delta a[j1,j2])^2),
	    }  
            where
	    \deqn{\Delta a_{j_1,j_2} = a_{j_1,j_2} - a_{j_1+1,j_2} -
	      a_{j_1, j_2+1} + a_{j_1+1,j_2+1}.}{%
              Delta a[j1,j2] = a[j1,j2] - a[j1+1,j2] -
	      a[j1,j2+1] + a[j1+1,j2+1].
	    }
	    Parameter \eqn{\lambda}{lambda} is a common smoothing hyperparameter.

        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{twelve.neighbors}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
           not (yet) implemented	
      }	
      \item{Ldistrib}{character specifying the prior distribution of the
        smoothing hyperparameters
	\eqn{\lambda_m, m=1,\dots,q}{lambda[m], m=1,...,q}
	(precision parameters of the Markov random fields in each
	margin) or of a common hyperparameter \eqn{\lambda}{lambda}
	It can be one of the following.
	
        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{fixed}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
            smoothing hyperparameters \eqn{\lambda}{lambda} are fixed to their
	    initial values and are not updated.
	
        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{gamma}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
            each of the smoothing hyperparameters
	    \eqn{\lambda_1,\ldots,\lambda_q}{lambda[1],\ldots,lambda[q]}
	    is assumed to follow a Gamma prior with the shapes specified
	    further by the item \code{Lshape} and the rates (inverse
	    scales) specified further by the item \code{LinvScale}.

            This is also a default choice when \code{Ldistrib} is not specified.
	
        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{sduniform}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
            square root of the inversion of each
	    smoothing hyperparameter, i.e.,
	    \eqn{\sqrt{\lambda_1^{-1}},\ldots,\sqrt{\lambda_q^{-1}}}{%
	      sqrt(lambda[1]^{-1}),\ldots,sqrt(lambda[q]^{-1})}
	    is assumed to follow apriori a uniform distribution on the
	    intervals \eqn{(0,\,S_m^{\lambda}).}{(0, S[m,lambda]).}
	    Values of
	    \eqn{S_1^{\lambda},\ldots,S_q^{\lambda}}{S[1,lambda],\ldots,S[q,lambda]}
	    are specified further by the item \code{Lupper}.
      }
      \item{Lequal}{logical indicating whether all smoothing
	hyperparameters should be kept equal.

	It is set to \code{FALSE} if not explicitely specified and
	\code{neighbor.system} is \code{uniCAR}. It is always
	\code{TRUE} when \code{neighbor.system} is different from \code{uniCAR}.
      }
      \item{Lshape}{number or vector with shape parameters for the \bold{gamma}
	priors of the smoothing hyperparameters \eqn{\lambda}{lambda}.

	If it is a single number and there is more than one smoothing
	hyperparameter \eqn{\lambda}{lambda} in the model it 
	may be recycled.
      }
      \item{LinvScale}{number or vector with rate (inverse scale) parameters for the \bold{gamma}
	priors of the smoothing hyperparameters \eqn{\lambda}{lambda}.

	If it is a single number and there is more than one smoothing
	hyperparameter \eqn{\lambda}{lambda} in the model it 
	may be recycled.
      }
      \item{Lupper}{number or vector with upper limits for the uniform
	distribution on \eqn{\sqrt{\lambda^{-1}}}{sqrt(lambda^{-1})}
	parameters when \code{Ldistrib} is \bold{sduniform}.

	If it is a single number and there is more than one smoothing
	hyperparameter \eqn{\lambda}{lambda} in the model it 
	may be recycled.	
      }
      \item{Aident}{character specifying in which way the transformed
	G-spline weights (\eqn{a}{a} coefficients) are identified. It
	can be one of the following.

        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{mean}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
            with this option, the \eqn{a}{a} coefficients
	    are forced to sum up to zero and have a zero mean.

            \bold{Note} that this option usually causes problems during
	    MCMC, especially with bivariate G-splines. The reason is
	    that if there are many almost zero weights, they lead to
	    many negative \eqn{a}{a} coefficients and to satisfy the
	    zero mean constrain, there must be some \eqn{a}'s which are
	    highly positive. When exponentiating them to get weights
	    \eqn{w}, an overflow occur.
	
        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{reference}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
            with this option, one of the \eqn{a}{a}
	    coefficients in each margin is chosen as the reference one
	    and is always equal to zero. Index of the reference
	    coefficient is specified by the item \code{Areference} (see
	    below).

	    This is also a default choice when \code{Aident} is not specified.
      }	
      \item{Areference}{vector or number (it is recycled) which
	  specifies the index of the reference \eqn{a}{a} coefficient in
	  each margin in the case \code{Aident} is equal to
	  \bold{reference}. For the \eqn{m}{m}-th margin, it must be an
	  integer between \eqn{-K_m,\dots,K_m}{-K[m],...,K[m]}.

	  To avoid numerical problems, the index of the reference
	  \eqn{a}{a} coefficient may change during the MCMC.
      }
      \item{AtypeUpdate}{character specifying in which way the
	transformed G-spline weights (\eqn{a}{a} coefficients) are
	updated. It can be one of the following.

        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{slice}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
            slice sampler of Neal (2003).
	
        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{ars.quantile}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
            adaptive rejection sampling of Gilks and
	    Wild (1992) with starting abscissae being quantiles of the
	    envelop at the previous iteration.
	
        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{ars.mode}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}
             adaptive rejection sampling of Gilks and Wild
	    (1992) with starting abscissae being the mode plus/minus 3
	    times estimated standard deviation of the full conditional
	    distribution.
	
        \eqn{\mbox{\hspace{\textwidth}}}{\cr}
        \dQuote{\code{block}}\eqn{\mbox{\hspace{\textwidth}}}{\cr}	
            all \eqn{a}{a} coefficients are updated in
	    1 block using the Metropolis-Hastings algorithm. This is
	    only available for univariate G-splines.

        \eqn{\mbox{\hspace{\textwidth}}}{\cr}	    
        Default is \code{slice}.
      }	
    }	
  }  
  \item{init.fixed}{optional vector with initial values for fixed
    effects, supplied in the same order as described above in the
    argument \code{prior.fixed}.

    If not given, initials are determined from the maximum-likelihood
    fit to the model where random effects are replaced by corresponding
    fixed effects.
  }
  \item{init.random}{optional list with initial values for the random
    effects and parameters determining their distribution. It can have
    the following components.
    \describe{
      \item{b}{vector or matrix with initial values of cluster specific
	random effects. Number of rows of the matrix or the length of
	the vector must be equal to
	\code{length(unique(cluster))}. Columns of the matrix correspond
	to the random effects in the same order as described above for
	argument \code{prior.random}.

	If not given, initials are taken to be equal over clusters and
	equal to the corresponding fixed effects from the
	maximum-likelihood fit to the model with fixed effects only.
      }	
      \item{mean}{vector giving the initial values of the means of
	random effects.

	If not given, initials are taken to be equal to the
	corresponding fixed effects from the maximum-likelihood fit to
	the model with fixed effects only.
      }
      \item{var}{matrix giving the initial value for the covariance
	matrix of the random effects when \code{drandom} is
	\bold{normal}.

	Vector giving the initial values of marginal overall variances
	of the random effects when \code{drandom} is \bold{gspline}.
      }	
    }
  }
  \item{init.gspline}{optional list with initial values related to the
    G-spline distribution of random effects (if \code{drandom} is equal
    to \bold{gspline}). It can have the following components.
    \describe{
      \item{lambda}{vector with initial values of the smoothing
	hyperparameters (precision parameters of the Markov random
	field). If not fixed and not given, the initials are sampled
	from the prior distribution.}
      \item{weights}{for \bold{univariate} G-splines: vector with
	initial weights.

	For \bold{bivariate} G-splines: matrix with initial weights.

	If the initial weights do not sum up to 1 they are re-scaled.
      }
      \item{alloc}{vector or matrix with initial component allocations
	for the individual random effects. For \bold{univariate} random
	effects this should be a vector with numbers from
	\eqn{\{-K,\ldots,K\}}{{-K,\ldots,K}}.
        For \bold{bivariate} random effects this should be a matrix with
	two columns where in the first column numbers from
	\eqn{\{-K_1,\ldots,K_1\}}{{-K1,\ldots,K1}} appear and in the
        second column numbers from
	\eqn{\{-K_2,\ldots,K_2\}}{{-K2,\ldots,K2}} appear.
      }
    }  
  }
  \item{nsimul}{list indicating the length of the MCMC. It should have
    the following components.
    \describe{
      \item{niter}{total number of the MCMC iterations after discarding
	the thinned values}
      \item{nthin}{thinning of the sample}
      \item{nburn}{length of the burn-in period}
      \item{nwrite}{frequency with which the iteration count
	changes. Further, during the burn-in, only every \code{nwrite}th
        sampled value is stored on the disk}      
    }  
  }
  \item{store}{list indicating which chains (out of these not stored by default) should be compulsory
    stored. The list has the logical components with the following
    names.
    \describe{
      \item{prob}{if \code{TRUE} values of individual predictive
	probabilities are stored.
      }
      \item{b}{if \code{TRUE} values of cluster specific random effects
	are stored.
      }
      \item{alloc}{if \code{TRUE} values of allocation indicators are stored.
      }
      \item{acoef}{if \code{TRUE} and distribution of random effects is
	given as a \bold{bivariate} G-spline values of log-G-spline
	weights (\eqn{a}{a} coefficients) are stored for all components.
      }
    }      
  }
  \item{dir}{character string specifying the directory in which the
    sampled values are stored.}  
  \item{precision}{precision with which the sampled values are written
    in files.}
}
\value{
  The function returns a complete list of parameters of the prior
  distribution and initial values.

  The main task of this function is to sample from the posterior
  distribution using MCMC. Sampled values are stored in various files
  which are described below.
}
\section{Files created}{
  \describe{
    \item{iteration.sim}{one column labeled \code{iteration} with
      indeces of MCMC iterations to which the stored sampled values
      correspond.
    }
    \item{betaF.sim}{
      sampled values of the fixed effects
      \eqn{\beta=(\beta^{(1)'},\ldots,\beta^{(C)'},\beta^{*'})'.}{beta=(beta[1]',\ldots,beta[C]',beta[*])'.}

      \bold{Note} that in models with \bold{G-spline} distributed random
      effects which are not hierarchically centered, the average effect
      of the covariates involved in the random effects (needed for
      inference) is obtained as a
      sum of the corresponding \eqn{\beta}{beta} coefficient and a
      scaled mean of the G-spline. \eqn{\beta}{beta} coefficients
      adjusted in this way are stored in the file \file{betaRadj.sim}
      (see below).
    }
    \item{betaR.sim}{sampled values of the location parameters
      \eqn{\alpha=(\alpha^{(1)'},\ldots,\alpha^{(C)'},\alpha^{*'})'}{alpha=(alpha[1]',\ldots,alpha[C]',alpha[*])'}
      of the random effects when the \bold{hierarchical centering} was
      used.

      \bold{Note} that in models with \bold{G-spline} distributed random
      effects which are hierarchically centered, the average effect
      of the covariates involved in the random effects (needed for
      inference) is obtained as a
      sum of the corresponding \eqn{\alpha}{alpha} coefficient and a mean
      of the G-spline.  \eqn{\alpha}{alpha} coefficients adjusted in
      this way are stored in the file \file{betaRadj.sim} (see below).      
    }
    \item{varR.sim}{variance components of the random effects. Format of
      the file depends on the assumed distribution of the random effects.\cr

      \describe{
      \item{\bold{Normal random effects}}{
      Let \eqn{q} be the dimension of the random effects. The first
      \eqn{0.5q(q+1)} columns of \file{varR.sim} contain a lower
      triangle (in column major order) of the matrix
      \eqn{\mathbf{D}_b}{D[b]},
      the second \eqn{0.5q(q+1)} columns of \file{varR.sim} contain a lower
      triangle of the matrix \eqn{\mathbf{D}_b^{-1}.}{D[b]^(-1).}
      }
      
      \item{\bold{Univariate G-spline random effects}}{
      The first column of \file{varR.sim} contains the G-spline variance
      parameter \eqn{\tau^2}{tau^2}, the second column its inverse.
      }
      
      \item{\bold{Bivariate G-spline random effects}}{
      The first two columns of \file{varR.sim} contain the G-spline variance
      parameters \eqn{\tau_1^2, \tau_2^2}{tau[1]^2, tau[2]^2}, the
      second two columns their inverse.
      }
      }    
    }      
    \item{loglik.sim}{sampled values of the log-likelihood (conditioned
      by the values of random effects).
    }  
    \item{probability.sim}{
      sampled values of category probabilities for each observations.
      
      Created only if \code{store$prob} is \code{TRUE}.
    }  
    \item{b.sim}{
      sampled values of individual random effects.
      
      Stores complete chains only if \code{store$b} is \code{TRUE}.
    }
  }
}
\section{Files created for models with G-spline distributed
  random effects}{
  \describe{
    \item{gspline.sim}{information concerning the fixed parameters of
      the G-spline which includes: dimension \eqn{q} of the G-spline, numbers of
      knots on each side of the reference knot for each margin
      (\eqn{K_1,\ldots,K_q}{K[1],\ldots,K[q]}), basis standard
      deviations
      \eqn{\sigma_1,\ldots,\sigma_q}{sigma[1],\ldots,sigma[q]} for each
      margin and knots
      \eqn{\mu_{1,-K_1},\ldots\mu_{1,K_1},}{mu[1,-K[1]],\ldots,mu[1,K[1]],}
      \eqn{\ldots,}
      \eqn{\mu_{q,-K_q},\ldots\mu_{q,K_q}}{mu[q,-K[q]],\ldots,mu[q,K[q]]} for each margin.      
    }  
    \item{weight.sim}{this file is created only for \bold{bivariate}
      G-splines and stores the weights \eqn{w} of the G-spline which are
      higher than a certain threshold value. That is, the weights that are
      numerically equal to zero are not recorded here. The link between
      the weights and G-spline components is provided by the file \file{knotInd.sim}.
    }  
    \item{knotInd.sim}{this file is created only for \bold{bivariate}
      G-splines. In its first column, it stores the number of G-spline
      components for which the weights are recorded on a corresponding
      row of the file \file{weight.sim}. Subsequently, it stores indeces
      of the G-spline components for which the weights are given in the
      file \file{weight.sim}. The indeces are stored as single indeces
      on the scale
      \eqn{0,\ldots,(2K_1+1)(2K_2+1)-1}{0,\ldots,(2K[1]+1)*(2K[2]+1)-1}
      such that
      index \eqn{0} corresponds to the component \eqn{(-K_1,-K_2)}{(K[1],K[2])},
      index \eqn{1} corresponds to the component \eqn{(-K_1+1,-K_2)}{(K[1]+1,K[2])}, \ldots,
      index \eqn{K_1-1}{K[1]-1} corresponds to the component
      \eqn{(K_1,-K_2)}{(K[1],-K[2])}, etc.      
    }
    \item{logweight.sim}{sampled values of the transformed G-spline
      weights \eqn{a}.

      For \bold{univariate} G-spline, this file
      stores always complete chains, irrespective of the value of
      \code{store$acoef}. For \bold{bivariate} G-splines,this file
      stores complete chains only if \code{store$acoef} is \code{TRUE}.      
    }  
    \item{gmoment.sim}{first two moments of the unshifted and unscaled
      G-spline at each iteration.

      For \bold{univariate} G-spline the column labeled \samp{gmean} is
      equal to
      \deqn{\mbox{gmean}=\sum_{j=-K}^K w_j\mu_j}{gmean=sum[j=-K]^K w[j]*mu[j]} and the
      column labeled \samp{gvar} is equal to
      \deqn{\mbox{gvar}=\sum_{j=-K}^K w_j(\mu_j - \mbox{gmean})^2 + \sigma^2.}{%
	gmean=sum[j=-K]^K w[j]*(mu[j]-gmean)^2 + sigma^2.}
      
      See
      Komárek and Lesaffre, E. (2008)
      for formulas that apply in the \bold{bivariate} case.

      Values stored here are the values of
      \eqn{\beta_1^*,\ldots,\beta_q^*}{beta[1]^*,\ldots,beta[q]^*}
      and \eqn{d_{1,1}^*, d_{2,1}^*,\ldots,d_{q,q}^*}{%
               d[1,1]^*, d[2,1]^*, \ldots, d[q,q]^*}
      as defined in Komárek and Lesaffre (2008).      
    }
    \item{betaRadj.sim}{sampled values of the average (overall) effects
      of the random effects. See notes under the files \file{betaF.sim}
      and \file{betaR.sim} above.

      Values stored here are the values of
      \eqn{\gamma_1,\ldots,\gamma_q}{gamma[1],\ldots,gamma[q]}
      as defined in Komárek and Lesaffre (2008).      
    }      
    \item{varRadj.sim}{sampled components of the variance-covariance
      matrix of the random-effects.

      Values stored here are the values of
      \eqn{d_{1,1}, d_{2,1},\ldots,d_{q,q}}{%
           d[1,1], d[2,1], \ldots, d[q,q]}
      as defined in Komárek and Lesaffre (2008).      
    }      
    \item{alloc.sim}{sampled values of the component allocations
      (in Komárek and Lesaffre (2008)
      denoted by \eqn{r_i}{r[i]}) for individual random effects.

      For \bold{univariate} G-spline, the allocations are stored on the
      scale \eqn{-K,\ldots,K}.

      For \bold{bivariate} G-spline, the allocation are stored
      as single indeces on the scale
      \eqn{0,\ldots,(2K_1+1)(2K_2+1)-1}{0,\ldots,(2K[1]+1)*(2K[2]+1)-1}
      where the link between the single and double indeces is the same
      as in the file \file{knotInd.sim}.
      
      Stores complete chains only if \code{store$alloc} is \code{TRUE}.      
    }
    \item{lambda.sim}{sampled values of smoothing hyperparameter(s)
      \eqn{\lambda}{lambda}.
    }      
  }  
}
\references{
  Agresti, A. (2002).
  \emph{Categorical Data Analysis. Second edition}.
  Hoboken: John Wiley \& Sons.

  Gelfand, A. E., Sahu, S. K., and Carlin, B. P. (1995).
  Efficient parametrisations for normal linear mixed models.
  \emph{Biometrika}, \bold{82}, 479--488.

  Gilks, W. R. and Wild, P. (1992).
  Adaptive rejection sampling for Gibbs sampling.
  \emph{Applied Statistics,} \bold{41}, 337--348.

  Neal, R. M. (2003).
  Slice sampling (with Discussion).
  \emph{The Annals of Statistics,} \bold{31}, 705--767.
  
  Komárek, A. and Lesaffre, E. (2008).
  Generalized linear mixed model with a penalized Gaussian mixture as a
  random-effects distribution.
  \emph{Computational Statistics and Data Analysis}, \bold{52}, 3441--3458.
  
  Molenberghs, G. and Verbeke, G. (2005).
  \emph{Models for Discrete Longitudinal Data}.
  New York: Springer Science+Business Media.
}  
\seealso{
  \code{\link{cumlogit}}, \code{\link{logpoissonRE}}, \code{\link{glm}}, \code{\link[MASS]{polr}}.
}  
\author{
  Arnošt Komárek \email{arnost.komarek[AT]mff.cuni.cz}
}  
\examples{
### See ex-Toenail.pdf and ex-Toenail.R
### available in the documentation
### to the package
}
\keyword{models}


