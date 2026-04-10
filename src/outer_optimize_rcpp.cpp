// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
using namespace Rcpp;
using namespace arma;

// --- small helpers -----------------------------------------------------------

static inline rowvec make_shift_row(std::size_t p, std::size_t numParsS,
                                    const vec& thetaS) {
  // row d = c(rep(0, p - numParsS), thetaS)
  rowvec d(p, fill::zeros);
  d.subvec(p - numParsS, p - 1) = thetaS.t();
  return d;
}

static inline vec kron_vec(const vec& a, const vec& b) {
  return kron(a, b);
}

static inline vec subvec_by_R_index(const vec& x, const IntegerVector& idx1based) {
  // R is 1-based; Armadillo is 0-based
  vec out(idx1based.size());
  for (int i = 0; i < idx1based.size(); ++i) out[i] = x[idx1based[i] - 1];
  return out;
}

static inline vec clamp_and_normalize_simplex(const vec& x, double eps = 1e-12) {
  vec y = x;
  y.transform([&](double v) { return std::isfinite(v) ? std::max(v, eps) : eps; });
  double s = accu(y);
  if (!(s > 0.0) || !std::isfinite(s)) {
    y.fill(1.0 / static_cast<double>(y.n_elem));
  } else {
    y /= s;
  }
  return y;
}

static inline vec mirror_ascent_step(const vec& x, const vec& grad, double step) {
  const double eps = 1e-12;
  vec xPos = clamp_and_normalize_simplex(x, eps);
  vec z = step * grad;
  z -= z.max(); // numerical stability
  vec y = xPos % exp(z);
  return clamp_and_normalize_simplex(y, eps);
}

static inline double simplex_grad_resid_inf(const vec& grad) {
  if (grad.n_elem == 0) return 0.0;
  double gbar = mean(grad);
  return abs(grad - gbar).max();
}

static inline double rel_scalar_step(double newVal, double oldVal) {
  return std::abs(newVal - oldVal) / std::max(1.0, std::abs(oldVal));
}

static inline double rel_vec_step_l1(const vec& newVal, const vec& oldVal) {
  return accu(abs(newVal - oldVal)) / std::max(1.0, accu(abs(oldVal)));
}

// [[Rcpp::export]]
void print_list2(List x) {
  Function printFunc("print");
  printFunc(x);
}

// [[Rcpp::export]]
List outer_optimize_rcpp(const arma::mat& scrsAug,
                         const arma::vec& betaInit,
                         const arma::vec& gammaInit,
                         double aInit,
                         int LTest,
                         int JTest,
                         int numParsS,
                         const IntegerVector& SepInd, // 1-based indices from R
                         int mOtr,
                         int LSmax,
                         int mInr,
                         double tolInr,
                         double tolOtr,
                         bool verb = false) {

  Rcpp::Environment melt_ns = Rcpp::Environment::namespace_env("melt");
  Rcpp::Function el_control = melt_ns["el_control"];
  Rcpp::Function el_mean    = melt_ns["el_mean"];

  const std::size_t p = scrsAug.n_cols;
  const double epsSimplex = 1e-12;
  const double epsAlpha = 1e-12;
  const int convNeedStep = 5;
  const int objWin = 25;
  const double zMin = -700.0;
  const double zMax = 700.0;

  // Convergence tolerances
  const double tolObjWin = 10.0 * tolOtr;      // rolling objective-change tolerance
  const double tolGradStrict = tolOtr;         // strict gradient threshold
  const double tolGradObj = 1e4 * tolOtr;   // looser gradient threshold for obj-based stop
  const double tolStep = 10.0 * tolOtr;        // step stagnation tolerance

  // Initial state (force interior for mirror descent and positivity for alpha)
  vec betaCur = clamp_and_normalize_simplex(betaInit, epsSimplex);
  vec gammaCur = clamp_and_normalize_simplex(gammaInit, epsSimplex);
  double aCur = std::max(aInit, epsAlpha);
  double zCur = std::log(aCur);

  Rcpp::S4 ctrl = el_control(_["maxit_l"] = mInr, _["tol_l"] = tolInr);

  vec thetaSCur = aCur * kron_vec(betaCur, gammaCur);
  rowvec d0 = make_shift_row(p, numParsS, thetaSCur);
  Rcpp::S4 elCur = el_mean(_["x"] = scrsAug, _["par"] = d0, _["control"] = ctrl);
  double logLikCur = as<double>(elCur.slot("logl"));
  List optimCur = elCur.slot("optim");
  vec lambdaCur = as<vec>(optimCur["lambda"]);
  vec lambdaSepCur = subvec_by_R_index(lambdaCur, SepInd);
  Rcpp::S4 elNew = elCur;

  NumericVector objTrace(mOtr + 1);
  objTrace[0] = logLikCur;
  NumericVector gradTrace(mOtr + 1);
  gradTrace[0] = NA_REAL;  // no gradient computed yet

  bool LSFail = false;
  bool conv = false;
  bool convObj = false;
  bool convGrad = false;
  bool convStep = false;
  double maxGrad = NA_REAL;
  double objWinChange = NA_REAL;
  double maxStep = NA_REAL;
  int stepStreak = 0;
  int niter = 0;

  while ((niter < mOtr) && !LSFail && !conv) {
    if (verb) Rcpp::Rcout << "Iteration " << (niter + 1) << "\n";

    // Gradients at current iterate; used for simultaneous updates.
    mat K_gamma = kron(betaCur.t(), eye(JTest, JTest));
    mat K_beta = kron(eye(LTest, LTest), gammaCur.t());

    vec gradGamma = aCur * (K_gamma * lambdaSepCur);                 // ascent gradient wrt gamma
    vec gradBeta  = aCur * (K_beta  * lambdaSepCur);                 // ascent gradient wrt beta
    double gradA  = as_scalar(kron_vec(betaCur, gammaCur).t() * lambdaSepCur); // ascent grad wrt alpha
    double gradZ  = aCur * gradA;                                    // chain rule: a = exp(z)

    double gradGammaRes = simplex_grad_resid_inf(gradGamma);
    double gradBetaRes  = simplex_grad_resid_inf(gradBeta);
    double gradZRes     = std::abs(gradZ);
    maxGrad = std::max(gradZRes, std::max(gradBetaRes, gradGammaRes));
    gradTrace[niter] = maxGrad;

    // Simultaneous backtracking line search: accept the first improving step.
    bool accepted = false;
    double step = 1.0;

    vec betaNew = betaCur;
    vec gammaNew = gammaCur;
    double aNew = aCur;
    double zNew = zCur;
    double logLikNew = logLikCur;
    vec thetaSNew = thetaSCur;
    vec lambdaSepNew = lambdaSepCur;
    Rcpp::S4 elCand = elCur;

    for (int nRed = 0; nRed < LSmax; ++nRed) {
      zNew = std::min(zMax, std::max(zMin, zCur + step * gradZ));
      aNew = std::exp(zNew);
      betaNew = mirror_ascent_step(betaCur, gradBeta, step);
      gammaNew = mirror_ascent_step(gammaCur, gradGamma, step);

      vec thetaSTmp = aNew * kron_vec(betaNew, gammaNew);
      rowvec d = make_shift_row(p, numParsS, thetaSTmp);
      mat gTmp = scrsAug.each_row() - d;

      Rcpp::S4 ctrlLS = el_control(_["maxit_l"] = mInr, _["tol_l"] = tolInr);
      Rcpp::S4 elTmp = el_mean(_["x"] = scrsAug, _["par"] = d, _["control"] = ctrlLS);
      double logl = as<double>(elTmp.slot("logl"));

      if (logl > logLikCur) {
        accepted = true;
        logLikNew = logl;
        thetaSNew = thetaSTmp;
        elCand = elTmp;
        List opt = elTmp.slot("optim");
        vec lambdaAll = as<vec>(opt["lambda"]);
        lambdaSepNew = subvec_by_R_index(lambdaAll, SepInd);
        break;
      }

      step *= 0.5;
    }

    if (!accepted) {
      // No improving step found.
      // Treat this as stagnation support, but only allow convergence if the gradient is small enough.
      convObj = true;
      convStep = true;
      convGrad = (maxGrad < tolGradStrict);

      if (maxGrad < tolGradStrict) {
        ++stepStreak;
      } else {
        stepStreak = 0;
      }

      bool convObjStop = (maxGrad < tolGradObj);
      conv = convObjStop || (stepStreak >= convNeedStep);
      LSFail = !conv;
      objTrace[niter + 1] = logLikCur;
      gradTrace[niter + 1] = maxGrad;
      ++niter;
      break;
    }

    // Rolling-window objective change, evaluated using the candidate accepted objective.
    if ((niter + 1) >= objWin) {
      double oldObj = objTrace[niter + 1 - objWin];
      objWinChange = std::abs(logLikNew - oldObj) / std::max(1.0, std::abs(oldObj));
    } else {
      objWinChange = NA_REAL;
    }
    double stepGamma = rel_vec_step_l1(gammaNew, gammaCur);
    double stepBeta = rel_vec_step_l1(betaNew, betaCur);
    double stepZ = rel_scalar_step(zNew, zCur);
    maxStep = std::max(stepZ, std::max(stepBeta, stepGamma));

    convObj = ((niter + 1) >= objWin) && R_finite(objWinChange) && (objWinChange < tolObjWin);
    convGrad = (maxGrad < tolGradStrict);
    convStep = (maxStep < tolStep);

    // Step-based convergence: require the stricter gradient threshold.
    if (convStep && (maxGrad < tolGradStrict)) {
      ++stepStreak;
    } else {
      stepStreak = 0;
    }

    bool convObjStop = convObj && (maxGrad < tolGradObj);
    conv = convObjStop || (stepStreak >= convNeedStep);

    // Accept current iterate.
    betaCur = betaNew;
    gammaCur = gammaNew;
    aCur = aNew;
    zCur = zNew;
    thetaSCur = thetaSNew;
    lambdaSepCur = lambdaSepNew;
    logLikCur = logLikNew;
    elNew = elCand;

    objTrace[niter + 1] = logLikCur;
    gradTrace[niter + 1] = maxGrad;

    if (verb) {
      Rcpp::Rcout << "  logLik=" << logLikCur
                  << ", objWinChange=" << objWinChange
                  << ", maxGrad=" << maxGrad
                  << ", maxStep=" << maxStep
                  << ", stepAccepted=" << step
                  << ", stepStreak=" << stepStreak << "\n";
    }

    ++niter;
  }

  // Trim objective trace to initial + completed iterations.
  NumericVector objTraceOut(niter + 1);
  for (int i = 0; i <= niter; ++i) objTraceOut[i] = objTrace[i];
  NumericVector gradTraceOut(niter + 1);
  for (int i = 0; i <= niter; ++i) gradTraceOut[i] = gradTrace[i];

  return List::create(
    _["elNew"] = elNew,
    _["beta"] = betaCur,
    _["gamma"] = gammaCur,
    _["a"] = aCur,
    _["thetaS"] = thetaSCur,
    _["logLik"] = logLikCur,
    _["objTrace"] = objTraceOut,
    _["gradTrace"] = gradTraceOut,
    _["lambdaSep"] = lambdaSepCur,
    _["niter"] = niter,
    _["LSFail"] = LSFail,
    _["converged"] = conv,
    _["convObj"] = convObj,
    _["convGrad"] = convGrad,
    _["convStep"] = convStep,
    _["stepStreak"] = stepStreak,
    _["maxGrad"] = maxGrad,
    _["objWinChange"] = objWinChange,
    _["maxStep"] = maxStep
  );
}
