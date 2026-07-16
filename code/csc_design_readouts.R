## =====================================================================================
## Causal Synergy Curve -- shared design construction + posterior g-computation readouts.
## Used by both the full-Bayes simulation and the applied NHANES harness.
##
## Estimating surface (locked 2026-07-09):
##   eta = alpha0 + W gamma + sum_j beta_j D_j + f(m) + sum_g kappa_g Gate_g + sum_T theta_T phi_T
##   baseline block Xb = [W, D, m, m^2, Gate_g]   (weakly-informative Normal prior)
##   coalition block Xc = [ prod_{j in T} D_j : 2 <= |T| <= order ]   (regularized horseshoe)
##
## All Synergy-Curve functionals are posterior g-computation functionals of the closure surface
##   C(S,alpha) = E[Y] - E[Y^{d_{S,alpha}}]
## computed by rebuilding the design under the shifted domain profile for every posterior draw.
## Regime is read on the LINK scale (structure); magnitude on the RESPONSE scale.
## =====================================================================================

inv_cloglog <- function(x) 1 - exp(-exp(pmin(x, 30)))
.wcm <- function(X, w) if (is.null(w)) colMeans(X) else as.numeric(crossprod(w / sum(w), X))

## ---- design construction -------------------------------------------------------------
## D: n x p domain-gap matrix (binary 0/1 unfavorable, or standardized continuous gap).
## W: n x q covariate matrix (numeric / already dummy-coded); may have 0 columns.
csc_build_design <- function(D, W, order = 2) {
  D <- as.matrix(D); n <- nrow(D); p <- ncol(D)
  if (is.null(W)) W <- matrix(0, n, 0)
  W <- as.matrix(W)
  m <- rowSums(D)
  gates <- sapply(seq_len(p), function(g) D[, g] * (m - D[, g]) / max(1, p - 1))
  Xb_raw <- cbind(W, D, m, m^2, gates)
  coal <- if (order >= 2) {
    unlist(lapply(2:order, function(r) if (r <= p) combn(p, r, simplify = FALSE) else NULL),
           recursive = FALSE)
  } else {
    list()
  }
  Xc_raw <- if (length(coal))
    sapply(coal, function(Tt) apply(D[, Tt, drop = FALSE], 1, prod)) else matrix(0, n, 0)
  Xc_raw <- matrix(as.numeric(Xc_raw), nrow = n)
  std <- function(X) {
    if (ncol(X) == 0) return(list(center = numeric(0), scale = numeric(0)))
    ce <- colMeans(X); sc <- apply(X, 2, sd); sc[!is.finite(sc) | sc < 1e-8] <- 1
    list(center = ce, scale = sc)
  }
  sb <- std(Xb_raw); sc <- std(Xc_raw)
  list(p = p, q = ncol(W), order = order, coal = coal,
       Xb = sweep(sweep(Xb_raw, 2, sb$center, "-"), 2, sb$scale, "/"),
       Xc = if (ncol(Xc_raw)) sweep(sweep(Xc_raw, 2, sc$center, "-"), 2, sc$scale, "/") else Xc_raw,
       sb = sb, sc = sc, W = W)
}

## Rebuild the standardized design for an arbitrary domain profile Dp (n x p), reusing the
## standardizers estimated on the observed data (so shifted profiles are on the same scale).
csc_std_design <- function(Dp, W, meta) {
  n <- nrow(Dp); p <- meta$p
  m <- rowSums(Dp)
  gates <- sapply(seq_len(p), function(g) Dp[, g] * (m - Dp[, g]) / max(1, p - 1))
  Xb_raw <- cbind(W, Dp, m, m^2, gates)
  Xc_raw <- if (length(meta$coal))
    matrix(as.numeric(sapply(meta$coal, function(Tt) apply(Dp[, Tt, drop = FALSE], 1, prod))), nrow = n)
  else matrix(0, n, 0)
  Xb <- sweep(sweep(Xb_raw, 2, meta$sb$center, "-"), 2, meta$sb$scale, "/")
  Xc <- if (ncol(Xc_raw)) sweep(sweep(Xc_raw, 2, meta$sc$center, "-"), 2, meta$sc$scale, "/") else Xc_raw
  list(Xb = Xb, Xc = Xc)
}

## Posterior linear predictor matrix (n x M) for a domain profile, given posterior draws.
## post: list(alpha0 = M-vector, beta_b = M x Pb, theta = M x Pc).
csc_eta <- function(Dp, W, meta, post) {
  X <- csc_std_design(Dp, W, meta)
  eta <- X$Xb %*% t(post$beta_b)                       # n x M
  if (ncol(X$Xc)) eta <- eta + X$Xc %*% t(post$theta)
  sweep(eta, 2, post$alpha0, "+")
}

## ---- closure enumeration -------------------------------------------------------------
## For every subset S (2^p masks), compute the closure draw-vector on link and response scales.
## alpha: intervention intensity. For binary domains alpha<1 is a STOCHASTIC flip-to-favorable
## (Monte-Carlo averaged over mc_flips); alpha=1 is the deterministic full shift.
## resp_fun(eta, m) -> response scale (m = draw index; for Weibull uses shape[m] & horizon).
csc_closures <- function(D, W, meta, post, alpha = 1, mc_flips = 16,
                         resp_fun, w = NULL) {
  n <- nrow(D); p <- meta$p; M <- length(post$alpha0)
  eta0 <- csc_eta(D, W, meta, post)                    # n x M, natural
  r0   <- resp_fun(eta0)                               # n x M response
  masks <- 0:(2^p - 1)
  clo_link <- matrix(0, length(masks), M)              # rows indexed by mask (mask=0 -> 0)
  clo_resp <- matrix(0, length(masks), M)
  Umats <- if (alpha < 1) lapply(seq_len(mc_flips), function(.) matrix(runif(n * p), n, p)) else NULL
  for (mi in seq_along(masks)) {
    mask <- masks[mi]; if (mask == 0) next
    S <- which(bitwAnd(mask, 2^(0:(p - 1))) > 0)
    if (alpha >= 1) {
      Dp <- D; Dp[, S] <- 0
      etaS <- csc_eta(Dp, W, meta, post)
      clo_link[mi, ] <- .wcm(eta0 - etaS, w)
      clo_resp[mi, ] <- .wcm(r0 - resp_fun(etaS), w)
    } else {
      accL <- 0; accR <- 0
      for (U in Umats) {
        Dp <- D
        for (j in S) Dp[U[, j] < alpha, j] <- 0
        etaS <- csc_eta(Dp, W, meta, post)
        accL <- accL + .wcm(eta0 - etaS, w)
        accR <- accR + .wcm(r0 - resp_fun(etaS), w)
      }
      clo_link[mi, ] <- accL / mc_flips
      clo_resp[mi, ] <- accR / mc_flips
    }
  }
  list(link = clo_link, resp = clo_resp, masks = masks, p = p, M = M)
}

## ---- functionals from a closure surface ----------------------------------------------
## clo: length-(2^p) numeric vector (one posterior draw, or the truth). Returns SC, ESR,
## Shapley phi, gates B/E/G, and the endpoint pairwise second-difference vector.
.mask_of <- function(S, p) sum(2^(S - 1))
csc_functionals <- function(clo, p, u_ref = NULL) {
  C <- function(S) clo[.mask_of(S, p) + 1]
  full <- .mask_of(seq_len(p), p) + 1
  ## Synergy Curve
  SC <- c(0, sapply(1:p, function(k)
    mean(apply(combn(p, k), 2, function(S) clo[.mask_of(S, p) + 1]))))
  ESR <- if (SC[p + 1] != 0) (SC[p + 1] - p * SC[2]) / abs(SC[p + 1]) else 0
  ## Shapley (exact, p <= ~10)
  phi <- numeric(p)
  for (j in seq_len(p)) {
    oth <- setdiff(seq_len(p), j)
    for (r in 0:(p - 1)) {
      w <- factorial(r) * factorial(p - r - 1) / factorial(p)
      combs <- if (r == 0) list(integer(0)) else combn(oth, r, simplify = FALSE)
      for (cc in combs) {
        base <- if (length(cc)) .mask_of(cc, p) + 1 else 1
        withj <- .mask_of(sort(c(cc, j)), p) + 1
        phi[j] <- phi[j] + w * (clo[withj] - clo[base])
      }
    }
  }
  ## endpoint pairwise second differences (S = empty)
  pd <- c()
  for (j in 1:(p - 1)) for (k in (j + 1):p)
    pd <- c(pd, clo[.mask_of(c(j, k), p) + 1] - clo[.mask_of(j, p) + 1] - clo[.mask_of(k, p) + 1])
  list(SC = SC, ESR = ESR, phi = phi, full = clo[full], pair_delta = pd)
}

## anchored gate contrast (link scale). u_ref: per-domain unfavorable reference (mean gap among gap>0).
csc_gates <- function(D, W, meta, post, u_ref, w = NULL) {
  p <- meta$p; M <- length(post$alpha0)
  Bg <- matrix(0, p, M); Eg <- matrix(0, p, M)
  for (g in seq_len(p)) {
    oth <- setdiff(seq_len(p), g)
    Bu <- D; Bu[, g] <- u_ref[g]; Bu2 <- Bu; Bu2[, oth] <- 0
    Bg[g, ] <- .wcm(csc_eta(Bu, W, meta, post) - csc_eta(Bu2, W, meta, post), w)
    Ef <- D; Ef[, g] <- 0; Ef2 <- Ef; Ef2[, oth] <- 0
    Eg[g, ] <- .wcm(csc_eta(Ef, W, meta, post) - csc_eta(Ef2, W, meta, post), w)
  }
  G <- (Eg - Bg) / (abs(Eg) + abs(Bg) + 1e-9)
  list(B = Bg, E = Eg, G = G)
}

## regime label from an endpoint pairwise second-difference vector + region of practical equivalence.
csc_regime <- function(pair_delta, rope = 0.05) {
  pos <- any(pair_delta >  rope); neg <- any(pair_delta < -rope)
  if (!pos && !neg) "additive" else if (pos && neg) "mixed" else if (pos) "complementary" else "redundant"
}
