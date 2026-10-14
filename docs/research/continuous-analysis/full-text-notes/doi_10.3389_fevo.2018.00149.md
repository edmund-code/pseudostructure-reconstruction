# doi_10.3389_fevo.2018.00149 - Simpson 2018, Frontiers in Ecology and Evolution 6:149: Modelling Palaeoecological Time Series Using Generalised Additive Models
Read status: FULL TEXT READ (publisher PDF, cached dl/doi_10.3389_fevo.2018.00149.txt; title matches)
Pages/sections read: all 21 pages: abstract, sections 1-5 (introduction, regression models, GAM bases, smoothness
selection, fitting, confidence/simultaneous intervals, periods of change, autocorrelation, GP smooths, adaptive
smoothing, age-model uncertainty, multivariate data, conclusions), all figure legends (Fig. 1-15) and Appendix 1
(simultaneous intervals). The Supplementary Material (an annotated R code document) is not in the file and was not
read; it reproduces the analyses rather than adding methods.

## What the paper actually does (from the full text)
- Model. A tutorial on GAMs for ONE irregularly spaced time series at a time: y_t = b0 + f(x_t) + e_t, f a penalised
  spline (cubic regression or thin plate regression spline), smoothness chosen by REML (recommended over GCV), basis
  size checked with gam.check (Sections 3-4). Extensions: CAR(1) residual autocorrelation via mgcv::gamm (Sec. 4.1,
  4.5), observation weights for heteroscedasticity (Sec. 4.2), Gaussian-process smooths with Matern or power
  exponential correlation (Sec. 4.6), adaptive smoothers with a smoothly varying lambda (Sec. 4.7). GLM families
  (Poisson, negative binomial, beta) are mentioned as possible (Sec. 3); both examples are Gaussian.
- Uncertainty. Bayesian credible intervals with across-the-function coverage (Nychka; Marra & Wood); posterior
  simulation of coefficients from MVN to draw plausible trends (Fig. 8); simultaneous intervals by the Ruppert et al.
  (2003) simulation: draw from N(0, Vb), take the 1 - alpha quantile of the maximum absolute standardised deviation
  over a grid (Sec. 4.3, Appendix 1). Simultaneous intervals were about 50-75% wider than across-the-function ones in
  the examples (critical values 3.07 and 3.43).
- Periods of change. First derivative by finite differences of the fitted trend; pointwise and simultaneous intervals
  for the derivative; "periods of significant change" are where the simultaneous interval on the first derivative
  excludes zero (Sec. 4.4, Fig. 10).
- Coordinate uncertainty (Sec. 4.8, Fig. 15). The time axis is itself estimated (an age-depth model). The author fits
  a monotone age-depth model (scam), draws 100 age models from its posterior, refits the trend GAM on each, then
  takes 50 posterior draws per fit (5,000 trends). The band widens by about 50% and the timing of the decline is
  reported as an interval (~1935-1945). This propagates uncertainty in the axis into the trend and into the timing of
  change.
- Null / calibration. No permutation or replicate-based null. Inference is model-based (REML, Bayesian covariance).
  The author warns that ignoring residual dependence gives narrower intervals on the trend and its derivative and so
  "anti-conservative identification of periods of change" (Sec. 4.5).
- Replicates. None. Each example is a single series (48 points; ~3,000-year record). No replicate random effects,
  no GAMM with sample-level terms; gamm() is used only for the CAR(1) residual process.
- Validation. Two worked examples only (Small Water d15N, Braya-So alkenones); no simulation study.
- Plots. Basis-function plots (Fig. 4); effect of lambda (Fig. 5); trend with across-the-function bands (Fig. 6);
  GCV/REML profiles (Fig. 7); posterior draws (Fig. 8); simultaneous vs across-the-function bands (Fig. 9); first
  derivative with simultaneous band (Fig. 10); estimated CAR(1) correlation (Fig. 11); GP smooths (Fig. 13); adaptive
  vs TPRS vs GP (Fig. 14); age-model uncertainty propagation (Fig. 15).
- Software. R with mgcv 1.8.24 (and scam for the monotone age model); annotated R code in the supplement. The paper
  does not mention gratia.
- Limitations the author states. A single smoothness parameter cannot represent step-like change: "If we anticipate
  abrupt change or step-like responses ... the method is ill-suited" (Sec. 4.7); adaptive smoothers help only partly
  and need much data. A wiggly trend and residual autocorrelation can be unidentifiable (Sec. 4.5). Multivariate
  (many-species) GAMs are not directly available (Sec. 4.9).
- Limitations for OUR design. (1) No replicates: bands from a GAM fitted to all units of a specimen would treat units
  as independent, which is exactly the anti-conservative case the author warns about; bands must come from specimen
  variation. (2) A penalised smooth with one lambda smears a step into a ramp, so a derivative-based "change density"
  will be biased toward spread-out (continuous) change at sharp boundaries. (3) One curve at a time; simultaneous
  bands control error along one curve, not across thousands of genes.

## Review check
- continuous-analysis-literature.md:158 - "GAM derivatives with simultaneous bands (Simpson 2018, doi ...)" -
  VERDICT: correct.
- continuous-analysis-literature.md:161-163 (Rec 5 "What it is": count genes whose derivative band excludes zero; ask
  whether change concentrates at the transitions) - VERDICT: misses an important point. Simpson states that a GAM
  with one smoothness parameter is ill-suited to step-like change, and that ignoring dependence makes periods of
  change anti-conservative. So a change-density track from per-gene smooths will (a) spread sharp boundary changes
  over neighbouring positions, biasing toward "continuous", and (b) over-call if bands come from unit-level fits.
  Replacement (add to Rec 5): "Derive bands from specimen-level variation, not from unit-level fits. Penalised smooths
  blur steps; compare the track with an explicit step model (Rec 1) or use adaptive smoothing before concluding that
  change is spread along the PT."
- a-gene-dynamics.md:51 - "Simpson takes the first derivative of a GAM smooth by posterior simulation, with
  simultaneous intervals. 'Periods of significant change' are where the simultaneous interval on the derivative
  excludes zero." - VERDICT: mostly correct; the derivative itself is computed by finite differences, and the
  simultaneous critical value comes from simulation (Ruppert et al. 2003). Replacement: "Simpson computes the first
  derivative of a GAM smooth by finite differences, with a simultaneous interval from simulation of the coefficient
  covariance. 'Periods of significant change' are where that interval excludes zero."
  "Replicates? Simpson: yes, via GAMM random effects" - VERDICT: wrong. The paper fits single series and has no
  replicates; gamm() is used only for CAR(1) residual autocorrelation. Replacement: "Simpson: no (single series;
  residual autocorrelation via CAR(1) only)."
  Software "R mgcv/gratia" - VERDICT: partly unsupported; the paper uses mgcv (and scam); gratia is not mentioned.
  Replacement: "R mgcv (scam for the monotone age model); annotated R code in the supplement."
  "controls error across the whole curve" - VERDICT: correct (for one curve).
- a-gene-dynamics.md:79-84 (idea 3: derivative maps and change-density track, Simpson as source) - VERDICT: misses an
  important point; same two caveats as Rec 5 above (smoothing blurs steps; bands from units are anti-conservative).
- a-gene-dynamics.md:74-75 (idea 2: intervals "by posterior simulation of the NB-GAM coefficients") - VERDICT: correct
  as a technique (Sec. 4.3); but those intervals reflect only within-fit uncertainty, not specimen variation.
- a-gene-dynamics.md:123 - "Propagating coordinate uncertainty into onset and peak intervals ... I found no method
  that carries it into switch- or peak-position intervals." - VERDICT: misses an important point. Simpson 2018
  (Sec. 4.8, Fig. 15) already does this for an estimated time axis: draw axes from the age-model posterior, refit the
  GAM on each, pool posterior draws, and report the timing of change as an interval. It is not genomics, but it is
  the same problem. Replacement: "Outside genomics, Simpson 2018 (Sec. 4.8) propagates uncertainty in an estimated
  axis (an age-depth model) by refitting the GAM over 100 draws of the axis and pooling posterior draws; the timing of
  change becomes an interval. We found no genomics method that does this for onset or peak positions."
- b-continuous-pathways.md:122 and c-spatial-gradients.md:83 (same gap: no method propagates axis uncertainty) -
  VERDICT: misses an important point; cite Simpson 2018 Sec. 4.8 as the non-genomic precedent and recipe.

## Fit to our design (2-4 lines, revised after reading the full text)
A clear recipe for derivative-based "where does it change" with simultaneous bands, and a ready recipe for carrying
axis uncertainty into trend and timing (refit over draws of the axis). It has no replicates and assumes smooth change,
so for us the bands must come from specimen-level variation, and step-like boundaries need an explicit step or
adaptive model rather than one penalised smooth.
