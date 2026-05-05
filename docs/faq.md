---
layout: default
title: FAQs 
next_page: "LineageandLegacy.html"
---

<style>
.accordion {
  background: #f1f1f1;
  cursor: pointer;
  padding: 12px;
  border: none;
  width: 100%;
  text-align: left;
  font-weight: bold;
  margin-top: 5px;
}

.panel {
  display: none;
  padding: 10px;
  background: white;
  border: 1px solid #ddd;
}
</style>

<details>
  <summary><strong>Is a non-converging Gundam fit a sign of a software problem?</strong></summary>
  <p>Not usually. Gundam uses MINUIT to minimize the negative log-likelihood (NLL) function and find the best-fit parameters. For this process to work reliably, the NLL surface must be smooth, well-behaved, and contain a clear global minimum. Attributes related to the model, parameterization or input uncertainties/correlations can create a likelihood surface that is difficult to optimize. Thus, non-convergence is often a diagnostic signal that the model or inputs need refinement, rather than a failure of Gundam.</p>
    <summary><strong>How does Gundam handle spline extrapolation?</strong></summary>
  <p>Gundam features multiple spline interpolation methods, and they have different approaches for when the spline extends beyond the boundaries. A 'not-a-knot' spline performs cubic extrapolation by continuing the cubic polynomial defined by the first/last two points of the dataset. This is to maintain agreement with the splines generated using the ROOT's TSpline3 class. While a 'Catmull-Rom' spline extrapolates linearly beyond the defined knots.</p>
</details>
