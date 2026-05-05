---
layout: default
title: FAQs 
next_page: "LineageandLegacy.html"
---

<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>FAQ Accordion</title>

  <style>
    body {
      font-family: Arial, sans-serif;
      max-width: 800px;
      margin: 40px auto;
      padding: 0 20px;
    }

    h2 {
      margin-top: 40px;
    }

    .accordion {
      background: #f5f5f5;
      cursor: pointer;
      padding: 14px;
      width: 100%;
      border: none;
      text-align: left;
      font-size: 16px;
      font-weight: bold;
      margin-top: 8px;
      border-radius: 6px;
      transition: background 0.2s ease;
    }

    .accordion:hover {
      background: #e8e8e8;
    }

    .accordion.active {
      background: #e0f0ff;
    }

    .panel {
      display: none;
      padding: 12px 14px;
      background: #ffffff;
      border: 1px solid #ddd;
      border-top: none;
      border-radius: 0 0 6px 6px;
    }

    .panel p {
      margin: 0;
    }
  </style>
</head>

<body>

<h1>FAQs</h1>

<h2>Gundam Fit & Spline Behavior</h2>

<button class="accordion">Is a non-converging Gundam fit a sign of a software problem?</button>
<div class="panel">
  <p>
    Not usually. Gundam uses MINUIT to minimize the negative log-likelihood (NLL)
    function and find the best-fit parameters. For this process to work reliably,
    the NLL surface must be smooth, well-behaved, and contain a clear global minimum.
    Attributes related to the model, parameterization, or input
    uncertainties/correlations can create a likelihood surface that is difficult
    to optimize. Thus, non-convergence is often a diagnostic signal that the model
    or inputs need refinement, rather than a failure of Gundam.
  </p>
</div>

<button class="accordion">How does Gundam handle spline extrapolation?</button>
<div class="panel">
  <p>
    Gundam features multiple spline interpolation methods, and they have different
    approaches for when the spline extends beyond the boundaries. A
    <strong>not-a-knot</strong> spline performs cubic extrapolation by continuing
    the cubic polynomial defined by the first/last two points of the dataset. This
    maintains agreement with splines generated using ROOT's <code>TSpline3</code>
    class. A <strong>Catmull-Rom</strong> spline extrapolates linearly beyond the
    defined knots.
  </p>
</div>

  <script>
    const accordions = document.querySelectorAll(".accordion");

    accordions.forEach((btn) => {
      btn.addEventListener("click", function () {

        // Close all panels
        accordions.forEach((item) => {
          item.classList.remove("active");
          item.nextElementSibling.style.display = "none";
        });

        // Toggle current
        const panel = this.nextElementSibling;

        if (panel.style.display === "block") {
          panel.style.display = "none";
        } else {
          this.classList.add("active");
          panel.style.display = "block";
        }
      });
    });
  </script>

</body>
</html>
