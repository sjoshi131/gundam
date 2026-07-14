## Core Architecture

GUNDAM is a configuration-driven fitting framework. A fit is mainly defined through YAML or JSON configuration files, while the input events are usually read from ROOT files. The configuration specifies which datasets are loaded, how events are selected and grouped into analysis samples, which parameters are allowed to vary, how systematic effects are applied to the Monte Carlo prediction, and how the prediction is compared with data.

A typical fit begins with the datasets defined in the configuration. GUNDAM reads the input ROOT files and makes the event-level information available to the analysis. The sample definitions then determine how those events are selected, weighted, and binned into histograms. These histograms form the data distributions and the Monte Carlo prediction used in the fit.

Fit parameters are managed separately from the event definitions. They may represent normalization factors, systematic uncertainties, spline parameters, oscillation parameters, or other quantities defined by the analysis. During the fit, the fitter or sampler changes the parameter values. These changes are applied to the Monte Carlo prediction through the Propagator, which updates the relevant event weights or systematic responses and rebuilds the prediction histograms.

After the prediction is updated, the Likelihood Interface compares it with the data. The resulting likelihood value is used by the fitter or sampler to choose the next parameter point. This process is repeated until the fit reaches a best-fit result or the requested sampling procedure is complete.

The main components discussed below are the Propagator, the Parameters Manager, Sample and Event Management, and Dataset Loading. Other components, such as the Likelihood Interface and the Fitter Engine or sampler, connect these parts to the likelihood calculation and fitting procedure.

### Propagator

The Propagator applies the current fit parameter values to the Monte Carlo prediction. It is the main component that connects the parameter configuration to the samples and histograms used in the likelihood calculation.

When the fitter or sampler changes parameter values, the Propagator updates the relevant Monte Carlo events according to the configured dials and other parameter responses. A normalization parameter typically multiplies the weights of all affected events by the same factor. By contrast, a shape systematic changes event weights according to event properties, which can modify the relative contents of the histogram bins. The Propagator relies on consistent connections between samples, parameters, systematic responses, and event weights. Defining a parameter is not enough by itself; the parameter must be connected to the intended effect on the Monte Carlo prediction. If this connection is missing or incorrect, the parameter will not influence the likelihood in the expected way.

After the Propagator updates the prediction histograms, those histograms are passed to the Likelihood Interface and compared with data. The resulting likelihood value is then used by the fitter or sampler to continue the fit.
