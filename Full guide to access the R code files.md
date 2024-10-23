# Overview of work
Fractional derivatives (FDs) are a generalization of classical integer-order derivatives and there are many ways to define a FD. This results in the existence of various FD definitions in the literature. FDs offer a powerful tool to model complex systems in various fields such as environmental sciences, physics, and biology. However, choosing a suitable FD definition is not always straightforward, as different FD definitions exhibit unique properties that can significantly impact the analysis results. Despite their theoretical significance, there is a gap in the literature regarding the data-driven selection of the most appropriate FD definition. This study addresses this gap by exploring the efficacy of six FD definitions based on simulated and real data, aiming to identify the optimal definition that aligns with the underlying process.

We have considered a simple fractional order initial value logistic growth model
$$
\left( {}_{0}^{C}D_t^\alpha \{N(t)\} \right) = r N(t) \left( 1 - \frac{N(t)}{K} \right), \quad r > 0, \, \alpha \in (0, 1]
$$
the fractional logistic growth model under six FD definitions namely, Caputo, Caputo-Fabrizio (CF), Atangana–Baleanu-Caputo (ABC), Riemann-Liouville (RL), RL with exponential decay (CF-RL), and Grünwald-Letnikov (GL). The mathematical expressions for these definitions are as follows;

<div align="center">
  <img src="equation_plot.png" alt="Flowchart for code file access" />
</div>
The motivation behind this work is 

<div align="center">
  <img src="Code%20file%20access%20guide%20flowchart.jpg" alt="Flowchart for code file access" />
  <p><em>Figure 1: This is the flowchart for accessing the code files.</em></p>
</div>
