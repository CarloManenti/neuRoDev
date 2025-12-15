neuRoDev is an R package that provides analysis tools to investigate the 
integrated transcriptomic references of cortex development, described in 
[Zonca&, Bot&, and Davila-Velderrain, 2025].

All the datasets needed as input are stored in the [Figshare database].
Functions in the package can install the main reference networks automatically 
(see [Network exploration]).

# Analysis tools

neuRoDev contains functions to perform:

- Reference networks exploration ([Network exploration])
- Expression enrichment analyses ([Analysis tools])
- Bulk RNAseq datasets mapping ([Mapping bulkRNAseq])
- Single-cell RNAseq datasets mapping ([Mapping scRNAseq])

In the [Tutorial], we provide an extensive description and the code to perform
all the mentioned analyses.

We also provide an interactive Shiny app to visualize eTraces on the [corticogenesis], [neurogenesis], and [gliogenesis] networks.

# Basic installation

To install `neuRoDev` from GitHub:

```{r}
install.packages("devtools")

devtools::install_github("https://github.com/davilavelderrainlab/neuRoDev", dependencies = TRUE)
```

# Bug report

Please use the [issues] to submit bug reports.

# Reference

If you use `neuRoDev` in your work, please cite

> **NeuRoDev resolves lifelong temporal and cellular variation in human cortical gene expression**
>
> Asia Zonca&, Erik Bot& & José Davila-Velderrain
>
> _Journal_ Date. doi: [doi](https://github.com/davilavelderrainlab/neuRoDev).

[Zonca&, Bot&, and Davila-Velderrain, 2025]: https://github.com/davilavelderrainlab/neuRoDev
[Figshare database]: https://github.com/davilavelderrainlab/neuRoDev
[Network exploration]: https://github.com/davilavelderrainlab/neuRoDev
[Analysis tools]: https://github.com/davilavelderrainlab/neuRoDev
[Mapping scRNAseq]: https://github.com/davilavelderrainlab/neuRoDev
[Mapping bulkRNAseq]: https://github.com/davilavelderrainlab/neuRoDev
[corticogenesis]: https://erikbot.shinyapps.io/etraceshinycortico/
[neurogenesis]: https://erikbot.shinyapps.io/etraceshinyneuro/
[gliogenesis]: https://erikbot.shinyapps.io/etraceshinyglio/
[Tutorial]: https://github.com/davilavelderrainlab/neuRoDev
[article]: https://github.com/davilavelderrainlab/neuRoDev
[issues]: https://github.com/davilavelderrainlab/neuRoDev/issues
