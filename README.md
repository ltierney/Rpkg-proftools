<!-- badges: start -->
[![Travis build status](https://travis-ci.org/ltierney/Rpkg-proftools.svg?branch=master)](https://travis-ci.org/ltierney/Rpkg-proftools)
<!-- badges: end -->

# PROFILE OUTPUT PROCESSING TOOLS FOR R

This package provides some simple tools for examining Rprof output
and, in particular, extracting and viewing call graph information.
Call graph information, including which direct calls where observed
and how much time was spent in these calls, can be very useful in
identifying performance bottlenecks.

One important caution: because of lazy evaluation a nested call
`f(g(x))` will appear on the profile call stack as if `g` had been
called by `f` or one of `f`'s callees, because it is the point at
which the value of `g(x)` is first needed that triggers the
evaluation.


## EXPORTED FUNCTIONS

The package exports these functions:

- `readProfileData` reads the data in the file produced by `Rprof`
  into a data structure used by the other functions in the package.
  The format of the data structure is subject to change.

- `flatProfile` is similar to `summaryRprof`.  It returns either a
   matrix with output analogous to `gprof`'s flat profile or a matrix
   like the `by.total` component returned by `summaryRprof`; which is
   returned depends on the value of an optional second argument.

- `printProfileCallGraph` produces a printed representation of the
  call graph.  It is analogous to the call graph produced by `gprof`
  with a few minor changes.  Reading the `gprof` manual section on the
  call graph should help understanding this output.  The output is
  similar enough to gprof `output` for the `cgprof`
  (http://mvertes.free.fr/) script to be able to produce a call graph
  via Graphviz.

- `profileCallGraph2Dot` prints out a Graphviz `.dot` file
  representing the profile graph.  Times spent in calls can be mapped
  to node and edge colors.  The resulting files can then be viewed
  with the Graphviz command line tools.

- `plotProfileCallGraph` uses the `graph` and `Rgraphviz` packages to
  produce call graph visualizations within R.  You will need to
  install these packages to use this function.

- Additional summary functions: `funSummary`, `callSummary`,
  `pathSummary`, `srcSummary`, and `hotPaths`.

- Additional functions: `filterProfileData`, `flameGraph`, `calleeTreeMap`
	`annotateSource`, and `profileExpr`.


## EXPORTED VARIABLES

The package also exports two variables:

- `plain.style`
- `google.style`

These are style specifications to be used with the call graph display
functions `plotProfileCallGraph` and `profileCallGraph2Dot`.


## A SIMPLE EXAMPLE

Collect profile information  for the examples for `glm`:

``` r
Rprof("glm.out")
example(glm)
Rprof()
pd <- readProfileData("glm.out")
```

Obtain flat profile information:

``` r
flatProfile(pd)
flatProfile(pd, FALSE)
```

Obtain hot paths information:
          
``` r
hotPaths(pd, maxdepth = 10)
```

Summaries can be obtained in a similar way:

``` r
funSummary(pd)
callSummary(pd)
pathSummary(pd)
```

Obtain a printed call graph on the standard output:

``` r
printProfileCallGraph(pd)
```

If you have the cgprof script and the Graphviz command line tools
available on a UNIX-like system, then you can save the printed graph
to a file,

``` r
printProfileCallGraph(pd, "glm.graph")
```

and either use

``` shell
cgprof -TX glm.graph
```

to display the graph in the interactive graph viewer `dotty`, or use

``` shell
cgprof -Tps glm.graph > glm.ps
gv glm.ps
```

to create a PostScript version of the call graph and display it with
`gv`.

Instead of using the printed graph and `cgprof` you can create a
Graphviz `.dot` file representation of the call graph with

``` r
profileCallGraph2Dot(pd, filename = "glm.dot", score = "total")
```

and view the graph interactively with `dotty` using

``` shell
dotty glm.dot
```

or as a postscript file with

``` shell
dot -Tps glm.dot > glm.ps
gv glm.ps
```

You can also write the profile data to a `callgrind` file to use with 
`kcachegrind` or `qcachegrind`

``` r
writeCallgrindFile(pd, file = "Rprof.cg")
```
          
If you have the packages `graph` and `Rgraphviz` from Bioconductor
installed, then you can view the call graph within R using

``` r
plotProfileCallGraph(pd, score = "total")

```

Both `plotProfileCallGraph` and `profileCallGraph2Dot` accept many
parameters for adjusting features of the display. You can specify
these parameters individually or with a single style parameter.  For
example,

``` r
plotProfileCallGraph(pd, style = google.style)
```

displays the call graph in a style similar to the one used by the
`pprof` tool in the Google Performance Tools suite.

Similarly, you can plot a flame graph and callee tree map using

``` r
flameGraph(pd)
calleeTreeMap(pd)
```

Finally, you can filter the profile data by selecting or dropping 
certain functions. For example, 

``` r
filteredPD <- filterProfileData(pd, select = "anova", focus = TRUE)
```

Now you can use `filteredPD` in you calls to summaries functions or 
plots, for example

``` r
hotPaths(filteredPD, maxdepth = 10)
flameGraph(filteredPD)
```


## OPEN ISSUES

My intention was to handle cycles roughly the same way that `gprof`
does.  I am not completely sure that I have managed to do this; I am
also not completely sure this is the best approach.

The graphs produced by `cgprof` and by `plotProfileGraph` and friends
when `mergeEdges` is false differ a bit.  I think this is due to the
heuristics of `cgprof` not handling cycle entries ideally and that the
`plotProfileGraph` graphs are actually closer to what is wanted.  When
`mergeEdges` is true the resulting graphs are DAGs, which simplifies
interpretation, but at the cost of lumping all cycle members together.

`gprof` provides options for pruning graph printouts by omitting
specified nodes.  It may be useful to allow this here as well.
