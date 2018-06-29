\name{plot.Brq}
\alias{plot.Brq}
\title{
Produce the histogram, autocorrelation and trace plots for the objects of \code{Brq} function}
\description{
\code{plot.Brq} is an S3 method that produces histogram, autocorrelation and trace plots for the objects of \code{Brq} function}
\usage{\method{plot}{Brq}(x, plottype=c("hist", "trace", "ACF", "traceACF", "histACF","tracehist",
 "traceACFhist"),Dim=1,breaks=30,lwd=1,col1=0,col2=1,col3=1,col4=1, ...)}
\arguments{
  \item{x}{an output object of the \code{Brq} function, i.e. an S3 object of class \code{Brq}.}
  \item{plottype}{should be "hist", "trace",  "ACF", "traceACF", "histACF", "tracehist" or "traceACFhist". }
  \item{Dim}{number of parameters that should be plotted}
  \item{breaks}{a single number giving the number of cells for the histogram,}
  \item{lwd}{line width relative to the default (default=1). 2 is twice as wide.}
  \item{col1}{color for the histograms }
  \item{col2}{color for the density of the histograms}
  \item{col3}{color for the autocorrelation}
  \item{col4}{color for the trace plots}
  \item{...}{ additional arguments that are passed to the generic \code{plot} function }  
}
