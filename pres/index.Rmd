---
title       : Causal Spatial Models
subtitle    : Joint work with Noel Cressie
author      : Andrew Zammit Mangion
job         : Statistical Computing Scientist, University of Wollongong
framework   : io2012        # {io2012, html5slides, shower, dzslides, ...}
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : tomorrow      # 
widgets     : [mathjax]            # {mathjax, quiz, bootstrap}
mode        : selfcontained # {standalone, draft}
knit        : slidify::knit2slides
---

## Introduction
type:exclaim

$$
\newcommand{\Deltab} {\boldsymbol{\Delta}}
\newcommand{\intd} {\mathrm{d}}                 
\newcommand{\Bmat} {\textbf{B}}                 
\newcommand{\Cmat} {\textbf{C}}                 
\newcommand{\cmat} {\textbf{c}}                 
\newcommand{\Imat} {\textbf{I}}                 
\newcommand{\bvec} {\textbf{b}}                 
\newcommand{\svec} {\textbf{s}}                 
\newcommand{\uvec} {\textbf{u}}                 
\newcommand{\omegab} {\boldsymbol {\omega}}     
\newcommand{\s}{\mathbf{s}}                     
\newcommand{\h}{\mathbf{h}}                     
\renewcommand{\b}{\mathbf{b}}                   
\newcommand{\z}{\mathbf{z}}                     
\renewcommand{\v}{\mathbf{v}}                   
\renewcommand{\u}{\mathbf{u}}                   
\newcommand{\w}{\mathbf{w}}                     
\renewcommand{\d}{\mathrm{d}}                   
\newcommand{\Z}{\mathbf{Z}}                     
\renewcommand{\x}{\mathbf{x}}                   
\newcommand{\Y}{\mathbf{Y}}                     
\newcommand{\Yvec}{\mathbf{Y}}                  
\newcommand{\Zvec}{\mathbf{Z}}                  
\newcommand{\epsilonb}{\boldsymbol{\varepsilon}}
\newcommand{\bI}{\mathbf{I}}                    
\newcommand{\bB}{\mathbf{B}}                    
\newcommand{\bbeta}{\boldsymbol{\beta}}         
\newcommand{\bzero}{\boldsymbol{0}}             
\newcommand{\bSigma}{\bm{\Sigma}}               
\newcommand{\E}{E}                              
\newcommand{\cov}{\mathrm{cov}}                 
\newcommand{\var}{\mathrm{var}}                 
\newcommand{\tr}{\mathrm{tr}}                   
\newcommand{\diag}{\mathrm{diag}}               
\newcommand{\vect}{\mathrm{vec}}                
\newcommand{\Gau}{\mathrm{Gau}}                 
\newcommand{\RR}{\mathbb{R}}     
$$

- Univariate spatial models

- Multivariate spatial models
  
  - Two or more interacting variates
  - We can learn about by one variate by observing the other



## Examples



$$V_t(S_t) = \max_{x_t \in \chi_t} \left(C(S_t, x_t) + 
            \gamma \sum_{s^{\prime} \in \mathcal{S}} \mathbb{P}(s^{\prime} | S_t^n, x_t) V_{t+1}^{n-1} s^{\prime} \right)$$
            
$$
\newcommand{\Yvec}{\mathbf{Y}}
C_{2|1} = \int_D \Yvec b(s,v)C_{11}(v,w)dwdu
$$

$$ \Yvec = 2 $$

1. Edit YAML front matter
2. Write using R Markdown
3. Use an empty line followed by three dashes to separate slides!






