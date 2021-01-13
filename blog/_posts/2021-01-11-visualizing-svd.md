---
layout: post
title: Visualizing the Singular Value Decomposition (SVD)
truncated_preview: true
excerpt_separator: <!--more-->
usemathjax: true
---

<div class="message">
    Recently I published a paper utilizing the Singular Value
    Decomposition (SVD), a well-known matrix decomposition method.
    With my first blog post on my site, I wish to explain how to
    intuitively visualize the results the singular value decomposition
    provides using a simple example.
</div>

<!--more-->

## Prerequisites
This post will require some understanding of *Julia*, *matplotlib* and Linear Algebra.  
I will try to explain most terms I deem are not fundamental knowledge in those three categories.

## The Singular Value Decomposition
Many might be more familiar with the Principal Component Analysis (PCA), which is built on the strength of the SVD.  
The SVD is a matrix decomposition method, which break single matrices down into a product of matrices, which offer advantages in a range of problems.  
Three matrices are produced with the SVD, **U**, **&Sigma;**, and **V<sup>T</sup>** according to the following formula:  
$$\mathbf{A} = \mathbf{U\Sigma V^T}$$
where $\mathbf{A}$ is the original matrix.

## Two Apple Vendors

## Visualizing the SVD
