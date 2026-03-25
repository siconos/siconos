# Notes Loris

## LCP

A Linear Complementariy Problem is the following problem:
$$
w = Mz + q, \, 0 \leq w \perp z \geq 0
$$

where the unknown is $z \in \mathbb{R}^n$ ($M \in \mathbb{R}^{n \times n}$ and $q \in \mathbb{R}^n$ define the problem).
### Solvers implemented

#### Basic lcp_pgs

What I could do:

- Clever complete error computation: I saw that the error computation `lcp_compute_error` takes as long as the inner loop (one iteration). First $z$ was updated:

$$
z^{n+1}_i = - \frac{1}{m_{i,i}} \left( q_i + \sum_{j=1}^{i-1} m_{i,j} z^{n+1}_j + \sum_{j=i+1}^n m_{i,j} z^{n}_j \right)
$$

then the error was computed: 

$$
w^{n+1} = Mz^{n+1} + q, \varepsilon = \left( \sum_{i=1}^n (z^{n+1}_i - \max(0, z^{n+1}_i - w^{n+1}_i))^2 \right)^{1/2}
$$

We can see that when we compute $z^{n+1}_i$ we actually compute a part of $w^{n+1}_i$. Let's define:
$$
L^n_i = q_i + \sum_{j=1}^{i-1} m_{i,j} z^{n+1}_j \text{ and } R^n_i = \sum_{j=i+1}^n m_{i,j} z^{n}_j
$$

Then:
$$
w^{n+1}_i = q_i + \sum_{j=1}^n m_{i,j} z^{n+1}_j = L^n_i + m_{i,i} z^{n+1}_i + R^{n+1}_i
$$

This means that when we compute $z^{n+1}_i$ we can store $L^n_i$ and when we compute $z^{n+2}_i$ we can compute $R^{n+1}_i$ so that we get $w^{n+1}_i$ without any additional computation.

#### "Parallel dense gauss-seidel algorithm on many-core processors"

From the following paper: "Hadrien Courtecuisse and Jérémie Allard. 2009. Parallel dense gauss-seidel algorithm on many-core processors. In 2009 11th IEEE International Conference on High Performance Computing and Communications. IEEE, 139–147 https://dl.acm.org/doi/10.1109/HPCC.2009.51"

The matrix $M$ is divided into $n_t$ blocks $B_i$ of size $s_i$. The update looks like:
$$
\forall i \in B_m, z^{n+1}_i = \max\left(0, - \frac{1}{m_{i,i}} \left( q_i + \sum_{j \in B_1}^{B_{m-1}} m_{i,j} z^{n+1}_j + \sum_{j \in B_m, j \neq i}^{B_{n_t}} m_{i,j} z^{n}_j \right) \right)
$$

#### Parallel Gauss-Seidel with Graph Coloring

Reference: https://erkaman.github.io/posts/gauss_seidel_graph_coloring.html

This algorithm first identifies "independent" lines in $M$, such that we can update several components of $z$ in parallel without causing the algorithm to do more iterations. Two lines $i$ and $j$ of $M$ are independent if $M_{i,j} = 0$. Groups of independent lines are formed using a graph coloring algorithm, which defines $p$ color sets $C_0, ..., C_{p-1}$ such that if two indices $i, j$ belong to a set $C_k$ then $M_{i, j} = 0$. Let $k$ be a color. Let
$$
\mathcal{L}_k = \bigcup_{a = 0}^{k - 1} C_a \text{ and } \mathcal{R}_k = \bigcup_{a = k + 1}^{p - 1} C_a
$$

be the sets of indices whose color index is lower / higher than $k$. Then, for each color $k$, we can do in parallel for the lines $i$ belonging to this color set:

$$
\forall i \in C_k, z^{n+1}_i = - \frac{1}{m_{i,i}} \left( q_i + \sum_{j \in \mathcal{L}_k } m_{i,j} z^{n+1}_j + \sum_{j \in \mathcal{R}_k} m_{i,j} z^{n}_j \right)
$$

It is still possible to compute the error during the iterations, but it's a bit more complicated. Indeed, we can no longer separate a line between what's before and after the diagonal to compute the two sums (as we did before), we now have to know which index goes to which sum (i.e. which component of $z$ has already been updated).

In order to do this, we first create an array of size $n$ called $\texttt{inv\_permutation}$, such that:
$$
\texttt{inv\_permutation[0, ..., n\_0 - 1]} = C_0 \\
\texttt{inv\_permutation[n\_0, ..., n\_0 + n\_1 - 1]} = C_1 \\
\texttt{inv\_permutation[n\_0 + n\_1, ..., n\_0 + n\_1 + n\_2 - 1]} = C_2 \\
...
$$

so that when going through $\texttt{inv\_permutation}$, we go through the indices of each color, in color order. Let us define the sum of sizes:
$$
\forall k \in \llbracket 0, p \rrbracket, s_k = \sum_{j = 0}^{k - 1} n_j
$$

Then, when suppose we are dealing with line $i$, whose color is $k$. We go through this line, and we have to find,
for each column $j$, if $m_{i,j} z_j$ goes to the left sum or the right sum. So we have to find if $j$ belongs to
$\mathcal{L}_k$ or $\mathcal{R}_k$. 

$$
\begin{align*}
j \in \mathcal{L}_k &\Longleftrightarrow j \in  \texttt{inv\_permutation[0, ..., s\_k - 1]} \\
&\Longleftrightarrow \exists i \in \llbracket 0, s_k \llbracket, \texttt{inv\_permutation[i]} = j \\
&\Longleftrightarrow \exists i \in \llbracket 0, s_k \llbracket, i = \texttt{permutation[j]} \\
&\Longleftrightarrow \texttt{permutation[j]} < s_k \\
\end{align*}
$$

where $\texttt{permutation}$ is the inverse permutation of $\texttt{inv\_permutation}$, which satisfies $\texttt{permutation[inv\_permutation[i]] = i}$. We can do the same for $\mathcal{R}_k$:

$$
\begin{align*}
j \in \mathcal{R}_k &\Longleftrightarrow j \in  \texttt{inv\_permutation[s\_\{k+1\}, ..., s\_p - 1]} \\
&\Longleftrightarrow \exists i \in \llbracket s_{k+1}, s_p \llbracket, \texttt{inv\_permutation[i]} = j \\
&\Longleftrightarrow \exists i \in \llbracket s_{k+1}, s_p \llbracket, i = \texttt{permutation[j]} \\
&\Longleftrightarrow \texttt{permutation[j]} \geq s_{k+1} \\
\end{align*}
$$

So we can easily compute the two sums by knowing the $\texttt{permutation}$ array. This is implemented in $\texttt{NM\_row\_prod\_graph}$.

#### Note

I have another version of this method, in which I permutate all the objects $z$, $w$, $M$ before starting the itreations, so that I don't have to deal with permutations inside the solver. The problem was that the code is a bit more complicated, and it takes more memory if we do not permutate arrays in place, without being much faster. The only advantage was that it seems better in the dense case, since we would not have to test $\texttt{permutation[j]}$ when computing the sums, we could directly use $\texttt{cblas\_ddot}$. But this algorithm is not meant to be used on dense matrices which is why I removed this version.

#### Remark on clever error computation

In $\texttt{lcp\_pgs\_parallel}$ and in $\texttt{lcp\_pgs\_graph\_permut}$ I found a way to compute the complementarity error efficiently. I will explain why these solvers always do 2 more iterations compared to the non-optimized versions $\texttt{lcp\_pgs\_multiline}$ (code in my test directory) and $\texttt{lcp\_pgs\_graph}$.

At iteration $k$, these solvers know $z_k$, $w_{k-1}$ and $\varepsilon_{k-1}$. Suppose the non-optimized solver stopped at iteration $N$, meaning $\varepsilon_N < \text{tol}$. The optimized algorithms need one more iteration to compute this error, so at iteration $N + 1$ they will know $z_{N+1}$, $w_N$ and $\varepsilon_N$. But we now need to synchronize all variables, i.e. we need to compute $w_{N+1}$ and $\varepsilon_{N+1}$, which is why we do one last iteration without updating $z$.

Two things:
- I could put the last iteration (the one to synchronize) out of the parallel loop, to make the code more readable. We would just call the regular $\texttt{lcp\_compute\_error}$ to compute $w_{N+1}$ and $\varepsilon_{N+1}$ from $z_{N+1}$.
- I heard Vincent say that $z$ does not matter (not sure though), only $w$ does. If it's the case, then we don't need that last synchronizing iteration.

### Bugs

#### [CORRECTED] CSparse matrix initialization in parallel

If the NumericsMatrix M is sparse, NM_block_prod will do:
```c
CSparseMatrix* S = NULL;
if (A->storageType == NM_SPARSE) {
    if (A->matrix2->origin == NSM_CSR) {
        S = NM_csr(A);
    } 
    else {
        S = NM_csc_trans(A);
    }
}
```

before computing the products. `lcp_pgs_parallel` calls `NM_bloc_prod` in parallel, so the above code gets executed multiple times in parallel, causing "double free or corruption" errors.
To prevent this, I execute these lines a first time before the parallel section, so that the sparse matrix is initialized without any issue.

#### [CORRECTED] Tests with one thread

$\texttt{lcp\_pgs\_parallel}$ was not passing a test (#106) when running with one thread. This is because when there is one thread, this method is Jacobi method, which does not converge with the test provided:
$$
M = \begin{bmatrix}
1 & 1 \\ 1 & 1
\end{bmatrix},
q = \begin{bmatrix}
-1 \\ -1
\end{bmatrix},
z_0 = \begin{bmatrix}
0 \\ 0
\end{bmatrix}
$$

Indeed, this will result in:
$$
z_1 = \begin{bmatrix}
1 \\ 1
\end{bmatrix},
z_2 = \begin{bmatrix}
0 \\ 0
\end{bmatrix},
z_3 = \begin{bmatrix}
1 \\ 1
\end{bmatrix},
\text{etc...}
$$

Now when $\texttt{lcp\_pgs\_parallel}$ is called with only one thread it defaults to $\texttt{lcp\_pgs}$.

### TO DO

- Performance on dahu

- Sparse block matrices

How to adapt my solvers to sparse block matrices? 

- Non-symmetric matrices

Check if it works on non-symmetric matrices. I'm not sure aobut the creation of a Petsc matrix from a NumericsMatrix (CSR format etc...)

## Friction Contact Problem

A Friction Contact problem is defined by:
- Dimension of the contact space: $d$
- Number of contacts: $n_c$
- Matrix: $M \in \mathbb{R}^{n \times n}$ where $n = d n_c$
- Vector: $q \in \mathbb{R}^n$
- Vector of friction coefficients: $\mu \in \mathbb{R}^{n_c}$

Let $r \in \mathbb{R}^n$ be the reaction.

### `fc2d_nsgs`

#### Questions

- why is there a "if" to compute the error?
- in the case of a SBM matrix, will all the blocks be 2x2? Because if it's the case, I might be able to speedup the construction of the PETSC matrix
- in LCP NM_row_prod_graph, I ignroe the 0 in the matrix $M$. Should I do it for FC2D ? Is it really faster ?
- if we have a large number of independent rows (i.e. big color sets), is there a way to use BLAS to do a parallel matrix product? How to deal with the columns to ignore?
- can I at least do a matrix vector product instead of two dto products in NM_row_prod_no_diag2_parallel

In 2D, $d = 2$.

Two important steps: `fc2d_nsgs_buildLocalProblem` and `fc2d_nsgs_local_solve`

Let $0 \leq c \leq n_c - 1$ be a contact number.

#### First step: `fc2d_nsgs_buildLocalProblem`

Creates a dense LCP of size (2, 2) with:
$$M_{\text{Loc}} = \begin{pmatrix} M_{2c, 2c} & M_{2c, 2c+1} \\ M_{2c+1, 2c} & M_{2c+1, 2c+1} \end{pmatrix}$$

$$q_{\text{Loc}} = \begin{pmatrix} q_{2c} \\ q_{2c+1} \end{pmatrix} + \begin{pmatrix} M_{2c, 0} & \cdots & M_{2c, 2c-1} & M_{2c, 2c+2} & \cdots & M_{2c, n - 1}  \\ M_{2c+1, 0} & \cdots & M_{2c+1, 2c-1} & M_{2c+1, 2c+2} & \cdots & M_{2c+1, n - 1} \end{pmatrix} \begin{pmatrix} r_0 \\ \vdots \\ r_{2c-1} \\ r_{2c+2} \\ \vdots \\ r_{n-1} \end{pmatrix}$$

(same as `NM_row_prod_nodiag1x1` in `lcp_pgs`)

#### Second step: `fc2d_nsgs_local_solve`

Parameters: $M_{\text{Loc}}, |M_{\text{Loc}}|, q_{\text{Loc}}, \mu_c$

Output: $r_{\text{Loc}} = \begin{pmatrix} r_{2c} \\ r_{2c+1} \end{pmatrix}$

Solving:
$$
\begin{align*}
u_{\text{Loc}} = M_{\text{Loc}} r_{\text{Loc}} + q_{\text{Loc}} \\
0 \leq u_{\text{Loc},N} \perp r_{\text{Loc}, N} \geq 0 \\
\end{align*}
$$

### fc2d_nsgs_dense

Trying to understand what `fc2d_nsgs_dense` does:
- two outputs, `reaction` and `velocity`, possibly like $z$ and $w$ from LCP?
- perhaps reaction is $r$ and velocity is $u$ ???

`vec` is the dense matrix $M$. `vec[i + j * M->size0]` $= M_{i,j}$

First step, loop on all contacts. For each contact $i$:
- $$\text{avn} = \sum_{j=0}^{2i - 1} \text{reaction}_j \text{vec}_{nj + 2i} = \sum_{j=0}^{2i - 1} r_j M_{2i, j}$$
- $$\text{avt} = \sum_{j=0}^{2i - 1} \text{reaction}_j \text{vec}_{nj + 2i + 1} = \sum_{j=0}^{2i - 1} r_j M_{2i + 1, j}$$
- $$\text{apn} = \sum_{k=2i+2}^{n - 1} \text{reaction}_k \text{vec}_{nk + 2i} = \sum_{k=2i+2}^{n - 1} r_k M_{2i, k}$$
- $$\text{apt} = \sum_{k=2i+2}^{n - 1} \text{reaction}_k \text{vec}_{nk + 2i + 1} = \sum_{k=2i+2}^{n - 1} r_k M_{2i+1, k}$$

Then:
$$
\text{zn} = -q_{2i} - \text{avn} - \text{apn} = -q_{2i} - \sum_{j=0, j \neq 2i, 2i+1}^{n-1} r_j M_{2i,j}
$$
$$
\text{zt} = -q_{2i+1} - \text{avt} - \text{apt} = -q_{2i+1} - \sum_{j=0, j \neq 2i, 2i+1}^{n-1} r_j M_{2i+1,j}
$$

So far its like a manual `NM_row_prod_no_diag1x1`, but two times, one for normal and one for tangential.

Then I think its the local solver. First if $\text{zn} \leq 0$ the local solver only does:
$$
r_{2i} = 0, r_{2i+1} = 0, u_{2i} = - \text{zn}, u_{2i+1} = - \text{zt}
$$

Else, the big part. First compute the determinant of the block:
$$
\text{det} = \begin{vmatrix} M_{2i,2i} & M_{2i, 2i+1} \\ M_{2i+1,2i} & M_{2i+1,2i+1} \end{vmatrix}
$$

If its too close to $0$ then partial pivoting is used to solve the system described below.

Else, do:
$$
r_{2i} = \frac{\text{zn} \cdot M_{2i+1,2i+1} - \text{zt} \cdot M_{2i, 2i+1}}{\text{det}} 
$$

$$
r_{2i+1} = \frac{-\text{zn} \cdot M_{2i+1,2i} + \text{zt} \cdot M_{2i, 2i}}{\text{det}} 
$$

So basically we just updated two components of $r$ by solving the 2-by-2 system:
$$
\begin{pmatrix} M_{2i,2i} & M_{2i, 2i+1} \\ M_{2i+1,2i} & M_{2i+1,2i+1} \end{pmatrix} \begin{pmatrix} r_{2i} \\ r_{2i+1} \end{pmatrix} = \begin{pmatrix} \text{zn} \\ \text{zt} \end{pmatrix}
$$

Then, check the following condition:
$$
r_{2i} \geq 0 \text{ and } \|r_{2i+1}\| - \mu_i r_{2i} \leq 0
$$

If it is true, then don't do anything. The velocities $u_{2i}$ and $u_{2i+1}$ stay at $0$.

Else compute:
$$
g_+ = M_{2i,2i} + \mu_i M_{2i, 2i+1}
$$

And:
$$
u_{2i+1} = - \text{zt} + \frac{\text{zn}}{g_+} \cdot (M_{2i,2i+1} + \mu_i M_{2i+1,2i+1})
$$
$$
r_{2i} = \frac{\text{zn}}{g_+}, r_{2i+1} = \mu_i r_{2i}
$$

Then, if $r_{2i} \geq 0$ and $u_{2i+1} \leq 0$ don't do anything. Else:
$$
g_- = M_{2i,2i} - \mu_i M_{2i, 2i+1}
$$

And:
$$
u_{2i+1} = - \text{zt} + \frac{\text{zn}}{g_-} \cdot (M_{2i,2i+1} - \mu_i M_{2i+1,2i+1})
$$
$$
r_{2i} = \frac{\text{zn}}{g_-}, r_{2i+1} = - \mu_i r_{2i}
$$