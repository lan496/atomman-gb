# Formulation

## Intersection of lattices

We define a lattice for basis {math}`\mathbf{B} = (\mathbf{b}_{1}, \dots, \mathbf{b}_{n})` as
```{math}
    \mathcal{L}(\mathbf{B})
        := \left\{
            \sum_{i=1}^{n} n_{i} \mathbf{b}_{i}
            \mid
            n_{i} \in \mathbb{Z} (\forall i)
            \right\}.
```
For lattice basis {math}`\mathbf{B}`, we define the dual basis as
```{math}
  \mathbf{D}
    &:= \mathbf{B} (\mathbf{B}^{\top} \mathbf{B})^{-1} =: (\mathbf{d}_{1}, \dots, \mathbf{d}_{n}) \\
  \mathbf{D}^{\top} \mathbf{B} &= \mathbf{I}.
```
We write the dual lattice of {math}`\mathcal{L}(\mathbf{B})` as
```{math}
  \mathcal{L}(\mathbf{B})^{\ast} := \mathcal{L}(\mathbf{D}).
```

The intersection of two lattices can be obtained as follows [^footnote]:
```{math}
  \mathcal{L}(\mathbf{B}) \cap \mathcal{L}(\mathbf{B}')
    = \mathcal{L}([\mathbf{D} \mid \mathbf{D}'])^{\ast}
      \quad
      (\mathcal{L}(\mathbf{D}) := \mathcal{L}(\mathbf{B})^{\ast}, \mathcal{L}(\mathbf{D}') := \mathcal{L}(\mathbf{B}')^{\ast})
```

[^footnote]: <https://cseweb.ucsd.edu/classes/wi10/cse206a/lec2.pdf>

## Coincidence-site lattice (CSL)

### Set up

Consider two lattices, {math}`L_{\mathbf{A}}` and {math}`L_{\mathbf{A}'}`, composed of lattice vectors {math}`\mathbf{A} = (\mathbf{a}_{1}, \mathbf{a}_{2}, \mathbf{a}_{3})` and {math}`\mathbf{A}' = (\mathbf{a}_{1}', \mathbf{a}_{2}', \mathbf{a}_{3}')`,
```{math}
  L_{\mathbf{A}}
    := \left\{
        \mathbf{An}
        \mid
        \mathbf{n} \in \mathbb{Z}^{3}
        \right\}.
```
The two lattices are related by a rotation matrix {math}`\mathbf{R}`,
```{math}
  \mathbf{A}' = \mathbf{R} \mathbf{A}.
```

### Rodrigue's formula
```{math}
  \mathbf{R} &= \mathbf{I} + \mathbf{K} \sin \theta + \mathbf{K}^{2} (1 - \cos \theta) \\
  \mathbf{K}
    &= \begin{pmatrix}
          0 & -k_{z} & k_{y} \\
          k_{z} & 0 & -k_{x} \\
          -k_{y} & k_{x} & 0
       \end{pmatrix},
```
where {math}`\mathbf{k}` is a unit vector along the rotation axis.
```{math}
  \overline{\mathbf{R}}
    &:= \mathbf{A}^{-1} \mathbf{RA} \\
  \mathbf{A}'
    &= \mathbf{RA} = \mathbf{A} \overline{\mathbf{R}}
```

Rodrigues vector
```{math}
    \mathbf{\rho}^{R}
        &:= \tan \frac{\theta}{2} (k_{z} \, k_{y} \, k_{z})^{\top} \\
    (1 + (\rho^{R})^{2}) R_{ij}
        &= (1 - (\rho^{R})^{2}) \delta_{ij}
        + 2\left( \rho^{R}_{i} \rho^{R}_{j} + \sum_{k} \epsilon_{ijk} \rho^{R}_{k} \right).
```

### CSL

When the intersection of the two lattice forms a lattice, we call it as the the coincidence-site lattice (CSL),
```{math}
  \mathrm{CSL}(L_{\mathbf{A}}, L_{\mathbf{A}'})
    := L_{\mathbf{A}} \cap L_{\mathbf{A}'}
```
We denote the ratio of volumes of the primitive cell of CSL and the two lattices as {math}`\Sigma \, (\geq 1)`.

The complete pattern-shift lattice (DSCL) is defined as a maximal sublattice of {math}`L_{\mathbf{A}}` such that also contains {math}`L_{\mathbf{A}'}`.

Grimmer showed the following {cite}`Grimmer:a28679` [^footnote2]:
- (Theorem 1) If the CSL is not empty, the rotation matrix {math}`\overline{\mathbf{R}}` is rational, that is, in {math}`\mathbb{Q}^{3 \times 3}`.
- (Theorem 2) {math}`\Sigma` is the least positive integer such that {math}`\Sigma \overline{\mathbf{R}}` and {math}`\Sigma \overline{\mathbf{R}}^{-1}` are integer matrices.

[^footnote2]: These are simplified cases from {cite}`Grimmer:a28679`.

Theorem 2 is shown with a Smith normal form (SNF).
Let {math}`\mathbf{E}` be a 3x3 integer matrix and its elements are coprime, and its SNF be {math}`\mathrm{diag}(e_{1}, e_{2}, e_{3})`.
Here, we write GCD of {math}`i \times i` minors of {math}`\mathbf{E}` as {math}`d_{i}`, we obtain {math}`e_{i} = d_{i}/d_{i-1}` ({math}`d_{0} := 1`).
Since the elements of {math}`\mathbf{E}` are coprime, we get {math}`d_{1} = 1`.
Also, using the fact that elements of {math}`|\mathbf{E}| \mathbf{E}^{-1}` are 2x2 minors of {math}`\mathbf{E}`, we get {math}`d_{2}` as GCD of the elements of {math}`|\mathbf{E}| \mathbf{E}^{-1}`.
Thus, we obtain {math}`e_{1}=1`, {math}`e_{2} = ` GCD of elements of {math}`|\mathbf{E}| \mathbf{E}^{-1}`, and {math}`e_{3} = |\mathbf{E}| / e_{2}`.

## {math}`\Sigma` values for cubic system {cite}`Grimmer:a23030`

For CSL of cubic lattices, we can take Rodrigues vector as
```{math}
    \mathbf{\rho}^{R} = \frac{n}{m} [hkl]
```
where {math}`m` and {math}`n` are positive integers and {math}`\mathrm{GCD}(m, n) = \mathrm{GCD}(h, l, k) = 1`

For a given rotation axis {math}`\hat{\mathbf{\rho}} \propto [hkl]`, the procedure to enumerate rotation angle for CSLs are as follows:
1. Enumerate {math}`\Sigma` up to some maximal value with
    ```{math}
        \alpha \Sigma = m^{2} + (h^{2} + k^{2} + l^{2}) n^{2},
    ```
    where {math}`\alpha = 1, 2, 4`.
2. Calculate {math}`\theta` for each {math}`\Sigma` as
    ```{math}
        \tan \frac{\theta}{2} = \frac{n}{m} (h^{2} + k^{2} + l^{2})^{1/2}
    ```

## References

```{bibliography}
:filter: docname in docnames
```
