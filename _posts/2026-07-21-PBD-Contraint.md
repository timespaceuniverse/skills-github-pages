---
layout: post
title: 'Constraint of Position-Based Dynamics (PBD)'
mathjax: true
author: "helloworld1943"
date: 2026-07-21
tags: [PBD,XPBD,Constraint,MATH,PHYSICS]
---

<!-- raw before LLM format , lazy work

We give a much more detailed anylsis of typical constraints in (X)PBD.
PBD Paper: https://matthias-research.github.io/pages/publications/posBasedDyn.pdf
XPBD Paper : https://matthias-research.github.io/pages/publications/XPBD.pdf


Lets still start from the basic assumption of PBD :

with C(x) , we could define the internal force easily :
the internal lambda_force is try to make the system rebalanced (C(x)=0) from the current position which is 
like hooks' law.  C(Xcurrent + delta x) = - 1/k * Lambda_force  (Hook's law : x = -1/k * F ,negtive sign for restoring force)
Left side use Taylor expansion C(Xcurrent)+ Transpose(gradient(C))* delta X  = 1/k * Lambda_force , 
delta X = M^-1 * gradient(C) * Lambda_force * delta t * delta t   ,  here Lambda_force is a scalar,
gradient(C) * Lambda_force the force vector is in the direction of the gradient which is reasonable.


1. Distance constraint

C(x1,x2) = ||x1 -x2|| - d 

from this
delta X = M^-1 * gradient(C) * Lambda_force * delta t * delta t 

we get 
delta x1 = 1/m1 * dC(x)/dx1  * Lambda_force * delta t * delta t 
delta x2 = 1/m2 * dC(x)/dx2  * Lambda_force * delta t * delta t 

dC(x)/dx1 = - dC(x)/dx2 , so delta x1 / delta x2 =m2/m1 , 

we can see that with this kind of definition of constraint, we we really define 
is that with two particles , while keep their distance , when distance change , mass bigger move smaller
that is what our distance constraint means. 

When you put a group of particles with adjecent particle having such distance constraint , the whole 
object will behave like a soft body object. No matter what you do with the stiffness k parameter the whole object
will still adjust itself with shape changing as you could image that with different k the object kind of shakes differently .


2. Shape constraint

If we want to simulate rigid body , what we want is at each step the whole shape keeps same.
So how do we define this kind of constraint ?

we define a group of particles (Xi ,Mi) which forms a whole shape , at rest we have the mass center Xc = sum(Xi_rest*Mi)/sum(Mi)

then after each particle moving we want to keep the whole shape the same as the reset shape at any moment.

if you want to keep the whole shape want you can have for the spatial change of the whole shape is  a translate T and a rotation R 
X_i_new_target = R*(X_i_rest-Xc)+T

The new position X_i_new after moving may be at chaotic places , the question is how do we adjust each of those X_i_new to X_i_new_target.
that is to solve the best value for the R and T given X_i_new.

Intuitively T is 3d vector which is easy to get , the new Xc_new= sum(Xi_new*Mi)/sum(Mi) 
T= Xc_new - Xc

Intuitively We define The loss(R) = sum[ ||R*(X_i_rest-Xc)+T - X_i_new||^2 * Mi ] 

Ai=X_i_rest-Xc , Bi=T - X_i_new

that is we want to minimize loss(R) = sum[ Mi*||R*Ai+Bi||^2 ]

loss(R) = sum [ Mi * Transpose(R*Ai+Bi) x (R*Ai+Bi)  ]
        = sum [ Mi * ( Transpose(Ai)*Transpose(R)*R*Ai + Transpose(Ai)*Transpose(R)*Bi + Transpose(Bi)*R*Ai + Transpose(Bi)*Bi ) ]
    
Transpose(R)*R = I ,  
as Transpose(Ai)*Transpose(R)*Bi is scalar, Transpose(Ai)*Transpose(R)*Bi = Transpose(Transpose(Ai)*Transpose(R)*Bi)=  Transpose(Bi)*R*Ai

loss(R) = sum [ Mi * ( Transpose(Ai)*Ai + 2*Transpose(Bi)*R*Ai + Transpose(Bi)*Bi ) ]

that is to minimize the loss(R_star)= sum[Mi*Transpose(Bi)*R_star*Ai ]
                               = Trace ( R_star * sum[ Mi*Ai*Transpose(Bi) ] ) 

define sum[ Mi*Ai*Transpose(Bi) ] = transpose(A)

loss(R_star)= Trace(R_star*transpose(A))

polar decompose A = R*S ,  
here if det(A)!=0 , then we always get a unique decompostion. 
We need to make sure actually det(A) != 0 
In the simulation if we meet det(A)=0 we just fallback to the previous Rotation.
For det(A)<0 , you need to do some work to handle it . 

loss(R_star)= Trace(R_star*transpose(S)*transpose(R)) = Trace(R_star*s*transpose(R)) = Trace(transpose(R)*R_star*S)

here matrix transpose(R)*R_star = G is also a rotation matrix, S is sysmetric matrix which can be diagnalized S= Q * Sigma * Transpose(Q)

loss(R_star) = Trace( G * Q * Sigma * Transpose(Q)) = Trace( Transpose(Q) * G * Q * Sigma ) , 

here Q is orthornormal matrix , G is orthornormal matrix (R_star and transpose(R) are both orthornormal matrix ) , 
so Transpose(Q) * G * Q is orthornormal matrix which means each colum sum to 1 then it easily to see that 
when Transpose(Q) * G * Q = -I will make the loss(R_star) minimized.

so that is Transpose(Q) * transpose(R)*R_star * Q = -I  => transpose(R) * R_star = -I => R_star = - R


svd decompose A = Usvd*Sigma*Transpose(Vsvd) ,as Usvd , Vsvd are othornormal , Transpose(Vsvd)*Vsvd = I 
so A =   Usvd*Transpose(Vsvd)*Vsvd*Sigma*Transpose(Vsvd) 
     =  (Usvd*Transpose(Vsvd)) * (Vsvd*Sigma*Transpose(Vsvd))

so Rpolar = (Usvd*Transpose(Vsvd)) and S = (Vsvd*Sigma*Transpose(Vsvd))

Rather than using SVD decomposition which is hard we use Newton-Schulz alogrithm to calculate the Rpolar which is fast and easy to implement.

Rk+1 = 1/2*(Rk+Transpose(Rk^-1))  , where R0 = A , 

Rk= Usvd*Sigma_k*Transpose(Vsvd) ,

Transpose(Rk^-1)=Usvd*Sigma_k^-1*Transpose(Vsvd) 

Rk+1 = Usvd * [1/2 * (Sigma_k+Sigma_k^-1)] * Transpose(Vsvd)

as we can see from k to k+1 only the middel term changes , actually it changes toward I

so after several times 6-9 times Rk will be just Usvd  * Transpose(Vsvd) which is what we need



!-->

# Mathematical Analysis of Geometric Constraints in Position-Based Dynamics (PBD) and XPBD

**References:**
* Müller, M., Heidelberger, B., Hennix, M., & Ratcliff, J. (2007). *Position Based Dynamics*. [[Paper Link](https://github.io)]
* Macklin, M., Müller, M., & Chentanez, N. (2016). *XPBD: Position-Based Simulation of Compliant Constrained Dynamics*. [[Paper Link](https://github.io)]

---

## 1. Theoretical Foundation and Governing Equations

Position-Based Dynamics (PBD) and its compliant extension (XPBD) omit the velocity layer during constraint projection to directly manipulate particle positions. This formulation eliminates the overshooting instabilities typical of explicit force-based integration schemes.

Let $\mathbf{x} \in \mathbb{R}^{3N}$ denote the system position vector. A given constraint is defined by a scalar or vector function $C(\mathbf{x}) = 0$. To establish a physical correspondence with classical mechanics, the internal constraint force $\mathbf{f}_{int}$ is modeled via a generalized restoring behavior analogous to Hooke's Law:

$$C(\mathbf{x} + \Delta\mathbf{x}) = -\frac{1}{k} \lambda$$

where $k$ denotes the structural stiffness and $\lambda$ represents the constraint force magnitude (Lagrange multiplier). Applying a first-order Taylor expansion to the left-hand side yields:

$$C(\mathbf{x}) + \nabla C(\mathbf{x})\,\Delta\mathbf{x} \approx -\frac{1}{k}\lambda$$

In a mass-weighted system, the position correction $\Delta\mathbf{x}$ must follow the constraint gradient direction to conserve linear and angular momentum:

$$\Delta\mathbf{x} = \mathbf{M}^{-1} \nabla C(\mathbf{x})^T \lambda \,\Delta t^2$$

where $\mathbf{M}^{-1}$ is the inverse mass matrix and $\Delta t$ is the simulation time step. This guarantees that internal forces act purely along the norm of the constraint manifold.

---

## 2. Distance Constraints

For a system containing two particles with positions $\mathbf{x}_1, \mathbf{x}_2$ and masses $m_1, m_2$, the distance constraint enforcing a rest length $d$ is written as:

$$C(\mathbf{x}_1, \mathbf{x}_2) = \Vert{}\mathbf{x}_1 - \mathbf{x}_2\Vert{} - d$$

The respective gradients with respect to individual particle coordinates are:

$$\nabla_{\mathbf{x}_1} C = \frac{\mathbf{x}_1 - \mathbf{x}_2}{\Vert{}\mathbf{x}_1 - \mathbf{x}_2\Vert{}} = \mathbf{n}, \quad \nabla_{\mathbf{x}_2} C = -\frac{\mathbf{x}_1 - \mathbf{x}_2}{\Vert{}\mathbf{x}_1 - \mathbf{x}_2\Vert{}} = -\mathbf{n}$$

Substituting these gradients into the position correction equation yields:

$$\Delta\mathbf{x}_1 = \frac{1}{m_1} \mathbf{n} \cdot \lambda \Delta t^2, \quad \Delta\mathbf{x}_2 = -\frac{1}{m_2} \mathbf{n} \cdot \lambda \Delta t^2$$

Taking the ratio of the displacements reveals a fundamental property of the solver:

$$\frac{\Vert{}\Delta\mathbf{x}_1\Vert{}}{\Vert{}\Delta\mathbf{x}_2\Vert{}} = \frac{m_2}{m_1}$$

This confirms that position updates are inversely proportional to particle masses, ensuring momentum conservation ($m_1\Delta\mathbf{x}_1 + m_2\Delta\mathbf{x}_2 = \mathbf{0}$). When an object is discretized into a network of these distance constraints, it models a stable deformable continuum whose material stiffness and dynamic oscillations can be modulated through $k$ (or compliance $\alpha = \frac{1}{k}$ in XPBD).

---

## 3. Shape Matching Constraints

To simulate rigid or highly stiff bodies without solving ill-conditioned global systems, shape matching maps a chaotic particle configuration back to a reference shape via an optimal rigid transformation.

Given a set of current particle positions $\mathbf{x}_i$ and reference positions $\mathbf{x}_i^0$ with masses $m_i$, we define the respective centers of mass as:

$$\mathbf{x}_c = \frac{\sum m_i \mathbf{x}_i}{\sum m_i}, \quad \mathbf{x}_c^0 = \frac{\sum m_i \mathbf{x}_i^0}{\sum m_i}$$

Let the relative rest coordinates be $\mathbf{a}_i = \mathbf{x}_i^0 - \mathbf{x}_c^0$. The goal is to compute an optimal rotation matrix $\mathbf{R}$ and translation vector $\mathbf{T}$ that minimizes the mass-weighted least-squares error:

$$\min_{\mathbf{R}, \mathbf{T}} \sum_i m_i \Vert{} \mathbf{R}\mathbf{a}_i + \mathbf{T} - (\mathbf{x}_i - \mathbf{x}_c) \Vert{}^2$$

The optimal translation aligns the centers of mass directly:

$$\mathbf{T} = \mathbf{x}_c - \mathbf{x}_c^0$$

Setting $\mathbf{b}_i = \mathbf{x}_i - \mathbf{x}_c$, the remaining objective reduces to finding the rotation $\mathbf{R}$ that minimizes:

$$\mathcal{L}(\mathbf{R}) = \sum_i m_i \Vert{} \mathbf{R}\mathbf{a}_i - \mathbf{b}_i \Vert{}^2$$

Expanding the quadratic norm yields:

$$\mathcal{L}(\mathbf{R}) = \sum_i m_i \left( \mathbf{a}_i^T \mathbf{R}^T \mathbf{R} \mathbf{a}_i - 2\mathbf{b}_i^T \mathbf{R} \mathbf{a}_i + \mathbf{b}_i^T \mathbf{b}_i \right)$$

Since $\mathbf{R}^T \mathbf{R} = \mathbf{I}$, the first and third terms are independent of rotation. Minimizing $\mathcal{L}(\mathbf{R})$ is mathematically equivalent to maximizing the cross term:

$$\max_{\mathbf{R}} \sum_i m_i \mathbf{b}_i^T \mathbf{R} \mathbf{a}_i = \max_{\mathbf{R}} \text{Tr}\left( \mathbf{R} \sum_i m_i \mathbf{a}_i \mathbf{b}_i^T \right)$$

We define the linear covariance matrix $\mathbf{A}$ as:

$$\mathbf{A} = \sum_i m_i \mathbf{b}_i \mathbf{a}_i^T \implies \mathbf{A}^T = \sum_i m_i \mathbf{a}_i \mathbf{b}_i^T$$

Thus, the objective function becomes $\max_{\mathbf{R}} \text{Tr}(\mathbf{R}\mathbf{A}^T)$.

### 3.1. Matrix Decomposition and Maximization Proof

Using Polar Decomposition, $\mathbf{A}$ can be factored into an orthogonal rotation matrix $\mathbf{R}_p$ and a symmetric positive semi-definite stretch matrix $\mathbf{S}$:

$$\mathbf{A} = \mathbf{R}_p \mathbf{S} \implies \mathbf{A}^T = \mathbf{S} \mathbf{R}_p^T$$

Substituting this into the trace maximization function:

$$\text{Tr}(\mathbf{R}\mathbf{A}^T) = \text{Tr}(\mathbf{R} \mathbf{S} \mathbf{R}_p^T) = \text{Tr}(\mathbf{R}_p^T \mathbf{R} \mathbf{S})$$

Let $\mathbf{G} = \mathbf{R}_p^T \mathbf{R}$. Since both $\mathbf{R}_p$ and $\mathbf{R}$ are proper rotation matrices, $\mathbf{G}$ is also an orthonormal rotation matrix ($\mathbf{G} \in \text{SO}(3)$). The stretch matrix $\mathbf{S}$ can be diagonalized using an orthogonal basis $\mathbf{Q}$:

$$\mathbf{S} = \mathbf{Q} \mathbf{\Sigma} \mathbf{Q}^T$$

where $\mathbf{\Sigma} = \text{diag}(\sigma_1, \sigma_2, \sigma_3)$ contains the non-negative eigenvalues. The trace equation expands to:

$$\text{Tr}(\mathbf{G}\mathbf{S}) = \text{Tr}(\mathbf{G}\mathbf{Q}\mathbf{\Sigma}\mathbf{Q}^T) = \text{Tr}(\mathbf{Q}^T\mathbf{G}\mathbf{Q}\mathbf{\Sigma})$$

Let $\mathbf{H} = \mathbf{Q}^T\mathbf{G}\mathbf{Q}$. Because $\mathbf{Q}$ and $\mathbf{G}$ are orthogonal, $\mathbf{H}$ is also orthogonal, meaning its diagonal components satisfy $H_{ii} \le 1$. The trace objective simplifies to a linear combination of its singular values:

$$\text{Tr}(\mathbf{H}\mathbf{\Sigma}) = \sum_{k=1}^3 H_{kk}\sigma_k$$

Since $\sigma_k \ge 0$, this sum reaches its global maximum if and only if $H_{kk} = 1$, implying $\mathbf{H} = \mathbf{I}$. 

$$\mathbf{H} = \mathbf{I} \implies \mathbf{Q}^T\mathbf{G}\mathbf{Q} = \mathbf{I} \implies \mathbf{G} = \mathbf{I}$$

Substituting back into our definition of $\mathbf{G}$:

$$\mathbf{R}_p^T \mathbf{R} = \mathbf{I} \implies \mathbf{R} = \mathbf{R}_p$$

The optimal rotation matrix is exactly the rotational component of the polar decomposition of the covariance matrix $\mathbf{A}$. This can alternatively be evaluated via Singular Value Decomposition (SVD), where $\mathbf{A} = \mathbf{U}\mathbf{\Sigma}\mathbf{V}^T$, yielding $\mathbf{R} = \mathbf{U}\mathbf{V}^T$. If $\det(\mathbf{A}) < 0$, a reflection must be handled by flipping the sign of the column corresponding to the smallest singular value.

### 3.2. Numerical Evaluation via Newton-Schulz Iteration

To avoid the high computational cost of explicit SVD on real-time architectures, the polar rotation component $\mathbf{R}$ can be efficiently approximated using the **Newton-Schulz** iterative algorithm:

$$\mathbf{R}_{k+1} = \frac{1}{2} \left( \mathbf{R}_k + (\mathbf{R}_k^{-1})^T \right)$$

with the initial state set to the scaled covariance matrix $\mathbf{R}_0 = \mathbf{A}$.

#### Proof of Convergence:
Using the SVD properties of the initial matrix $\mathbf{A} = \mathbf{U}\mathbf{\Sigma}_0\mathbf{V}^T$:

$$\mathbf{R}_k = \mathbf{U}\mathbf{\Sigma}_k\mathbf{V}^T \implies (\mathbf{R}_k^{-1})^T = \mathbf{U}\mathbf{\Sigma}_k^{-1}\mathbf{V}^T$$

Substituting these terms back into the recurrence relation:

$$\mathbf{U}\mathbf{\Sigma}_{k+1}\mathbf{V}^T = \frac{1}{2}\left( \mathbf{U}\mathbf{\Sigma}_k\mathbf{V}^T + \mathbf{U}\mathbf{\Sigma}_k^{-1}\mathbf{V}^T \right)$$

$$\mathbf{\Sigma}_{k+1} = \frac{1}{2}\left( \mathbf{\Sigma}_k + \mathbf{\Sigma}_k^{-1} \right)$$

As $k \to \infty$, each singular value independent component drives quadratically toward identity ($\sigma \to 1$). Consequently, $\mathbf{\Sigma}_k \to \mathbf{I}$, yielding the proper orthogonal rotation matrix:

$$\lim_{k \to \infty} \mathbf{R}_k = \mathbf{U}\mathbf{I}\mathbf{V}^T = \mathbf{U}\mathbf{V}^T$$

Typically, 6 to 9 iterations are sufficient to guarantee high numerical precision for real-time physics frameworks.

