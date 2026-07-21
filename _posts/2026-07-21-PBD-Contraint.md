---
layout: post
title: 'Mathematical Analysis of Geometric Constraints in Position-Based Dynamics (PBD) and XPBD'
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
Left side use Taylor expansion C(Xcurrent)+ Transpose(gradient(C))* delta X  = - 1/k * Lambda_force , 
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

we define a group of particles (Xi ,Mi) which forms a whole shape , at rest we have the mass center Xc_rest = sum(Xi_rest*Mi)/sum(Mi)

then after each particle moving we want to keep the whole shape the same as the reset shape at any moment.

if you want to keep the whole shape want you can have for the spatial change of the whole shape is  a translate T and a rotation R 
X_i_new_target = R*(X_i_rest-Xc_rest)+T

The new position X_i_new after moving may be at chaotic places , the question is how do we adjust each of those X_i_new to X_i_new_target.
that is to solve the best value for the R and T given X_i_new.

Intuitively T is 3d vector which is easy to get , the new Xc_new= sum(Xi_new*Mi)/sum(Mi) 
T= Xc_new - Xc_rest

Intuitively We define The loss(R) = sum[ ||R*(X_i_rest-Xc_rest)+ T - X_i_new||^2 * Mi ] 

Ai = X_i_rest - Xc_rest , Bi= X_i_new - T

that is we want to minimize loss(R) = sum[ Mi*||R*Ai-Bi||^2 ]

loss(R) = sum [ Mi * Transpose(R*Ai+Bi) x (R*Ai+Bi)  ]
        = sum [ Mi * ( Transpose(Ai)*Transpose(R)*R*Ai - Transpose(Ai)*Transpose(R)*Bi - Transpose(Bi)*R*Ai + Transpose(Bi)*Bi ) ]
    
Transpose(R)*R = I ,  
as Transpose(Ai)*Transpose(R)*Bi is scalar, Transpose(Ai)*Transpose(R)*Bi = Transpose(Transpose(Ai)*Transpose(R)*Bi)=  Transpose(Bi)*R*Ai

loss(R) = sum [ Mi * ( Transpose(Ai)*Ai - 2*Transpose(Bi)*R*Ai + Transpose(Bi)*Bi ) ]

that is to find the R_star to maxmize the sum[Mi*Transpose(Bi)*R_star*Ai ]
                               = Trace ( R_star * sum[ Mi*Ai*Transpose(Bi) ] ) 

define sum[ Mi*Ai*Transpose(Bi) ] = transpose(A)

to maxmize Trace(R_star*transpose(A))

polar decompose A = R*S ,  
here if det(A)!=0 , then we always get a unique decompostion. 
We need to make sure actually det(A) != 0 
In the simulation if we meet det(A)=0 we just fallback to the previous Rotation.
For det(A)<0 , you need to do some work to handle it . 

Trace(R_star*transpose(A)) = Trace(R_star*transpose(S)*transpose(R)) = Trace(R_star*s*transpose(R)) = Trace(transpose(R)*R_star*S)

here matrix transpose(R)*R_star = G is also a rotation matrix, S is sysmetric matrix which can be diagnalized S= Q * Sigma * Transpose(Q)

so maxmize Trace( G * Q * Sigma * Transpose(Q)) = Trace( Transpose(Q) * G * Q * Sigma ) , 

here Q is orthornormal matrix , G is orthornormal matrix (R_star and transpose(R) are both orthornormal matrix ) , 
so Transpose(Q) * G * Q is orthornormal matrix which means each colum sum to 1 then it easily to see that 
when Transpose(Q) * G * Q = I will make it maxmized.

so that is Transpose(Q) * transpose(R)*R_star * Q = I  => transpose(R) * R_star = I => R_star =  R


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

**References:**
* Müller, M., Heidelberger, B., Hennix, M., & Ratcliff, J. (2007). *Position Based Dynamics*. [[Paper Link](https://github.io)]
* Macklin, M., Müller, M., & Chentanez, N. (2016). *XPBD: Position-Based Simulation of Compliant Constrained Dynamics*. [[Paper Link](https://github.io)]

---
# Detailed Analysis of Constraints in PBD and XPBD

## 1. PBD Basic Assumption

Starting from a constraint function $C(\mathbf{x})$, the internal restoring force acts to rebalance the system toward $C(\mathbf{x}) = 0$. Following Hooke's Law ($\mathbf{x} = -\frac{1}{k}\mathbf{F}$):

$$C(\mathbf{x}_{\text{current}} + \Delta \mathbf{x}) = -\frac{1}{k} \lambda_{\text{force}}$$

Applying a first-order Taylor expansion yields:

$$C(\mathbf{x}_{\text{current}}) + \nabla C(\mathbf{x}) \cdot \Delta \mathbf{x} = - \frac{1}{k} \lambda_{\text{force}}$$

The position update is defined along the constraint gradient vector:

$$\Delta \mathbf{x} = \mathbf{M}^{-1} \nabla C(\mathbf{x}) \lambda_{\text{force}} \Delta t^2$$

Here, $\lambda_{\text{force}}$ is a scalar, and having the force vector in the direction of the gradient $\nabla C(\mathbf{x}) \lambda_{\text{force}}$ is physically reasonable.

---

## 2. Distance Constraint

For two particles, the constraint is:

$$C(\mathbf{x}_1, \mathbf{x}_2) = \|\mathbf{x}_1 - \mathbf{x}_2\| - d$$

Using the position update equation:

$$\Delta \mathbf{x}_1 = \frac{1}{m_1} \frac{\partial C(\mathbf{x})}{\partial \mathbf{x}_1} \lambda_{\text{force}} \Delta t^2$$

$$\Delta \mathbf{x}_2 = \frac{1}{m_2} \frac{\partial C(\mathbf{x})}{\partial \mathbf{x}_2} \lambda_{\text{force}} \Delta t^2$$

Since $\frac{\partial C(\mathbf{x})}{\partial \mathbf{x}_1} = -\frac{\partial C(\mathbf{x})}{\partial \mathbf{x}_2}$, the displacement ratio scales inversely with mass:

$$\frac{\|\Delta \mathbf{x}_1\|}{\|\Delta \mathbf{x}_2\|} = \frac{m_2}{m_1}$$

We can see that with this kind of definition of constraint, what we really define is that with two particles, while keeping their distance, when the distance changes, the bigger mass moves smaller. That is what our distance constraint means. 

When you put a group of particles with adjacent particles having such distance constraints, the whole object will behave like a soft body object. No matter what you do with the stiffness $k$ parameter, the whole object will still adjust itself with shape changing, as you could imagine that with different $k$ the object kind of shakes differently.

---

## 3. Shape Constraint

If we want to simulate a rigid body, what we want is at each step the whole shape keeps the same. So how do we define this kind of constraint?

We define a group of particles $(\mathbf{X}_i, m_i)$ which forms a whole shape. At rest, we have the mass center:

$$\mathbf{X}_{c\_\text{rest}} = \frac{\sum \mathbf{X}_{i\_\text{rest}} m_i}{\sum m_i}$$

Then, after each particle moving, we want to keep the whole shape the same as the rest shape at any moment. If you want to keep the whole shape, what you can have for the spatial change of the whole shape is a translation $\mathbf{T}$ and a rotation $\mathbf{R}$:

$$\mathbf{X}_{i\_\text{new\_target}} = \mathbf{R}(\mathbf{X}_{i\_\text{rest}} - \mathbf{X}_{c\_\text{rest}}) + \mathbf{T}$$

The new position $\mathbf{X}_{i\_\text{new}}$ after moving may be at chaotic places. The question is how do we adjust each of those $\mathbf{X}_{i\_\text{new}}$ to $\mathbf{X}_{i\_\text{new\_target}}$. That is to solve the best value for $\mathbf{R}$ and $\mathbf{T}$ given $\mathbf{X}_{i\_\text{new}}$.

Intuitively, $\mathbf{T}$ is a 3D vector which is easy to get from the new center of mass $\mathbf{X}_{c\_\text{new}} = \frac{\sum \mathbf{X}_{i\_\text{new}} m_i}{\sum m_i}$:

$$\mathbf{T} = \mathbf{X}_{c\_\text{new}} - \mathbf{X}_{c\_\text{rest}}$$

Intuitively, we define the loss function as:

$$\mathcal{L}(\mathbf{R}) = \sum m_i \|\mathbf{R}(\mathbf{X}_{i\_\text{rest}} - \mathbf{X}_{c\_\text{rest}}) + \mathbf{T} - \mathbf{X}_{i\_\text{new}}\|^2$$

Letting relative coordinates be exactly defined as:

$$\mathbf{A}_i = \mathbf{X}_{i\_\text{rest}} - \mathbf{X}_{c\_\text{rest}}, \quad \mathbf{B}_i = \mathbf{X}_{i\_\text{new}} - \mathbf{T}$$

That is, we want to minimize:

$$\mathcal{L}(\mathbf{R}) = \sum m_i \|\mathbf{R}\mathbf{A}_i - \mathbf{B}_i\|^2$$

Expanding the expression:

$$\mathcal{L}(\mathbf{R}) = \sum m_i \left( (\mathbf{R}\mathbf{A}_i - \mathbf{B}_i)^T (\mathbf{R}\mathbf{A}_i - \mathbf{B}_i) \right)$$

$$\mathcal{L}(\mathbf{R}) = \sum m_i \left( \mathbf{A}_i^T \mathbf{R}^T \mathbf{R} \mathbf{A}_i - \mathbf{A}_i^T \mathbf{R}^T \mathbf{B}_i - \mathbf{B}_i^T \mathbf{R} \mathbf{A}_i + \mathbf{B}_i^T \mathbf{B}_i \right)$$

Since $\mathbf{R}^T \mathbf{R} = \mathbf{I}$, and because $\mathbf{A}_i^T \mathbf{R}^T \mathbf{B}_i$ is a scalar, we have $\mathbf{A}_i^T \mathbf{R}^T \mathbf{B}_i = (\mathbf{A}_i^T \mathbf{R}^T \mathbf{B}_i)^T = \mathbf{B}_i^T \mathbf{R} \mathbf{A}_i$. This simplifies to:

$$\mathcal{L}(\mathbf{R}) = \sum m_i \left( \mathbf{A}_i^T \mathbf{A}_i - 2\mathbf{B}_i^T \mathbf{R} \mathbf{A}_i + \mathbf{B}_i^T \mathbf{B}_i \right)$$

Thus, minimizing $\mathcal{L}(\mathbf{R})$ is to find the $\mathbf{R}^*$ to maximize:

$$\sum m_i \mathbf{B}_i^T \mathbf{R}^* \mathbf{A}_i = \text{Tr}\left(\mathbf{R}^* \sum m_i \mathbf{A}_i \mathbf{B}_i^T\right)$$

Defining $\mathbf{A}^T = \sum m_i \mathbf{A}_i \mathbf{B}_i^T$, we maximize $\text{Tr}(\mathbf{R}^* \mathbf{A}^T)$.

Using Polar Decomposition ($\mathbf{A} = \mathbf{R}\mathbf{S}$):

$$\text{Tr}(\mathbf{R}^* \mathbf{A}^T) = \text{Tr}(\mathbf{R}^* (\mathbf{R}\mathbf{S})^T) = \text{Tr}(\mathbf{R}^* \mathbf{S}^T \mathbf{R}^T) = \text{Tr}(\mathbf{R}^T \mathbf{R}^* \mathbf{S})$$

Here, if $\det(\mathbf{A}) \neq 0$, then we always get a unique decomposition. We need to make sure actually $\det(\mathbf{A}) \neq 0$. In the simulation, if we meet $\det(\mathbf{A}) = 0$, we just fallback to the previous rotation. For $\det(\mathbf{A}) < 0$, you need to do some work to handle it.

The matrix $\mathbf{G} = \mathbf{R}^T \mathbf{R}^*$ is also a rotation matrix. Since $\mathbf{S}$ is a symmetric matrix, it can be diagonalized as $\mathbf{S} = \mathbf{Q} \mathbf{\Sigma} \mathbf{Q}^T$:

$$\max \text{Tr}(\mathbf{G} \mathbf{Q} \mathbf{\Sigma} \mathbf{Q}^T) = \max \text{Tr}(\mathbf{Q}^T \mathbf{G} \mathbf{Q} \mathbf{\Sigma})$$

Here, $\mathbf{Q}$ is an orthonormal matrix and $\mathbf{G}$ is an orthonormal matrix, so $\mathbf{Q}^T \mathbf{G} \mathbf{Q}$ is an orthonormal matrix which means each column sum squares to 1. Then it is easy to see that setting $\mathbf{Q}^T \mathbf{G} \mathbf{Q} = \mathbf{I}$ will make it maximized. 

That is $\mathbf{Q}^T \mathbf{R}^T \mathbf{R}^* \mathbf{Q} = \mathbf{I} \implies \mathbf{R}^T \mathbf{R}^* = \mathbf{I} \implies \mathbf{R}^* = \mathbf{R}$.

Alternatively, using SVD Decomposition ($\mathbf{A} = \mathbf{U}_{\text{svd}} \mathbf{\Sigma} \mathbf{V}_{\text{svd}}^T$). As $\mathbf{U}_{\text{svd}}$ and $\mathbf{V}_{\text{svd}}$ are orthonormal, $\mathbf{V}_{\text{svd}}^T \mathbf{V}_{\text{svd}} = \mathbf{I}$:

$$\mathbf{A} = \mathbf{U}_{\text{svd}} \mathbf{V}_{\text{svd}}^T \mathbf{V}_{\text{svd}} \mathbf{\Sigma} \mathbf{V}_{\text{svd}}^T = (\mathbf{U}_{\text{svd}} \mathbf{V}_{\text{svd}}^T) (\mathbf{V}_{\text{svd}} \mathbf{\Sigma} \mathbf{V}_{\text{svd}}^T)$$

So $\mathbf{R}_{\text{polar}} = \mathbf{U}_{\text{svd}} \mathbf{V}_{\text{svd}}^T$ and $\mathbf{S} = \mathbf{V}_{\text{svd}} \mathbf{\Sigma} \mathbf{V}_{\text{svd}}^T$.

### Newton-Schulz Algorithm

Rather than using SVD decomposition, which can be difficult to compute, we use the Newton-Schulz algorithm to calculate $\mathbf{R}_{\text{polar}}$ because it is fast and easy to implement:

$$\mathbf{R}_{k+1} = \frac{1}{2} \left( \mathbf{R}_k + (\mathbf{R}_k^{-1})^T \right), \quad \text{where } \mathbf{R}_0 = \mathbf{A}$$

By substituting the SVD form into the iteration, we can analyze its behavior:

$$\mathbf{R}_k = \mathbf{U}_{\text{svd}} \mathbf{\Sigma}_k \mathbf{V}_{\text{svd}}^T$$

$$(\mathbf{R}_k^{-1})^T = \mathbf{U}_{\text{svd}} \mathbf{\Sigma}_k^{-1} \mathbf{V}_{\text{svd}}^T$$

$$\mathbf{R}_{k+1} = \mathbf{U}_{\text{svd}} \left[ \frac{1}{2} \left( \mathbf{\Sigma}_k + \mathbf{\Sigma}_k^{-1} \right) \right] \mathbf{V}_{\text{svd}}^T$$

As we can see, from $k$ to $k+1$ only the middle term changes; specifically, it drives toward $\mathbf{I}$. 

Therefore, after a few iterations (typically 6–9 times), $\mathbf{R}_k$ will converge directly to $\mathbf{U}_{\text{svd}} \mathbf{V}_{\text{svd}}^T$, which is precisely the optimal rotation matrix we need.
