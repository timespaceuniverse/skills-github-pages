---
layout: post
title: "A Rigorous Mathematical Formulation of PBD and XPBD"
mathjax: true
---

<!-- raw before LLM format , lazy work

We give a much more logic clear way to represent the math of PBD and XPBD in this blog.
PBD Paper: https://matthias-research.github.io/pages/publications/posBasedDyn.pdf
XPBD Paper : https://matthias-research.github.io/pages/publications/XPBD.pdf

ABBREviate: PBD (Position based dynamics) and XPBD (Extended PBD) can simulate world physics(rigid body ,soft body, fluid)
using a unified particles physics (Newton's Law no quantum physics involved). 
XPBD framework can be in the future integrated into AI as a basic building block to construct a world model.
The math in the XPBD paper is fussy, we here give a much clear step by step deduction.


the math of (X)PBD is easy :
Given a n of partiles with initial position pi(x) ,massi , 
the newton's moment law gives: Delta Position = Delta T * V = Delta T * F / M  
Simulation : Delta T , Mass are knowns in computer simulation environment , Force is not known.
Force can be split to two categroies : 1. active force (external) 2. inactive force .(internal)
Examples: Active force like people push someobject which can be defined as sytem external force .
Inactive force like friction force , like spring internal force.  
In Reality it is the way how we could define it easily . For External force which can be directly 
defined and for internal force we don't define it directly but deduced from Constraint(x). 
There are serval types of constraint which used to describe different type of object (fluid,rigid,soft)
So for the given many paritiles we can define many constraints with different constraints types.

Example :  
          given two particles in one dimension space . (x1,m1) , (x2,m2) , a distance constraint 
          can be defined as C(x1,x2)= ||x1-x2|| - L , you can think of it as we always want the two 
          particles keep a distance L , if they are too close C < 0 , too far C> 0 , distant L when C=0.
          So we can always change x1 , x2 to let it get balanced C=0 .  Imagine a group of particles with 
          each has a distance constraint to close particles . This kind of structure could represent 
          a softbody. For other type of constraints(rigid body ,fluid ) i will not explain here , 
          you can refer to the paper for details. 

with C(x) , we could define the internal force easily :
the internal lambda_force is try to make the system rebalanced (C(x)=0) from the current position which is 
like hooks' law.  C(Xcurrent + delta x) = - 1/k * Lambda_force  (Hook's law : x = -1/k * F ,negtive sign for restoring force)
Left side use Taylor expansion C(Xcurrent)+ Transpose(gradient(C))* delta X  = 1/k * Lambda_force , 
delta X = M^-1 * gradient(C) * Lambda_force * delta t * delta t   ,  here Lambda_force is a scalar,
gradient(C) * Lambda_force the force vector is in the direction of the gradient which is reasonable.
So put together replace delta x gives :
C(x) + Transpose(gradient(C)) * (M^-1 * gradient(C) * Lambda_force * delta t * delta t ) =  - 1/k * Lambda_force ,
for the paper in original XPBD they define lambda = lambda_force * delta t * delta t , and 1/k = alpha  , 1/(k * delta t * delta t) = alpha~
so we can write :
C(x) + Transpose(gradient(C)) * (M^-1 * gradient(C) * Lambda_force * delta t * delta t ) =  - 1/(k * delta t * delta t) * Lambda_force *(delta t * delta t)
which is :
C(x) + Transpose(gradient(C)) * (M^-1 * gradient(C) * Lambda ) =  - alpha~ * Lambda , 

In the Simulation the Lambda is calculating in K steps loop , for each loop the Position x get updated and 
the Lambda is updated accordingly. The actually logic is you get deltaX as you solve constraint Cj(X) one by one ,
After you have the delta x you are going to solve a delta lambda 

for the Step K you have :
C(Xk+deltax) = C(Xk) + Transpose ( gradient C(x)) *(M^-1* gradient C(Xk)*lambda k *delta t ^ 2) 

solve the lambda k , solve the delta x , get X k+1 = X k 

for the Step K+1 you have : 

C(Xk+deltax + small delta x ) = C(Xk) + Transpose ( Gradient C(Xk))*(M^-1 *Gradient C(Xk)(Lamdba k + small delta Lambda) * delta t ^ 2)

= C(Xk+1) +  Transpose ( Gradient C(Xk))*(M^-1 *Gradient C(Xk)(small delta Lambda) * delta t ^ 2)  =  - alpha * ( Lambda k + small delta Lambda ) * delta t ^ 2

Intersting here :

You can assume : Gradient C(Xk) = Gradient C(Xk+1) or you stick to C(X0) for better performance.

from the last equation :
M^-1 is just a diagnal matrix , redefine the variables with delta t ^2 abbsorbed we finally have (same with the paper) : 

small delta Lambda  =  (-C(Xk+1)-Lambda k * ~alpha)/( Transpose ( Gradient C(Xk)) * M^-1 *Gradient C(Xk) + ~alpha)

--> 

# A Rigorous Mathematical Formulation of Position-Based Dynamics (PBD) and Extended PBD (XPBD)

## 1. Introduction
This document presents a structured, mathematically transparent derivation of Position-Based Dynamics (PBD) and Extended Position-Based Dynamics (XPBD). Grounded in Newtonian mechanics, these frameworks establish a unified particle-based representation capable of simulating macro-scale physical entities, including rigid bodies, soft bodies, and fluids. Beyond traditional computer graphics animations, the algorithmic robustness of the XPBD framework positions it as a foundational inductive bias for physics-informed artificial intelligence and the construction of predictive world models. 

While the original literature introduces algebraic complexities that mask the underlying intuition, this text delivers a clear, step-by-step mathematical deduction.

* **Primary References:**
  * Müller et al., [Position-Based Dynamics (2007)](https://matthias-research.github.io/pages/publications/posBasedDyn.pdf)
  * Macklin et al., [XPBD: Extended Position-Based Dynamics (2016)](https://matthias-research.github.io/pages/publications/XPBD.pdf)

---

## 2. Kinematic Foundations & Force Classification
Consider a discrete dynamical system consisting of $n$ particles. Each particle $i$ is parameterized by its position $\mathbf{x}_i \in \mathbb{R}^3$ and elemental mass $m_i$. Let $\mathbf{M}$ denote the diagonal mass matrix. Under a standard discrete-time Newtonian formulation with a temporal step size $\Delta t$, the position update vector $\Delta \mathbf{x}$ maps directly to the generalized forces acting on the system:

$$\Delta \mathbf{x} = \Delta t \, \mathbf{v} = \Delta t^2 \mathbf{M}^{-1} \mathbf{F}$$

In a discrete computer simulation environment, variables such as $\Delta t$ and $\mathbf{M}$ are static and known. However, the generalized forces $\mathbf{F}$ must be resolved. We decompose $\mathbf{F}$ into two orthogonal categories:
1. **External (Active) Forces ($\mathbf{F}_{ext}$):** Directly defined user or environmental inputs (e.g., gravity, manual perturbations).
2. **Internal (Inactive/Constraint) Forces ($\mathbf{F}_{int}$):** Restoring forces that preserve structural topology. Instead of defining these forces explicitly via stiff, numerically volatile differential equations, they are implicitly deduced from a system of geometric constraints.

---

## 3. Geometric Constraints & Structural Representation
Internal mechanics are governed by a set of $m$ scalar constraint functions:

$$C_j(\mathbf{x}) : \mathbb{R}^{3n} \rightarrow \mathbb{R}, \quad j \in \{1, \dots, m\}$$

The system achieves equilibrium when all constraints are satisfied, such that $C_j(\mathbf{x}) = 0$. 

### Example: One-Dimensional Distance Constraint
For a minimal two-particle system with coordinates $\mathbf{x}_1, \mathbf{x}_2$ and a target rest length $L$, the distance constraint is expressed as:

$$C(\mathbf{x}_1, \mathbf{x}_2) = \|\mathbf{x}_1 - \mathbf{x}_2\| - L$$

* $C(\mathbf{x}) < 0 \implies$ Compression regime.
* $C(\mathbf{x}) > 0 \implies$ Tension regime.
* $C(\mathbf{x}) = 0 \implies$ Equilibrium state.

By coupling particles globally via localized distance networks, complex macroscopic structures like deformable soft bodies can be represented seamlessly.

---

## 4. Constraint-Force Derivation via Hooke's Law
To drive a deformed system at state $\mathbf{x}$ back to its equilibrium manifold ($C(\mathbf{x} + \Delta \mathbf{x}) = 0$), we introduce an internal generalized constraint force vector along the constraint gradient, scaled by a scalar magnitude $\lambda_{force}$. Analogous to a generalized continuous Hooke's law:

$$C(\mathbf{x} + \Delta \mathbf{x}) = -\frac{1}{k} \lambda_{force}$$

where $k$ represents physical stiffness, and the negative sign denotes its restoring characteristic. 

Applying a first-order Taylor expansion to the left-hand side yields:

$$C(\mathbf{x}) + \nabla C(\mathbf{x})^T \Delta \mathbf{x} = -\frac{1}{k} \lambda_{force}$$

The displacement $\Delta \mathbf{x}$ induced by the constraint force over time step $\Delta t$ is defined as:

$$\Delta \mathbf{x} = \mathbf{M}^{-1} \nabla C(\mathbf{x}) \lambda_{force} \Delta t^2$$

Substituting this displacement equation back into the Taylor expansion results in:

$$C(\mathbf{x}) + \nabla C(\mathbf{x})^T \left( \mathbf{M}^{-1} \nabla C(\mathbf{x}) \lambda_{force} \Delta t^2 \right) = -\frac{1}{k} \lambda_{force}$$

### XPBD Parameter Scaling
To eliminate time-step dependencies, XPBD re-parameterizes the stiffness $k$ and force magnitude $\lambda_{force}$. We define the total Lagrange multiplier $\lambda$ and compliance $\alpha, \tilde{\alpha}$ as:

$$\lambda = \lambda_{force} \Delta t^2, \quad \alpha = \frac{1}{k}, \quad \tilde{\alpha} = \frac{\alpha}{\Delta t^2} = \frac{1}{k \Delta t^2}$$

Multiplying the right side by $\frac{\Delta t^2}{\Delta t^2}$ yields:

$$C(\mathbf{x}) + \nabla C(\mathbf{x})^T \left( \mathbf{M}^{-1} \nabla C(\mathbf{x}) \lambda \right) = -\left(\frac{1}{k \Delta t^2}\right) (\lambda_{force} \Delta t^2)$$

This simplifies to the fundamental continuous XPBD constraint equation:

$$C(\mathbf{x}) + \nabla C(\mathbf{x})^T \mathbf{M}^{-1} \nabla C(\mathbf{x}) \lambda = -\tilde{\alpha} \lambda$$

---

## 5. Iterative Solver & Localized Update Rules
In practical computer simulations, constraints are solved sequentially in an inner loop across $K$ sub-steps. Let $k$ denote the current sub-step index within the solver sequence. Position adjustments $\Delta \mathbf{x}$ are determined by solving individual constraints $C_j(\mathbf{x})$ one by one. After computing the base displacement, an incremental "small delta" change for the Lagrange multiplier ($\delta \lambda$) is calculated to advance to the next step.

For sub-step $k$, the constraint evaluation under a structural correction satisfies:

$$C(\mathbf{x}_k + \Delta \mathbf{x}) = C(\mathbf{x}_k) + \nabla C(\mathbf{x}_k)^T \left( \mathbf{M}^{-1} \nabla C(\mathbf{x}_k) \lambda_k \Delta t^2 \right)$$

Solving for $\lambda_k$ and $\Delta \mathbf{x}$ yields the advanced intermediate position state:

$$\mathbf{x}_{k+1} = \mathbf{x}_k + \Delta \mathbf{x}$$

For the subsequent sub-step $k+1$, introducing a small incremental position change $\delta \mathbf{x}$ alongside its corresponding incremental multiplier update $\delta \lambda$ yields:

$$C(\mathbf{x}_k + \Delta \mathbf{x} + \delta \mathbf{x}) = C(\mathbf{x}_k) + \nabla C(\mathbf{x}_k)^T \left( \mathbf{M}^{-1} \nabla C(\mathbf{x}_k) (\lambda_k + \delta \lambda) \Delta t^2 \right)$$

$$= C(\mathbf{x}_{k+1}) + \nabla C(\mathbf{x}_k)^T \left( \mathbf{M}^{-1} \nabla C(\mathbf{x}_k) \delta \lambda \Delta t^2 \right) = -\alpha (\lambda_k + \delta \lambda) \Delta t^2$$

> **Computational Note:** For optimized runtime efficiency, it is standard practice to assume that the constraint gradient is invariant across localized updates ($\nabla C(\mathbf{x}_k) \approx \nabla C(\mathbf{x}_{k+1})$), or to stick strictly to the initial configuration state $\nabla C(\mathbf{x}_0)$.

Given that $\mathbf{M}^{-1}$ is a diagonal matrix, we absorb the $\Delta t^2$ factor into our normalized compliance parameter ($\tilde{\alpha} = \frac{\alpha}{\Delta t^2}$) to match original literature conventions. Isolating the small incremental multiplier $\delta \lambda$ from the matching sides of the final equations reveals the definitive XPBD system update rule:

$$\delta \lambda = \frac{-C(\mathbf{x}_{k+1}) - \lambda_k \tilde{\alpha}}{\nabla C(\mathbf{x}_k)^T \mathbf{M}^{-1} \nabla C(\mathbf{x}_k) + \tilde{\alpha}}$$






