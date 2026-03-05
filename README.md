# Modelling the Oral Epithlium in PhysiCell

## Table of contents
- [Introduction](#introduction)
- [Presentation of the Oral Epithelium](#presentation-of-the-oral-epithelium)
  - [Structure and composition](#structure-and-composition)
  - [Properties](#properties)
- [Key Mechanisms in Epithelium Homeostasis](#key-mechanisms-in-epithelium-homeostasis)
  - [Division model](#division-model)
  - [Treadmill Movement](#treadmill-movement)
  - [Differentiation](#differentiation)
  - [Death Regulation](#death-regulation)
- [PhysiCell Overview](#physicell-overview)
  - [A Quick Guide to the PhysiUniverse](#a-quick-guide-to-the-physiuniverse)
  - [Agents and the interactions](#agents-and-the-interactions)
    - [What is a PhysiCell agent?](#what-is-a-physicell-agent)
    - [How do agents interact with other agents?](#how-do-agents-interact-with-other-agents)
    - [How do agents evolve in time?](#how-do-agents-evolve-in-time)
  - [PhysiCell in practice](#physicell-in-practice)
- [Running Epithelium Simulations](#running-epithelium-simulations)
  - [Analysis Pipeline](#analysis-pipeline)
    - [Initialization](#initialization)
    - [Parameters Update](#parameters-update)
    - [Running](#running)
    - [Post Processing](#post-processing)
- [Sensibility Analysis](#sensibility-analysis)
  - [The SALib package](#the-salib-package)
  - [Sobol Analysis: input/output](#sobol-analysis-inputoutput)
    - [In practice the scripts](#in-practice-the-scripts)
    - [Experiments](#experiments)
    - [Results](#results)
- [Parametrization](#parametrization)
  - [Gradient Descent](#gradient-descent)
    - [PhysiCOOL Package](#physicool-package)
    - [Workflow and output](#workflow-and-output)
    - [Modifications to PhysiCell](#modifications-to-physicell)
  - [Bayesian Optimization](#bayesian-optimization)
- [Perspectives](#perspectives)
  - [Implementing Lineage Tracing](#implementing-lineage-tracing)
  - [Implementing Dilution Experiments](#implementing-dilution-experiments)

---

## Introduction
Head and Neck Squamous cell carcinoma (HNSCC) was the 7th most prevalent cancer in 2020, it leads to approximately 450,000 deaths per year. The oral cavity is the most frequently impacted anatomical region. This carcinoma are known as Oral Squamous Cell Carcinoma (OSCC). It is the deadliest among HNSCC. Moreover, it leads to severe disfugurements and functional impairments. Understanding the multiple factors at play in OSCC oncogenesis therefore appears essentail. This would enable better prognostic predictions as well as early detection, which is key in improving patients outcome and quality of life. [these yannick]
A lot remains unknown about OSCC carcinogenesis. **We aim to decipher the conditions of OSCC establishment and maintenance.**

We assume that mechanical effects as well as chemical signal (through contact and diffusion) from differentiated cells are necessary to ensure the treadmill behavior of the epithelium observed in vivo. To this end, we built a Agent Based Modelling (ABM) framework of an healthy epithelium, using PhysiCell<sup><a href="#ref16">16</a></sup>, to integrate spatial considerations. PhysiCell is a physics-based multicellular simulator. We first attempt to identify minimal conditions for homestasis emergence in a PhysiCell simulation. We will try to determine if we can build a model encompassing the key biological propeties epithelium. Through sensibility analysis, we will then establish which parameters influence epithelium stability. Finally, we will try to asses the existence of a parameter space enabling epithelium stability and coherent with dynamics of oral epithelium in vivo and in vitro. 

The following [poster](figures/poster.png) summarizes the project.

## Presentation of the Oral Epithelium

The oral epithelium is complex stratified tissue formed of multiple layers of different type of cells. Despite, renewing very fast (~15 days for total renewal) its topology is highly conserved in time. This robustness is essential as this epithelium is the first barrier to physical and chemical insults from mastication, food, and microorganisms. It remains unclear how homeostasis is established and maintain in time. Pathways for epithelial stem cell division regulation are still under investigation. Understanding epithelium maintenance and renewal is essential to understand the mechanisms of cancer maintenance.
![Schematic representation of the oral epithelium](figures/epi_schema.png) *Shematic representation of oral mucose structure from Int. J. Mol. Sci. 2021, 22, 7821. [https://doi.org/10.3390/ijms22157821](https://doi.org/10.3390/ijms22157821)*

### Structure and composition
In keratinized oral mucosa, the epithelium stratifies into four layers (stratum basale, stratum spinosum, stratum granulosum, and stratum corneum), while on the other hand, it consists of three layers (stratum basale, stratum filamentosum, and stratum distendum) in non-keratinized tissue. (Mu, 2024). In this model, we will focus on non-keratinized epithelhium. Moreover, we will focus on a simplified model of this epithelium. We will consider 2 layers in the epithelium (2 cells type) rather than 3, dividing epithelium basal cell and differentiated epithelium cells. It is also important to note that we will not take into account morphological deformations of cells that elongates as they differentiate. Indeed, PhysiCell only allows for spherical cell agents. This difference in morphology will certainly affect the dynamics of homeostasis, however we assume this difference should not prevent the emergence of homeostasis and robustness under similar cues.<br>
Underneath the epithelium, there is a layer of connective tissue called the lamina propria composed of blood vessels, nerves, fibroblasts, macrophages, mast cells, and inflammatory cells fibers all immersed in an amorphous substance formed by proteoglycans and glycoproteins.[citation Brizuela] We can distinguish two layers in the connective tissue. On the top we find thin collangen fiber irregularly oriented, connecting with th eepithelium. Below, we find thicker collagen fiber oriented parallel to the surface.[citation Brizuela]<br>

### Key Mechanisms in Epithelium Homeostasis

The oral epithelial cells are frequently replaced by cell division, around each 14 to 21 days. This is because the oral cavity is constantly exposed to high functional demands, which necessitate frequent turnover. The replenishment process starts in the basal layer and then a differentiation and migration process. Proliferation drives the delamination of nearby cells through a density dependent mechanisms. This coordination of behaviors is thought to maintain stem cell numbers and local density over time, allowing constitutive stem cell divisions to be compensated by the later exit of neighboring cells via delamination.[mesa citation] Tissue homeostasis requires differentiation and desquamation at the epithelial surface to be matched by cell division. Many factors, including aging and disease, can alter this balance so that an epithelium may become thicker (hyperplastic) or thinner (atrophic) than normal.[squier citation]

#### Division and differentiation
To maintain a constant number of proliferating cells, on average each cell division must generate one daughter that will go on to divide and one that will differentiate after first exiting the cell cycle. However, the nature  of the dividing cell population was subjected to controversy.
Lineage tracing has ruled out older deterministic models of a proliferative hierarchy of asymmetrically dividing stem cells generating ‘transit amplifying’ cells that undergo a fixed number of divisions prior to differentiation.
Tracking of labeled cells in transgenic mice revealed that the most likely division model was the **single progenitor division hypothesis**: All dividing keratinocytes are functionally equivalent and generate dividing and differentiating daughters with equal probability.<br>

Additionally, **the orientation of these divisions are not random, they occur predominantly in the basal layer**. Cells that exits the basal layer eventually differentiate. This stratification occurs trough a mechanical process, cells are squizzed out from the basal layer due to other division occuring in the basal layer. The switch between division and differentiation is modulated by the presence of multiple factors, such as extracellular calcium[squier citation]. Cells in the basal layer are attached by integrin-containing focal adhesions, and differentiation involves migration with a loss of integrin expression and an increase in cadherin-mediated adhesion via close intercellular junctions or desmosomes. **Therefore, here we assume that stratification is the main driver of differentiation,** even though basal layer cell can differentiate at a fixed rate.

#### Death Regulation
We assume that cell death occurs predominantly through **anoikis, and apoptosis**. Anoikis is a subset of apoptosis triggered by inadequate or inappropriate cell–matrix contacts. It maintains the correct cell number of high-turnover epithelial tissues. [Frisch citation]The cell–ECM and cell–cell adhesion is mediated by cell surface adhesion molecules. Integrins are major adhesion molecules that mediate cell–ECM contact. These molecules can sense the mechanical forces arising from ECM and convert the stimuli to downstream signals modulating cell viability. Furthermore, integrins regulate the activity of many growth factor receptors (GFR). On the other hand, GFR
engagement also cross-talk with the integrin-activated survival signaling. Cooperative function between integrins and GFR is believed to be necessary for proper cell survival and tissue homeostasis [1]. In the absence of attachment to ECM, cells undergo an intrinsically programmed cell death.[zhong citation]
Moreover, the number of cell in the epithelium is regulated through differentiation signal, differentiated cells favor stem cell differentiation, therefore in absence of a sufficient number of differentiated cells, division of basal cell will be favored.

## PhysiCell Overview

PhysiCell is an open source physics-based cell simulator for 3-D multicellular systems. Interactions between agents are modeled as set of potential functions of simplified mechanical interactions.[PhysiCell citation].

### A Quick Guide to the PhysiUniverse
This is only a brief introduction to the PhysiCell framework. For details and tutorials please refer to the official documentation and tutorials:

- Installing PhysiCell:
  - [Linux](https://github.com/physicell-training/ws2023/blob/main/setup/physicell_setup_poweruser_linux_v20250227a.pdf) 
  - [Windows](https://github.com/physicell-training/ws2023/blob/main/setup/PhysiCell_ws2023_Windows_setup.pdf)
  - [MacOS](https://github.com/physicell-training/ws2023/blob/main/setup/PhysiCell_ws2023_macOS_setup.pdf)

- Tutorial:
  - [2023 Workshop](https://github.com/physicell-training/ws2023/tree/main/sessions)
- PhysiCell Addons and tools:
  - [PhysiMeSS](https://github.com/PhysiMeSS/PhysiMeSS) : Modelling the ExtraCellular Matrix (ECM) as agents
  - [PhysiCool](https://physicool.readthedocs.io/en/latest/) : Module to perform gradient descent
  - [PhysiS&S](https://github.com/smilies-polito/PhysiSandS) : Module to start and resume a simulation from a saved state
  - [UQPhysiCell](https://uq-physicell.readthedocs.io/en/latest/) : Sensitivity analysis and Bayesian optimization 
  - [PhysiCell Studio](https://github.com/PhysiCell-Tools/PhysiCell-Studio) : GUI for PhysiCell
  - [PhysiPKPD](https://github.com/drbergman/PhysiPKPD): A pharmacokinetics and pharmacodynamics module for PhysiCell
  - [BIWT](https://github.com/drbergman/BIWT-Paper): a bioinformatics walkthrough for embedding spatial multiomics in agent-based models for virtual cells

- [Non-exhaustive Bibliography of tools and addons of PhysiCell but also papers using the PhysiCell framework](figures/physicell_bibli.bib)

### Framework

PhysiCell is designed to for building multicellular simulations that investigate the relationship between (diffusional) substrate limitations, multicellular biochemical communication, and essential phenotypic processes.<br> 
It uses a latticefree, physics-based approach. It provides optimized, biologically realistic functions for key cell behaviors, including: cell cycling (multiple models for in vitro and in vivo-focused simulations), cell death (apoptosis and necrosis), volume regulation (fluid and solid biomass; nuclear and cytoplasmic sub-volumes), motility, and cell-cell mechanical interactions. It is fully coupled to a fast multi substrate diffusion code (BioFVM) [cite BioFVM] that solves for vectors of diffusing substrates, so that users can tie cell phenotype to many diffusing signals.<br> 
Diffusive biotransport occurs at relatively fast time scales (on the order of 0.1 min or faster) compared to cell mechanics (* 1 min) and cell processes (* 10 to 100 min or slower). PhysiCell takes advantage of this by using three separate time step sizes (Δtdiff, Δtmech, and Δtcells). In particular, the cell phenotypes and arrangement (operating on slow time scales) can be treated as quasi-static when advancing the solution to the biotransport PDEs, so BioFVM can be called without modification with the cell arrangements fixed. PhysiCell provides default time step sizes that should suffice for typical applications in cancer biology and tissue engineering. 

![Schematic representation of PhysiCell functionning](figures/physicell_schema.png) *Figure from [PhysiBoSS 2.0 paper citation] shows a schematic representation of an agent-based model of a 3D multicellular system in a microenvironment defined by the domain divided into fixed volumes together with examples of different intracellular models including signalling network, metabolism and cell cycle. b. depicts examples of the different simulators and time scales for a multi-scale model including, the Δtdiff time scale where diffusion, uptake and secretion processes are updated; Δtmec where the mechanics (movement and physical interactions) are updated; Δtcell in which cell processes such as volume, cell cycle and death models are updated; and Δtreg the regulatory time scale in which Boolean models are updated.*

### Agents and Interactions

#### More details about PhysiCell agents
<p align="center">
  <img src="figures/cell_rules_1.png" alt="Agent in PhysiCell" style="width:50%;" />
</p>

*Schematic representation of a Cell agent in PhysiCell (figure from [Cell rule paper citation])* 

#### Agents mechanics

As defined in PhysiCell and a previous model presenting agent-based cell mechanics , we consider that the unknown cell morphology can be approximated by a spherical cell of equivalent volume. Cells are able to adhere or be repelled by other agents in a radius Ra. Each cell is attributed the position of its center, a velocity and a radius, that can evolve in time. To account for cell deformation, they are able to partially overlap with other agents. The user can set the deformability ability for each cell type. Cells can move at a user defined speed, and the direction of migration depends on chemotactics signals and stochastic brownian movement. Cells’ velocity is modified upon interactions. We consider inertia negligible, as it has been observed experimentally for cells. Hence, we make the assumption:
$m_i v̇_i ≈ 0 (1)$

Thus, once an agent is no longer subjected to forces, its motion ceases in the order of a timestep.

#### How do agents move in time?
At each timestep, we determine the position of an agent by computing its velocity. To solve for each agent's velocity, we use Newton's second law of motion:

The inertialess assumption gives us: $m_i * dv_i/dt ≈ 0$

Newton's second law for agent motion becomes:

$m_i * dv_i/dt = Σ(F_{cca}^{ij} + F_{ccr}^{ij}) + F_{loc}^i + F_{drag}^i$

where the first sum represents cell-cell interactions (adhesion and repulsion forces) and the second sum represents cell-fibre interactions.

We take *i* a cell among the *N(t)* agents at time *t*, with a velocity **$v_i$** and a mass $*m_i*$.  $F_{cca}^{ij}$ and **$F_{ccr}^{ij}$** are respectively the force of adhesion and repulsion on *i* exerted by a cell agent *j* in proximity of *i*. For details about cell-cell interactions, we refer you to Macklin & al. DCIS model in which they were defined. **$F_{loc}^i$** corresponds to the motility force of cell's *i*. **$F_{drag}^i$** represents the drag of the microenvironment, that we can describe as **$F_{drag}^i = -ν_iv_i. ν_i$** was not explicitly described in PhysiCell. The user is rather expected to adjust the different cell's mechanic parameter to account for it.

Given the equations above and the inertialess assumption, we obtain:

$v_i = (1/ν_i) * [Σ(F_{cca}^{ij} + F_{ccr}^{ij}) + F_{loc}^i]$

The agent's position is then updated using the second-order Adam's Bashforth discretization.

#### How do the cells' phenotype evolve in time?
Phenotype updates, cell cycle progression, death, and movement.

### PhysiCell in practice

I made modifications to the PhysiCell 1.14.2 (see [changes.md](changes.md)) to integrate the PhysiS&S module, prevent crashes at runtime, constrain cell types in the simulated space, and modified XML parsing to facilitate user interaction. I also created a number of custom functions to initialize simulations and to choose the orientation of division. However, **any project working for the published version of PhysiCell 1.14.2 should compile and run in this version as well** and produce the same results.

You can look at this presentation, from a 2023 PhysiCell workshop, to understand the structure of the PhysiCell repository, it is up to date with the 1.14.2 version (no major changes): [Working with
PhysiCell Projects](https://github.com/physicell-training/ws2023/blob/main/sessions/session_01/PhysiCell_ws2023_Session01.pdf)

Moreover, you can look at this file: [Running epithelium simulations](tutorial/epithelium_simulation.md)

## Running Epithelium Simulations

### The different user projects
User projects store different configurations. The cell types and the environment remain the cell but the cells position, the cell rule files and some parameters differ. Fpr every case, we model four different types of agent:
- conjonctif: These represent the connective tisue below the epithelium.
- membrane: These represent the top layer of the connective tissue below the epthelium. It is primarly composed of collagen fibers.
- epi_basal: These are the epithelium stem cell. They are attached to the agents representing the basement membrane, they are the only cells dividing. Divisions occur preferentially parallel to the basement mebrane.
- epi_inter: These represent the differentiated cells, they lose there ability to divide. They can undergo apoptosis. 

Additionally, whe have a div_inhib substrate that can be secreted by the epi_inter cells.<br>

test_start_and_stop models a growing epithelium, permeability only the membrane agents. Tthe other projects simulate full epithelium, with a setting file to generate a four layer epithelium and save the data when the epithelium reaches the desired size and an other setting file to run the simulation from he saved data.


### Analysing simulations
![Pipeline Scheme](figures/scheme_pipe.png)

### In practice the scripts

To run such a pipeline we define a number of scripts defned in `scripts_sensitivity_analysis`. We have a script with all general functions, a script for plotting univariate hill functions. Then, we have a number of functions and associated run scripts, repectively used to define specific metrics to each sensitivity analysis and to execute the said analysis.<br>
One sensitivity analysis run can be parallelized. We execute multiple instances of one project for different sets of parameters. This is accomplished with a bash file `launch_run.sh`. This bash file enable to run different instances of the same exe file usingonly one PhysiCell repository. It workks by copying files that will be modified in temporary folders and then regrouping the results in one output folder.

## Sensitivity Analysis
For one simulation in PhysiCell, more than one hundred parameters needs to be set. The number of parameters directly depends from the complexity of the model, specifically the number of cell types. Here we will focus on only 4 cell types, which corresponds to about 120 parameters. For each parameter, we use at start default PhysiCell values, those are set from litterature or from experiments with PhysiCell. We therefore need to drastically reduce model dimensionality to be able to explore parameter space and establish conditions of healthy epithelium stability. To do so, we will perform Global Sensitivity Analysis (GSA). We establish a pipeline following the above scheme for sensitivity analysis. Sensitivity analysis is “the study of how the uncertainty in the output of a mathematical model or system (numerical or otherwise) can be apportioned to
different sources of uncertainty in its inputs.”  The sensitivity of each input is often represented by a numeric value, called the ***sensitivity index***.

We will use the SaLib python package to perform the sensibility analysis.[SaLib citation] This package supports several methods for sensitivity analysis. We will use the Sobol method, a variance based approach [Sobol citation, Saltelli and Annoni 2010 citation]. It is a common method, enablling to explore nth order interactions.This method relies on the variance decomposition of the model output and attribution to a parameter or a set of parameter of the model. The relative influence of the different parameters on the output is summarized by a set of indices. As a output of the Sobol method in SaLib we obtain per parameter or group of parameter: 

1. **First-order indices (S1):** measure the contribution to the output variance of a single
model input.
2. **Second-order indices (S2):** measure the contribution to the output variance caused by
the interaction of two model inputs.
3. **Total-order index (ST):** measure the contribution to the output variance caused by
a model input, including both its first-order effects and all higher-order interactions.

Along with each indices, ST, S1, and S2 we have a corresponding confidence intervals, by default with a confidence level of 95%.

#### Brief overview of Sobol indices computation:
![Schematic Representation of the Sobol Method](figures/sobol_schema.png) (Schematic representation from [Weerasinghe] ).
The method implemented in SaLib is not the original Sobol algorithm published in 1990 [Sobol Citation] but rather an improved more recent version of the algorithm (Saltelli 2010)
1. We define a space of inputs, parameters or group of parameters and their bounds. We consider the model as a black box, written above as the function F, therefore GSA works with any form model from ODE to ABM.
2. Then a quasi-random sampling method is used to obtain an independent uniformly distributed set of inputs within the hypercube. This enables us to write the model output as:

$$Y = F_0 + \sum_{i=1}^{d} F_i(x_i) + \sum_{i<j}^{d} F_{ij}(x_i, x_j) + \cdots + F_{1,2,\dots,d}(x_1, x_2, \dots, x_N)$$

From this equation we can derive the variance of the output:

$$V(Y) = \sum V_i + \sum_i \sum_{i<j} V_{ij} + \cdots + V_{1,2,\dots,N}$$

and also compute the partial variance which for the first order is written as:

$$V_i = V_{x_i}[E_{x_{\sim i}}(x_i)]$$

where $x_{\sim i}$ indicates all the variables except $x_i$.

This gives the Sobol index of the first order:

$$S_i = \frac{V_i}{\text{Var}(Y)}$$

Sobol indices of higher order are of the same form (partial variance over global variance). The form of these indices displays the importance of evaluating independent parameters or groups of parameters.

### Experiments
We are going to run a number of sensitivity analysis, to identify the parameters essential for the epithelium properties we aim to recreate:
  - Permeability: Cells (like immune cells) should be able to move through the different layers of the epithelium without compromising the overall stratification.
  - Structure: As it grows ( is submitted to other kind of mechanical stimuli), the epithelium should preserve its four distinct layers
  - Localized division: Division should not occur anywhere in the tissue, they should occur predominantly in the basal layer
  -Localized death: Death should occur predominantly at the highest position of the epithlium in a shedding process.
  - Size stability: Balancing of division and death should result in an epithelium of a stable thickness.

To evaluate these different properties, we need to define a number of metrics. Then, we can compute the Sobol indices for these metrics.
 - Structure: 
   - Conservation of membrane agents' neighbors: The basement membrane should be able to deform but not rupture. As it is modelled like a layer of two agent thickness, we expect the membrane agents to not reorganized drastically in space relatively to one another. 
   - Non-mixing of cell types: Conjonctives cell should not be neighbors with upper layers cell types.
   - Basal Cell positions: Basal cells are expected to be in the neighborhood of membrane cells

Those metrics return boolean value. They are evaluated in growing epithelium. Division exhert forces on the surrounding cells which lead to loss of stratification.
- Permeability:
  - Speed of crossing: The time a cell attracted to a chemotactic signal takes to cross the basemement membrane.
The structure metrics also need to be computed. Permeability is not relevant in a non-stratified tissue
- Localized division: 
  - Position of the division events
- Localized death:
  - Median position of cell death events
  - Number of apoptosis events
- Stability:
  - Cells count
  - Thickness of the epithelium when simulation ends
  - Mean growth rate

Metrics returning boolean values are not ideal for sensitivity analysis with Sobol. It could be worth identifying addtional metrics.

We won't run sobol analysis on all possible parameter. For each epithelium property, we use information from the litterature to choose a set of potentially relevant parameters. A file summarizing all the parameters tested, their bounds can be found [here](scripts_sensibility_analysis/parameters_output_description.md).

#### For structure: 
We ran one set of seven parameters. From litterature [citation], we expect adhesion structure between the different cell types to be essential in stratification maintenance. We evaluate the importance of **cell adhesion affinities** between the different cell types. As described in the supplemental information of the PhysiCell artcle [citation PhysiCell], cell-cell adhesive forces are described as:
$$F^{ij}_{cca} = −C_{cca}A_iA_j \nabla \phi_{n_{cca},R_{i,A}+R_{j,A}} (x_j − x_i)$$
where $\phi _{n,R_a}(r)$ is a potential function for adhesive interactions dependant of distance between cells, $R_a$ is the maximum adhesion distance
Cell adhesion affinity is used to compute the effective adhesion coefficient. It controls how much cells forming an adhesion stay close to one an other.  Moreover, it controls the probability of attachment of cells.For each cell $j$ in the neighbors list, it forms an attachment with probability: 
  
$$\textrm{Prob attach } i \textrm{ to cell } j = \textrm{adhesion affinity}_j \cdot \textrm{attachment rate}_i \cdot \Delta t$$
    
The attachment is only formed if both cell $i$ and $j$ have not exceeded their maximum number of attachments.

We fix the other adhesion parameters, such as maximum adhesive interaction distance, adhesion strength.

We evaluate these parameters in a growing epithelium:
<p style="display:flex; gap:2%; align-items:flex-start;">
  <img src="figures/loss_integrity.gif" alt="Loss integrity gif" style="width:30%; height:auto;" />
  <img src="figures/integrity.gif" alt="Integrity preserved gif" style="width:30%; height:auto;" />
</p>

We need to simulate a growing epithlium that keeps its layers, so that we can generate a stable epithlium to start simulation from. We use as a metric the number of membrane breaks, a break event corresponds to a loss of neighborhood for a membrane cell.<br>
**Results:**<br>

**Total-order indices (ST)**

| Parameter | ST | ST_conf |
|---|---:|---:|
| adhesion affinity epi_basal epi_basal | 0.053483 | 0.043232 |
| adhesion affinity epi_basal epi_inter | 0.043062 | 0.028244 |
| adhesion affinity epi_inter epi_basal | 0.043012 | 0.025796 |
| adhesion affinity epi_basal membrane | 0.172296 | 0.109545 |
| adhesion affinity membrane epi_basal | 0.153518 | 0.080007 |
| adhesion affinity membrane membrane | 0.740994 | 0.279493 |
| adhesion affinity membrane conjonctif | 0.094185 | 0.059679 |
| adhesion affinity conjonctif membrane | 0.095021 | 0.051743 |

**First-order indices (S1)**

| Parameter | S1 | S1_conf |
|---|---:|---:|
| adhesion affinity epi_basal epi_basal | 0.024983 | 0.045464 |
| adhesion affinity epi_basal epi_inter | -0.026464 | 0.026157 |
| adhesion affinity epi_inter epi_basal | -0.021823 | 0.032449 |
| adhesion affinity epi_basal membrane | 0.037489 | 0.058834 |
| adhesion affinity membrane epi_basal | -0.007186 | 0.058452 |
| adhesion affinity membrane membrane | 0.585717 | 0.224913 |
| adhesion affinity membrane conjonctif | 0.031189 | 0.054848 |
| adhesion affinity conjonctif membrane | 0.062607 | 0.051264 |

**Second-order indices (S2):**

| Parameter pair | S2 | S2_conf |
|---|---:|---:|
| (adhesion affinity epi_basal epi_basal, adhesion affinity epi_basal epi_inter) | 0.026376 | 0.045729 |
| (adhesion affinity epi_basal epi_basal, adhesion affinity epi_inter epi_basal) | 0.038370 | 0.052539 |
| (adhesion affinity epi_basal epi_basal, adhesion affinity epi_basal membrane) | 0.034536 | 0.053524 |
| (adhesion affinity epi_basal epi_basal, adhesion affinity membrane epi_basal) | 0.021657 | 0.054782 |
| (adhesion affinity epi_basal epi_basal, adhesion affinity membrane membrane) | -0.062009 | 0.076652 |
| (adhesion affinity epi_basal epi_basal, adhesion affinity membrane conjonctif) | 0.013398 | 0.051786 |
| (adhesion affinity epi_basal epi_basal, adhesion affinity conjonctif membrane) | 0.010449 | 0.053976 |
| (adhesion affinity epi_basal epi_inter, adhesion affinity epi_inter epi_basal) | 0.045803 | 0.049770 |
| (adhesion affinity epi_basal epi_inter, adhesion affinity epi_basal membrane) | 0.040690 | 0.057924 |
| (adhesion affinity epi_basal epi_inter, adhesion affinity membrane epi_basal) | 0.043837 | 0.048685 |
| (adhesion affinity epi_basal epi_inter, adhesion affinity membrane membrane) | 0.020143 | 0.086420 |
| (adhesion affinity epi_basal epi_inter, adhesion affinity membrane conjonctif) | 0.031252 | 0.048591 |
| (adhesion affinity epi_basal epi_inter, adhesion affinity conjonctif membrane) | 0.030662 | 0.054175 |
| (adhesion affinity epi_inter epi_basal, adhesion affinity epi_basal membrane) | 0.044334 | 0.050421 |
| (adhesion affinity epi_inter epi_basal, adhesion affinity membrane epi_basal) | 0.054166 | 0.065728 |
| ( adhesion affinity epi_inter epi_basal,  adhesion affinity membrane membrane) | -0.004233 | 0.084227 |
| ( adhesion affinity epi_inter epi_basal,  adhesion affinity membrane conjonctif) | 0.040303 | 0.050912 |
| ( adhesion affinity epi_inter epi_basal,  adhesion affinity conjonctif membrane) | 0.054067 | 0.059146 |
| ( adhesion affinity epi_basal membrane,  adhesion affinity membrane epi_basal) | 0.018064 | 0.106036 |
| ( adhesion affinity epi_basal membrane,  adhesion affinity membrane membrane) | 0.001941 | 0.140764 |
| ( adhesion affinity epi_basal membrane,  adhesion affinity membrane conjonctif) | -0.005040 | 0.100469 |
| ( adhesion affinity epi_basal membrane,  adhesion affinity conjonctif membrane) | -0.006023 | 0.102180 |
| ( adhesion affinity membrane epi_basal,  adhesion affinity membrane membrane) | 0.235555 | 0.331385 |
| ( adhesion affinity membrane epi_basal,  adhesion affinity membrane conjonctif) | 0.116889 | 0.158891 |
| ( adhesion affinity membrane epi_basal,  adhesion affinity conjonctif membrane) | 0.143828 | 0.197991 |
| ( adhesion affinity membrane membrane,  adhesion affinity membrane conjonctif) | 0.139392 | 0.403463 |
| ( adhesion affinity membrane membrane,  adhesion affinity conjonctif membrane) | 0.067131 | 0.408451 |
| ( adhesion affinity membrane conjonctif,  adhesion affinity conjonctif membrane) | -0.031389 | 0.086157 |

**Adhesion affinity between membrane cells is the most dominant parameter for layer maintenance.**

#### For permeability:
We look at the movement of a cell attracted by a chemotactic signal through the membrane alone. We evaluate membrane cell adhesion affinity, attracted cell adhesion affinity, migration bias and attracted cell's volume.

<p style="display:flex; gap:2%; align-items:flex-start;">
  <img src="figures/no_crossing.gif" alt="No crossing gif" style="width:30%; height:auto;" />
  <img src="figures/crossing.gif" alt="Crossing gif" style="width:30%; height:auto;" />
</p>

**Total-order indices (ST):**

| Parameter | ST | ST_conf |
|---|---:|---:|
| adhesion_affinity membrane membrane | 0.312576 | 0.117755 |
| adhesion_affinity membrane attracted | 0.250061 | 0.109815 |
| adhesion_affinity attracted membrane | 0.375092 | 0.120881 |
| migration_bias attracted | 0.922100 | 0.166867 |
| volume total attracted | 0.328205 | 0.128530 |

**First-order indices (S1):**

| Parameter | S1 | S1_conf |
|---|---:|---:|
| adhesion_affinity membrane membrane | -0.036142 | 0.151516 |
| adhesion_affinity membrane attracted | -0.043468 | 0.131743 |
| adhesion_affinity attracted membrane | -0.038584 | 0.154831 |
| migration_bias attracted | 0.798291 | 0.187767 |
| volume total attracted | -0.024176 | 0.142589 |

**Second-order indices (S2):**

| Parameter pair | S2 | S2_conf |
|---|---:|---:|
| (membrane membrane, membrane attracted) | 0.103297 | 0.196148 |
| (membrane membrane, attracted membrane) | 0.103297 | 0.194847 |
| (membrane membrane, migration_bias attracted) | -0.084249 | 0.184272 |
| (membrane membrane, volume total attracted) | 0.072039 | 0.198889 |
| (membrane attracted, attracted membrane) | 0.072283 | 0.198822 |
| (membrane attracted, migration_bias attracted) | -0.115263 | 0.196659 |
| (membrane attracted, volume total attracted) | 0.072283 | 0.179976 |
| (attracted membrane, migration_bias attracted) | -0.136996 | 0.191302 |
| (attracte membrane, volume total attracted) | 0.113065 | 0.190951 |
| (migration_bias attracted, volume total attracted) | -0.026618 | 0.221700 |

The only parameter which seem to have an influence is the migration bias. Mechanical properties don't come into play and hinder the attracted cells movement when considering only the membrane. This is because cells are in "the void" per say. There is no drag exherted by the environment only by the other cells.

**We do the same in a complete epithelium**. We look therefore at more adhesion affinities parameter, volume and speed of the migrating cell, and additionally repulsion strength coefficient, and attachment rate.
We define groups for the equivalent parameters:

- Group_m1: adhesion_affinity membrane membrane" 
- Group_m2: adhesion affinity membrane attracted, adhesion affinity attracted membrane
- Group_m3: adhesion_affinity membrane conjonctif, adhesion_affinity conjonctif membrane
- Group_c1: adhesion_affinity conjonctif conjonctif
- Group_sattr: speed attracted
- Group_vattr: volume total attracted
- Group_m_att: attachment_rate membrane
- Group_c_att: attachment_rate conjonctif
- Group_at_rep: cell_cell_repulsion_strength attracted 
- Group_m_rep: cell_cell_repulsion_strength membrane

**Total-order indices (ST):**

| Parameter | ST | ST_conf |
|---|---:|---:|
| Group_m1 | 0.921811 | 0.325106 |
| Group_m2 | 0.638683 | 0.217890 |
| Group_m3 | 0.928395 | 0.388862 |
| Group_c1 | 0.974486 | 0.298880 |
| Group_sattr | 1.060082 | 0.520448 |
| Group_vattr | 0.651852 | 0.274817 |
| Group_m_att | 0.941564 | 0.359918 |
| Group_c_att | 0.783539 | 0.294527 |
| Group_at_rep | 0.671605 | 0.392607 |
| Group_m_rep | 0.816461 | 0.211981 |

**First-order indices (S1):**

| Parameter | S1 | S1_conf |
|---|---:|---:|
| Group_m1 | 0.018201 | 0.069454 |
| Group_m2 | 0.049677 | 0.048739 |
| Group_m3 | 0.098284 | 0.134930 |
| Group_c1 | 0.060812 | 0.134495 |
| Group_sattr | -0.060704 | 0.286030 |
| Group_vattr | 0.020235 | 0.044743 |
| Group_m_att | 0.050641 | 0.073517 |
| Group_c_att | 0.019271 | 0.050406 |
| Group_at_rep | 0.022269 | 0.079479 |
| Group_m_rep | 0.059848 | 0.070277 |

**Second-order indices (S2):**

| Parameter pair | S2 | S2_conf |
|---|---:|---:|
| (Group_m1, Group_m2) | -0.073981 | 0.094013 |
| (Group_m1, Group_m3) | -0.021306 | 0.091544 |
| (Group_m1, Group_c1) | -0.060812 | 0.117922 |
| (Group_m1, Group_sattr) | -0.021306 | 0.194771 |
| (Group_m1, Group_vattr) | -0.047643 | 0.089879 |
| (Group_m1, Group_m_att) | -0.060812 | 0.107161 |
| (Group_m1, Group_c_att) | -0.034474 | 0.084447 |
| (Group_m1, Group_at_rep) | -0.021306 | 0.138918 |
| (Group_m1, Group_m_rep) | -0.073981 | 0.088766 |
| (Group_m2, Group_m3) | -0.060812 | 0.079232 |
| (Group_m2, Group_c1) | 0.031369 | 0.212730 |
| (Group_m2, Group_sattr) | -0.060812 | 0.104622 |
| (Group_m2, Group_vattr) | -0.034474 | 0.093946 |
| (Group_m2, Group_m_att) | -0.087149 | 0.101227 |
| (Group_m2, Group_c_att) | -0.008137 | 0.133139 |
| (Group_m2, Group_at_rep) | -0.073981 | 0.083105 |
| (Group_m2, Group_m_rep) | -0.047643 | 0.082777 |
| (Group_m3, Group_c1) | -0.064988 | 0.238622 |
| (Group_m3, Group_sattr) | 0.264230 | 0.374499 |
| (Group_m3, Group_vattr) | -0.012313 | 0.153667 |
| (Group_m3, Group_m_att) | -0.078157 | 0.125413 |
| (Group_m3, Group_c_att) | 0.027193 | 0.187336 |
| (Group_m3, Group_at_rep) | -0.038650 | 0.169929 |
| (Group_m3, Group_m_rep) | -0.091325 | 0.169122 |
| (Group_c1, Group_sattr) | 0.078905 | 0.363409 |
| (Group_c1, Group_vattr) | -0.000108 | 0.112562 |
| (Group_c1, Group_m_att) | 0.013061 | 0.125262 |
| (Group_c1, Group_c_att) | 0.065736 | 0.168625 |
| (Group_c1, Group_at_rep) | 0.052567 | 0.189132 |
| (Group_c1, Group_m_rep) | -0.052782 | 0.122943 |
| (Group_sattr, Group_vattr) | 0.280504 | 0.381323 |
| (Group_sattr, Group_m_att) | 0.267335 | 0.451056 |
| (Group_sattr, Group_c_att) | 0.267335 | 0.363536 |
| (Group_sattr, Group_at_rep) | 0.333179 | 0.434566 |
| (Group_sattr, Group_m_rep) | 0.254167 | 0.347841 |
| (Group_vattr, Group_m_att) | 0.058670 | 0.189983 |
| (Group_vattr, Group_c_att) | 0.005995 | 0.175915 |
| (Group_vattr, Group_at_rep) | -0.033511 | 0.111973 |
| (Group_vattr, Group_m_rep) | -0.046680 | 0.113257 |
| (Group_m_att, Group_c_att) | -0.044753 | 0.139299 |
| (Group_m_att, Group_at_rep) | -0.031584 | 0.201606 |
| (Group_m_att, Group_m_rep) | -0.044753 | 0.139299 |
| (Group_c_att, Group_at_rep) | 0.011027 | 0.104925 |
| (Group_c_att, Group_m_rep) | -0.028479 | 0.089803 |
| (Group_at_rep, Group_m_rep) | -0.042611 | 0.104401 |

#### For growth evaluation:
We want to evaluate how the substrate div_inhib secreted by the epi_inter cells control the growth rate of the epithelium. We look at the property of div_inhib and its secretion.

**Total-order indices (ST):**

| Parameter | ST | ST_conf |
|---|---:|---:|
| secretion_rate epi_inter div_inhib | 0.989686 | 1.123561 |
| diffusion_coefficient div_inhib | 0.531357 | 0.828716 |
| decay_rate div_inhib | 0.727531 | 1.091096 |

**First-order indices (S1):**

| Parameter | S1 | S1_conf |
|---|---:|---:|
| secretion_rate;epi_inter div_inhib | 0.475714 | 0.772000 |
| diffusion_coefficient div_inhib | -0.000118 | 0.049567 |
| decay_rate div_inhib | -0.051153 | 0.145835 |

**Second-order indices (S2):**

| Parameter pair | S2 | S2_conf |
|---|---:|---:|
| (secretion_rate epi_inter div_inhib, diffusion_coefficient div_inhib) | -0.522589 | 0.753387 |
| (secretion_rate epi_inter div_inhib, decay_rate div_inhib) | -0.442425 | 0.793593 |
| (diffusion_coefficient div_inhib, decay_rate div_inhib) | 0.101211 | 0.306338 |

It appears that the secretion rate of div_inhib migth be the most relevant parameters to growth control. We additionaly analyse the rule controlling the effect of div_inhib on the division rate of the epi basal cell.
The parameter (cell_rules;1;5) corresponds to the base value of the hill function that describes how div_inhib decrease epi_basal cell division rate. 

**Total-order indices (ST):**

| Parameter | ST | ST_conf |
|---|---:|---:|
| secretion_rat epi_inter div_inhib | 0.005418582 | 0.009935 |
| diffusion_coefficient div_inhib | 0.00006502501 | 0.000148 |
| decay_rate div_inhib | 0.0000005238454 | 0.000001 |
| cell_rules;1;5 | 0.8973821 | 0.200514 |

**First-order indices (S1):**

| Parameter | S1 | S1_conf |
|---|---:|---:|
| secretion_rate epi_inter div_inhib | 0.003938 | 0.009427 |
| diffusion_coefficient div_inhib | 0.000420 | 0.000960 |
| decay_rate div_inhib | -0.000038 | 0.000086 |
| cell_rules;1;5 | 0.903839 | 0.395694 |

**Second-order indices (S2):**

| Parameter pair | S2 | S2_conf |
|---|---:|---:|
| (secretion_rate epi_inter div_inhib, diffusion_coefficient div_inhib) | 0.000075 | 0.007935 |
| (secretion_rate epi_inter div_inhib, decay_rate div_inhib) | 0.000075 | 0.007935 |
| (secretion_rate epi_inter div_inhib, cell_rules;1;5) | -0.002563 | 0.009682 |
| (diffusion_coefficient div_inhib, decay_rate div_inhib) | -0.000420 | 0.000960 |
| (diffusion_coefficient div_inhib, cell_rules;1;5) | -0.000420 | 0.000960 |
| (decay_rate div_inhib, cell_rules;1;5) | 0.000038 | 0.000086 |

#### For apoptsis evaluation:

We evaluate the effect of the rule controlling the death rate in function of the number of attachment on the median position of apoptosis. 

**Total-order indices (ST):**

| Parameter | ST | ST_conf |
|---|---:|---:|
| death_rate epi_inter | 0.400636 | 0.213056 |
| cell_rule;2;5 | 0.937232 | 0.469965 |
| cell_rule;2;6 | 0.729218 | 0.349677 |

**First-order indices (S1):**

| Parameter | S1 | S1_conf |
|---|---:|---:|
| death_rate epi_inter | -0.004370 | 0.203810 |
| cell_rule;2;5 | 0.050811 | 0.354527 |
| cell_rule;2;6 | 0.244015 | 0.258628 |

**Second-order indices (S2):**

| Parameter pair | S2 | S2_conf |
|---|---:|---:|
| (death_rate epi_inter, cell_rule;2;5) | -0.081238 | 0.383975 |
| (death_rate epi_inter, cell_rule;2;6) | -0.261255 | 0.296031 |
| (cell_rule;2;5, cell_rule;2;6) | -0.395250 | 0.589492 |

It is not possible to draw conclusion from the sensitivity analysis. We should consider the three as relevant in the position of death events.


### Perpectives 
If Sobol method was chosen as a first approach, because it is commonly used in analysis for GSA in biological models, it could be worth investigating other approaches. If it appears that different methods converge on the identification of the dominant parameters, they differ in their computational cost and quantitative abilities. [Crusenberry citation]
For our simulations, we are limited by computational time. Simulation runtime are in the order of the minutes to the tens of minutes. Thus, depending on the number of parameters evaluated it can take tens of hours to run all the simulations necessary to compute sobol output. Even then, we can get very poor confidence interval, making it impossible to draw definite conclusion on the most dominant parameter and not allowing any interpretation of the 2nd order indices. **There might be GSA methods requiring less runs to obtain better first order results.**

## Perspectives
The sensitivity analysis is only the first step of the identification of a parameter space enabling the modelling of a stable epithelium.
### Exploring parameter space

The PhysiCOOL module [PhysiCool citation] allows users to create a "black-box model" with three main components:

  - A function that updates the PhysiCell configuration file with new input parameters values;
  - The PhysiCell model;
  - A function that reads the model outputs and computes the desired output metric.

This allows to easily configure and execute PhysiCell simulations to evaluate the effect of the different parameters on the model output. Moreover, PhysiCOOL implements a multilevel parameter sweep class that is aimed at identifying the parameters that best fit a target data set. The parameter sweep considers two PhysiCell parameters, and the user should provide an initial value for each of them. At each level, MultiLevelSweep creates a search grid based on these two values, the number of points per direction and the percentage per direction. These values should be configured by the user.

However, this module was published in 2021. Since, some python used as been deprecated and needs to be updated. Moreover, some features were not developped in PhysiCOOL. It lacks the possibility to update cell rules and not only the xml settings file. The pipeline doesn't allow for the loading of saved data with start and stop

### Bayesian Optimization

Bayesian optimization is used on problems of the form max $x ∈ X f ( x ) {\textstyle \max _{x\in X}f(x)}$, with $X$ being the set of all possible parameters x {\textstyle x}. Bayesian optimization is useful for problems where f ( x ) {\textstyle f(x)}, with f the objective function, is difficult to evaluate due to its computational cost, which is the case here. The Bayesian strategy is to treat it as a random function and place a prior over it. The prior captures beliefs about the behavior of the function. After gathering the function evaluations, which are treated as data, the prior is updated to form the posterior distribution over the objective function. The posterior distribution, in turn, is used to construct an acquisition function that determines the next query point. <br>
The UQPhysiCell module provides a Bayesian optimization framework for the calibration of PhysiCell models. The framework is designed to efficiently find optimal parameter configurations that minimize the discrepancy between model predictions and observed experimental data.

Bayesian optimization requires to define feasible bounds for each parameter, define quantities of interests (QoIs) that can be computed on the simulated data and the biological data. The goal is to determine the Pareto-optimal set of parameter. To do so we also need to define distance metrics to measure the discrepancy between model and observed data, and as well a fitness function, to compute a fitnesss value to be maximized.<br>
These are implemented in the UQ-PhysiCell model, apart from specific QoIs. The most important thing to define is the set of biological data to establish the comparison. 

This module is still being maintained and updated, however it will need to be modified to adapt our simulation as it is not compatible as is. The following modifications are to be made:
- Enabling cell rule modification
- Modification of the pipeline to load saved data at the start of the simulation
- Adding the possibility to modify the saved data
- Adapting the reading of the output file (indexing problems might arise)


## Bibliography

[1] Bai, Yuchen, Jarryd Boath, Gabrielle R. White, Uluvitike G. I. U. Kariyawasam, Camile S. Farah, et Charbel Darido. 2021. « The Balance between Differentiation and Terminal Differentiation Maintains Oral Epithelial Homeostasis ». Cancers 13 (20): 5123. https://doi.org/10.3390/cancers13205123.<br>
[2] Barkley, Dalia, Reuben Moncada, Maayan Pour, et al. 2022. « Cancer Cell States Recur across Tumor Types and Form Specific Interactions with the Tumor Microenvironment ». Nature Genetics 54 (8): 1192‑201. https://doi.org/10.1038/s41588-022-01141-9.<br>
[3] Bergman, Daniel R, Jeanette Johnson, Marwa Naji, et Maxwell Booth. s. d. BIWT: A Bioinformatics Walkthrough for Embedding Spatial Multiomics in Agent-Based Models for Virtual Cells.<br>
[4] Blanpain, Cédric, William E. Lowry, H. Amalia Pasolli, et Elaine Fuchs. 2006. « Canonical Notch Signaling Functions as a Commitment Switch in the Epidermal Lineage ». Genes & Development 20 (21): 3022‑35. https://doi.org/10.1101/gad.1477606.<br>
[5] Brizuela, Melina, et Ryan Winters. 2025. « Histology, Oral Mucosa ». In StatPearls. StatPearls Publishing. http://www.ncbi.nlm.nih.gov/books/NBK572115/.<br>
[6] Byrd, Kevin M., Natalie C. Piehl, Jeet H. Patel, et al. 2019. « Heterogeneity within stratified epithelial stem cell populations maintains the oral mucosa in response to physiological stress ». Cell stem cell 25 (6): 814-829.e6. https://doi.org/10.1016/j.stem.2019.11.005.<br>
[7] Chamoli, Ambika, Abhishek S. Gosavi, Urjita P. Shirwadkar, et al. 2021. « Overview of oral cavity squamous cell carcinoma: Risk factors, mechanisms, and diagnostics ». Oral Oncology 121 (octobre): 105451. https://doi.org/10.1016/j.oraloncology.2021.105451.<br>
[8] Cruchley, Alan T., et Lesley Ann Bergmeier. 2018. « Structure and Functions of the Oral Mucosa ». In Oral Mucosa in Health and Disease, édité par Lesley Ann Bergmeier. Springer International Publishing. https://doi.org/10.1007/978-3-319-56065-6_1.<br>
[9] Crusenberry, Cody & Sobey, Adam & Termaath, Stephanie. (2023). Evaluation of global sensitivity analysis methods for computational structural mechanics problems. Data-Centric Engineering. 4. 10.1017/dce.2023.23.<br>
[9] Damen, Mareike, Lisa Wirtz, Ekaterina Soroka, et al. 2021. « High Proliferation and Delamination during Skin Epidermal Stratification ». Nature Communications 12 (1): 3227. https://doi.org/10.1038/s41467-021-23386-4.<br>
[10] Dawson, D. V., D. R. Drake, J. R. Hill, K. A. Brogden, C. L. Fischer, et P. W. Wertz. 2013. « Organization, barrier function and antimicrobial lipids of the oral mucosa ». International Journal of Cosmetic Science 35 (3): 220‑23. https://doi.org/10.1111/ics.12038.<br>
[11] Dunn, Gavin P., Lloyd J. Old, et Robert D. Schreiber. 2004. « The Three Es of Cancer Immunoediting ». Annual Review of Immunology 22 (1): 329‑60. https://doi.org/10.1146/annurev.immunol.22.012703.104803.<br>
[12] Dvoretzki, Svetlana. s. d. Intéractions Cellules-Cellules dans le cas du cancer de la cavité buccale.<br>
[13] Fischer, Matthias M., Hanspeter Herzel, et Nils Blüthgen. 2022. « Mathematical Modelling Identifies Conditions for Maintaining and Escaping Feedback Control in the Intestinal Epithelium ». Scientific Reports 12 (1): 5569. https://doi.org/10.1038/s41598-022-09202-z.<br>
[14] Frisch, Steven M, et Robert A Screaton. 2001. « Anoikis mechanisms ». Current Opinion in Cell Biology 13 (5): 555‑62. https://doi.org/10.1016/S0955-0674(00)00251-9.<br>
[15] Ghaffarizadeh, Ahmadreza, Samuel H. Friedman, et Paul Macklin. 2015. « BioFVM: an efficient, parallelized diffusive transport solver for 3-D biological simulations ». Bioinformatics 32 (8): 1256‑58. https://doi.org/10.1093/bioinformatics/btv730.<br>
<a id="ref16"></a>
[16] Ghaffarizadeh, Ahmadreza, Randy Heiland, Samuel H. Friedman, Shannon M. Mumenthaler, et Paul Macklin. 2018. « PhysiCell: An Open Source Physics-Based Cell Simulator for 3-D Multicellular Systems ». PLOS Computational Biology 14 (2): e1005991. https://doi.org/10.1371/journal.pcbi.1005991.
[17] Gonçalves, Inês G., David A. Hormuth Ii, Sandhya Prabhakaran, et al. 2023. « PhysiCOOL: A Generalized Framework for Model Calibration and Optimization Of modeLing Projects ». Gigabyte 2023 (février): 1‑11. https://doi.org/10.46471/gigabyte.77.<br>
[18] Groeger, Sabine E., et Joerg Meyle. 2015. « Epithelial Barrier and Oral Bacterial Infection ». Periodontology 2000 69 (1): 46‑67. https://doi.org/10.1111/prd.12094.<br>
[19] Groeger, Sabine, et Joerg Meyle. 2019. « Oral Mucosal Epithelial Cells ». Frontiers in Immunology 10 (février): 208. https://doi.org/10.3389/fimmu.2019.00208.
Herms, Albert, David Fernandez-Antoran, Maria P. Alcolea, et al. 2024.<br> 
[20] Herms, Albert, David Fernandez-Antoran, Maria P. Alcolea, et al. 2024. « Self-Sustaining Long-Term 3D Epithelioid Cultures Reveal Drivers of Clonal Expansion in Esophageal Epithelium ». Nature Genetics 56 (10): 2158‑73. https://doi.org/10.1038/s41588-024-01875-8.
[21] Iwanaga, T., Usher, W., & Herman, J. (2022). Toward SALib 2.0: Advancing the accessibility and interpretability of global sensitivity analyses. Socio-Environmental Systems Modelling, 4, 18155. doi:10.18174/sesmo.18155<br>
[21] Johnson, Jeanette A. I., Daniel R. Bergman, Heber L. Rocha, et al. 2025. « Human interpretable grammar encodes multicellular systems biology models to democratize virtual cell laboratories ». Cell 188 (17): 4711-4733.e37. https://doi.org/10.1016/j.cell.2025.06.048.<br>
[22] Jones, Kyle B., Sachiko Furukawa, Pauline Marangoni, et al. 2019. « Quantitative Clonal Analysis and Single-Cell Transcriptomics Reveal Division Kinetics, Hierarchy, and Fate of Oral Epithelial Progenitor Cells ». Cell Stem Cell 24 (1): 183-192.e8. https://doi.org/10.1016/j.stem.2018.10.015.<br>
[23] Jones, Kyle B., et Ophir D. Klein. 2013. « Oral Epithelial Stem Cells in Tissue Maintenance and Disease: The First Steps in a Long Journey ». International Journal of Oral Science 5 (3): 121‑29. https://doi.org/10.1038/ijos.2013.46.<br>
[24] Kitsukawa, Yoshiaki, Chonji Fukumoto, Toshiki Hyodo, et al. 2024. « Difference between Keratinized- and Non-Keratinized-Originating Epithelium in the Process of Immune Escape of Oral Squamous Cell Carcinoma ». International Journal of Molecular Sciences 25 (7): 3821. https://doi.org/10.3390/ijms25073821.<br>
[25] Koster, Maranke I., et Dennis R. Roop. 2007. « Mechanisms Regulating Epithelial Stratification ». Annual Review of Cell and Developmental Biology 23: 93‑113. https://doi.org/10.1146/annurev.cellbio.23.090506.123357.<br>
[26] Lefort, Karine, et G. Paolo Dotto. 2004. « Notch signaling in the integrated control of keratinocyte growth/differentiation and tumor suppression ». Seminars in Cancer Biology, Notch Signaling and Cancer, vol. 14 (5): 374‑86. https://doi.org/10.1016/j.semcancer.2004.04.017.<br>
[27] Liu, Chunzi, et Gerald G. Fuller. 2023. « Air-Liquid Interface Induced Epithelial Delamination ». Prépublication, Biophysics, août 16. https://doi.org/10.1101/2023.08.14.553291.<br>
[28] Liu, Chunzi, et Gerald G. Fuller. 2025. « Delamination of epithelia induced by air–liquid interfaces ». Molecular Biology of the Cell 36 (8): ar90. https://doi.org/10.1091/mbc.E24-11-0500.<br>
[29] Mackenzie, Ian C., et Murray W. Hill. 1981. « Maintenance of Regionally Specific Patterns of Cell Proliferation and Differentiation in Transplanted Skin and Oral Mucosa ». Cell and Tissue Research 219 (3): 597‑607. https://doi.org/10.1007/BF00209997.<br>
[30] Mesa, Kailin R., Kyogo Kawaguchi, Katie Cockburn, et al. 2018. « Homeostatic Epidermal Stem Cell Self-Renewal Is Driven by Local Differentiation ». Cell Stem Cell 23 (5): 677-686.e4. https://doi.org/10.1016/j.stem.2018.09.005.<br>
[31] Moriyama, Mariko, André-Dante Durham, Hiroyuki Moriyama, et al. 2008. « Multiple Roles of Notch Signaling in the Regulation of Epidermal Development ». Developmental Cell 14 (4): 594‑604. https://doi.org/10.1016/j.devcel.2008.01.017.<br>
[32] Mu, Xindi, Mitsuaki Ono, Ha Thi Thu Nguyen, et al. 2024. « Exploring the Regulators of Keratinization: Role of BMP-2 in Oral Mucosa ». Cells 13 (10): 807. https://doi.org/10.3390/cells13100807.<br>
[33] Ngo, Phuong A., Markus F. Neurath, et Rocío López-Posadas. 2022. « Impact of Epithelial Cell Shedding on Intestinal Homeostasis ». International Journal of Molecular Sciences 23 (8): 4160. https://doi.org/10.3390/ijms23084160.<br>
[34] Papagerakis, Silvana, Giuseppe Pannone, Li Zheng, et al. 2014. « Oral epithelial stem cells – implications in normal development and cancer metastasis ». Experimental cell research 325 (2): 111‑29. https://doi.org/10.1016/j.yexcr.2014.04.021.
Pereira, Diana, et Inês Sequeira. 2021. « A Scarless Healing Tale: Comparing Homeostasis and Wound Healing of Oral Mucosa With Skin and Oesophagus ». Frontiers in Cell and Developmental Biology 9 (juillet): 682143. https://doi.org/10.3389/fcell.2021.682143.<br>
[35] Piedrafita, Gabriel, Vasiliki Kostiou, Agnieszka Wabik, et al. 2020. « A Single-Progenitor Model as the Unifying Paradigm of Epidermal and Esophageal Epithelial Maintenance in Mice ». Nature Communications 11 (1): 1429. https://doi.org/10.1038/s41467-020-15258-0.<br>
[36] Presland, Richard B., et Richard J. Jurevic. 2002. « Making Sense of the Epithelial Barrier: What Molecular Biology and Genetics Tell Us about the Functions of Oral Mucosal and Epidermal Tissues ». Journal of Dental Education 66 (4): 564‑74.<br>
[37] Recka, Nicole, Andrean Simons, Robert A. Cornell, et Eric Van Otterloo. 2024. « Epidermal Loss of PRMT5 Leads to the Emergence of an Atypical Basal Keratinocyte-like Cell Population and Defective Skin Stratification ». Prépublication, bioRxiv, novembre 8. https://doi.org/10.1101/2024.11.08.620904.<br>
[38] Saitoh, Masao, Takuya Shirakihara, Akira Fukasawa, et al. 2013. « Basolateral BMP Signaling in Polarized Epithelial Cells ». PLOS ONE 8 (5): e62659. https://doi.org/10.1371/journal.pone.0062659.<br>
[39] Sakamoto, Kei, Takuma Fujii, Hiroshi Kawachi, et al. 2012. « Reduction of NOTCH1 Expression Pertains to Maturation Abnormalities of Keratinocytes in Squamous Neoplasms ». Laboratory Investigation; a Journal of Technical Methods and Pathology 92 (5): 688‑702. https://doi.org/10.1038/labinvest.2012.9.<br>
[40] Saltelli, Andrea. 2002. « Making best use of model evaluations to compute sensitivity indices ». Computer Physics Communications 145 (2): 280‑97. https://doi.org/10.1016/S0010-4655(02)00280-1.
[41] Saltelli, Andrea, Paola Annoni, Ivano Azzini, Francesca Campolongo, Marco Ratto, et Stefano Tarantola. 2010. « Variance based sensitivity analysis of model output. Design and estimator for the total sensitivity index ». Computer Physics Communications 181 (2): 259‑70. https://doi.org/10.1016/j.cpc.2009.09.018.
[40] Schreiber, Robert D., Lloyd J. Old, et Mark J. Smyth. 2011. « Cancer Immunoediting: Integrating Immunity’s Roles in Cancer Suppression and Promotion ». Science 331 (6024): 1565‑70. https://doi.org/10.1126/science.1203486.<br>
[41] Schroeder, H. E. 1981. Differentiation of Human Oral Stratified Epithelia. Karger Medical and Scientific Publishers.<br>
[42] Smeriglio, Riccardo, Roberta Bardini, Alessandro Savino, et Stefano Di Carlo. 2025. « Start & Stop: A PhysiCell and PhysiBoSS 2.0 Add-on for Interactive Simulation Control ». BMC Bioinformatics 26 (1): 158. https://doi.org/10.1186/s12859-025-06144-x.
Squier, Christopher A., et Mary J. Kremer. 2001. « Biology of Oral Mucosa and Esophagus ». JNCI Monographs 2001 (29): 7‑15. https://doi.org/10.1093/oxfordjournals.jncimonographs.a003443.<br>
[43] Su, Dan, Tadkamol Krongbaramee, Samuel Swearson, et al. s. d. « Irx1 mechanisms for oral epithelial basal stem cell plasticity during reepithelialization after injury ». JCI Insight 10 (1): e179815. https://doi.org/10.1172/jci.insight.179815.
Tadokoro, Tomomi, Xia Gao, Charles C. Hong, Danielle Hotten, et Brigid L. M. Hogan. 2016. « BMP signaling and cellular dynamics during regeneration of airway epithelium from basal progenitors ». Development 143 (5): 764‑73. https://doi.org/10.1242/dev.126656.<br>
[44] Tan, Yunhan, Zhihan Wang, Mengtong Xu, et al. 2023. « Oral Squamous Cell Carcinomas: State of the Field and Emerging Directions ». International Journal of Oral Science 15 (1): 44. https://doi.org/10.1038/s41368-023-00249-w.
« The maintenance of an oral epithelial barrier | 10.1016/j.lfs.2019.04.029_Sci-hub ». s. d. Consulté le 1 octobre 2025. https://www.pismin.com/10.1016/j.lfs.2019.04.029.<br>
[45] Trepat, Xavier, Michael R. Wasserman, Thomas E. Angelini, et al. 2009. « Physical Forces during Collective Cell Migration ». Nature Physics 5 (6): 426‑30. https://doi.org/10.1038/nphys1269.<br>
[46] Van Liedekerke, P., M. M. Palm, N. Jagiella, et D. Drasdo. 2015. « Simulating Tissue Mechanics with Agent-Based Models: Concepts, Perspectives and Some Novel Results ». Computational Particle Mechanics 2 (4): 401‑44. https://doi.org/10.1007/s40571-015-0082-3.<br>
[47] Wang, Sha-Sha, Ya-Ling Tang, Xin Pang, Min Zheng, Ya-Jie Tang, et Xin-Hua Liang. 2019. « The Maintenance of an Oral Epithelial Barrier ». Life Sciences 227 (juin): 129‑36. https://doi.org/10.1016/j.lfs.2019.04.029.<br>
[48] Watt, Fiona M, Soline Estrach, et Carrie A Ambler. 2008. « Epidermal Notch signalling: differentiation, cancer and adhesion ». Current Opinion in Cell Biology, Cell regulation, vol. 20 (2): 171‑79. https://doi.org/10.1016/j.ceb.2008.01.010<br>
[49] Weerasinghe, Gihan & Kannan, Ramaseshan & Bandara, Samila. (2025). Making global sensitivity analysis feasible using neural network surrogates. Data-Centric Engineering. 6. 10.1017/dce.2025.10029.<br>
[49] Zhong, Xiaoling, et Frederick J. Rescorla. 2012. « Cell surface adhesion molecules and adhesion-initiated signaling: Understanding of anoikis resistance mechanisms and therapeutic opportunities ». Cellular Signalling 24 (2): 393‑401. https://doi.org/10.1016/j.cellsig.2011.10.005.<br>

