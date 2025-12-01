# ARCADE: A New Paradigm for a Fully Automated and Interpretable Design of Complex Catalysts

Jianzhuo Wu¹,#, Xue Jia²,#,\*, Yongfeng Guo¹,#, Honglin Li¹, Yu Li¹, Yuhang Wang²,  
Jialu Li¹, Guangsheng Liu¹, Yufei Ding⁴, Hao Li²,\* and Wan-Lu Li¹,³,\*

¹ Aiiso Yufeng Li Family Department of Chemical and Nano Engineering, University of California, San Diego, La Jolla, CA, USA  
² Advanced Institute for Materials Research (WPI-AIMR), Tohoku University, Sendai, Japan  
³ Program in Materials Science and Engineering, University of California, San Diego, La Jolla, CA, USA  
⁴ Department of Computer Science and Engineering, University of California, San Diego, La Jolla, CA, USA  

\*Correspondence: li.hao.b8@tohoku.ac.jp; wal019@ucsd.edu  

\# These authors contributed equally.

---

## Abstract

The rational design of high-performance catalysts is pivotal for modern industry but is often hindered by the time-consuming nature of trial-and-error methods. While computational and data-driven approaches have accelerated materials discovery, significant challenges remain, especially for complex systems like multielement nanoparticles. These challenges include immense computational costs, difficulties in modeling stability and heterogeneous active sites, and a lack of interpretability in "black-box" AI models. To address these multifaceted challenges, we introduce ARCADE (Automated Rational CAtalyst DEsign), a generalizable and transferable framework for the automated and interpretable design of complex catalysts. Our framework integrates a high-throughput workflow that synergizes realistic structure generation, active-site screening, state-of-the-art pH-dependent microkinetic modeling, and short-range order (SRO) guided interpretability, establishing an end-to-end design pipeline that leverages public databases. Using the challenging design of multielement nanoparticle catalysts for CO₂ reduction reaction (CO₂RR) as a representative system, we demonstrate that ARCADE can systematically navigate the design process, from evaluating material stability to identifying site-specific activity and selectivity. Our results show that this methodology substantially accelerate the development cycle for complex catalysts, confirming its feasibility and broad applicability. Crucially, this work establishes a new, fully integrated paradigm for the rational design and discovery of complex materials for challenging chemical reactions.

---

## Introduction

Electrocatalysts play a pivotal role in driving key energy conversion reactions, enabling the efficient transformation of abundant molecules (e.g., H₂O, CO₂, and N₂) into value-added products such as hydrogen, synthetic fuels, and ammonia.¹ Over the past two decades, theory-guided catalyst discovery, such as the development of adsorption-energy scaling relations²¹ and the seminal d-band theory², has provided fundamental insights into heterogeneous catalysis. More recently, advances in data-driven artificial intelligence, including the use of machine learning (ML) to predict electrocatalyst properties and to accelerate the understanding of dynamic catalytic mechanisms through Machine Learning Potentials (MLPs), have further transformed catalyst research and optimization.¹⁹˒²⁰ However, despite these methodological advances, rationally designing complex catalytic systems, particularly multielement nanoparticles, remains limited due to their inherent structural and compositional complexity.⁶˒¹⁷

Multielement nanoparticles, typically composed of three to five or more principal elements, have recently emerged as promising candidates for heterogeneous catalysis owing to their abundant active sites, tunable surface chemistry, and entropy-stabilized structure.⁵⁰˒⁶⁰ Nevertheless, several formidable challenges remain in guiding their rational design as electrocatalysts: (i) the structural complexity driven by random multielement mixing make it difficult to identify thermodynamically stable atomic-packing structure models;²² (ii) the irregular morphology and complex surface structure of nanoparticles make the manual placement of adsorbates difficult;¹⁵ (iii) the lack of interpretable structure–performance relationships that link the atomic structure of multielement nanoparticles with specific catalytic reaction schemes hinders rational catalyst design.⁶¹ To overcome these challenges, several theoretical and ML approaches have been proposed;¹⁹ however, their applicability remains limited due to difficulties in accurately capturing structural disorder and chemical complexity.

Traditional nanoparticle modeling relies on Wulff construction and global structure search. Wulff construction offers a simple way to rapidly predict equilibrium shapes from density functional theory (DFT)⁶⁹˒⁷⁰-derived surface energies but neglect nonequilibrium growth, solvent or ligand effects, and multielement reconstruction.³⁻⁵ Global optimization packages such as the USPEX or CALYPSO identify low-energy configurations, yet the search space grows combinatorially, and energy evaluations are costly.⁷˒⁸ More importantly, they mainly utilize symmetry to extract structural features, which is often applied to crystals rather than nanoparticles. In multielement nanoparticles, thermodynamics drives intrinsically nonrandom mixing, which manifests as short-range order (SRO), and the advanced characterization studies have directly visualized SRO and confirmed its significant impact on structural defects and mechanical properties.⁹ Explicitly accounting for SRO is therefore essential for generating realistic multielement nanoparticle structural models. On the other hand, ML interatomic potentials achieve near-DFT accuracy at molecular dynamics speed, enabling large-scale nanosecond simulations that capture realistic structural ensembles and dynamic reconstruction beyond the limits of Wulff or global optimization methods.¹¹˒¹² For nanoparticle modeling, developing a stable structure-search methodology that combines SRO-based structural feature representation with MLP-accelerated optimization is highly promising.

Adsorption-site determination commonly begins with symmetry-based enumeration on ideal slabs (e.g., the Pymatgen’s AdsorbateSiteFinder¹³, CatKit/ASAP¹⁵), which is efficient but presumes crystallinity and becomes intractable as coverage increases or disorder grows. Global co-optimization frameworks that relax adsorbate and surface concurrently using surrogate MLPs or advanced potential energy surface (PES) explorers (Global Optimization using Bayesian Ensemble Estimation (GOFEE)¹⁶, Solid-State Surface Walk (SSW)¹⁸, San Diego Global Minimum Search (SDGMS)⁷⁴) better escape local minima, yet they are computationally demanding and require careful uncertainty management on heterogeneous, low-symmetry nanoparticles. Recent high-throughput studies also highlight the combinatorial explosion when enumerating multi-site coverages and thus resort to stochastic sampling/filters²³˒²⁴ underscoring the need for automated, geometry-aware active-site search beyond ideal surfaces to handle edges, kinks, and disordered nanoparticle facets. Automated and comprehensive active-site search frameworks are thus required to explore possible adsorption configurations and evaluate their adsorbate binding energies. 

Regarding interpretable structure–activity relationships, descriptor-based catalytic volcano models are widely used.²⁵˒²⁶ But the conventional use of a fixed correction term (0.059 × pH) for pH-dependent energetics is inadequate to describe activity on the reversible hydrogen electrode (RHE) scale.⁶⁸ Incorporating pH-dependent thermodynamics/kinetics into descriptor-based analyses yields pH-dependent microkinetic volcano models that more faithfully capture activity evolution under realistic operating conditions and the RHE scale.²⁷⁻³⁰ Additionally, the beneficial effects of SRO on catalytic properties, including selectivity and activity, are increasingly being recognized.³¹⁻³³ Notably, even though enhancements in catalytic performance related to SRO have been observed, establishing more systematic links between SRO features and catalytic performance is crucial for a comprehensive, mechanism-level understanding of complex multielement nanoparticle electrocatalysts. Ultimately, leveraging SRO-performance relationships to understand catalytic performance would enable the rational and interpretable design of catalysts.

To address these challenges, herein, we propose a new paradigm for automated and interpretable catalyst design, termed ARCADE (Automated Rational CAtalyst DEsign). This generalizable framework encompasses thermodynamically stable nanoparticle generation by introducing SRO-based structural feature, automated active-site identification and binding energy evaluation, selectivity analysis through comparisons of Gibbs free energies among competing reactions, and activity analysis via advanced pH-dependent microkinetic modeling. SRO-derived descriptors are further employed to construct structure-performance relationships. Notably, the entire workflow is implemented within the Digital Catalysis Platform (DigCat)³⁴, where data generation, model construction, and theoretical analyses can be executed automatically. The resulting datasets and mechanistic insights are continuously fed back into the DigCat platform’s database, establishing a closed-loop ecosystem that accelerates the rational design and discovery of next-generation electrocatalysts.

In this work, the electrochemical CO₂ reduction reaction (CO₂RR) was selected as a representative system to demonstrate the efficiency of ARCADE. Among the possible CO₂RR products, formic acid/formate (HCOOH) is particularly attractive due to its high energy density, ease of storage, and wide applicability in chemical, pharmaceutical, and energy-related industries.⁷⁵˒⁷⁶ Guided by initial compositional screening on the DigCat platform, two representative nanoparticle systems, medium-entropy La-Sn-Bi and high-entropy La-Sn-Bi-Co-Ni, were investigated. The most thermodynamically stable structures were rapidly determined via a Monte Carlo (MC)³⁸-based SRO-constrained conditional generative protocol, coupled with MLP-accelerated global search algorithm. After systematically searching all possible adsorption sites, the adsorption energies of \*OCHO, \*COOH, and \*H intermediates were calculated. Both compositions exhibited selectivity toward HCOOH formation and reached the peak of the theoretical volcano model, indicating high catalytic activity superior to that of the corresponding elemental nanoparticles. The constructed SRO-performance score heatmaps successfully decoupled the intertwined dependencies of stability, selectivity, and activity on local structural features, revealing clear structure-function relationships in compositionally complex catalysts. Our analysis reveals that surface Bi enrichment and bulk Co/Ni enrichment govern thermodynamic stability, while local environments enriched in La and Sn are the primary drivers for selectivity and activity, respectively, with specific atomic pairwise patterns further elucidating the origins of these metrics. Ultimately, this work establishes a transformative paradigm for decoding “interpretable complexity,” offering a systematic route to map intricate atomic-level disorder onto predictable catalytic function, thereby paving the way for the rational and efficient design of complex catalysts.

---

## Results

### Overall Workflow of ARCADE

The ARCADE workflow establishes a fully data-driven framework for nanoalloy catalyst discovery, integrating structure generation, active-site identification, selectivity and activity analysis, and structure-performance relationship evaluation. All components of this workflow have been implemented within the DigCat platform, enabling a closed-loop discovery framework. The workflow is shown in Fig. 1. First, we deploy our upgraded ApolloX (Automatic Prediction by generative mOdel for Large-scaLe Optimization of X-composition materials)¹⁰ engine to generate thermodynamically stable multicomponent nanoparticles based on user-defined compositions and parent phases, ensuring realistic structural configurations (Fig. 1a). By integrating MLPs with Particle Swarm Optimization (PSO)³⁹ and MC³⁸ sampling, ApolloX explicitly incorporates realistic SRO, enabling the construction of physically reasonable atomic configurations for multielement nanoparticles.

Next, the workflow automatically identifies all possible adsorption sites based on Voronoi tessellation⁴⁰˒⁴¹ and constructs absorbed structures utilizing the Atomic Simulation Environment (ASE)¹⁴, followed by the structural optimization of key reaction intermediates to determine their binding energies (Fig. 1b). Then, the local adsorption environment of each intermediate can be characterized using SRO descriptors. By identifying similar local environments among different intermediates, the corresponding adsorption free energies (e.g., of \*COOH, \*OCHO, and \*H) can be compared to reveal the preferred reaction pathways and determine product selectivity (Fig. 1c).

Subsequently, pH-dependent microkinetic volcano modeling implemented on DigCat translates these energies into predicted turnover frequencies (TOFs) (Fig. 1d), explicitly capturing the variation of catalytic performance across different pH environments (detailed theoretical derivations are provided in the Supporting Information).

Finally, by aggregating data across diverse compositions and structures, the workflow integrates SRO descriptors, specifically the Local Density Deviation (LDD) and Warren-Cowley (W-C) parameters³⁵⁻³⁷ (formally defined in the Methods Section), with property descriptors such as energies, adsorption free energy differences, and turnover frequencies (TOFs). By calculating a contribution score for specific SRO patterns, ARCADE constructs interpretable structure-performance maps that directly correlate local atomic motifs and elemental mixing patterns with catalytic performance, enabling rational and high-throughput catalyst design (Fig. 1e). Note that while this interpretability module is presented here as the final integration step, its specific analyses are applied iteratively throughout the study to decode stability, selectivity, and activity in turn.

**Fig. 1. ARCADE workflow for data-driven nanoalloy catalysis.**  
*(a)* Nanoparticle generation: Given only candidate compositions and a parent phase specified via our Digital Catalysis Platform (DigCat), the ApolloX multicomponent amorphous-structure modeling engine samples thermodynamically stable nanoalloy structures and returns atomistic nanoparticles. *(b)* Adsorption-site searching: On each nanoparticle, possible adsorption motifs of key intermediates are automatically enumerated and relaxed to obtain binding energies. *(c)* Selectivity analysis: Computed Gibbs free energies for competing intermediates (e.g., vs., , etc.) are cross-compared to map selectivity trends. *(d)* Theoretical activity: pH-dependent microkinetic modeling (implemented on the DigCat) translates adsorption energetics into predicted turnover frequencies and pH-dependent microkinetic volcano relationships. *(e)* Structure-property relationships: Aggregating results across composition/structure space enables interpretable maps that link local atomic motifs and elemental mixing patterns to stability, selectivity, and activity.

> **Note**: Replace this with an actual image link when figures are ready, e.g.  
> `![Fig. 1. ARCADE workflow for data-driven nanoalloy catalysis.](fig1.png)`

---

### Workflow of Structure Generation and Linking Thermodynamic Stability to Short-Range Order

The first step in catalyst design is to identify thermodynamically stable structures for potential catalysts of the target reaction. To achieve this, we integrate initial compositional determination based on the DigCat platform and stable structure generation implemented by the ApolloX. DigCat hosts comprehensive literature and structural data corresponding to diverse chemical reactions and target products, facilitating the initial identification of candidate compositions and structural parent phases. After determining the catalyst components and the structural parent phase by DigCat, we deploy an upgraded version of ApolloX to generate stable structures for complex materials. The main workflow (Fig. 2a) of ApolloX operates as an iterative optimization loop. We treat SRO as a key structural feature, transforming the highly challenging structure-energy iterative optimization into a practical and feasible structure-SRO-energy iterative optimization. SRO is described using the W-C parameter³⁵⁻³⁷, where i and j represent different types of elements. On one hand, a MC scheme³⁸ is adopted to probe elemental mixing and local motif rearrangements, tightly linking SRO to structural configurations. This serves as a Generating Algorithm, enabling the generation of structures corresponding to a given SRO. On the other hand, a PSO algorithm³⁹ is employed to tune the SRO-energy relationship. Ultimately, this approach allows for the efficient search and generation of thermodynamically stable structures. The Generating Algorithm (Fig. 2b) is a hybrid, two-stage process. The first stage, which corresponds to the "global search" mentioned in the text, is a rapid greedy local search (Fig. 2c). It performs random atomic swaps and only accepts changes that reduce the SRO error (defined as Δ), quickly producing a set of low-error candidate structures. This candidate pool is then passed to the second stage, or "detailed search," which employs a Metropolis criterion. This Metropolis search allows for a controlled exploration of the configurational space (by probabilistically accepting some higher-error swaps), enabling the algorithm to escape local minima and converge more precisely on the target α values.

Focusing on the electrochemical CO₂RR to formate as a representative system, we first utilized the DigCat platform to perform a rigorous compositional screening to define the precise catalytic systems for our structural generation and subsequent benchmarking. From 8,439 retrieved literature records on CO₂RR electrocatalysts yielding formate, 30.95% employed metal or alloy catalysts. A comprehensive summary of Faradaic efficiency data reported over the past decade is provided in Table S1, while the statistical distribution of constituent elements involved is presented in Fig. S1. Analysis of these entries revealed that Bi and Sn are the most frequently reported metal components⁴²⁻⁴⁵. It is also worth noting that recent studies indicate that La incorporation can markedly enhance catalytic performance.⁴⁶ Consequently, we targeted La-Sn-Bi as our medium-entropy alloy (MEA) system. Furthermore, guided by elements commonly incorporated into high-entropy alloys⁴⁷⁻⁴⁹, we expanded this composition with Co and Ni to design a La-Sn-Bi-Co-Ni high-entropy alloy (HEA) system. Since most structural entries in the DigCat platform correspond to face-centered cubic (FCC) parent lattices, we adopted FCC lattices to construct our nanoparticle prototypes. Specifically, for both the ternary and quinary systems, we first constructed sufficiently large equimolar FCC supercells. By selecting the geometric center of the structure and applying varying cutoff radii, we generated spherical nanoparticle parent phases with five distinct sizes ranging from 2 to 8 nm, as illustrated in Fig. S2.

These generated nanoparticle models served as specific benchmarks to evaluate the performance of the updated generation module within the ApolloX engine. Focusing on the complex quinary system, the capability to generate structures based on target SRO is detailed in Fig. 2d. The algorithm demonstrates exceptional scalability, producing ~5000 structures for a ~5000-atom system in under 50 seconds, while maintaining near-constant generation time per structure even as system size increases (Fig. 2d, left). Furthermore, it shows high fidelity: a statistical analysis of 2000 generated structures reveals a tight distribution of SRO root mean square errors (RMSE) centered around an average of 0.035 (Fig. 2d, middle). A representative comparison for a structure with this average error confirms the close match between the target and generated α values (Fig. 2d, right). Additionally, we conducted performance tests on ternary catalysts across varying sizes. As shown in Fig. S3, the module exhibits similarly outstanding performance, demonstrating the broad applicability and robustness of this generation function.

Although the MC-based generation method and MLPs are capable of handling significantly larger structures, we prioritized computational feasibility for the subsequent high-precision DFT calculations of adsorption energies. Therefore, we selected the 2 nm-sized models, specifically La₂₃Sn₂₃Bi₂₂ and La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃ as the final catalyst compositions. The specific atomic ratios in these models represent the natural outcome of cutting an equimolar FCC supercell into a nanoparticle shape. Subsequently, we applied the complete ApolloX workflow for structural generation. For each composition, we performed 15 rounds of PSO iterative cycles, with each cycle producing 100 structures. We utilized Crystal Hamiltonian Graph Neural Network (CHGNet)⁶³ as the surrogate MLP, and its performance benchmark (Fig. S4) confirms that it meets the required accuracy for this study. The relationship between the minimum energy of all structures and the iteration count for ternary and quinary catalysts is presented in Fig. S5 and Fig. S6, respectively. The results indicate that the total minimum energy gradually decreases and stabilizes as iterations progress, signifying a successful convergence toward the global minimum. To further validate this optimization capability, we compared the energy distribution of the 1,500 structures generated by ApolloX against a baseline of randomly shuffled structures (Fig. S7). The results demonstrate that ApolloX identifies a significantly higher proportion of structures within the low-energy regime and achieves lower absolute minimum energies than the random method, further underscoring the effectiveness of our approach.

**Fig. 2. The updated ApolloX workflow for structure generation and its performance.**  
*(a)* Main workflow of ApolloX for structural modeling, including the model’s iterative generation process. *(b)* Details of the hybrid Generating algorithm. This algorithm consists of two stages: First, two successive steps of a greedy local search (defined in c) are applied to quickly generate a set of low-error lattices; second, these structures are further refined using a Metropolis-based search to precisely match the target SRO parameters (α). *(c)* The greedy local search routine used in the first stage of (b). It performs random atom swaps and only accepts changes that reduce the error, enabling rapid initial optimization. *(d)* Performance benchmarking of quinary catalysts with varying sizes and quantities. (Left) Computational time as a function of the number of structures (for different numbers of atoms) on an NVIDIA 4090 GPU. (Middle) Statistical distribution of the SRO root mean square error (RMSE) for 2000 generated structures, showing an average RMSE of 0.035. (Right) A representative comparison between target α values (Target) and the generated values (Generated) for a structure with the average RMSE of 0.035.

> Placeholder:  
> `![Fig. 2. The updated ApolloX workflow for structure generation and its performance.](fig2.png)`

Leveraging the extensive structural ensemble generated by the ApolloX, we next systematically investigated how SRO within the nanoparticles governs their thermodynamic stability. To disentangle surface and interior effects, each nanoparticle was partitioned into a surface shell and a bulk shell (Fig. 3a–b). The partition was defined using two concentric spheres centered at the geometric center of all catalyst atoms; the outer radius equals the maximum atomic distance from the center, and the inner radius is half of this maximum.

We characterized SRO using two descriptors: (i) the W-C parameter αᵢⱼ, which captures aggregation vs dispersion between atomic pairs, and (ii) the LDD parameter δⱼ, which measures local enrichment or depletion of element j. Based on their mathematical definitions, αᵢⱼ can be directly interpreted as describing the distribution of j atoms around i atoms, while δⱼ represents the distribution of the j element itself. Therefore, LDD parameter allows for a quick and convenient understanding of the individual contributions of different elements, while W-C parameter provides more detailed insights into the short-range ordering and structural characteristics. Notably, the two descriptors exhibit both consistency and complementarity. Their progressive relationship reflects a shared description of local chemical environments, while the enrichment-depletion trends captured by δⱼ and the aggregation-dispersion trends represented by αᵢⱼ highlight distinct structural features. In this sense, the two descriptors offer complementary perspectives on short-range ordering.

**Fig. 3. Statistical contribution of short-range order to thermodynamic stability and schematic structures.**  
*(a–b)* Schematics of La₂₃Sn₂₃Bi₂₂ (a) and La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃ (b) indicating inner-outer shells and the surface–bulk partition. *(c–d)* Surface and bulk Warren-Cowley parameter α score maps for La₂₃Sn₂₃Bi₂₂ (c) and La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃ (d). The horizontal bar below reports surface Local Density Deviation parameter δ contributions. The color scale is S, where S = ΔE⁺ - ΔE⁻.

For each composition, we evaluated SRO separately in the bulk and surface shells across 1,500 structures. For each descriptor (α or δ), structures were partitioned by the sign of the descriptor. We then computed the energy contrast ΔE as a score to quantify the stabilizing influence of that SRO motif, where ΔE = E⁺ - E⁻. The results for all descriptors were combined and visualized as score maps (Fig. 3c–d). Negative values (blue) indicate that α > 0 (or δ > 0) is associated with lower energy, i.e., aggregation of the pair or local enrichment of the element stabilizes the structure; positive values (red) indicate the opposite tendency, favoring dispersion of the corresponding pair or depletion of the element.

Applying this framework yields the following trends. For the La₂₃Sn₂₃Bi₂₂, at the surface, local composition tends to be enriched in Bi and depleted in Sn and La. Pairwise ordering in both surface and bulk tends to favor aggregation of Bi–La and Sn–La pairs, while disfavoring Sn–Bi and Bi–Bi pairs through dispersion. For the La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃, at the surface, local composition tends to be enriched in Bi and La and depleted in Ni and Co. Surface ordering favors aggregation of Bi–La and disfavors Bi–Bi, Ni–Bi, and Sn–Bi. In the bulk area, ordering tends to favor aggregation among Co- and Ni-containing pairs (Co–Co, Ni–Co, Sn–Co, La–Ni, Sn–Ni), while disfavoring Co–Bi and Sn–Bi. Overall, Bi enrichment in the surface region and Co/Ni enrichment in the bulk both enhance the structural stability of the nanoparticle catalyst. La₂₃Sn₂₃Bi₂₂ exhibits a strong consistency between bulk and surface in the aggregation/dispersion tendencies of atomic pairs, whereas La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃ shows pronounced differences between bulk and surface due to the bulk aggregation propensity of Co and Ni, highlighting the significant impact of multicomponent high-entropy effects on structure. Specifically, driven by the differences in surface energies and cohesive energies among constituent elements, high-entropy nanoparticles often exhibit non-uniform elemental distribution (e.g., surface segregation), resulting in distinct SRO between the surface and the bulk.⁹˒⁵⁰ These observations both substantiate the rationale for selecting Co and Ni as components of the HEA catalyst and reveal that the role of Bi extends beyond its commonly noted influence on activity⁵¹⁻⁵⁵ to include a substantial contribution to stability.

---

### Automated Adsorption Sites Search and Absorbate Binding Energy Calculations

After structure generation, ARCADE systematically integrates adsorption-site identification and energy evaluation for the most stable structure. We employ an autonomous, geometry-aware search protocol based on Voronoi tessellation⁴⁰˒⁴¹ to distinguish surface adsorption sites from interior atoms (Fig. 4). The Voronoi cell of each atom is constructed from the convex hull of its Voronoi vertices; interior atoms, being fully surrounded, have finite and relatively small cell volumes, whereas surface atoms lack neighbors in certain directions and thus exhibit much larger (sometimes unbounded) volumes (see Fig. S8 for the spatial partitioning illustration). We therefore apply a volume threshold Vₜₕ: atoms with Voronoi cell volumes exceeding Vₜₕ are classified as surface sites, while the remainder are fixed as interior atoms to mimic bulk constraints. This geometry-aware rule is robust to disorder and does not require predefining facets.

For each identified surface site, the workflow derives configuration parameters through two general geometric criteria. First, for single-site adsorption, the adsorption direction is defined by the radial vector extending from the nanoparticle center through the target surface atom, with the initial bond length estimated based on the sum of the covalent radii of the adsorption site atom and the binding atom of the intermediate. Second, to establish setups for multi-site adsorption, the workflow identifies neighboring surface atom pairs and constructs the initial geometry by aligning the adsorbate's binding atoms along the radial vectors of these adjacent sites—analogous to the single-site protocol. The adsorption height is then adjusted along these directions to reconcile the intermediate's internal bond distances with the spatial separation of the surface sites, ensuring a physically plausible starting configuration prior to relaxation. For the CO₂RR-to-formate pathway studied here, these protocols were applied to generate configurations for three key intermediates: \*COOH and \*H are treated as monodentate species using the single-site criteria (binding via the C atom and H atom, respectively), while \*OCHO is modeled as a bidentate species utilizing the bridge setup logic (binding via its two O atoms).

**Fig. 4. Schematic workflow for autonomous adsorption site identification.**  
Voronoi tessellation distinguishes surface sites from interior atoms based on a volume threshold. For identified surface sites, the workflow determines single-site adsorption geometries (direction and bond length) and identifies neighboring atom pairs for bridge setups, integrating these parameters to generate the final atomic configurations. Blue and green modules represent information derived from the surface and the integration of interior atoms, respectively. Note: "Nano." in the flowchart is an abbreviation for nanoparticle.

Based on the ARCADE screening results, the most stable candidates, La₂₃Sn₂₃Bi₂₂ and La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃, were identified for further electrocatalytic analysis. Following the procedure illustrated in Fig. 4, the outermost atoms were determined and designated as adsorption sites, from which the adsorption structures of \*COOH, \*OCHO, and \*H were constructed. For La₂₃Sn₂₃Bi₂₂ nanoparticles, the outermost layer contains 53 atoms, yielding 53 possible adsorption structures each for \*COOH and \*H. For \*OCHO, which requires diatomic adsorption sites, 145 possible configurations were obtained. In the case of LaSnBiCoNi nanoparticles, the outermost layer comprises 41 atoms, corresponding to 41 adsorption structures for \*COOH and \*H, and 112 for \*OCHO. Adsorption energies for all configurations were subsequently calculated. Specifically, the calculated adsorption energies for \*COOH, \*H, and \*OCHO on La₂₃Sn₂₃Bi₂₂ are detailed in Tables S2–S4, respectively, while the corresponding data for La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃ are provided in Tables S5–S7. The resulting distributions of these adsorption energies are presented in Fig. S9. Detailed computational methods are provided in the Methods Section. Representative examples of adsorption configurations before and after optimization for different compositions and intermediates are illustrated in Fig. S10. As observed in these examples, although adsorbates were initially placed on predefined sites, they were allowed to relax freely during DFT geometry optimizations and, where energetically favorable, migrated to lower-energy positions. For La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃, the optimized sites are: \*COOH on Sn, Bi, Ni, Co, or La; \*OCHO on diatomic sites composed of any two among Sn, Bi, Ni, Co, and La; and \*H on Sn, Bi, Ni, Co, La, or bridge sites. For LaSnBi, the optimized sites are: \*COOH on Sn, Bi, La, or bridge; \*OCHO on diatomic sites composed of any two among Sn, Bi, and La; and \*H on Sn, Bi, La, or bridge sites.

---

### Selectivity Analysis and the Relationship with SRO

Starting with initial compositional guesses, followed by generating stable structures and performing adsorption site searches, we have established the structural modeling foundation for catalyst design. The structure-energy data obtained during structure generation have been used for catalyst stability analysis, while the adsorption site and adsorption free energy data for different intermediates from the adsorption site search will be utilized to evaluate the selectivity among competing pathways and the activity of the preferred reaction pathway.

To assess the selectivity of CO₂RR, we compared the Gibbs free energies of \*OCHO, \*COOH, and \*H adsorbed on the same site³⁰˒⁵⁶˒⁵⁷. However, as noted above, the adsorption sites for these intermediates are not strictly identical: \*COOH and \*H typically adsorb on a single site, whereas \*OCHO involves two adjacent atoms, and adsorbates may migrate during optimization. To address this inconsistency, we proposed using SRO to identify adsorption environments with high similarity, thereby enabling a fair comparison of their energies. We matched structures using: (i) the nearest adsorption-site atom and (ii) SRO descriptors characterizing the local environment around the adsorbate. Specifically, the “nearest adsorption-site atom” is defined as the catalyst atom participating in the pair with the shortest relative distance to any adsorbate atom, where the relative distance is the Euclidean separation divided by the sum of the two atoms’ covalent radii. To delineate a chemically meaningful local environment, we adopt an “effective interaction range” per adsorbate atom equal to the average of the sum of the van der Waals radii of the adsorbate and catalyst atoms, approximately 4.5 Å. The union of catalyst atoms lying within this range for all atoms of the adsorbate defines the adsorbate’s local chemical environment. Within this local environment, we computed: (i) the W-C parameter αᵢⱼ and (ii) the LDD parameter δⱼ to quantify SRO. For each relaxed adsorption structure, the set of αᵢⱼ and δⱼ values within the local environment are concatenated into a single feature vector **s** used as a label. Structures are then matched if they share the same nearest adsorption-site atom and if the cosine similarity between their local-environment SRO vectors is ≥ 0.9. We chose the 0.9 threshold to accommodate minor differences in adsorbate orientation and in the range of the local environment, as well as small catalyst distortions induced during relaxation; this tolerance balances robustness with specificity. A schematic of the matching protocol for medium/high-entropy catalysts is shown in Fig. 5a.

**Fig. 5. Catalytic selectivity analysis and the role of short-range order.**  
*(a)* Matching of intermediates configurations by local chemical environment near adsorbates, which is symbolized by the nearest adsorption site and the local SRO vector **s**. *(b–e)* Selectivity analysis of \*OCHO vs \*COOH and \*OCHO vs \*H for La₂₃Sn₂₃Bi₂₂ (b–c) and La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃ (d–e). *(f–i)* Score heat maps showing the contribution of short-range order on selectivity, quantified as S\_sel, where S\_sel = ⟨ΔE⟩₊ - ⟨ΔE⟩₋. *(f)* and *(h)* report the Local Density Deviation parameter δ contributions for La₂₃Sn₂₃Bi₂₂ and La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃, respectively; *(g)* and *(i)* report the Warren-Cowley parameter α contributions for La₂₃Sn₂₃Bi₂₂ and La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃, respectively.

For elemental catalysts, structures of a given adsorbate are matched using a geometric criterion: we computed the adsorbate position vector as the displacement from the catalyst geometric center to the adsorbate geometric center, and we required a cosine similarity ≥ 0.9 between these position vectors to establish a match.

For La₂₃Sn₂₃Bi₂₂, the matching results and SRO cosine similarity data for \*OCHO vs. \*H and \*OCHO vs. \*COOH are provided in Tables S8 and S9, respectively. For La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃, the corresponding data for \*OCHO vs. \*H and \*OCHO vs. \*COOH are detailed in Tables S10 and S11, respectively. After matching structures corresponding to each adsorbate, catalytic selectivity was analyzed as shown in Fig. 5b–e and Fig. S11. We quantify selectivity by comparing the binding energies of \*OCHO to those of the other two adsorbates. The difference in binding energy between adsorbates, ΔΔE, indicates both the direction and magnitude of selectivity: a positive and larger ΔΔE implies that the reaction favors the \*OCHO pathway, i.e., higher selectivity. Fig. 5b–c reports the \*OCHO vs. \*COOH and \*OCHO vs. \*H selectivity for La₂₃Sn₂₃Bi₂₂, and Fig. 5d–e provides the corresponding results for La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃, with the nearest adsorption-site atom annotated for each structure. All the results indicate the high CO₂RR-to-HCOOH selectivity for the nanoparticles. Fig. S11 shows the \*OCHO vs. \*COOH and \*OCHO vs. \*H selectivity for elemental catalysts. Comparing MEA/HEA catalysts with elemental catalysts, the MEA/HEA data points cluster preferentially toward the upper-left of the diagonal, indicating superior selectivity relative to the elemental references.

To explicitly assess how SRO affects catalytic selectivity, we group structures by the sign of each descriptor (α and δ), and then compare the difference between the groups’ averages of ΔΔE. Because the extremes of ΔΔE (its maximum and minimum) are most consequential—setting the upper bound of selectivity and indicating whether any non-OCHO-favorable cases occur—we assign linearly varying weights to data points:

> (weighting scheme described in Supporting Information)

where ΔΔE\_max and ΔΔE\_min are the maximum and minimum of ΔΔE, respectively; ΔΔE\_center represents the center of the ΔΔE range; and η controls the weighting strength, set to 0.8 in this work. The difference between the groups’ weighted means is then:

> (score definition as described in the Supporting Information)

where ⟨·⟩ denotes the weighted average. For each SRO descriptor, we compute the weighted-mean difference in ΔΔE. The results for all descriptors were combined and visualized as score heat maps (Fig. 5f–i). Fig. 5f–g shows the impacts of the LDD parameters δⱼ and the W-C parameters αᵢⱼ on La₂₃Sn₂₃Bi₂₂ selectivity; Fig. 5h–i shows the corresponding impacts for La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃. Positive values (dark blue) indicate that the δ > 0 group has larger ΔΔE, i.e., stronger preference for the \*OCHO pathway; negative values (light yellow) indicate that the δ < 0 group is more selective. Larger absolute values are shown by colors closer to the ends of the scale, indicating a stronger influence of the SRO descriptor on selectivity.

Synthesizing the site-specific comparisons (Fig. 5b–e) with the SRO–selectivity analysis (Fig. 5f–i) reveals some key trends governing the catalytic selectivity. For the La₂₃Sn₂₃Bi₂₂ catalyst, enhanced selectivity is associated with adsorption sites where La is the nearest atomic neighbor, a local environment enriched in La and Sn but depleted in Bi, and the aggregation of Sn–Sn and Sn–La pairs alongside the dispersion of Bi–Sn and Bi–Bi pairs. A similar pattern emerges in the La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃ system, where proximity to a La atom enhances selectivity while proximity to a Co atom suppresses it. In this system, a local environment enriched in La and Ni and depleted in Co is favorable, as is the aggregation of Co–La, Ni–La, Co–Ni, and Ni–Ni pairs coupled with the dispersion of Co–Bi and Co–Co pairs. Critically, a comparison of the two systems highlights a unifying principle: the presence of La, both as the nearest neighbor to the adsorption site and as an enriching element in the local chemical environment, is a key determinant of selectivity. In other words, beyond its experimentally reported role in La-doped SnO₂, where La enhances reaction selectivity via pinning effects and water activation⁴⁶, we also observe a pronounced influence of La on selectivity in our MEA/HEA catalysts. These findings further underscore the substantial application potential of La elements in CO₂RR catalysis.

---

### Activity Analysis and the Relationship with SRO

Structures with different adsorbates but similar local chemical environments are matched through the cosine similarity of the short-range order vector **s**, and the \*OCHO favorable sites are identified by comparing ΔG\_OCHO with ΔG\_H and ΔG\_COOH. Subsequently, the TOFs are predicted using our recently developed pH-dependent microkinetic modeling for CO₂RR under the RHE scale,³⁰˒⁷¹ by comprehensively considering the electric field effects,⁷² potential of zero-charges under explicit solvent conditions,⁷³ and solvation corrections via the H-bonding. Notably, this pH-dependent modeling method has excellent agreement with experimental results, in particular for those with CO₂RR-to-formate³⁰˒⁷¹ and CO₂RR-to-CO.⁷² Details of the microkinetic modeling methods, the considered elementary steps, and all relevant parameters can be found in the Supplementary Information. These pH-dependent microkinetic volcano models are shown in Fig. 6a–c for La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃, La₂₃Sn₂₃Bi₂₂, and the corresponding elemental catalysts, respectively. The activity of M/HEA catalysts is clearly higher than that of elemental catalysts, as several points are located close to the peak of the volcano models for both La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃ and La₂₃Sn₂₃Bi₂₂, while no such points appear for the single elemental catalysts. However, the number of constituent elements does not monotonically determine catalytic activity. Both La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃ and La₂₃Sn₂₃Bi₂₂ exhibit adsorption sites near the volcano peak, suggesting that they each possess high-activity sites for CO₂RR-to-formate.

**Fig. 6. Catalytic activity analysis and the role of short-range order.**  
*(a–c)* pH-dependent microkinetic models for La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃ (a), La₂₃Sn₂₃Bi₂₂ (b) and elementary substances catalysts (c) at an applied potential of -1.0 V vs. RHE. *(d–g)* Heat maps showing the contribution of short-range order on activity, quantified as S\_act, where S\_act = ⟨A⟩₊ - ⟨A⟩₋. *(d)* and *(e)* report the Local Density Deviation parameter δ contributions for La₂₃Sn₂₃Bi₂₂ and La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃, respectively; *(f)* and *(g)* report the Warren-Cowley parameter α contributions for La₂₃Sn₂₃Bi₂₂ and La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃, respectively. *(i)* Structures with the optimal active adsorption sites for La₂₃Sn₂₃Bi₂₂ and La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃ under different pH conditions.

To further investigate the influence of SRO on catalytic performance, we define an activity descriptor as A = -log₁₀(TOF), where TOF is measured in s⁻¹. A larger value of this descriptor indicates poorer activity of the adsorption site. Structures with different adsorption sites of La₂₃Sn₂₃Bi₂₂ or La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃ under different pH conditions are divided into two groups based on each SRO descriptor α or δ. The average difference in activity descriptors, ΔA, is calculated between the two groups as the score of the descriptor. The results for all descriptors were combined and reported as the score heatmaps in Fig. 6d–g, where d–e denote the effect of δ for La₂₃Sn₂₃Bi₂₂ and La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃, respectively, and f–g denote the effect of α for La₂₃Sn₂₃Bi₂₂ and La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃, respectively. Negative values (dark blue) indicate that the α > 0 (or δ > 0) group exhibits higher catalytic activity, while positive values (light yellow) indicate stronger activity for the α < 0 (or δ < 0) group. Larger absolute values are shown by colors closer to the ends of the scale, indicating a stronger influence of the SRO descriptor on activity. Fig. 6h depicts the adsorption structures corresponding to the highest activity for the two catalyst compositions across three pH levels. As La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃ exhibits the same optimal active site at pH = 1 and 7, five representative adsorption configurations are presented in total.

The SRO–activity analysis presented in Fig. 6d–g reveals that the structure-performance relationships for the MEA/HEA catalysts are approximately consistent across different pH conditions. For the La₂₃Sn₂₃Bi₂₂ system, enhanced catalytic activity correlates with a local adsorbate environment enriched in Bi and Sn and depleted in La. This is structurally accompanied by the aggregation of Sn–Bi and Bi–Bi pairs and the dispersion of La–La and La–Sn pairs. Similarly, for the La₁₄Sn₁₄Bi₁₄Co₁₃Ni₁₃ catalyst, higher activity is linked to a local environment enriched in Sn but depleted in La and Bi. The corresponding favorable structural motifs include the aggregation of Co–Ni, Sn–Sn, and Bi–Ni pairs, along with the dispersion of Co–La, Sn–Ni, Sn–Bi, and Ni–La pairs. A comparison between the two systems indicates that the SRO–activity relationships share key similarities, as catalytic performance in both is enhanced by a local environment that is enriched in Sn and depleted in La. The activity-enhancing role of Sn in both catalyst systems is readily understood, as Sn-mediated promotion of CO₂RR-to-formate has been extensively documented⁵⁸˒⁵⁹. However, Bi, although a common catalyst component, exhibits a pronounced promotional effect only in the MEA catalyst. In the HEA catalyst, its impact on activity is markedly weaker than that of Sn and can even be slightly detrimental. This phenomenon reflects the unique complexity of high-entropy alloys (e.g., the cocktail effect and complex site distributions), where the catalytic behavior of individual elements is strongly modulated by their distinct local environments and interaction with other elements, thus warranting careful discrimination.⁶˒⁶⁰ In addition, La shows a modest suppressive effect on activity, which may partly explain the historically limited exploration of La in CO₂RR-to-formate catalysts; nevertheless, its strong enhancement of selectivity leads to an overall positive contribution to catalytic performance.

---

## Discussion

In this work, we have introduced ARCADE, a hierarchical and fully integrated computational framework designed to bridge the gap between the realistic structural complexity of high-entropy alloys and the rational design of their catalytic performance. This paradigm is built upon four synergistic pillars: (i) the updated SRO-constrained ApolloX engine for generating thermodynamically stable nanoparticles; (ii) an autonomous, geometry-aware search algorithm for identifying comprehensive adsorption ensembles; (iii) state-of-the-art pH-dependent microkinetic modeling for the precise evaluation of theoretical activity under realistic operating conditions; and (iv) crucially, a quantitative interpretability methodology utilizing SRO–performance score heatmaps. This final component successfully decouples the intricate relationships between local SRO environments and multiple performance metrics within complex catalytic ensembles.

Applying this paradigm to the design of MEA and HEA nanoparticle catalysts for CO₂RR-to-formate, we utilized these SRO-guided score heatmaps to quantitatively decode the distinct roles of constituent elements and atomic pairs. Our analysis reveals general design principles where different performance metrics are governed by distinct structural features. Catalyst stability is primarily governed by specific segregation behaviors driven by high-entropy effects: the enrichment of Bi at the surface and Co/Ni in the bulk. Regarding this, trends are consistent across both systems, where the aggregation of Bi–La pairs at the surface and Co–Co/Co–Ni pairs in the core align perfectly with the observed elemental enrichment. Thus, while LDD analysis provides a clear understanding of individual elemental contributions, W-C parameters further offer granular structural insights, such as Bi's specific tendency to aggregate with La rather than itself.

Distinct from the consistent trends observed for stability, the influence of SRO on catalytic function (selectivity and activity) exhibits specific structural dependencies between the two catalysts. Based on LDD analysis, La and Sn are identified as the primary drivers for selectivity and activity, respectively. Specifically, high selectivity correlates with La-enriched local environments (favored by Sn–La aggregation in MEA and Co–La/Ni–Ni pairs in HEA), while high activity is promoted by Sn-enriched but La-depleted environments (supported by Sn–Bi aggregation in MEA and Co–Ni/Sn–Sn aggregation in HEA). These results demonstrate that while Sn and La act as the primary functional drivers for CO₂RR-to-formate, their performance is critically modulated by specific pairwise interactions defined by the surrounding matrix. Collectively, these findings underscore the necessity of our dual-descriptor approach: without LDD parameters, the individual contributions of each element would remain unclear; without W-C parameters, the detailed characteristics of pairwise interactions governing specific reaction pathways could not be revealed. This comprehensive interpretability enables the optimization of catalysts by systematically tuning both local elemental enrichment and atomic pairwise patterns.

Finally, this computational workflow is fully integrated with the DigCat platform, establishing a closed-loop ecosystem for catalyst design. By leveraging DigCat for initial composition screening and parent structure determination, followed by our automated generation-and-analysis cycle, we successfully derived quantitative relationships between structural features and multiple performance metrics. Crucially, the newly generated structures and insights are deposited back into the DigCat database. Consequently, subsequent optimization becomes more targeted: starting from known structural prototypes, one can generate models using SRO motifs associated with superior performance and proceed directly to analysis. Iterating this loop can markedly accelerate interpretable catalyst design.

Beyond the specific success demonstrated here, the ARCADE framework represents a significant leap in digital materials discovery. By treating structural complexity as a programmable design parameter rather than a hindrance, this protocol offers a generalized methodology applicable to diverse electrochemical processes operating under varying pH environments. This approach effectively closes the loop between atomic-level disorder and multiple catalytic performance metrics, providing a clear decoding of how local chemical environments dictate overall system behavior. The ARCADE framework thus provides a powerful and transferable paradigm for the discovery of complex catalyst materials for relatively complicated catalytic reactions, moving beyond black-box screening toward the rational, SRO-informed design of next-generation catalysts. By integrating these capabilities into the DigCat platform, we empower the research community to navigate the vast compositional space of high-entropy materials with unprecedented efficiency and interpretability.

---

## Methods

### Warren-Cowley parameter

A widely-used descriptor to capture the chemical short-range order of structure is the Warren-Cowley (W-C) parameter, which is defined as

> (formal definition given in Supporting Information)

where i and j denote atomic species, Pᵢⱼ is the conditional probability of encountering an atom of type j at a certain distance cutoff (or a given neighbor shell) from a central atom of type i, and cⱼ is the overall concentration of element j in the system.

When αᵢⱼ = 0, the distribution of j around i follows the global composition, i.e., random mixing. Values αᵢⱼ < 0 indicate that j occurs less frequently near i than expected from randomness, i.e., j is depleted around i, i–j atom pairs tend to disperse; while αᵢⱼ > 0 signifies an overrepresentation of j near i, i.e., j is enriched around i, i–j atom pairs tend to aggregate. The deviation of αᵢⱼ from zero quantitatively characterizes the aggregation or dispersion propensity of an i–j pair at separation r.

### Local Density Deviation parameter

The adsorption process on catalysts is strongly related to the local chemical environments near the adsorbates, especially the compositions of the local environments. To more directly and transparently capture short-range order in local elemental concentrations, we introduce, by analogy to the Warren–Cowley parameter, a new descriptor termed the Local Density Deviation (LDD) parameter. It is defined as

> (formal definition given in Supporting Information)

where j denotes the atomic species, cⱼ^local is the local concentration of element j within the chosen local region, and cⱼ is the overall concentration of element j in the system.

When δⱼ = 0, the local density of element j is equal to its overall concentration. Positive values (δⱼ > 0) indicate that element j is locally depleted relative to its global concentration, whereas negative values (δⱼ < 0) signify local enrichment of j.

Although the LDD parameter shows certain similarity to the W-C parameter in form, their meanings are obviously non-equivalent. The LDD parameters focus on the enrichment or depletion of elements, while W-C parameters show more about aggregation or dispersion relationship between atom pairs. So LDD offers practical advantages—direct interpretability, low dimensionality, and robustness, making it well suited for SRO analysis for M/HEA catalysts.

### Density functional theory (DFT) calculation method

The Vienna Ab initio Simulation Package (VASP)⁶² was employed to perform density functional theory (DFT) calculations. The interaction between ions and electrons was described using the projector augmented-wave (PAW) method. Exchange-correlation effects were treated with the revised Perdew-Burke-Ernzerhof (RPBE) functional within the generalized gradient approximation (GGA)⁶⁴˒⁶⁵, a choice informed by prior benchmarking comparisons between experimental data and theoretical predictions.⁶⁶˒⁶⁷ The DFT-D3 method is employed to correct the van der Waals interactions of the system. A plane-wave kinetic energy cutoff of 500 eV was applied for the expansion of the Kohn–Sham wave functions. Electronic self-consistency was converged to 10⁻⁵ eV, and structural relaxations were performed until the forces were below 0.05 eV Å⁻¹. For Brillouin zone integration, a Γ-centered Monkhorst-Pack grid was employed, ensuring that the product of k-points and basis vector length (k × a) exceeded 20 Å. Additionally, a vacuum spacing of 15 Å was applied perpendicular to the surface to prevent interactions between periodic images.

---

## Acknowledgments

The calculations were partly performed using the San Diego Supercomputer Center (SDSC) Expanse at UC San Diego through allocation MAT240028 from the Advanced Cyberinfrastructure Coordination Ecosystem: Services & Support (ACCESS) program. Additionally, H. Li acknowledges the Center for Computational Materials Science, Institute for Materials Research, Tohoku University for the use of MASAMUNE-IMR (No. 202412-SCKXX0211) and the Institute for Solid State Physics (ISSP) at the University of Tokyo for the use of their supercomputers. W.-L. Li gratefully acknowledges support from the ACS Petroleum Research Fund under Doctoral New Investigator Grant 69037-DNI10. W.-L. Li and Y. Ding acknowledge the Jacobs School of Engineering Early Career Faculty Development Award, University of California San Diego. W.-L. Li also acknowledges support of the Hellman Fellowship. 

---

## Data availability

Adsorption energy data, as well as data used for selectivity and activity analysis, are available in the paper. All structural data have been uploaded to DigCat (www.digcat.org). Additional data supporting this project are available from the corresponding author upon reasonable request.

---

## Code availability

The entire workflow can be implemented through DigCat (www.digcat.org). The code for ApolloX is also available on GitHub (https://github.com/FNC001/ApolloX).

---

## References

1. Seh, Zhi Wei, et al. Combining theory and experiment in electrocatalysis: Insights into materials design. *Science* **355**.6321 (2017): eaad4998.  

2. Hammer, Bjørk, and Jens K. Nørskov. Why gold is the noblest of all the metals. *Nature* **376**.6537 (1995): 238–240.  

3. Barmparis, Georgios D., et al. Nanoparticle shapes by using Wulff constructions and first-principles calculations. *Beilstein Journal of Nanotechnology* **6**.1 (2015): 361–368.  

4. Ringe, E., Richard P. Van Duyne, and L. D. Marks. Wulff construction for alloy nanoparticles. *Nano Letters* **11**.8 (2011): 3399–3403.  

5. Alsunni, Yousef A., and Charles B. Musgrave. Effect of applied potential on metal surfaces: Surface energy, Wulff shape and charge distribution. *Applied Surface Science* **610** (2023): 155147.  

6. Xin, Yue, et al. High-entropy alloys as a platform for catalysis: progress, challenges, and opportunities. *ACS Catalysis* **10**.19 (2020): 11280–11306.  

7. Lyakhov, Andriy O., et al. New developments in evolutionary structure prediction algorithm USPEX. *Computer Physics Communications* **184**.4 (2013): 1172–1182.  

8. Wang, Hui, et al. CALYPSO structure prediction method and its wide application. *Computational Materials Science* **112** (2016): 406–415.  

9. Moniri, Saman, et al. Three-dimensional atomic structure and local chemical order of medium-and high-entropy nanoalloys. *Nature* **624**.7992 (2023): 564–569.  

10. Li, Honglin, et al. Conditional Generative Modeling for Amorphous Multi-Element Materials. *arXiv preprint* arXiv:2503.07043 (2025).  

11. Zeni, Claudio, et al. Data-driven simulation and characterisation of gold nanoparticle melting. *Nature Communications* **12**.1 (2021): 6056.  

12. Ju, Shin-Pon, et al. Tailoring mechanical performance in bulk nanoparticle-structured ZnO and Al₂O₃: Insights from deep learning potential molecular dynamics simulations. *Materials Today Communications* **42** (2025): 111161.  

13. Ong, Shyue Ping, et al. Python Materials Genomics (pymatgen): A robust, open-source python library for materials analysis. *Computational Materials Science* **68** (2013): 314–319.  

14. Larsen, Ask Hjorth, et al. The atomic simulation environment—a Python library for working with atoms. *Journal of Physics: Condensed Matter* **29**.27 (2017): 273002.  

15. Boes, Jacob R., et al. Graph theory approach to high-throughput surface adsorption structure generation. *The Journal of Physical Chemistry A* **123**.11 (2019): 2281–2285.  

16. Bisbo, Malthe K., and Bjørk Hammer. Efficient global structure optimization with a machine-learned surrogate model. *Physical Review Letters* **124**.8 (2020): 086102.  

17. Ma, Huan, et al. Machine learning predicts atomistic structures of multielement solid surfaces for heterogeneous catalysts in variable environments. *The Innovation* **5**.2 (2024).  

18. Shang, Cheng, and Zhi-Pan Liu. Stochastic surface walking method for structure prediction and pathway searching. *Journal of Chemical Theory and Computation* **9**.3 (2013): 1838–1845.  

19. Jia, Xue, et al. Advancing electrocatalyst discovery through the lens of data science: State of the art and perspectives. *Journal of Catalysis* (2025): 116162.  

20. Chen, Yuanzheng, et al. Data-Driven Strategies for Designing Multicomponent Molten Catalysts to Accelerate the Industrialization of Methane Pyrolysis. *ACS Catalysis* **15** (2025): 11003–11012.  

21. Zhao, Zhi-Jian, et al. Theory-guided design of catalytic materials using scaling relationships and reactivity descriptors. *Nature Reviews Materials* **4**.12 (2019): 792–804.  

22. Odetola, Peter Ifeolu, et al. Exploring high entropy alloys: a review on thermodynamic design and computational modeling strategies for advanced materials applications. *Heliyon* **10**.22 (2024).  

23. Tran, Kevin, and Zachary W. Ulissi. Active learning across intermetallics to guide discovery of electrocatalysts for CO₂ reduction and H₂ evolution. *Nature Catalysis* **1**.9 (2018): 696–703.  

24. Wang, Qi, and Yonggang Yao. Harnessing machine learning for high-entropy alloy catalysis: a focus on adsorption energy prediction. *npj Computational Materials* **11**.1 (2025): 91.  

25. Nørskov, Jens Kehlet, et al. Origin of the overpotential for oxygen reduction at a fuel-cell cathode. *The Journal of Physical Chemistry B* **108**.46 (2004): 17886–17892.  

26. Nørskov, Jens K., et al. Density functional theory in surface chemistry and catalysis. *Proceedings of the National Academy of Sciences* **108**.3 (2011): 937–943.  

27. Liu, Heng, et al. Reversible hydrogen electrode (RHE) scale dependent surface pourbaix diagram at different pH. *Langmuir* **40**.14 (2024): 7632–7638.  

28. Zhang, Di, et al. Unraveling the pH-dependent oxygen reduction performance on single-atom catalysts: from single-to dual-sabatier optima. *Journal of the American Chemical Society* **146**.5 (2024): 3210–3219.  

29. Ye, Songbo, et al. Decoding pH-dependent electrocatalysis through electric field models and microkinetic volcanoes. *Journal of Materials Chemistry A* (2025).  

30. Wang, Yuhang, et al. Divergent Activity Shifts of Tin‐Based Catalysts for Electrochemical CO₂ Reduction: pH‐Dependent Behavior of Single‐Atom Versus Polyatomic Structures. *Angewandte Chemie International Edition* **64**.8 (2025): e202418228.  

31. Chen, Guanzhen, et al. A long‐range disordered RuO₂ catalyst for highly efficient acidic oxygen evolution electrocatalysis. *Angewandte Chemie International Edition* **63**.50 (2024): e202411603.  

32. Yang, Yiyuan, et al. Chemical short-range order in multi-principal element alloy with ordering effects on water electrolysis performance. *Materials Today* **72** (2024): 97–108.  

33. Wang, Xuefeng, et al. RuO₂ with short‐range ordered tantalum single atoms for enhanced acidic oxygen evolution reaction. *Advanced Energy Materials* **15**.6 (2025): 2403388.  

34. Zhang, Di, and Hao Li. Digital catalysis platform (DigCat): a gateway to big data and AI-powered innovations in catalysis. *ChemRxiv* (2024).  

35. Cowley, John M. An approximate theory of order in alloys. *Physical Review* **77**.5 (1950): 669.  

36. Cowley, J. M. Short-and long-range order parameters in disordered solid solutions. *Physical Review* **120**.5 (1960): 1648.  

37. Cowley, J. M. Short-range order and long-range order parameters. *Physical Review* **138**.5A (1965): A1384.  

38. Metropolis, Nicholas, and Stanislaw Ulam. The Monte Carlo method. *Journal of the American Statistical Association* **44**.247 (1949): 335–341.  

39. Kennedy, James, and Russell Eberhart. Particle swarm optimization. *Proceedings of ICNN'95-International Conference on Neural Networks.* Vol. 4. IEEE, 1995.  

40. Aurenhammer, Franz. Voronoi diagrams—a survey of a fundamental geometric data structure. *ACM Computing Surveys* **23**.3 (1991): 345–405.  

41. Okabe, Atsuyuki, et al. *Spatial Tessellations: Concepts and Applications of Voronoi Diagrams.* (2009).  

42. Du, Yadong, et al. Engineering efficient Self-Supporting SnO₂/Carbon nanofiber electrode for electrochemical CO₂ reduction. *Chemical Engineering Science* (2025): 122147.  

43. Shi, Kaige, et al. Effect of Ammonium-Based Cations on CO₂ Electroreduction. *ACS Catalysis* **15**.5 (2025): 3647–3659.  

44. Zhu, Jiaye, et al. Vanadium Oxide Clusters Mediated Bismuth‐Tin Alloy for Accelerated Dynamics of Electrocatalytic CO₂ Conversion. *Advanced Functional Materials* **35**.16 (2025): 2420177.  

45. Wang, Xiaowen, et al. Steering geometric reconstruction of bismuth with accelerated dynamics for CO₂ electroreduction. *Angewandte Chemie* **136**.34 (2024): e202407665.  

46. Wang, Yanlin, et al. Boosting Electrochemical CO₂ Reduction to Formate over La-Doped SnO₂ via Pinning Effect and Water Activation. *Journal of the American Chemical Society* (2025).  

47. Cantor, Brian, et al. Microstructural development in equiatomic multicomponent alloys. *Materials Science and Engineering: A* **375** (2004): 213–218.  

48. Yeh, J‐W., et al. Nanostructured high‐entropy alloys with multiple principal elements: novel alloy design concepts and outcomes. *Advanced Engineering Materials* **6**.5 (2004): 299–303.  

49. Dey, Gaurav R., et al. Chemical insights into the formation of colloidal high entropy alloy nanoparticles. *ACS Nano* **17**.6 (2023): 5943–5955.  

50. Yao, Yonggang, et al. High-entropy nanoparticles: Synthesis-structure-property relationships and data-driven discovery. *Science* **376**.6589 (2022): eabn3103.  

51. Li, Fang, et al. Highly stable two-dimensional bismuth metal-organic frameworks for efficient electrochemical reduction of CO₂. *Applied Catalysis B: Environmental* **277** (2020): 119241.  

52. Han, Na, et al. Ultrathin bismuth nanosheets from in situ topotactic transformation for selective electrocatalytic CO₂ reduction to formate. *Nature Communications* **9**.1 (2018): 1320.  

53. Sabouhanian, Negar, Jacek Lipkowski, and Aicheng Chen. Unveiling the potential of bismuth-based catalysts for electrochemical CO₂ reduction. *Industrial Chemistry & Materials* (2025).  

54. Yao, Zhengjie, Zhenjie Cheng, and Jiacheng Wang. Heteroatom‐Doped Bismuth‐Based Electrocatalysts for CO₂ Reduction to Formic Acid: Advancement and Perspective. *cMat* **2**.3 (2025): e70016.  

55. Luo, Yuqing, et al. Bismuth-Catalyzed Electrochemical Carbon Dioxide Reduction to Formic Acid: Material Innovation and Reactor Design. *Accounts of Materials Research* **6**.4 (2025): 462–472.  

56. Feaster, Jeremy T., et al. Understanding selectivity for the electrochemical reduction of carbon dioxide to formic acid and carbon monoxide on metal electrodes. *ACS Catalysis* **7**.7 (2017): 4822–4827.  

57. Kortlever, Ruud, et al. Catalysts and reaction pathways for the electrochemical reduction of carbon dioxide. *The Journal of Physical Chemistry Letters* **6**.20 (2015): 4073–4082.  

58. Hori, Yoshio. Electrochemical CO₂ reduction on metal electrodes. *Modern Aspects of Electrochemistry* (2008): 89–189.  

59. Feaster, Jeremy T., et al. Understanding selectivity for the electrochemical reduction of carbon dioxide to formic acid and carbon monoxide on metal electrodes. *ACS Catalysis* **7**.7 (2017): 4822–4827.  

60. Löffler, Tobias, et al. What makes high‐entropy alloys exceptional electrocatalysts? *Angewandte Chemie International Edition* **60**.52 (2021): 26894–26903.  

61. Batchelor, Thomas A.A., et al. High-entropy alloys as a discovery platform for electrocatalysis. *Joule* **3**.3 (2019): 834–845.  

62. Kresse, Georg, and Jürgen Furthmüller. Efficiency of ab-initio total energy calculations for metals and semiconductors using a plane-wave basis set. *Computational Materials Science* **6**.1 (1996): 15–50.  

63. Deng, Bowen, et al. CHGNet as a pretrained universal neural network potential for charge-informed atomistic modelling. *Nature Machine Intelligence* **5**.9 (2023): 1031–1041.  

64. Perdew, John P., Kieron Burke, and Matthias Ernzerhof. Generalized gradient approximation made simple. *Physical Review Letters* **77**.18 (1996): 3865.  

65. Hammer, B. H. L. B., Lars Bruno Hansen, and Jens Kehlet Nørskov. Improved adsorption energetics within density-functional theory using revised Perdew-Burke-Ernzerhof functionals. *Physical Review B* **59**.11 (1999): 7413.  

66. Araujo, Rafael B., et al. Adsorption energies on transition metal surfaces: towards an accurate and balanced description. *Nature Communications* **13**.1 (2022): 6853.  

67. Wellendorff, Jess, et al. A benchmark database for adsorption bond energies to transition metal surfaces and comparison to selected DFT functionals. *Surface Science* **640** (2015): 36–44.  

68. Kelly, Sara R., et al. Electric field effects in oxygen reduction kinetics: rationalizing pH dependence at the Pt(111), Au(111), and Au(100) electrodes. *The Journal of Physical Chemistry C* **124**.27 (2020): 14581–14591.  

69. Hohenberg, Pierre, and Walter Kohn. Inhomogeneous electron gas. *Physical Review* **136**.3B (1964): B864.  

70. Kohn, Walter, and Lu Jeu Sham. Self-consistent equations including exchange and correlation effects. *Physical Review* **140**.4A (1965): A1133.  

71. Wang, Yuhang, et al. Bridging Theory and Experiment: Machine Learning Potential‐Driven Insights into pH‐Dependent CO₂ Reduction on Sn‐Based Catalysts. *Advanced Functional Materials* **35**.36 (2025): e06314.  

72. Chu, Yue, et al. Data-driven discovery of single-atom catalysts for CO₂ reduction considering the pH-dependency at the reversible hydrogen electrode scale. *The Journal of Chemical Physics* **162**.17 (2025).  

73. Zhang, Di, and Hao Li. The potential of zero charge and solvation effects on single-atom M–N–C catalysts for oxygen electrocatalysis. *Journal of Materials Chemistry A* **12**.23 (2024): 13742–13750.  

74. Burkhardt, Jordan, Yinglu Jia, and Wan-Lu Li. Structure Search with the Strategic Escape Algorithm. *Journal of Chemical Theory and Computation* **21**.7 (2025): 3765–3773.  

75. Grasemann, Martin, and Gábor Laurenczy. Formic acid as a hydrogen source–recent developments and future trends. *Energy & Environmental Science* **5**.8 (2012): 8171–8181.  

76. Eppinger, Jorg, and Kuo-Wei Huang. Formic acid as a hydrogen energy carrier. *ACS Energy Letters* **2**.1 (2017): 188–195.  
