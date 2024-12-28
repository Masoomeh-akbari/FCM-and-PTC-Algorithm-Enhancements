# Fuzzy Cognitive Maps and Probabilistic Transitive Closure: Algorithm Enhancements

This repository contains Python implementations of advanced algorithms designed to **enhance the computation of probabilistic transitive closure (PTC)** for **bipolar weighted digraphs**. These enhancements, based on the Reduction-Recovery Algorithm, improve computational efficiency when analyzing indirect impacts in complex systems modeled by **Fuzzy Cognitive Maps (FCMs)**.

To provide context, let’s briefly explore the foundational concepts of FCMs and PTC, highlighting their significance in understanding complex systems.

---

## What Are Fuzzy Cognitive Maps (FCMs)?

Fuzzy Cognitive Maps (FCMs) are graphical models that capture and represent relationships among various factors in complex systems. Each relationship is characterized by its weight (a value in the range [0,1]) and its sign (positive or negative). These models are widely applied in diverse fields such as:
- **Education** and **economics** (e.g., modeling determinants of outcomes),
- **Communication networks**, 
- **Management** and **policy planning**.

In graph theory, FCMs are represented as **bipolar weighted digraphs**, where:
- **Vertices** represent the factors in the system.
- **Arcs** denote direct impacts, assigned weights and signs to quantify their nature and strength.

FCMs simplify the representation of interdependencies, but understanding their full implications often requires the analysis of indirect relationships—this is where probabilistic transitive closure comes into play.

---

## What is Probabilistic Transitive Closure (PTC)?

The **probabilistic transitive closure (PTC)** extends the insights provided by FCMs by computing the **indirect impacts** of one factor on another, traversing multiple pathways. In the PTC framework:
- The weight of an arc represents the **probability** of an indirect impact.
- The sign of the arc indicates whether the net effect is positive or negative.

The result of this computation is a new bipolar weighted digraph, where each arc corresponds to an indirect relationship. This analysis is particularly critical for exploring the net effects of complex interactions in systems modeled by FCMs.

While existing methods for computing the PTC, such as **ProbTC**, are effective, they can be computationally intensive for large or complex graphs. This repository introduces enhancements that optimize the computation process through the Reduction-Recovery Algorithm.

---

## Algorithms in This Repository

This repository provides two implementations of the **Reduction-Recovery Algorithm**, an enhancement designed to optimize the computation of PTC for bipolar weighted digraphs. The Reduction-Recovery Algorithm simplifies the input graph, making computations faster and more practical:

1. **Reduction:** Iteratively removes vertices with limited connectivity (indegree or outdegree at most one).
2. **Recovery:** Uses the PTC of the simplified graph to reconstruct the PTC of the original graph.

### Implementations:
- **SRR-PTC (Stack-based):** A memory-efficient implementation using a stack for managing reductions.
- **RRR-PTC (Recursion-based):** A simpler, recursion-based implementation leveraging Python's built-in capabilities.

These approaches enhance the computational process, particularly for reducible graphs, while maintaining comparable performance on irreducible graphs.

---

## Advantages of These Enhancements

Compared to the existing **ProbTC** software, these implementations offer:
- **Faster computation for reducible digraphs.** 
- **Comparable performance for irreducible digraphs.** 

---

## Repository Contents

- **[`SRR_PTC.py`](./SRR_PTC.py) and [`RRR_PTC.py`](./RRR_PTC.py):** Python implementations of the Reduction-Recovery Algorithm.
- **Master's Thesis (PDF):** A comprehensive reference with mathematical derivations, pseudocode, and a detailed comparison of algorithm performance.

For an in-depth understanding of the algorithms, their theoretical foundations, and performance evaluations, readers are encouraged to refer to the **Master's Thesis** included in this repository.
