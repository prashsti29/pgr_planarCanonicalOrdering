# Planar Canonical Ordering- Sample Data Analysis

Planar canonical ordering establishes an incremental vertex ordering for planar graphs by repeatedly selecting and removing vertices along the outer face while preserving planarity and connectivity. Reversing this process yields a canonical sequence that supports linear-time planar embedding and drawing algorithms.
The code builds an undirected graph from pgRouting-style edge data, runs Boyer-Myrvold planarity testing, makes the graph biconnected and maximally planar as prerequisites, then applies the canonical ordering to produce a structured vertex sequence with prior neighbor counts per node.

Output for the Planar Canonical Ordering for the generated graph:
<img width="917" height="516" alt="image" src="https://github.com/user-attachments/assets/c518d918-c5fe-4535-a165-5d6c8334b94f" />
