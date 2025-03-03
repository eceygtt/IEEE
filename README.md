# shiny for python
This Shiny application allows users to upload a `.zip` file containing single-cell RNA-seq data in `.h5ad` format and visualize the data using UMAP (Uniform Manifold Approximation and Projection). The app processes the data, performs dimensionality reduction, and generates a UMAP plot to help users explore cell clusters and their relationships.

Features
- Upload ZIP File: Users can upload a `.zip` file containing `.h5ad` data.
- UMAP Visualization: The app automatically processes the data and generates a UMAP plot.
- Interactive Interface: The Shiny interface makes it easy to upload data and view results.

Benefits
- User-Friendly: No coding required; simply upload your data and view the results.
- Quick Insights: Quickly visualize cell clusters and identify patterns in single-cell RNA-seq data.
- Customizable: The UMAP plot can be customized to highlight different cell features.

How to Use
1. Clone this repository.
2. Install the required Python packages:
   ```bash
   pip install shiny scanpy matplotlib and then
   run the app with "shiny run" command :)

<img width="430" alt="Ekran Resmi 2025-03-03 17 18 51" src="https://github.com/user-attachments/assets/03a1de76-c73d-400a-ae83-45355a6ab24e" />



