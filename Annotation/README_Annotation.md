# Cell Type Annotation for Bee Brain Spatial Data

A comprehensive toolkit for annotating cell types in enhanced spatial transcriptomics data using literature-curated marker genes.

## 🎯 Overview

This toolkit provides **best practices** for cell type annotation in spatial transcriptomics data, specifically designed for bee brain analysis but adaptable to other systems.

### Key Features

- **Marker Gene-Based Scoring**: Multiple scoring methods (mean, median, scanpy)
- **Confidence-Based Assignment**: Prevents over-confident assignments
- **Spatial Validation**: Ensures assignments make anatomical sense
- **Comprehensive Visualization**: UMAP, spatial, and statistical plots
- **Quality Control**: Built-in QC metrics and recommendations

## 📁 Files in this Directory

- `cell_type_annotation.py` - Main annotation toolkit class
- `annotation_workflow_example.py` - Complete example workflow
- `README.md` - This documentation file

## 🚀 Quick Start

### 1. Basic Usage

```python
from cell_type_annotation import BeebrainCellTypeAnnotator
import scanpy as sc

# Load your enhanced spatial data
adata = sc.read_h5ad("your_enhanced_data.h5ad")

# Initialize annotator with your marker genes
annotator = BeebrainCellTypeAnnotator(adata, marker_genes=your_markers)

# Check marker availability
annotator.check_marker_availability()

# Score and assign cell types
annotator.score_cell_types(method='mean_expression')
annotator.assign_cell_types(confidence_threshold=0.1)

# Visualize results
annotator.plot_comprehensive_annotation_results()

# Save annotated data
annotator.save_annotations("annotated_data.h5ad")
```

### 2. Run Complete Workflow

```bash
cd /path/to/your/p5_SvsF/code/annotation
python annotation_workflow_example.py
```

## 📚 Best Practices for Cell Type Annotation

### **1. Marker Gene Selection**

#### Literature Review Strategy:
- **Start broad, go specific**: Begin with pan-neuronal markers, then refine to subtypes
- **Multiple sources**: Cross-reference multiple papers and databases
- **Validation studies**: Prefer markers validated in multiple studies
- **Species consideration**: Prioritize bee/insect studies over Drosophila when possible

#### Quality Criteria:
- **Specificity**: Genes expressed in one cell type but not others
- **Sensitivity**: Expressed in most/all cells of that type
- **Reliability**: Consistent expression across conditions
- **Availability**: Present in your dataset (check with `check_marker_availability()`)

### **2. Scoring Methods**

| Method | Best For | Pros | Cons |
|--------|----------|------|------|
| `mean_expression` | Initial screening | Fast, interpretable | Sensitive to outliers |
| `median_expression` | Robust scoring | Outlier-resistant | Less sensitive |
| `scanpy_score` | Publication-quality | Normalized, comparable | More complex |

**Recommendation**: Start with `mean_expression`, validate with `scanpy_score`

### **3. Assignment Strategies**

#### Confidence Thresholds:
- **Conservative (0.2-0.5)**: Fewer assignments, higher confidence
- **Moderate (0.1-0.2)**: Balanced assignments
- **Liberal (0.05-0.1)**: More assignments, lower confidence

#### Quality Control:
```python
# Check unassigned percentage
unassigned_pct = (adata.obs['cell_type_assignment'] == 'Unassigned').sum() / len(adata) * 100

if unassigned_pct > 50:
    print("Consider lowering confidence threshold or adding more markers")
elif unassigned_pct < 5:
    print("Consider raising confidence threshold to avoid over-assignment")
```

### **4. Spatial Validation**

#### Anatomical Coherence:
- **Regional enrichment**: Cell types should cluster in expected brain regions
- **Spatial continuity**: Similar cell types should be spatially adjacent
- **Known anatomy**: Validate against established bee brain atlases

#### Validation Questions:
- Do Kenyon cells cluster in the mushroom body region?
- Are olfactory neurons enriched in antennal lobe areas?
- Do visual system markers appear in optic lobe regions?

### **5. Hierarchical Annotation Strategy**

```python
# Level 1: Broad categories
broad_markers = {
    'Neurons': ['elav', 'nSyb', 'Syt1'],
    'Glia': ['repo', 'Gs2', 'Eaat1'],
    'Hemocytes': ['Hml', 'srp', 'gcm']
}

# Level 2: Neuronal subtypes  
neuron_markers = {
    'Kenyon_Cells': ['mb247', 'FasII', 'OK107'],
    'Projection_Neurons': ['GH146', 'Mz19', 'NP225'],
    'Local_Neurons': ['Gad1', 'VGlut', 'ChAT']
}

# Level 3: Specific subtypes
KC_markers = {
    'sKC': ['mb247', 'FasII', 'arm', 'Lac'],
    'lKC': ['mb247', 'FasII', 'VGlut', 'ChAT'],
    'mKC': ['mb247', 'FasII', 'Synapsin', 'VGlut']
}
```

## 🔬 Scientific Workflow

### **Step 1: Data Quality Assessment**
```python
# Check data characteristics
annotator._validate_data()

# Assess marker availability
availability = annotator.check_marker_availability()
```

### **Step 2: Initial Annotation**
```python
# Score with multiple methods
annotator.score_cell_types(method='mean_expression', score_name_suffix='_mean')
annotator.score_cell_types(method='scanpy_score', score_name_suffix='_scanpy')

# Conservative assignment
annotator.assign_cell_types(confidence_threshold=0.2)
```

### **Step 3: Quality Control**
```python
# Generate comprehensive plots
annotator.plot_comprehensive_annotation_results()

# Check spatial coherence
# (Manual inspection of spatial plots)

# Review assignment summary
summary = annotator.get_annotation_summary()
print(summary)
```

### **Step 4: Refinement**
```python
# Adjust parameters based on QC
annotator.assign_cell_types(confidence_threshold=0.15)  # Adjusted threshold

# Re-evaluate results
annotator.plot_spatial_annotations()
```

### **Step 5: Validation**
- Compare with existing clustering results
- Validate against known bee brain anatomy
- Check for biological plausibility
- Cross-reference with literature

## 📊 Output Files

After running the annotation workflow, you'll get:

```
annotation_results/
├── score_distributions.png           # Distribution of cell type scores
├── umap_cell_types.png              # UMAP with cell type annotations
├── spatial_cell_types.png           # Spatial plot with annotations
├── comprehensive_annotation_results.png  # 12-panel summary figure
├── annotation_summary.csv           # Cell type statistics
├── bee_brain_markers_detailed.csv   # Marker gene template
└── f11_s31_col116_annotated.h5ad   # Annotated spatial data
```

## 🧬 Bee Brain Cell Types

### Major Categories:

1. **Neurons** (80-90% of brain cells)
   - Kenyon cells (mushroom body)
   - Projection neurons (olfactory)
   - Local neurons (various regions)
   - Motor neurons
   - Interneurons

2. **Glia** (5-15% of brain cells)
   - Astrocyte-like glia
   - Ensheathing glia
   - Cortex glia

3. **Hemocytes** (<5% of brain cells)
   - Immune cells
   - Phagocytic cells

### Specific Subtypes:

#### Mushroom Body (Learning & Memory):
- **sKC**: Small Kenyon cells (γ neurons)
- **lKC**: Large Kenyon cells (α/β neurons)  
- **mKC**: Medium Kenyon cells (α'/β' neurons)

#### Olfactory System:
- **OLC**: Olfactory lobe cells (antennal lobe)
- **PN**: Projection neurons (olfactory pathway)
- **LN**: Local neurons (lateral processing)

#### Visual System:
- **Photoreceptors**: R1-R8 cells
- **Lamina neurons**: L1-L5 cells
- **Medulla neurons**: Various types

## ⚠️ Important Considerations

### **Data-Specific Factors:**

1. **Gene Nomenclature**: Ensure consistent gene naming between your data and markers
2. **Expression Levels**: Log-transformed data works best (values 0-10 range)
3. **Gene Coverage**: Enhanced data should have 10,000+ genes for best results
4. **Spatial Resolution**: Cell-level resolution required for accurate annotation

### **Biological Considerations:**

1. **Developmental Stage**: Marker expression varies with bee age/caste
2. **Brain Region**: Not all markers are expressed in all brain regions
3. **Caste Differences**: Soldier vs. forager brains may have different cell type proportions
4. **Technical Artifacts**: Ensure assignments aren't driven by batch effects

### **Statistical Considerations:**

1. **Multiple Testing**: Consider correction for multiple cell type tests
2. **Confidence Intervals**: Report confidence ranges for assignments
3. **Validation Set**: Hold out data for independent validation
4. **Cross-Validation**: Test marker robustness across samples

## 🔧 Troubleshooting

### Common Issues:

**Problem**: High percentage of unassigned cells (>50%)
**Solutions**:
- Lower confidence threshold
- Add more marker genes for underrepresented cell types
- Check marker gene availability
- Verify data preprocessing

**Problem**: All cells assigned to one type
**Solutions**:
- Increase confidence threshold
- Check for batch effects in data
- Verify marker gene specificity
- Examine score distributions

**Problem**: Assignments don't match spatial anatomy
**Solutions**:
- Review marker gene literature
- Check for contamination in marker lists
- Validate preprocessing steps
- Consider hierarchical annotation

**Problem**: Low marker gene availability
**Solutions**:
- Update gene names/symbols
- Check case sensitivity
- Cross-reference with gene databases
- Consider ortholog mapping

## 📖 Literature Resources

### Key Databases:
- **FlyBase**: Drosophila gene information
- **UniProt**: Protein sequences and functions
- **NCBI Gene**: Gene annotations and orthologs
- **BeeBase**: Honey bee genomic resources

### Important Papers:
- Aso et al. (2014) - Mushroom body connectivity
- Couto et al. (2005) - Olfactory receptor neurons
- Wong et al. (2002) - Projection neuron anatomy
- Awasaki et al. (2008) - Glial cell types

## 🤝 Contributing

To improve this annotation toolkit:

1. **Add new marker genes** from recent literature
2. **Improve scoring methods** for better accuracy
3. **Add spatial validation metrics** for quality control
4. **Extend to other brain regions** beyond current coverage

## 📞 Support

For questions about cell type annotation:

1. **Check the troubleshooting section** above
2. **Review the example workflow** for proper usage
3. **Validate marker genes** against literature
4. **Cross-reference** with existing clustering results

Remember: **Cell type annotation is as much art as science** - use multiple approaches and always validate results against known biology! 