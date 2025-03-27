# T1D ML Project

Type 1 diabetes (T1D) is a chronic autoimmune disease in which the immune system destroys insulin-producing beta cells in the pancreas. The primary treatment is insulin injections, which do not prevent severe complications. A new approach aims to target autoreactive T cells, the specificity of which is defined by the T-cell receptor (TCR). In this project, we classify the TCR repertoires of patients with T1D and healthy individuals and search for T1D-associated TCR clonotypes (features that have the biggest impact on the model). The results can potentially be applied for diagnostics and new therapy development. The project implements ML methods for published and our own unpublished data of patients with T1D and healthy controls.

More details about the project you can find in the [report](report.pdf) and [presentation](presentation.pdf)

# Repository structure

The project has three parts:

**1. ML models trained on features selected by F-test or TCRnet.**

:file_folder: Feature tables for training and test dataset constructed on features selected by fisher test or TCRnet are stored in `feature_tables` folder. 

:clipboard: The notebook with classifiers based on these feature tables `ml_models_on_features_selected_by_fisher_or_TCRnet_30k_norm.ipynb`.

:herb: yml file to create conda env to execute the notebook `env_for_ml_with_classical_feature_engineering.yml`

**2. ML models trained on embeddings.**

:file_folder: TCR repertoires in form suitable for embeddings construction are stored in `data_30k_umi/repertoires`. Metadata for our files is `data_30k_umi/metadata`

:clipboard: Script for embeddings construction `esm_embed_script.py` 

:herb: yml file to create conda env to execute the script `env_for_ml_on_embeddings.yml`

:clipboard: Notebook for ML models on embeddings `ml_analysis_on_embeddings.ipynb`

:herb: yml file to create conda env to execute the notebook `env_for_ml_on_embeddings.yml`

**3. DeepRC model implementation**

:file_folder: `DeepRC-master` is a modified copy of [original DeepRC github repo](https://github.com/ml-jku/DeepRC). Our TCR reperotires in the form suitable for DeepRC are stored in `DeepRC-master/deeprc/datasets/example_dataset/repertoires`. Metadata for our files is `DeepRC-master/deeprc/datasets/example_dataset/metadata.tsv`. 

:clipboard: The notebook for our model is `DeepRC-master/deeprc/example_single_task_cnn.py`.

:herb: yml file to create conda env to execute the notebook `deeprc.yaml`

# How to execute the code

All requirements are stored in yml files. It is better to create new environment for each part of the projects by running

```
conda env create -f environment.yml
# environment activation
conda activate name_of_the_env
# than you can run the notebook of python script in created environment
```