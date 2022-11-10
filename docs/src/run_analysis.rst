Steps required to run an analysis with Slim-TPCA

=============================================================================================
STEP 1: Train final model using optimized hyperparameters and evaluate predictive performance
=============================================================================================

After optimizing the kernel type and optimal number of inducing points, final models are built and trained for each set of conditions that you want to compare (e.g. control vs treatment 1, control vs. treatment 2, etc.).
Following model training, the final predictive performance of the model can be evaluated using the held-out test partition of the dataset.

See :doc:`this notebook <./notebooks/Slim_TPCA_test>` for a detailed workflow.

==============================
