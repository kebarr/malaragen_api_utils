python generate_plot_pca_csv_multi_rows.py
python ../papermill_meta.py plot_pca.ipynb  --notebooks-dir ./notebooks   --output-dir generated_runners
python generated_runners/papermill_plot_pca.py --config perf_testing_notebooks/plot_pca_config_multi_row.csv --output-dir executed_runs

