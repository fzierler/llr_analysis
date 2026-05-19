include: "rules.smk"


metadata_file = "metadata/runs_su3.csv"
metadata = pd.read_csv(metadata_file)


rule all_su3:
    input:
        metadata_file,
        overview_plots=expand(
            "tmp/su3/plots/overview/{dataset}.pdf",
            dataset=list_datasets(metadata),
        ),
        an_trajectory_plots=expand(
            "tmp/su3/plots/an_trajectories/{dataset}.pdf",
            dataset=list_datasets(metadata),
        ),
        free_energy_plots=expand(
            "assets/su3/plots/free_energy/free_energy_{dataset}.pdf",
            dataset=list_datasets(metadata),
        ),
        an_volume_plots="assets/su3/plots/an_volume_Nt4.pdf",
        entropy_plots="assets/su3/plots/entropy_Nt4.pdf",
        critical_beta_volume_plots="assets/su3/plots/critical_beta_ratios_Nt4.pdf",
        double_gaussian_volume_plots="assets/su3/plots/plaquette_distribution_Nt4_volumes.pdf",
        double_gaussian_plots=expand(
            "assets/su3/plots/plaquette_distribution/pd_{dataset}_1:1.pdf",
            dataset=list_datasets(metadata),
        ),
        plot_binder_cumulant="assets/su3/plots/binder_cumulant_Nt4.pdf",
        plot_specific_heat="assets/su3/plots/specific_heat_Nt4.pdf",
        surface_tension_plot="assets/su3/plots/surface_tension_term.pdf",
        ensembles_combined="assets/su3/tables/runs.tex",
        file_out="assets/su3/tables/beta_critical.tex",
        out_cumulants="data_assets/su3/critical_beta_cumulants.csv",
        out_critical="data_assets/su3/critical_beta_1:1.csv",
        out_ratio_21="data_assets/su3/critical_beta_2:1.csv",
        out_ratio_12="data_assets/su3/critical_beta_1:2.csv",
        supercooling_csv="data_assets/su3/supercooling_param.csv",
        polyakov_loop_hist=expand(
            "assets/su3/plots/polyakov_loop_hist/polyakov_{dataset}.pdf",
            dataset=list_datasets_poly(metadata),
        ),
        polyakov_loop_vs_beta="data_assets/su3/polyakov/polyakov_loop_data.hdf5",
