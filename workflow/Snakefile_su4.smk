include: "rules.smk"


metadata = pd.read_csv("metadata/runs_su4.csv")


rule all:
    input:
        overview_plots=expand(
            "tmp/su4/plots/overview/{dataset}.pdf",
            dataset=list_datasets(metadata),
        ),
        an_trajectory_plots=expand(
            "tmp/su4/plots/an_trajectories/{dataset}.pdf",
            dataset=list_datasets(metadata),
        ),
        free_energy_plots=expand(
            "assets/su4/plots/free_energy/free_energy_{dataset}.pdf",
            dataset=list_datasets(metadata),
        ),
        an_volume_plots=expand("assets/su4/plots/an_volume_Nt{Nt}.pdf", Nt=[5, 6]),
        entropy_plots=expand("assets/su4/plots/entropy_Nt{Nt}.pdf", Nt=[5, 6]),
        critical_beta_volume_plots=expand(
            "assets/su4/plots/critical_beta_volumes_Nt{Nt}.pdf", Nt=[5, 6]
        ),
        double_gaussian_volume_plots=expand(
            "assets/su4/plots/plaquette_distribution_Nt{Nt}_volumes.pdf", Nt=[5, 6]
        ),
        double_gaussian_plots=expand(
            "assets/su4/plots/plaquette_distribution/pd_{dataset}_1:1.pdf",
            dataset=list_datasets(metadata),
        ),
        plot_binder_cumulant=expand(
            "assets/su4/plots/binder_cumulant_Nt{Nt}.pdf", Nt=[5, 6]
        ),
        plot_specific_heat=expand(
            "assets/su4/plots/specific_heat_Nt{Nt}.pdf", Nt=[5, 6]
        ),
        surface_tension_plot="assets/su4/plots/surface_tension_term.pdf",
        ensembles_combined="assets/su4/tables/runs.tex",
        file_out="assets/su4/tables/beta_critical.tex",
        out_cumulants="data_assets/su4/critical_beta_cumulants.csv",
        out_critical="data_assets/su4/critical_beta_1:1.csv",
        out_ratio_21="data_assets/su4/critical_beta_1:2.csv",
        out_ratio_12="data_assets/su4/critical_beta_2:1.csv",
