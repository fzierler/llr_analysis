include: "rules.smk"


metadata = pd.read_csv("metadata/runs_su4.csv")


rule all_su4:
    input:
        overview_plots=expand(
            "assets/su4/plots/overview/{dataset}.pdf",
            dataset=list_datasets(metadata),
        ),
        an_trajectory_plots=expand(
            "assets/su4/plots/an_trajectories/{dataset}.pdf",
            dataset=list_datasets(metadata),
        ),
        #free_energy_plots=expand(
        #    "assets/su4/plots/free_energy/free_energy_{dataset}.pdf",
        #    dataset=list_datasets(metadata),
        #),
        an_volume_plots=expand("assets/su4/plots/an_volume_Nt{Nt}.pdf", Nt=[5, 6, 7, 8]),
        entropy_plots=expand("assets/su4/plots/entropy_Nt{Nt}.pdf", Nt=[5, 6]),
        #critical_beta_ratio_plots=expand(
        #    "assets/su4/plots/critical_beta_ratios_Nt{Nt}.pdf", Nt=[5, 6]
        #),
        critical_beta_volume_plots=expand(
            "assets/su4/plots/beta_critical_Nt{Nt}.pdf", Nt=[5, 6, 7, 8]
        ),
        double_gaussian_volume_plots=expand(
            "assets/su4/plots/plaquette_distribution_Nt{Nt}_volumes.pdf", Nt=[5, 6]
        ),
        #double_gaussian_plots=expand(
        #    "assets/su4/plots/plaquette_distribution/pd_{dataset}_1:1.pdf",
        #    dataset=list_datasets(metadata),
        #),
        plot_binder_cumulant=expand(
            "assets/su4/plots/binder_cumulant_Nt{Nt}.pdf", Nt=[5, 6, 7, 8]
        ),
        plot_specific_heat=expand(
            "assets/su4/plots/specific_heat_Nt{Nt}.pdf", Nt=[5, 6, 7, 8]
        ),
        #surface_tension_plot="assets/su4/plots/surface_tension_term.pdf",
        #ensembles_combined="assets/su4/tables/runs.tex",
        #file_out="assets/su4/tables/beta_critical.tex",
        supercooling_csv="data_assets/su4/supercooling_param.csv",
        an_replicas=expand(
            "assets/su4/plots/an_replicas_Nt{Nt}_Ns{Ns}.pdf", zip, Nt=[5], Ns=[28]
        ),
        plot_traj="assets/su4/plots/energies/energy_trajectory_5x32_96replicas.pdf",
        plot_hist="assets/su4/plots/energies/energy_histogram_5x32_96replicas.pdf",

