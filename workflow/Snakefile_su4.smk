include: "rules.smk"


metadata = pd.read_csv("metadata/runs_su4.csv")
Nts = [5]

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
        free_energy_plots=expand(
            "assets/su4/plots/free_energy/free_energy_{dataset}.pdf",
            dataset=list_datasets(metadata),
        ),
        an_volume_plots=expand("assets/su4/plots/an_volume_Nt{Nt}.pdf", Nt=Nts),
        entropy_plots=expand("assets/su4/plots/entropy_Nt{Nt}.pdf", Nt=Nts),
        critical_beta_ratio_plots=expand(
            "assets/su4/plots/critical_beta_ratios_Nt{Nt}.pdf", Nt=Nts),
        critical_beta_volume_plots=expand(
            "assets/su4/plots/beta_critical_Nt{Nt}.pdf", Nt=Nts),
        double_gaussian_volume_plots=expand(
            "assets/su4/plots/plaquette_distribution_Nt{Nt}_volumes.pdf",Nt=Nts),
        plot_binder_cumulant=expand(
            "assets/su4/plots/binder_cumulant_Nt{Nt}.pdf", Nt=Nts),
        plot_specific_heat=expand(
            "assets/su4/plots/specific_heat_Nt{Nt}.pdf", Nt=Nts),
        supercooling_csv="data_assets/su4/supercooling_param.csv",
        an_replicas=expand(
            "assets/su4/plots/an_replicas_Nt{Nt}_Ns{Ns}.pdf", zip,
            Nt=[ 5, 5],
            Ns=[28,32],
        ),
        plot_traj_hist=expand(
            "assets/su4/plots/{name}/energy_{plot}_{dataset}.pdf",
            dataset=list_datasets(metadata),
            plot=("trajectory","histogram"),
            name=("energies","energies_extra_therm")
        ),

