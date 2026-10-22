include: "rules.smk"


metadata = pd.read_csv("metadata/runs_su4.csv")
Nts = [5]

rule all_su4:
    input:
        plots_per_run=expand([
            "assets/su4/plots/overview/{dataset}.pdf",
            "assets/su4/plots/overview_zoom/{dataset}.pdf",
            "assets/su4/plots/an_trajectories/{dataset}.pdf",
            "assets/su4/plots/free_energy/free_energy_{dataset}.pdf",
            ], dataset=list_datasets(metadata),
        ),
        plot_per_Nt=expand([
            "assets/su4/plots/an_volume_Nt{Nt}.pdf",
            "assets/su4/plots/entropy_Nt{Nt}.pdf",
            "assets/su4/plots/critical_beta_ratios_Nt{Nt}.pdf",
            "assets/su4/plots/beta_critical_Nt{Nt}.pdf",
            "assets/su4/plots/plaquette_distribution_Nt{Nt}_volumes.pdf",
            "assets/su4/plots/binder_cumulant_Nt{Nt}.pdf",
            "assets/su4/plots/specific_heat_Nt{Nt}.pdf",
            ], Nt=Nts),
        plot_traj_hist=expand(
            "assets/su4/plots/{name}/energy_{plot}_{dataset}.pdf",
            dataset=list_datasets(metadata),
            plot=("trajectory","histogram"),
            name=("energies","energies_extra_therm")
        ),
        supercooling_csv="data_assets/su4/supercooling_param.csv",
        an_replicas=(
            "assets/su4/plots/an_replicas_Nt5_Ns28.pdf",
            "assets/su4/plots/an_replicas_Nt5_Ns32.pdf",
            "assets/su4/plots/an_Nt5_Ns32_Nrep128.pdf",
        ),
        tables = (
            "assets/su4/tables/runs.tex",
            "assets/su4/tables/beta_critical.tex",
        ),

