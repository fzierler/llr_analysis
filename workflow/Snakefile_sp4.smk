include: "rules.smk"


metadata = pd.read_csv("metadata/runs_sp4.csv")


rule all_sp4:
    input:
        overview=expand(
            "tmp/sp4/plots/overview/{dataset}.pdf",
            dataset=list_datasets(metadata),
        ),
        an_trajectory=expand(
            "tmp/sp4/plots/an_trajectories/{dataset}.pdf",
            dataset=list_datasets(metadata),
        ),
        free_energy=expand(
            "assets/sp4/plots/free_energy/free_energy_{dataset}.pdf",
            dataset=list_datasets(metadata),
        ),
        an_volume=expand("assets/sp4/plots/an_volume_Nt{Nt}.pdf", Nt=[4, 5]),
        entropy=expand("assets/sp4/plots/entropy_Nt{Nt}.pdf", Nt=[4, 5]),
        critical_beta_volume=expand(
            "assets/sp4/plots/critical_beta_ratios_Nt{Nt}.pdf", Nt=[4, 5]
        ),
        double_gaussian_volume=expand(
            "assets/sp4/plots/plaquette_distribution_Nt{Nt}_volumes.pdf", Nt=[4, 5]
        ),
        double_gaussian=expand(
            "assets/sp4/plots/plaquette_distribution/pd_{dataset}_{ratio}.pdf",
            dataset=list_datasets(metadata),
            ratio=["1:1", "2:1", "1:2"],
        ),
        an_replicas=expand(
            "assets/sp4/plots/an_replicas_Nt{Nt}_Ns{Ns}.pdf",
            zip,
            Nt=[5, 5],
            Ns=[48, 56],
        ),
        polyakov_loop_hist=expand(
            "assets/sp4/plots/polyakov_loop_hist/polyakov_{dataset}.pdf",
            dataset=list_datasets_poly(metadata),
        ),
        polyakov_loop_compact=expand(
            "assets/sp4/plots/polyakov_loop_compact/polyakov_{dataset}.pdf",
            dataset=list_datasets_poly(metadata),
        ),
        polyakov_loop_vs_beta="data_assets/sp4/polyakov/polyakov_loop_data.hdf5",
        polyakov_loop_susceptibility=expand("assets/sp4/plots/polyakov_susceptibility_Nt{Nt}.pdf", Nt=[4, 5]),
        binder_cumulant=expand("assets/sp4/plots/binder_cumulant_Nt{Nt}.pdf", Nt=[4, 5]),
        specific_heat=expand("assets/sp4/plots/specific_heat_Nt{Nt}.pdf", Nt=[4, 5]),
        surface_tension="assets/sp4/plots/surface_tension_term.pdf",
        ensembles_combined="assets/sp4/tables/runs.tex",
        file_out="assets/sp4/tables/beta_critical.tex",
        out_cumulants="data_assets/sp4/critical_beta_cumulants.csv",
        out_critical11="data_assets/sp4/critical_beta_1:1.csv",
        out_critical12="data_assets/sp4/critical_beta_2:1.csv",
        out_critical21="data_assets/sp4/critical_beta_1:2.csv",
        supercooling_csv="data_assets/sp4/supercooling_param.csv",
        karsch_out="data_assets/sp4/scale_setting.csv",
        latent_heat_plot="assets/plots/sp4/latent_heat.pdf",

